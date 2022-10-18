#!/usr/bin/env python3

import argparse
import sys
import time
from pathlib import Path
import numpy as np
from yaml import load
from io import StringIO
from contextlib import contextmanager, redirect_stdout

import pyslow5

import pyguppy_client_lib
from pyguppy_client_lib.pyclient import PyGuppyClient
from pyguppy_client_lib import helper_functions

from ._version import __version__


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

@contextmanager
def start_guppy_server_and_client(args, server_args):
    """
    Thanks to Alex Payne for basis of this code to appropriately handle the server and client connections
    https://gist.github.com/alexomics/043bb120c74161e5b93e1b68fb00206c

    Starts server and connects client
    TODO: allow a connection to existing guppy server
    """
    server_args.extend(["--log_path", args.log,
                        "--config", args.config,
                        # "--port", args.port,
                        # "--max_queued_reads", args.max_queued_reads,
                        # "--chunk_size", args.chunk_size,
                        ])
    params = {"priority": PyGuppyClient.high_priority}
    
    if args.call_mods:
        params["move_and_trace_enabled"] = True
    
    # if args.align_ref:
    #     server_args.extend(["--align_ref", args.align_ref])
    #     params["align_ref"] = args.align_ref

    # This function has it's own prints that may want to be suppressed
    with redirect_stdout(StringIO()) as fh:
        server, port = helper_functions.run_server(server_args, bin_path=args.guppy_bin)

    if port == "ERROR":
        raise RuntimeError("Server couldn't be started")

    if port.startswith("ipc"):
        address = "{}".format(port)
    else:
        address = "localhost:{}".format(port)
    client = PyGuppyClient(address=address, config=args.config, throttle=3)



    sys.stderr.write("Setting params...\n")
    client.set_params(params)
    sys.stderr.write("Connecting...\n")
    try:
        with client:
            yield client
    finally:
        server.terminate()



# def start_guppy_server(args):
#     """
#     Start guppy server in separate process
#     """
#     server_args = ["--log_path", args.log,
#                    "--config", args.config,
#                    "--port", args.port,
#                    "--use_tcp",
#                    "-x", args.device,
#                    "--max_queued_reads", args.max_queued_reads,
#                    "--chunk_size", args.chunk_size]
#     # for i in guppy_server_args:
#     #     server_args.append(i)
#     ret = helper_functions.run_server(server_args, bin_path=args.guppy_bin)
#     return ret


def calibration(digitisation, range):
    """
    input:
        digitisation: float
        range: float
    output:
        scale: float
    """
    return range / digitisation


# def write_fastq(OUT, header, seq, qscore):
#     """
#     crude but effective fastq writter
#     """
#     OUT.write("{}\n".format(header))
#     OUT.write("{}\n".format(seq))
#     OUT.write("+\n")
#     OUT.write("{}\n".format(qscore))

def sam_header(OUT, sep='\t'):
    """
    Format a string sam header.
    This is taken from Bonito by Chris Seymour at ONT.
    https://github.com/nanoporetech/bonito/blob/master/bonito/io.py#L103
    """
    HD = sep.join([
        '@HD',
        'VN:1.5',
        'SO:unknown',
    ])
    PG1 = sep.join([
        '@PG',
        'ID:basecaller',
        'PN:guppy',
        'VN:%s' % pyguppy_client_lib.__version__,
    ])
    PG2 = sep.join([
        '@PG',
        'ID:wrapper',
        'PN:buttery-eel',
        'VN:%s' % __version__,
        'CL:buttery-eel %s' % ' '.join(sys.argv[1:]),
        'DS:guppy wrapper',
    ])
    OUT.write("{}\n".format(HD))
    OUT.write("{}\n".format(PG1))
    OUT.write("{}\n".format(PG2))


def write_output(OUT, read_id, header, seq, qscore, SAM_OUT, read_qscore, sam="", mods=False):
    
    if SAM_OUT:
        if mods:
            OUT.write("{}\n".format(sam))
        else:
            OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tNM:i:0\tqs:i:{}\n".format(read_id, seq, qscore, read_qscore))
    else:
        OUT.write("{}\n".format(header))
        OUT.write("{}\n".format(seq))
        OUT.write("+\n")
        OUT.write("{}\n".format(qscore))


def submit_read(client, read):
    """
    submit a read to the basecall server
    """
    skipped = ""
    read_id = read['read_id']
    # calculate scale
    scale = calibration(read['digitisation'], read['range'])
    result = False
    tries = 0
    while not result:
        result = client.pass_read(
                helper_functions.package_read(
                    read_id=read_id,
                    raw_data=np.frombuffer(read['signal'], np.int16),
                    daq_offset=read['offset'],
                    daq_scaling=scale,
                )
            )
        if tries > 1:
            time.sleep(1)
        tries += 1
        if tries >= 5:
            if not result:
                sys.stderr.write("Skipped a read: {}\n".format(read_id))
                skipped = read_id
                break
    return result, skipped


def get_reads(client, OUT, SAM_OUT, mods, read_counter, qscore_cutoff):
    """
    Get the reads back from the basecall server after being basecalled
    bcalled object contains 1 or more called reads, which contain various data
    """
    SPLIT_PASS = False
    if qscore_cutoff:
        SPLIT_PASS = True
        qs_cutoff = float(qscore_cutoff)
    done = 0
    while done < read_counter:
        bcalled = client.get_completed_reads()
        if not bcalled:
            time.sleep(0.2)
            continue
        else:
            for call in bcalled:
                sam_record = ""
                done += 1
                if len(call) != 1:
                    # possible split reads?
                    sys.stderr.write("Call is longer than 1: {}\n".format(len(call)))
                read_id = call[0]['metadata']['read_id']
                read_qscore = call[0]['metadata']['mean_qscore']
                int_read_qscore = int(read_qscore)
                # @read_id runid=bf... sampleid=NA12878_SRE read=476 ch=38 start_time=2020-10-26T19:58:23Z model_version_id=2021-05-17_dna_r9.4.1_minion_96_29d8704b
                # model_version_id = get_model_info(args.config, args.guppy_bin)
                header = "@{} model_version_id={} mean_qscore={}".format(call[0]['metadata']['read_id'], call[0]['metadata']['model_version_id'], int_read_qscore)
                sequence = call[0]['datasets']['sequence']
                qscore = call[0]['datasets']['qstring']
                # when calling mods, can just output sam_record value
                # otherwise, write_output will handle unaligned sam with no mods
                if mods:
                    try:
                        sam_record = call[0]['metadata']['alignment_sam_record']
                    except:
                        # TODO: add warning that mods model not being used, and exit
                        sam_record = ""
                        mods = False
                if SPLIT_PASS:
                    if read_qscore >= qs_cutoff:
                        # pass
                        out = OUT[0]
                    else:
                        # fail
                        out = OUT[1]
                else:
                    out = OUT
                write_output(out, read_id, header, sequence, qscore, SAM_OUT, int_read_qscore, sam=sam_record, mods=mods)
    done = 0

# How we get data out of the model files if they are not provided by the metadata output

# def get_model_info(config, guppy_bin):
#     config = os.path.join(guppy_bin,"../data/", config)
#     model = ""
#     with open(config, 'r') as f:
#         for line in f:
#             line = line.strip('\n')
#             line = line.replace(" ", "")
#             line = line.split("=")
#             if len(line) > 1:
#                 if line[0] == "model_file":
#                     model = line[1]
#                     break
#
#     if len(model) == 0:
#         logger.warning("could not deduce model for fastq writing, falling back to default, model_version_id=conf")
#         model_version_id = config
#     else:
#         model_json = os.path.join(guppy_bin,"../data/", model)
#         with open(model_json, 'r') as f:
#             jdata = json.load(f)
#             model_version_id = jdata["version"]["id"]
#
#     return model_version_id


def main():
    # ==========================================================================
    # Software ARGS
    # ==========================================================================
    """
    Example:

    buttery-eel --guppy_bin /install/ont-guppy-6.1.3/bin --use_tcp --chunk_size 200 \
    --max_queued_reads 1000 -x "cuda:all" --config dna_r9.4.1_450bps_fast.cfg --port 5558 \
    -i /Data/test.blow5 -o /Data/test.fastq

    """

    VERSION = __version__

    parser = MyParser(description="buttery-eel - wrapping guppy for file agnostic basecalling",
    epilog="Citation:...",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    # Args for the wrapper, and then probably best to just have free form args for guppy
    parser.add_argument("-i", "--input", required=True,
                        help="input blow5 file for basecalling")
    parser.add_argument("-o", "--output", required=True,
                        help="output .fastq or unaligned .sam file to write")
    parser.add_argument("-g", "--guppy_bin", type=Path, required=True,
                        help="path to ont_guppy/bin folder")
    parser.add_argument("--config", default="dna_r9.4.1_450bps_fast.cfg", required=True,
                        help="basecalling model config")
    parser.add_argument("--call_mods", action="store_true",
                        help="output MM/ML tags for methylation - will output sam - use with appropriate mod config")
    parser.add_argument("-q", "--qscore", type=int,
                        help="A mean q-score to split fastq/sam files into pass/fail output")
    parser.add_argument("--slow5_threads", type=int, default=4,
                        help="Number of threads to use reading slow5 file")
    parser.add_argument("--slow5_batchsize", type=int, default=4000,
                        help="Number of reads to process at a time reading slow5")
    # Disabling alignment because sam file headers are painful and frankly out of scope. Just use minimap2.
    # parser.add_argument("-a", "--align_ref",
    #                     help="reference .mmi file. will output sam. (build with: minimap2 -x map-ont -d ref.mmi ref.fa )")
    # parser.add_argument("--port", default="5558",
    #                     help="port to use between server/client")
    parser.add_argument("--log", default="buttery_guppy_logs",
                        help="guppy log folder path")
    # parser.add_argument("--max_queued_reads", default="2000",
    #                     help="Number of reads to send to guppy server queue")
    # parser.add_argument("--chunk_size", default="2000",
    #                     help="signal chunk size, lower this for lower VRAM GPUs")
    # parser.add_argument("-x", "--device", default="auto",
    #                     help="Specify GPU device: 'auto', or 'cuda:<device_id>'")
    parser.add_argument("-v", "--version", action='version', version="buttery-eel - wraping guppy for file agnostic basecalling version: {}".format(VERSION),
                        help="Prints version")
    # parser.add_argument("--debug", action="store_true",
    #                     help="Set logging to debug mode")

    # args = parser.parse_args()
    # This collects known and unknown args to parse the server config options
    args, other_server_args = parser.parse_known_args()


    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    sys.stderr.write("\n")
    sys.stderr.write("               ~  buttery-eel - SLOW5 Guppy Basecalling  ~\n")
    sys.stderr.write("==========================================================================\n  ARGS\n==========================================================================\n")
    sys.stderr.write("args:\n {}\n{}\n".format(args, other_server_args))

    # guppy_server_args = None
    # guppy_client_args = None

    SAM_OUT = False
    # check version, as only 6.3.0+ will support MM/ML tags correctly
    if args.call_mods:
        check = True
        check_major = 6
        check_minor = 3
        major, minor, patch = [int(i) for i in pyguppy_client_lib.__version__.split(".")]
        if major < check_major:
            check = False
        elif major == check_major:
            if minor < check_minor:
                check = False
        sys.stderr.write("\n")
        sys.stderr.write("MOD CALLING VERSION CHECK: >6.3.0? {}\n".format(check))
        sys.stderr.write("\n")
        if not check:
            sys.stderr.write("ERROR: Please use guppy and ont-pyguppy-client-lib version 6.3.0 or higher for modification calling\n")
            sys.stderr.write("\n")
            sys.exit(1)

    # ==========================================================================
    # Start guppy_basecall_server
    # ==========================================================================
    sys.stderr.write("\n\n")
    sys.stderr.write("==========================================================================\n  Starting Guppy Basecalling Server\n==========================================================================\n")
    with start_guppy_server_and_client(args, other_server_args) as client:
        print(client)
        sys.stderr.write("guppy_basecall_server started...\n")
        sys.stderr.write("\n")


        # ==========================================================================
        # Connect to server with guppy_basecall_client
        # ==========================================================================

        # TODO: add guppy_client_args
        sys.stderr.write("==========================================================================\n  Connecting to server\n==========================================================================\n")
        sys.stderr.write("Connection status: ")
        sys.stderr.write("{}".format(client.get_status()))
        # print(client.get_barcode_kits("127.0.0.1:{}".format(args.port), 10))
        # print(client.get_protocol_version())
        # print(client.get_server_information("127.0.0.1:{}".format(args.port), 10))
        # print(client.get_software_version())

        sys.stderr.write("\n\n")

        # ==========================================================================
        # Read signal file
        # ==========================================================================
        sys.stderr.write("==========================================================================\n  Files\n==========================================================================\n")
        sys.stderr.write("Reading from: {}\n".format(args.input))
        # sys.stderr.write("Writing to: {}\n".format(args.output))
        if args.call_mods or args.output.split(".")[-1]=="sam":
            SAM_OUT = True
            if args.qscore:
                file = args.output.split(".")
                # doing [-1:] rather than [-1] gives a list back
                name, ext = [".".join(file[:-1])], file[-1:]
                pass_file = ".".join(name + ["pass"] + ext)
                fail_file = ".".join(name + ["fail"] + ext)
                PASS = open(pass_file, 'w') 
                FAIL = open(fail_file, 'w')
                OUT = (PASS, FAIL)
                sam_header(PASS)
                sam_header(FAIL)
                sys.stderr.write("Writing to: {}\n".format(pass_file))
                sys.stderr.write("Writing to: {}\n".format(fail_file))
            else:
                OUT = open(args.output, 'w')
                sam_header(OUT)
                sys.stderr.write("Writing to: {}\n".format(args.output))
        else:
            # TODO: check output ends in .fastq
            # if args.output.split(".")[-1] not in ["fastq", "fq"]:
            #   some error!
            if args.qscore:
                file = args.output.split(".")
                # doing [-1:] rather than [-1] gives a list back
                name, ext = [".".join(file[:-1])], file[-1:]
                pass_file = ".".join(name + ["pass"] + ext)
                fail_file = ".".join(name + ["fail"] + ext)
                PASS = open(pass_file, 'w') 
                FAIL = open(fail_file, 'w')
                OUT = (PASS, FAIL)
                sys.stderr.write("Writing to: {}\n".format(pass_file))
                sys.stderr.write("Writing to: {}\n".format(fail_file))
            else:
                OUT = open(args.output, 'w')
                sys.stderr.write("Writing to: {}\n".format(args.output))
        
        s5 = pyslow5.Open(args.input, 'r')
        reads = s5.seq_reads_multi(threads=args.slow5_threads, batchsize=args.slow5_batchsize)
        sys.stderr.write("\n")

        # ==========================================================================
        # Process reads and send to basecall server
        # ==========================================================================
        sys.stderr.write("==========================================================================\n  Basecalling\n==========================================================================\n")
        sys.stderr.write("\n")

        total_reads = 0
        read_counter = 0
        skipped = []
        for read in reads:
            res, skip = submit_read(client, read)
            if not res:
                skipped.append(skip)
                continue
            else:
                read_counter += 1
                total_reads += 1
            if read_counter >= 1000:
                get_reads(client, OUT, SAM_OUT, args.call_mods, read_counter, args.qscore)
                read_counter = 0
            sys.stderr.write("\rprocessed reads: %d" % total_reads)
            sys.stderr.flush()

        # collect any last leftover reads
        if read_counter > 0:
            get_reads(client, OUT, SAM_OUT, args.call_mods, read_counter, args.qscore)
            read_counter = 0

        sys.stderr.write("\n\n")
        sys.stderr.write("Basecalling complete!\n\n")

        # ==========================================================================
        # Finish up, close files, disconnect client and terminate server
        # ==========================================================================
        sys.stderr.write("\n")
        sys.stderr.write("==========================================================================\n  Summary\n==========================================================================\n")
        sys.stderr.write("Processed {} reads\n".format(total_reads))
        sys.stderr.write("skipped {} reads\n".format(len(skipped)))
        sys.stderr.write("\n")
        # close file
        if type(OUT) == tuple:
            OUT[0].close()
            OUT[1].close()
        else:
            OUT.close()

    sys.stderr.write("==========================================================================\n  Cleanup\n==========================================================================\n")
    sys.stderr.write("Disconnecting client\n")
    sys.stderr.write("Disconnecting server\n")
    sys.stderr.write("Done\n")

if __name__ == '__main__':
    main()
