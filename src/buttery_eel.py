#!/usr/bin/env python3

import argparse
import sys
import os
import time
from pathlib import Path
import numpy as np
from yaml import load

import pyslow5

from pyguppy_client_lib.pyclient import PyGuppyClient
from pyguppy_client_lib import helper_functions

from ._version import __version__


total_guppy_poll_time = 0
total_fastq_write_time = 0
total_slow5_read_time = 0

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def start_guppy_server(args):
    """
    Start guppy server in separate process
    """
    server_args = ["--log_path", args.log,
                   "--config", args.config,
                   "--port", args.port,
                   "--use_tcp",
                   "-x", args.device,
                   "--max_queued_reads", args.max_queued_reads,
                   "--chunk_size", args.chunk_size]
    # for i in guppy_server_args:
    #     server_args.append(i)
    ret = helper_functions.run_server(server_args, bin_path=args.guppy_bin)
    return ret


def calibration(digitisation, range):
    """
    input:
        digitisation: float
        range: float
    output:
        scale: float
    """
    return range / digitisation


def write_fastq(fq, header, seq, qscore):
    global total_fastq_write_time
    t0 = time.time()
    fq.write("{}\n".format(header))
    fq.write("{}\n".format(seq))
    fq.write("+\n")
    fq.write("{}\n".format(qscore))
    total_fastq_write_time = total_fastq_write_time + (time.time()-t0)


def get_reads(client, fq, read_counter):
    global total_guppy_poll_time
    t0 = time.time()
    done = 0
    while done < read_counter:
        bcalled = client.get_completed_reads()
        if not bcalled:
            time.sleep(0.2)
            continue
        else:
            for call in bcalled:
                done += 1
                if len(call) != 1:
                    # possible split reads?
                    sys.stderr.write("Call is longer than 1: {}\n".format(len(call)))


                # @read_id runid=bf... sampleid=NA12878_SRE read=476 ch=38 start_time=2020-10-26T19:58:23Z model_version_id=2021-05-17_dna_r9.4.1_minion_96_29d8704b
                # model_version_id = get_model_info(args.config, args.guppy_bin)
                header = "@{} model_version_id={}".format(call[0]['metadata']['read_id'], call[0]['metadata']['model_version_id'])
                sequence = call[0]['datasets']['sequence']
                qscore = call[0]['datasets']['qstring']
                write_fastq(fq, header, sequence, qscore)
    done = 0
    total_guppy_poll_time = total_guppy_poll_time + (time.time()-t0)

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
    tt = time.time()
    # ==========================================================================
    # Software ARGS
    # ==========================================================================
    """
    Example:
    ./buttery-eel.py --guppy_bin /home/jamfer/Downloads/ont-guppy-6.1.3/bin \
    --max_queued_reads 1000 --port 5558 -x "auto" \
    -i ~/Data/bench/1_slow5/PAF25452_pass_bfdfd1d8_11.blow5 \
    -o ~/Data/bench/buttery_test/test.fastq
    """

    VERSION = __version__

    parser = MyParser(description="buttery-eel - wrapping guppy for file agnostic basecalling",
    epilog="Citation:...",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    # Args for the wrapper, and then probably best to just have free form args for guppy
    parser.add_argument("-i", "--input", required=True,
                        help="input blow5 file for basecalling")
    parser.add_argument("-o", "--output", required=True,
                        help="output .fastq file to write")
    parser.add_argument("--guppy_bin", type=Path, required=True,
                        help="path to ont_guppy/bin folder")
    parser.add_argument("--config", default="dna_r9.4.1_450bps_fast.cfg", required=True,
                        help="basecalling model config")
    parser.add_argument("--port", default="5558",
                        help="port to use between server/client")
    parser.add_argument("--log", default="buttery_guppy_logs",
                        help="guppy log folder path")
    parser.add_argument("--max_queued_reads", default="2000",
                        help="Number of reads to send to guppy server queue")
    parser.add_argument("--chunk_size", default="2000",
                        help="signal chunk size, lower this for lower VRAM GPUs")
    parser.add_argument("-x", "--device", default="auto",
                        help="Specify GPU device: 'auto', or 'cuda:<device_id>'")
    # parser.add_argument("--guppy_server_args",
    #                     help="config file containing any extra args to set on guppy_basecall_server")
    # parser.add_argument("--guppy_client_args",
    #                     help="config file containing any extra args to set on guppy_basecall_client")
    parser.add_argument("-v", "--version", action='version', version="buttery-eel - wraping guppy for file agnostic basecalling version: {}".format(VERSION),
                        help="Prints version")
    # parser.add_argument("--debug", action="store_true",
    #                     help="Set logging to debug mode")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    sys.stderr.write("\n")
    sys.stderr.write("               ~  buttery-eel - SLOW5 Guppy Basecalling  ~\n")
    sys.stderr.write("==========================================================================\n  ARGS\n==========================================================================\n")
    sys.stderr.write("args:\n {}\n".format(args))

    # guppy_server_args = None
    # guppy_client_args = None

    # if args.guppy_server_args:
    #     with open(args.guppy_server_args) as f:
    #         guppy_server_args = yaml.safe_load(f)
    #
    # if args.guppy_client_args:
    #     with open(args.guppy_client_args) as f:
    #         guppy_client_args = yaml.safe_load(f)

    # ==========================================================================
    # Start guppy_basecall_server
    # ==========================================================================
    total_start_server_time = 0
    sst = time.time()
    sys.stderr.write("\n\n")
    sys.stderr.write("==========================================================================\n  Starting Guppy Basecalling Server\n==========================================================================\n")
    serv = start_guppy_server(args)
    server = serv[0]
    p = serv[1]
    sys.stderr.write("guppy_basecall_server started...\n")
    sys.stderr.write("\n")
    total_start_server_time = total_start_server_time + (time.time() - sst)


    # ==========================================================================
    # Connect to server with guppy_basecall_client
    # ==========================================================================

    # TODO: add guppy_client_args
    total_client_connect_time = 0
    sst = time.time()
    sys.stderr.write("==========================================================================\n  Connecting to server\n==========================================================================\n")
    client = PyGuppyClient(
    "127.0.0.1:{}".format(args.port),
    "{}".format(args.config))

    sys.stderr.write("Setting params...\n")
    client.set_params({"priority": PyGuppyClient.high_priority})

    sys.stderr.write("Connecting...\n")
    client.connect()
    sys.stderr.write("Connected, testing connection...\n")
    sys.stderr.write("Connection status:\n")
    print(client.get_status())
    # print(client.get_barcode_kits("127.0.0.1:{}".format(args.port), 10))
    # print(client.get_protocol_version())
    # print(client.get_server_information("127.0.0.1:{}".format(args.port), 10))
    # print(client.get_software_version())
    total_client_connect_time = total_client_connect_time + (time.time() - sst)

    sys.stderr.write("\n")

    # ==========================================================================
    # Read signal file
    # ==========================================================================
    sys.stderr.write("==========================================================================\n  Files\n==========================================================================\n")
    sys.stderr.write("Reading from: {}\n".format(args.input))
    sys.stderr.write("Writing to: {}\n".format(args.output))
    fq = open(args.output, 'w')
    s5 = pyslow5.Open(args.input, 'r')
    reads = s5.seq_reads()
    sys.stderr.write("\n")

    # ==========================================================================
    # Process reads and send to basecall server
    # ==========================================================================
    sys.stderr.write("==========================================================================\n  Basecalling\n==========================================================================\n")
    sys.stderr.write("\n")

    get_raw = True
    total_reads = 0
    read_counter = 0
    done = 0
    skipped = []
    global total_slow5_read_time
    for read in reads:
        t0 = time.time()
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
                    skipped.append(read_id)
                    break

        if not result:
            continue

        else:
            read_counter += 1
            total_reads += 1
        total_slow5_read_time = total_slow5_read_time + (time.time() - t0)
        if read_counter >= 1000:
            get_reads(client, fq, read_counter)
            read_counter = 0
        sys.stderr.write("\rprocessed reads: %d" % total_reads)
        sys.stderr.flush()

    # collect any last leftover reads
    if read_counter > 0:
        get_reads(client, fq, read_counter)
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
    fq.close()
    s5.close()

    sys.stderr.write("==========================================================================\n  Cleanup\n==========================================================================\n")
    # terminate guppy_basecall_server
    sys.stderr.write("Disconnecting client\n")
    client.disconnect()
    sys.stderr.write("Disconnecting server\n")
    server.terminate()
    sys.stderr.write("Done\n")

    sys.stderr.write("\n")
    sys.stderr.write("DEBUG: Timing info:\n")
    sys.stderr.write("total_start_server_time: {}s\n".format(total_start_server_time))
    sys.stderr.write("total_client_connect_time: {}s\n".format(total_client_connect_time))
    sys.stderr.write("total_guppy_poll_time: {}s\n".format(total_guppy_poll_time))
    sys.stderr.write("total_fastq_write_time: {}s\n".format(total_fastq_write_time))
    sys.stderr.write("total_slow5_read_time: {}s\n".format(total_slow5_read_time))
    sys.stderr.write("Total script time: {}s\n".format(time.time() - tt))



if __name__ == '__main__':
    main()
