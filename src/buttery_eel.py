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
    
    if args.do_read_splitting:
        params["do_read_splitting"] = True
        params["min_score_read_splitting"] = args.min_score_read_splitting
    
    if args.detect_adapter:
        params["detect_adapter"] = True
        params["min_score_adapter"] = args.min_score_adapter
        
    if args.detect_mid_strand_adapter:
        params["detect_mid_strand_adapter"] = True

    if args.trim_adapters:
        params["trim_adapters"] = True

    
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
    client = PyGuppyClient(address=address, config=args.config, move_and_trace_enabled=args.moves_out)


    print("Setting params...")
    client.set_params(params)
    print("Connecting...")
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


def write_output(args, OUT, read_id, header, seq, qscore, SAM_OUT, read_qscore, sam="", mods=False, moves=False, move_table=None, model_stride=None, num_samples=None, trimmed_samples=None):
    '''
    TODO: use args rather than specific args
    add these fields:
    
    std::vector<std::string> Read::generate_read_tags() const {
    // GCC doesn't support <format> yet...
    std::vector<std::string> tags = {
        "qs:i:" + std::to_string(static_cast<int>(std::round(utils::mean_qscore_from_qstring(qstring)))),
        "ns:i:" + std::to_string(num_samples),
        "ts:i:" + std::to_string(num_trimmed_samples),
        "mx:i:" + std::to_string(attributes.mux),
        "ch:i:" + std::to_string(attributes.channel_number),
        "st:Z:" + attributes.start_time,
        "rn:i:" + std::to_string(attributes.read_number),
        "f5:Z:" + attributes.fast5_filename
    };

    return tags;
    '''
    
    if SAM_OUT:
        if mods:
            OUT.write("{}\n".format(sam))
        elif moves:
            m = move_table.tolist()
            move_str = ','.join(map(str, m))
            if args.do_read_splitting:
                OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tmv:B:c,{},{}\tNM:i:0\tqs:i:{}\n".format(read_id, seq, qscore, model_stride, move_str, read_qscore))
            else:
                # do ns and ts tags
                OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tmv:B:c,{},{}\tNM:i:0\tqs:i:{}\tns:i:{}\tts:i:{}\n".format(read_id, seq, qscore, model_stride, move_str, read_qscore, num_samples, trimmed_samples))
        else:
            if args.do_read_splitting:
                OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tNM:i:0\tqs:i:{}\n".format(read_id, seq, qscore, read_qscore))
            else:
                # do ns and ts tags
                OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tNM:i:0\tqs:i:{}\tns:i:{}\tts:i:{}\n".format(read_id, seq, qscore, read_qscore, num_samples, trimmed_samples))
    else:
        # write fastq
        OUT.write("{}\n".format(header))
        OUT.write("{}\n".format(seq))
        OUT.write("+\n")
        OUT.write("{}\n".format(qscore))

def write_summary(summary, data):
    """
    write summary file output
    """
    summary.write("{}\n".format(data))

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
            time.sleep(client.throttle)
        tries += 1
        if tries >= 1000:
            if not result:
                print("Skipped a read: {}".format(read_id))
                skipped = read_id
                break
    return result, skipped


def get_reads(args, client, OUT, SAM_OUT, SUMMARY, mods, moves, read_counter, qscore_cutoff, header_array, read_groups, aux_data, filename_slow5):
    """
    Get the reads back from the basecall server after being basecalled
    bcalled object contains 1 or more called reads, which contain various data
    """
    SPLIT_PASS = False
    move_table = None
    model_stride = None
    passes_filtering = "-"
    if qscore_cutoff:
        SPLIT_PASS = True
        qs_cutoff = float(qscore_cutoff)
    done = 0
    while done < read_counter:
        bcalled = client.get_completed_reads()
        if not bcalled:
            time.sleep(client.throttle)
            continue
        else:
            for calls in bcalled:
                sam_record = ""
                done += 1
                split_reads = False
                if len(calls) > 1:
                    split_reads = True
                for call in calls:
                    read_id = call['metadata']['read_id']
                    parent_read_id = read_id
                    if split_reads:
                        read_id = call['metadata']['strand_id']
                    read_qscore = call['metadata']['mean_qscore']
                    int_read_qscore = int(read_qscore)
                    # @read_id runid=bf... sampleid=NA12878_SRE read=476 ch=38 start_time=2020-10-26T19:58:23Z model_version_id=2021-05-17_dna_r9.4.1_minion_96_29d8704b
                    # model_version_id = get_model_info(args.config, args.guppy_bin)
                    header = "@{} parent_read_id=@{} model_version_id={} mean_qscore={}".format(read_id, parent_read_id, call['metadata']['model_version_id'], int_read_qscore)
                    sequence = call['datasets']['sequence']
                    qscore = call['datasets']['qstring']
                    # when calling mods, can just output sam_record value
                    # otherwise, write_output will handle unaligned sam with no mods
                    if moves:
                        move_table = call['datasets']['movement']
                        model_stride = call['metadata']['model_stride']
                    if mods:
                        try:
                            sam_record = call['metadata']['alignment_sam_record']
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
                    
                    if args.do_read_splitting:
                        num_samples = None
                        trimmed_samples = None
                    else:
                        raw_num_samples = len(call['datasets']['raw_data'])
                        trimmed_samples = call['metadata']['trimmed_samples']
                        trimmed_duration = call['metadata']['trimmed_duration']
                        num_samples = trimmed_duration + trimmed_samples
                        if num_samples != raw_num_samples:
                            print("WARNING: {} ns:i:{} != raw_num_samples:{}".format(read_id, num_samples, raw_num_samples))
                    
                    # create summary data
                    if SUMMARY is not None:
                        read_group = read_groups[parent_read_id]
                        minknow_events = call['metadata']['num_minknow_events']
                        duration = call['metadata']['duration']
                        num_events = call['metadata']['num_events']
                        median = round(call['metadata']['median'], 6)
                        med_abs_dev = round(call['metadata']['med_abs_dev'], 6)
                        pore_type = header_array[read_group]['pore_type']
                        experiment_id = header_array[read_group]['protocol_group_id']
                        run_id = header_array[read_group]["run_id"]
                        sample_id = header_array[read_group]["sample_id"]
                        strand_score_template = round(call['metadata']['call_score'], 6)
                        sequence_length = call['metadata']['sequence_length']
                        channel = aux_data[parent_read_id]['channel_number']
                        mux = aux_data[parent_read_id]['start_mux']
                        start_time = aux_data[parent_read_id]['start_time']
                        end_reason_val = aux_data[parent_read_id]['end_reason']
                        end_reason = aux_data[parent_read_id]['end_reason_labels'][end_reason_val]
                        output_name = out.name.split("/")[-1]
                        sum_out = "\t".join([str(i) for i in  [output_name, filename_slow5, parent_read_id, read_id, run_id, channel, mux, minknow_events,
                                start_time, duration, passes_filtering, "-", num_events, "-",
                                sequence_length, read_qscore, strand_score_template, median, med_abs_dev, pore_type,
                                experiment_id, sample_id, end_reason]])
                        
                        # write the data
                        write_summary(SUMMARY, sum_out)
                    
                    write_output(args, out, read_id, header, sequence, qscore, SAM_OUT, int_read_qscore, sam=sam_record, mods=mods, moves=moves, move_table=move_table, model_stride=model_stride, num_samples=num_samples, trimmed_samples=trimmed_samples)
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

    parser = MyParser(description="buttery-eel - wrapping guppy for SLOW5 basecalling",
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
    parser.add_argument("--guppy_batchsize", type=int, default=4000,
                        help="number of reads to send to guppy at a time.")
    parser.add_argument("--call_mods", action="store_true",
                        help="output MM/ML tags for methylation - will output sam - use with appropriate mod config")
    parser.add_argument("-q", "--qscore", type=int,
                        help="A mean q-score to split fastq/sam files into pass/fail output")
    parser.add_argument("--slow5_threads", type=int, default=4,
                        help="Number of threads to use reading slow5 file")
    parser.add_argument("--slow5_batchsize", type=int, default=4000,
                        help="Number of reads to process at a time reading slow5")
    parser.add_argument("--quiet", action="store_true",
                        help="Don't print progress")
    parser.add_argument("--moves_out", action="store_true",
                        help="output move table (sam format only)")
    parser.add_argument("--do_read_splitting", action="store_true",
                        help="Perform read splitting based on mid-strand adapter detection")
    parser.add_argument("--min_score_read_splitting", type=float, default=50.0,
                        help="Minimum mid-strand adapter score for reads to be split")
    parser.add_argument("--detect_adapter", action="store_true",
                        help="Enable detection of adapters at the front and rear of the sequence")
    parser.add_argument("--min_score_adapter", type=float, default=60.0,
                        help="Minimum score for a front or rear adapter to be classified. Default is 60.")
    parser.add_argument("--trim_adapters", action="store_true",
                        help="Flag indicating that adapters should be trimmed. Default is False.")
    parser.add_argument("--detect_mid_strand_adapter", action="store_true",
                        help="Flag indicating that read will be marked as unclassified if the adapter sequence appears within the strand itself. Default is False.")
    # Disabling alignment because sam file headers are painful and frankly out of scope. Just use minimap2.
    # parser.add_argument("-a", "--align_ref",
    #                     help="reference .mmi file. will output sam. (build with: minimap2 -x map-ont -d ref.mmi ref.fa )")
    # parser.add_argument("--port", default="5558",
    #                     help="port to use between server/client")
    parser.add_argument("--log", default="buttery_guppy_logs",
                        help="guppy log folder path")
    parser.add_argument("--seq_sum", action="store_true",
                        help="[Experimental] - Write out sequencing_summary.txt file")
    # parser.add_argument("--max_queued_reads", default="2000",
    #                     help="Number of reads to send to guppy server queue")
    # parser.add_argument("--chunk_size", default="2000",
    #                     help="signal chunk size, lower this for lower VRAM GPUs")
    # parser.add_argument("-x", "--device", default="auto",
    #                     help="Specify GPU device: 'auto', or 'cuda:<device_id>'")
    parser.add_argument("-v", "--version", action='version', version="buttery-eel - wraping guppy for SLOW5 basecalling version: {}".format(VERSION),
                        help="Prints version")
    # parser.add_argument("--debug", action="store_true",
    #                     help="Set logging to debug mode")

    # args = parser.parse_args()
    # This collects known and unknown args to parse the server config options
    args, other_server_args = parser.parse_known_args()


    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    print()
    print("               ~  buttery-eel - SLOW5 Guppy Basecalling  ~")
    print("==========================================================================\n  ARGS\n==========================================================================")
    print("args:\n {}\n{}".format(args, other_server_args))

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
        print()
        print("MOD CALLING VERSION CHECK: >6.3.0? {}".format(check))
        print()
        if not check:
            sys.stderr.write("ERROR: Please use guppy and ont-pyguppy-client-lib version 6.3.0 or higher for modification calling\n")
            sys.stderr.write("\n")
            sys.exit(1)

    # ==========================================================================
    # Start guppy_basecall_server
    # ==========================================================================
    print()
    print()
    print("==========================================================================\n  Starting Guppy Basecalling Server\n==========================================================================")
    with start_guppy_server_and_client(args, other_server_args) as client:
        print(client)
        print("guppy_basecall_server started...")
        print()


        # ==========================================================================
        # Connect to server with guppy_basecall_client
        # ==========================================================================

        # TODO: add guppy_client_args
        print("==========================================================================\n  Connecting to server\n==========================================================================")
        print("Connection status:")
        print("status: {}".format(client.get_status()))
        print("throttle: {}".format(client.throttle))
        # print(client.get_barcode_kits("127.0.0.1:{}".format(args.port), 10))
        # print(client.get_protocol_version())
        # print(client.get_server_information("127.0.0.1:{}".format(args.port), 10))
        # print(client.get_software_version())

        print()
        print()

        # ==========================================================================
        # Read signal file
        # ==========================================================================
        print("==========================================================================\n  Files\n==========================================================================")
        print("Reading from: {}".format(args.input))
        # print("Writing to: {}".format(args.output))
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
                print("Writing to: {}".format(pass_file))
                print("Writing to: {}".format(fail_file))
            else:
                OUT = open(args.output, 'w')
                sam_header(OUT)
                print("Writing to: {}".format(args.output))
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
                print("Writing to: {}".format(pass_file))
                print("Writing to: {}".format(fail_file))
            else:
                OUT = open(args.output, 'w')
                print("Writing to: {}".format(args.output))
        
        if args.seq_sum:
            if "/" in args.output:
                SUMMARY = open("{}/sequencing_summary.txt".format("/".join(args.output.split("/")[:-1])), "w")
                print("Writing summary file to: {}/sequencing_summary.txt".format("/".join(args.output.split("/")[:-1])))
            else:
                SUMMARY = open("./sequencing_summary.txt", "w")
                print("Writing summary file to: ./sequencing_summary.txt")

            SUMMARY_HEADER = "\t".join(["filename_fastq", "filename_slow5", "parent_read_id",
                                        "read_id", "run_id", "channel", "mux", "minknow_events", "start_time", "duration",
                                        "passes_filtering", "template_start", "num_events_template", "template_duration",
                                        "sequence_length_template", "mean_qscore_template", "strand_score_template",
                                        "median_template", "mad_template", "pore_type", "experiment_id", "sample_id", "end_reason"])
            write_summary(SUMMARY, SUMMARY_HEADER)
        
        else:
            SUMMARY = None
        
        filename_slow5 = args.input.split("/")[-1]
        s5 = pyslow5.Open(args.input, 'r')
        # reads = s5.seq_reads()
        if args.seq_sum:
            reads = s5.seq_reads_multi(aux='all', threads=args.slow5_threads, batchsize=args.slow5_batchsize)
        else:
            reads = s5.seq_reads_multi(threads=args.slow5_threads, batchsize=args.slow5_batchsize)

        print()

        # ==========================================================================
        # Process reads and send to basecall server
        # ==========================================================================
        print("==========================================================================\n  Basecalling\n==========================================================================")
        print()

        total_reads = 0
        read_counter = 0
        skipped = []
        header_array = {}
        aux_data = {}
        read_groups = {}
        for read in reads:
            if args.seq_sum:
                # get header once for each read group
                read_group = read["read_group"]
                if read_group not in header_array:
                    header_array[read_group] = s5.get_all_headers(read_group=read_group)
                # get aux data for ead read
                readID = read["read_id"]
                read_groups[readID] = read_group
                

                aux_data[readID] = {"channel_number": read["channel_number"], 
                                    "start_mux": read["start_mux"],
                                    "start_time": read["start_time"],
                                    "read_number": read["read_number"],
                                    "end_reason": read["end_reason"],
                                    "median_before": read["median_before"],
                                    "end_reason_labels": s5.get_aux_enum_labels('end_reason')
                                    }
            
            res, skip = submit_read(client, read)
            if not res:
                skipped.append(skip)
                continue
            else:
                read_counter += 1
                total_reads += 1
            if read_counter >= args.guppy_batchsize:
                get_reads(args, client, OUT, SAM_OUT, SUMMARY, args.call_mods, args.moves_out, read_counter, args.qscore, header_array, read_groups, aux_data, filename_slow5)
                read_counter = 0
                aux_data = {}
                read_groups = {}
            if not args.quiet:
                sys.stdout.write("\rprocessed reads: %d" % total_reads)
                sys.stdout.flush()

        # collect any last leftover reads
        if read_counter > 0:
            get_reads(args, client, OUT, SAM_OUT, SUMMARY, args.call_mods, args.moves_out, read_counter, args.qscore, header_array, read_groups, aux_data, filename_slow5)
            read_counter = 0

        print()
        print()
        print("Basecalling complete!")

        # ==========================================================================
        # Finish up, close files, disconnect client and terminate server
        # ==========================================================================
        print()
        print("==========================================================================\n  Summary\n==========================================================================")
        print("Processed {} reads".format(total_reads))
        print("skipped {} reads".format(len(skipped)))
        print()
        # close file
        if type(OUT) == tuple:
            OUT[0].close()
            OUT[1].close()
        else:
            OUT.close()
        if args.seq_sum:
            SUMMARY.close()

    print("==========================================================================\n  Cleanup\n==========================================================================")
    print("Disconnecting client")
    print("Disconnecting server")
    print("Done")

if __name__ == '__main__':
    main()
