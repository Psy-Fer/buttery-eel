#!/usr/bin/env python3

import argparse
import sys
import os
import multiprocessing as mp

try:
    import pybasecall_client_lib
except ImportError:
    # maybe i can do a version check or something? hard to do from this side of the software
    # print("Could not load pybasecall, trying for version earlier versions <=7.2.15 pyguppy lib")
    try:
        import pyguppy_client_lib
    except ImportError:
        print("Can't import pybasecall_client_lib or pyguppy_client_lib, please check environment and try again.")
        sys.exit(1)

from ._version import __version__
from .cli import get_args
from .reader import read_worker
from .writer import write_worker
from .basecaller import start_guppy_server_and_client, submit_read

# constants
total_reads = 0
div = 50
skipped = 0

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

    # get args from cli
    args, other_server_args, arg_error = get_args()

    # get version to set secret flags
    try:
        major, minor, patch = [int(i) for i in pybasecall_client_lib.__version__.split(".")]
    except:
        major, minor, patch = [int(i) for i in pyguppy_client_lib.__version__.split(".")]
    if major >= 7 and minor >= 3:
        above_7310_flag = True
    else:
        above_7310_flag = False


    # add super sneaky hidden flags the user can't interact with but makes global sharing easier
    extra_args = argparse.Namespace(
        above_7310=above_7310_flag, # is the version >= 7.3.* where the name and inputs change?
    )

    # now merge them. This will all get printed into the arg print below which also helps with troubleshooting
    args = argparse.Namespace(**vars(args), **vars(extra_args))

    if len(sys.argv) == 1:
        arg_error(sys.stderr)
        sys.exit(1)
    
    if args.slow5_batchsize > args.max_read_queue_size:
        print("slow5_batchsize > max_read_queue_size, please alter args so max_read_queue_size is the larger value")
        arg_error(sys.stderr)
        sys.exit(1)
    
    print()
    print("               ~  buttery-eel - SLOW5 Guppy/Dorado Basecalling  ~")
    print("                            version: {}".format(VERSION))
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
        try:
            major, minor, patch = [int(i) for i in pybasecall_client_lib.__version__.split(".")]
        except:
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
            print("ERROR: Please use guppy/dorado and ont-pyguppy-client-lib/pybasecall_client_lib version 6.3.0 or higher for modification calling")
            print()
            sys.exit(1)
    


    # ==========================================================================
    # Start guppy_basecall_server
    # ==========================================================================
    print("\n")
    print("==========================================================================\n  Starting Guppy/Dorado Basecalling Server\n==========================================================================")
    with start_guppy_server_and_client(args, other_server_args) as client_one:
        client, address, config, params = client_one
        print(client)
        print("guppy/dorado started...")
        print("basecaller version:", "{}".format(".".join([str(i) for i in client.get_software_version()])))
        print()


        # ==========================================================================
        # Connect to server with guppy_basecall_client
        # ==========================================================================

        # TODO: add guppy_client_args
        print("==========================================================================\n  Connecting to server\n==========================================================================")
        print("Connection status:")
        print("status: {}".format(client.get_status()))
        print("throttle: {}".format(client.throttle))
        # print("Client Basecalling config:")
        # print(client.get_basecalling_config())
        # print("Server Basecalling config:")
        # print(client.get_server_information("127.0.0.1:5000", 10))
        # print(client.get_barcode_kits("127.0.0.1:{}".format(args.port), 10))
        # print(client.get_protocol_version())
        # print(client.get_software_version())

        print("\n")

        # ==========================================================================
        # Read signal file
        # ==========================================================================
        print("==========================================================================\n  Files\n==========================================================================")
        print("Reading from: {}".format(args.input))

        print("Output: {}".format(args.output))
        if args.output.split(".")[-1] not in ["fastq", "sam"]:
            print("output file is not a fastq or sam file")
            arg_error(sys.stderr)
            sys.exit(1)

        # check that the output dir exists
        if "/" in args.output:
            # get everyting but the name of the file
            output_path = "/".join(args.output.split("/")[:-1])
            if not os.path.exists(output_path):
                # If it doesn't exist, create the directory
                print("{} does not exist, creating it".format(output_path))
                os.makedirs(output_path)
        
        if args.call_mods or args.output.split(".")[-1]=="sam":
            SAM_OUT = True
            if args.qscore:
                file = args.output.split(".")
                # doing [-1:] rather than [-1] gives a list back
                name, ext = [".".join(file[:-1])], file[-1:]
                pass_file = ".".join(name + ["pass"] + ext)
                fail_file = ".".join(name + ["fail"] + ext)
                OUT = {"pass": pass_file, "fail": fail_file}
                print("Writing to: {}".format(pass_file))
                print("Writing to: {}".format(fail_file))
            else:
                OUT = {"single": args.output}
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
                OUT = {"pass": pass_file, "fail": fail_file}
                print("Writing to: {}".format(pass_file))
                print("Writing to: {}".format(fail_file))
            else:
                OUT = {"single": args.output}
                print("Writing to: {}".format(args.output))
        
        print()

        # ==========================================================================
        # Process reads and send to basecall server
        # ==========================================================================
        print("==========================================================================\n  Basecalling\n==========================================================================")
        print()

        mp.set_start_method('spawn')
        input_queue = mp.JoinableQueue()
        result_queue = mp.JoinableQueue()
        processes = []
        reader = mp.Process(target=read_worker, args=(args, input_queue), name='read_worker')
        reader.start()
        out_writer = mp.Process(target=write_worker, args=(args, result_queue, OUT, SAM_OUT), name='write_worker')
        out_writer.start()
        for i in range(args.procs):
            basecall_worker = mp.Process(target=submit_read, args=(args, input_queue, result_queue, address, config, params, i), daemon=True, name='basecall_worker_{}'.format(i))
            basecall_worker.start()
            processes.append(basecall_worker)

        reader.join()
        for p in processes:
            p.join()
        result_queue.put(None)
        out_writer.join()
        print("\n")
        print("Basecalling complete!\n")

        # ==========================================================================
        # Finish up, close files, disconnect client and terminate server
        # ==========================================================================
        print("\n")
        # print("==========================================================================\n  Summary\n==========================================================================")
        # global total_reads
        # print("Processed {} reads\n".format(total_reads))
        # print("skipped {} reads\n".format(len(skipped)))
        # print("\n")

    print("==========================================================================\n  Cleanup\n==========================================================================")
    print("Disconnecting client")
    print("Disconnecting server")
    print("Done")

if __name__ == '__main__':
    main()
