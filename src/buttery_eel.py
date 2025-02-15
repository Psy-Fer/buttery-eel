#!/usr/bin/env python3

# region IMPORTS
import argparse
import sys
import os
import multiprocessing as mp
import platform
import time

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
from .reader import read_worker, duplex_read_worker, duplex_read_worker_single
from .writer import write_worker
from .basecaller import start_guppy_server_and_client, basecaller_proc

# region constants
# total_reads = 0
# div = 50
# skipped = 0

# How we get data out of the model files if they are not provided by the metadata output

# def get_model_info(config, basecaller_bin):
#     config = os.path.join(basecaller_bin,"../data/", config)
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
#         model_json = os.path.join(basecaller_bin,"../data/", model)
#         with open(model_json, 'r') as f:
#             jdata = json.load(f)
#             model_version_id = jdata["version"]["id"]
#
#     return model_version_id

# region main
def main():
    # ==========================================================================
    # Software ARGS
    # ==========================================================================
    """
    Example:

    buttery-eel --basecaller_bin /install/ont-guppy-6.1.3/bin --use_tcp --chunk_size 200 \
    --max_queued_reads 1000 -x "cuda:all" --config dna_r9.4.1_450bps_fast.cfg --port 5558 \
    -i /Data/test.blow5 -o /Data/test.fastq

    """
    # region version checks

    VERSION = __version__

    # get version to set secret flags
    above_7310_flag = False
    above_7412_flag = False
    above_768_flag = False
    try:
        major, minor, patch = [int(i) for i in pybasecall_client_lib.__version__.split(".")]
    except:
        major, minor, patch = [int(i) for i in pyguppy_client_lib.__version__.split(".")]
    if major >= 7:
        if minor >= 3:
            above_7310_flag = True
        if minor >= 4:
            above_7412_flag = True
        if minor >= 6:
            above_768_flag = True

    # get args from cli
    args, other_server_args, arg_error = get_args(above_7310_flag, above_7412_flag, above_768_flag)

    if len(sys.argv) == 1:
        arg_error(sys.stderr)
        sys.exit(1)
    
    if args.slow5_batchsize > args.max_read_queue_size:
        print("slow5_batchsize > max_read_queue_size, please alter args so max_read_queue_size is the larger value")
        arg_error(sys.stderr)
        sys.exit(1)
    
    # region start of pipeline
    print()
    print("               ~  buttery-eel - SLOW5 Guppy/Dorado Server Basecalling  ~")
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

    if args.resume is not None:
        if args.resume.split(".")[-1] not in ["fastq", "sam"]:
            print("ERROR: resume file {} is not a fastq or sam file".format(args.resume))
            arg_error(sys.stderr)
            sys.exit(1)
        elif not os.path.isfile(args.resume):
            print("ERROR: resume file {} is does not exist".format(args.resume))
            arg_error(sys.stderr)
            sys.exit(1)
        else:
            args.resume_run = True

    # ==========================================================================
    # region Start guppy_basecall_server
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
        # region connect client
        # TODO: add guppy_client_args
        print("==========================================================================\n  Connecting to server\n==========================================================================")
        print("Connection status:")
        print("status: {}".format(client.get_status()))
        print("throttle: {}".format(client.throttle))
        # print("Client Basecalling config:")
        # print(client.get_basecalling_config())
        bc_config = client.get_basecalling_config()[0]
        print(bc_config)
        # print("model: {}".format(bc_config["model_version_id"]))
        model_version_id = bc_config["model_version_id"]
        model_config_name = bc_config["config_name"]
        # print("Server Basecalling config:")
        # print("get_server_internal_state():", client.get_server_internal_state(address, 10))
        # print(client.get_server_information("127.0.0.1:5000", 10))
        # print(client.get_barcode_kits("127.0.0.1:{}".format(args.port), 10))
        # print(client.get_protocol_version())
        # print(client.get_software_version())

        print("\n")

        # ==========================================================================
        # Read signal file
        # ==========================================================================
        # region file handler
        print("==========================================================================\n  Files\n==========================================================================")
        print("Reading from: {}".format(args.input))
        
        print("Output: {}".format(args.output))
        if args.output.split(".")[-1] not in ["fastq", "sam"]:
            print("ERROR: output file is not a fastq or sam file")
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
        # region run
        print("==========================================================================\n  Basecalling\n==========================================================================")
        print()

        # check if model is for RNA. If so, print a warning about U/T and the flag --U2T
        if "rna" in model_version_id or "RNA" in model_version_id:
            # if args.output.split(".")[-1]=="sam":
            print("==========================================================================\n  RNA U/T warning \n==========================================================================")
            print("RNA model detected: {}/{}".format(bc_config["config_name"], model_version_id))
            if not args.U2T:
                print("By default, Uracil (U) will be written. To instead write Thymine (T), use the --U2T flag")
            else:
                print("--U2T flag has been enabled. Uracil (U) bases will be converted to Thymine (T) bases. If this was not intended, please remove the --U2T flag")
            print("\n")

        mp.set_start_method('spawn')

        if platform.system() == "Darwin":
            im = mp.Manager()
            rm = mp.Manager()
            sm = mp.Manager()
            input_queue = im.JoinableQueue()
            result_queue = rm.JoinableQueue()
            skip_queue = sm.JoinableQueue()
        else:
            input_queue = mp.JoinableQueue()
            result_queue = mp.JoinableQueue()
            skip_queue = mp.JoinableQueue()
        
        # track total samples for samples/s calculation
        total_samples = mp.Value('i', 0)
        sample_time_start = time.perf_counter()

        processes = []

        if args.duplex:
            if platform.system() == "Darwin":
                 print("MacOS not currently supported for duplex calling")
                 sys.exit(1)
            if args.single:
                print("Duplex mode active - a duplex model must be used to output duplex reads")
                print("Buttery-eel does not have checks for this, as the model names are in flux")
                print("SINGLE MODE ACTIVATED - FOR TESTING")
                print()
                duplex_pre_queue = mp.JoinableQueue()
                # create the same number of queues as there are worker processes so each has its own queue
                # queue_names = range(args.procs)
                # duplex_queues = {name: mp.JoinableQueue() for name in queue_names}
                duplex_queue = mp.JoinableQueue()
                reader = mp.Process(target=duplex_read_worker_single, args=(args, duplex_queue, duplex_pre_queue), name='duplex_read_worker_single')
                reader.start()
                out_writer = mp.Process(target=write_worker, args=(args, result_queue, OUT, SAM_OUT, model_version_id, model_config_name), name='write_worker')
                out_writer.start()
                # set up each worker to have a unique queue, so it only processes 1 channel at a time
                basecall_worker = mp.Process(target=basecaller_proc, args=(args, duplex_queue, result_queue, skip_queue, address, config, params, 0), daemon=True, name='basecall_worker_{}'.format(0))
                basecall_worker.start()
                processes.append(basecall_worker)

            else:
                print("Duplex mode active - a duplex model must be used to output duplex reads")
                print("Buttery-eel does not have checks for this, as the model names are in flux")
                print()
                duplex_pre_queue = mp.JoinableQueue()
                # create the same number of queues as there are worker processes so each has its own queue
                queue_names = range(args.procs)
                duplex_queues = {name: mp.JoinableQueue() for name in queue_names}
                reader = mp.Process(target=duplex_read_worker, args=(args, duplex_queues, duplex_pre_queue), name='duplex_read_worker')
                reader.start()
                out_writer = mp.Process(target=write_worker, args=(args, result_queue, OUT, SAM_OUT, model_version_id, model_config_name), name='write_worker')
                out_writer.start()
                # set up each worker to have a unique queue, so it only processes 1 channel at a time
                for name in queue_names:
                    basecall_worker = mp.Process(target=basecaller_proc, args=(args, duplex_queues[name], result_queue, skip_queue, address, config, params, name), daemon=True, name='basecall_worker_{}'.format(name))
                    basecall_worker.start()
                    processes.append(basecall_worker)
        else:
            reader = mp.Process(target=read_worker, args=(args, input_queue, total_samples), name='read_worker')
            reader.start()
            out_writer = mp.Process(target=write_worker, args=(args, result_queue, OUT, SAM_OUT, model_version_id, model_config_name), name='write_worker')
            out_writer.start()
            for i in range(args.procs):
                basecall_worker = mp.Process(target=basecaller_proc, args=(args, input_queue, result_queue, skip_queue, address, config, params, i), daemon=True, name='basecall_worker_{}'.format(i))
                basecall_worker.start()
                processes.append(basecall_worker)

        sample_time_start = time.perf_counter()

        # Anakin, the Process supervisor
        # Monitors all procs for a non-zero exit code. If found, it termintates all the children.
        # If all the exitcodes are 0, it breaks the while loop and continues to join() calls.
        while True:
            # print("reader exit code:", reader.exitcode)
            if reader.exitcode is not None:
                if reader.exitcode != 0:
                    print("ERROR: Reader process encountered an error. exitcode: ", reader.exitcode)
                    for child in mp.active_children():
                        child.terminate()
                    sys.exit(1)
            # print("writer exit code:", out_writer.exitcode)
            if out_writer.exitcode is not None:
                if out_writer.exitcode != 0:
                    print("ERROR: Writer process encountered an error. exitcode: ", out_writer.exitcode)
                    for child in mp.active_children():
                        child.terminate()
                    sys.exit(1)
            for p in processes:
                # print("proc exit code:", p.exitcode)
                if p.exitcode is not None:
                    if p.exitcode != 0:
                        print("ERROR: Worker client encountered an error. exitcode: ", p.exitcode)
                        for child in mp.active_children():
                            child.terminate()
                        sys.exit(1)
            if reader.exitcode == 0:
                p_sum = 0
                for p in processes:
                    if p.exitcode != 0:
                        p_sum += 1
                if p_sum == 0:
                    result_queue.put(None)
                    time.sleep(3)
                    if out_writer.exitcode == 0:
                        print("\n\nProc supervisor: all processes completed without detected error")
                        break
            time.sleep(5)

        sample_time_end = time.perf_counter()
        final_total_samples = 0
        with total_samples.get_lock():
            final_total_samples = total_samples.value

        total_time = sample_time_end - sample_time_start
        samples_per_sec = 0
        if final_total_samples > 0:
            samples_per_sec = float(final_total_samples) / float(total_time)
        
        #Basecalled @ Samples/s: 3.450401e+07
        print("\nBasecalled @ Samples/s:", "{0:.6e}".format(samples_per_sec))

        # Join() calls for all procs
        reader.join()
        if reader.exitcode != 0:
            print("ERROR: Reader process encountered an error. exitcode: ", reader.exitcode)
            for child in mp.active_children():
                child.terminate()
            sys.exit(1)
        for p in processes:
            p.join()
            if p.exitcode != 0:
                print("ERROR: Worker client encountered an error. exitcode: ", p.exitcode)
                for child in mp.active_children():
                    child.terminate()
                sys.exit(1)
        # result_queue.put(None)
        out_writer.join()
        if out_writer.exitcode != 0:
            print("ERROR: Writer process encountered an error. exitcode: ", out_writer.exitcode)
            for child in mp.active_children():
                child.terminate()
            sys.exit(1)
        

        if skip_queue.qsize() > 0:
            # print("1")
            skipped = 0
            skip_queue.put(None)
            if "/" in args.output:
                SKIPPED = open("{}/skipped_reads.txt".format("/".join(args.output.split("/")[:-1])), "w")
                print("Skipped reads detected, writing details to file: {}/skipped_reads.txt".format("/".join(args.output.split("/")[:-1])))
            else:
                SKIPPED = open("./skipped_reads.txt", "w")
                print("Skipped reads detected, writing details to file: ./skipped_reads.txt")

            SKIPPED.write("read_id\tstage\terror\n")
            # print("2")

            while True:
                read = skip_queue.get()
                if read is None:
                    break
                read_id, stage, error = read
                skipped += 1
                SKIPPED.write("{}\t{}\t{}\n".format(read_id, stage, error))

            # print("3")
            SKIPPED.close()
            print("Skipped reads total: {}".format(skipped))
        
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
