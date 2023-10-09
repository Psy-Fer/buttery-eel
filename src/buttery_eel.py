#!/usr/bin/env python3

import argparse
import sys
import os
import time
from pathlib import Path
import numpy as np
from io import StringIO
from contextlib import contextmanager, redirect_stdout
from itertools import chain
import multiprocessing as mp

import pyslow5

import pyguppy_client_lib
from pyguppy_client_lib.pyclient import PyGuppyClient
from pyguppy_client_lib import helper_functions

import cProfile, pstats, io

from ._version import __version__

# constants
total_reads = 0
div = 50
skipped = 0


total_reads = 0
div = 50
skipped = 0


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

    if args.moves_out:
        params["move_and_trace_enabled"] = True
    
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
    
    if args.barcode_kits:
        params["barcode_kits"] = args.barcode_kits
        params["min_score_barcode_front"] = args.min_score_barcode_front
        params["min_score_barcode_rear"] = args.min_score_barcode_rear
        params["min_score_barcode_mid"] = args.min_score_barcode_mid
        # docs are a bit wonky on this, enable_trim_barcodes vs barcode_trimming_enabled
        params["enable_trim_barcodes"] = args.enable_trim_barcodes
        params["require_barcodes_both_ends"] = args.require_barcodes_both_ends
        params["detect_mid_strand_barcodes"] = args.detect_mid_strand_barcodes

    # This function has it's own prints that may want to be suppressed
    with redirect_stdout(StringIO()) as fh:
        server, port = helper_functions.run_server(server_args, bin_path=args.guppy_bin)

    if port == "ERROR":
        raise RuntimeError("Server couldn't be started")

    if port.startswith("ipc"):
        address = "{}".format(port)
    else:
        address = "localhost:{}".format(port)
    client = PyGuppyClient(address=address, config=args.config)


    print("Setting params...")
    client.set_params(params)
    print("Connecting...")
    try:
        with client:
            yield [client, address, args.config, params]
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

def write_summary(summary, data):
    """
    write summary file output
    """
    summary.write("{}\n".format(data))

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
        'PN:ont basecaller',
        'VN:%s' % pyguppy_client_lib.__version__,
    ])
    PG2 = sep.join([
        '@PG',
        'ID:wrapper',
        'PN:buttery-eel',
        'VN:%s' % __version__,
        'CL:buttery-eel %s' % ' '.join(sys.argv[1:]),
        'DS:ont basecaller wrapper',
    ])
    OUT.write("{}\n".format(HD))
    OUT.write("{}\n".format(PG1))
    OUT.write("{}\n".format(PG2))

def read_worker(args, iq):
    '''
    single threaded worker to read slow5 (with multithreading)
    '''
    if args.profile:
        pr = cProfile.Profile()
        pr.enable()

    header_array = {}
    # is dir, so reading recursivley
    if os.path.isdir(args.input):
        # this adds a limit to how many reads it will load into memory so we
        # don't blow the ram up
        max_limit = int(args.max_read_queue_size / args.slow5_batchsize)
        for dirpath, _, files in os.walk(args.input):
            for sfile in files:
                if sfile.endswith(('.blow5', '.slow5')):
                    s5 = pyslow5.Open(os.path.join(dirpath, sfile), 'r')
                    reads = s5.seq_reads_multi(threads=args.slow5_threads, batchsize=args.slow5_batchsize, aux='all')
                    if args.seq_sum:
                        num_read_groups = s5.get_num_read_groups()
                        for read_group in range(num_read_groups):
                            header_array[read_group] = s5.get_all_headers(read_group=read_group)
                    batches = get_slow5_batch(args, s5, reads, size=args.slow5_batchsize, slow5_filename=sfile, header_array=header_array)
                    # put batches of reads onto the queue
                    for batch in chain(batches):
                        # print(iq.qsize())
                        if iq.qsize() < max_limit:
                            iq.put(batch)
                        else:
                            while iq.qsize() >= max_limit:
                                time.sleep(0.01)
                            iq.put(batch)
        
    else:
        s5 = pyslow5.Open(args.input, 'r')
        filename_slow5 = args.input.split("/")[-1]
        reads = s5.seq_reads_multi(threads=args.slow5_threads, batchsize=args.slow5_batchsize, aux='all')
        if args.seq_sum:
            num_read_groups = s5.get_num_read_groups()
            for read_group in range(num_read_groups):
                header_array[read_group] = s5.get_all_headers(read_group=read_group)
        batches = get_slow5_batch(args, s5, reads, size=args.slow5_batchsize, slow5_filename=filename_slow5, header_array=header_array)
        # this adds a limit to how many reads it will load into memory so we
        # don't blow the ram up
        max_limit = int(args.max_read_queue_size / args.slow5_batchsize)
        # put batches of reads onto the queue
        for batch in chain(batches):
            # print(iq.qsize())
            if iq.qsize() < max_limit:
                iq.put(batch)
            else:
                while iq.qsize() >= max_limit:
                    time.sleep(0.01)
                iq.put(batch)
    for _ in range(args.procs):
        iq.put(None)
    
    # if profiling, dump info into log files in current dir
    if args.profile:
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        with open("read_worker.log", 'w') as f:
            print(s.getvalue(), file=f)


def write_worker(args, q, files, SAM_OUT):
    '''
    single threaded worker to process results queue
    '''
    if args.profile:
        pr = cProfile.Profile()
        pr.enable()
    
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
    
    if args.barcode_kits:
        bc_files = {}
        if "/" in args.output:
            BARCODE_SUMMARY = open("{}/barcoding_summary.txt".format("/".join(args.output.split("/")[:-1])), "w")
            print("Writing summary file to: {}/barcoding_summary.txt".format("/".join(args.output.split("/")[:-1])))
        else:
            BARCODE_SUMMARY = open("./barcoding_summary.txt", "w")
            print("Writing summary file to: ./barcoding_summary.txt")
        BARCODE_SUMMARY_HEADER = "\t".join(["parent_read_id", "read_id", "barcode_arrangement", "barcode_full_arrangement", "barcode_kit", "barcode_variant", "barcode_score",
                                            "barcode_front_id", "barcode_front_score", "barcode_front_refseq", "barcode_front_foundseq", "barcode_front_foundseq_length",
                                            "barcode_front_begin_index", "barcode_rear_id", "barcode_rear_score", "barcode_rear_refseq", "barcode_rear_foundseq", "barcode_rear_foundseq_length",
                                            "barcode_rear_end_index"])
        write_summary(BARCODE_SUMMARY, BARCODE_SUMMARY_HEADER)


    if SAM_OUT:
        if args.qscore:
            PASS = open(files["pass"], 'w') 
            FAIL = open(files["fail"], 'w')
            sam_header(PASS)
            sam_header(FAIL)
            OUT = {"pass": PASS, "fail": FAIL}
        else:
            single = open(files["single"], 'w')
            sam_header(single)
            OUT = {"single": single}
    else:
        if args.qscore:
            PASS = open(files["pass"], 'w') 
            FAIL = open(files["fail"], 'w')
            OUT = {"pass": PASS, "fail": FAIL}
        else:
            single = open(files["single"], 'w')
            OUT = {"single": single}
    while True:
        bcalled_list = []
        bcalled_list = q.get()
        if bcalled_list is None:
            break
        for read in bcalled_list:
            fkey = "single"
            if args.qscore:
                if read["read_qscore"] >= float(args.qscore):
                    fkey = "pass"
                else:
                    fkey = "fail"
            # write sequencing_summary file
            if SUMMARY is not None:
                summary_str = read["sum_out"]
                sum_out = files[fkey] + "\t" + summary_str
                write_summary(SUMMARY, sum_out)
            
            if args.barcode_kits:
                # write barcode summary
                bc_summary_str = read["bc_sum_out"]
                write_summary(BARCODE_SUMMARY, bc_summary_str)
                
                # prep barcode writing
                barcode = read["barcode_arrangement"]
                barcode_name = barcode
                if fkey == "pass":
                    barcode_name = barcode + "_" + "pass"
                elif fkey == "fail":
                    barcode_name = barcode + "_" + "fail"
                
                # create file for new detected barcodes
                if barcode_name not in bc_files:
                    fff = args.output.split(".")
                    # doing [-1:] rather than [-1] gives a list back
                    name, ext = [".".join(fff[:-1])], fff[-1:]
                    # if just a single output
                    bcod_file = ".".join(name + [barcode] + ext)
                    # otherwise split on pass/fail
                    if fkey == "pass":
                        bcod_file = ".".join(name + ["pass"] + [barcode] + ext)
                    elif fkey == "fail":
                        bcod_file = ".".join(name + ["fail"] + [barcode] + ext)

                    bc_files[barcode_name] = open(bcod_file, 'w')
                    if SAM_OUT:
                        bc_writer = bc_files[barcode_name]
                        sam_header(bc_writer)

                # write the barcode split sam/fastq
                bc_writer = bc_files[barcode_name]
                if SAM_OUT:
                    if args.call_mods:
                        bc_writer.write("{}\tBC:Z:{}\n".format(read["sam_record"], barcode))
                    elif args.moves_out:
                        m = read["move_table"].tolist()
                        move_str = ','.join(map(str, m))
                        if args.do_read_splitting:
                            bc_writer.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tmv:B:c,{},{}\tNM:i:0\tqs:i:{}\tpi:Z:{}\tBC:Z:{}\n".format(read["read_id"], read["sequence"], read["qscore"], read["model_stride"], move_str, read["int_read_qscore"], read["parent_read_id"], barcode))
                        else:
                            # do ns and ts tags
                            bc_writer.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tmv:B:c,{},{}\tNM:i:0\tqs:i:{}\tns:i:{}\tts:i:{}\tBC:Z:{}\n".format(read["read_id"], read["sequence"], read["qscore"], read["model_stride"], move_str, read["int_read_qscore"], read["num_samples"], read["trimmed_samples"], barcode))
                    else:
                        if args.do_read_splitting:
                            bc_writer.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tNM:i:0\tqs:i:{}\tpi:Z:{}\tBC:Z:{}\n".format(read["read_id"], read["sequence"], read["qscore"], read["int_read_qscore"], read["parent_read_id"], barcode))
                        else:
                            # do ns and ts tags
                            bc_writer.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tNM:i:0\tqs:i:{}\tns:i:{}\tts:i:{}\tBC:Z:{}\n".format(read["read_id"], read["sequence"], read["qscore"], read["int_read_qscore"], read["num_samples"], read["trimmed_samples"], barcode))

                else:
                    bc_writer.write("{} barcode={}\n".format(read["header"], barcode))
                    bc_writer.write("{}\n".format(read["sequence"]))
                    bc_writer.write("+\n")
                    bc_writer.write("{}\n".format(read["qscore"]))

            write_output(args, read, OUT[fkey], SAM_OUT)
        q.task_done()
    
    if len(OUT.keys()) > 1:
        OUT["pass"].close()
        OUT["fail"].close()
    else:
        OUT["single"].close()
    if args.barcode_kits:
        for fffile in bc_files:
            bc_files[fffile].close()
    
    if args.profile:
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        with open("write_worker.log", 'w') as f:
            print(s.getvalue(), file=f)

def write_output(args, read, OUT, SAM_OUT):
    '''
    write the ouput to the file
    '''
    read_id = read["read_id"]
    if SAM_OUT:
        if args.call_mods:
            OUT.write("{}\n".format(read["sam_record"]))
        elif args.moves_out:
            m = read["move_table"].tolist()
            move_str = ','.join(map(str, m))
            if args.do_read_splitting:
                OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tmv:B:c,{},{}\tNM:i:0\tpi:Z:{}\tqs:i:{}\n".format(read_id, read["sequence"], read["qscore"], read["model_stride"], move_str, read["parent_read_id"], read["int_read_qscore"]))
            else:
                # do ns and ts tags
                OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tmv:B:c,{},{}\tNM:i:0\tqs:i:{}\tns:i:{}\tts:i:{}\n".format(read_id, read["sequence"], read["qscore"], read["model_stride"], move_str, read["int_read_qscore"], read["num_samples"], read["trimmed_samples"]))
        else:
            if args.do_read_splitting:
                OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tNM:i:0\tpi:Z:{}\tqs:i:{}\n".format(read_id, read["sequence"], read["qscore"], read["parent_read_id"], read["int_read_qscore"]))
            else:
                # do ns and ts tags
                OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tNM:i:0\tqs:i:{}\tns:i:{}\tts:i:{}\n".format(read_id, read["sequence"], read["qscore"], read["int_read_qscore"], read["num_samples"], read["trimmed_samples"]))
    else:
        # write fastq
        OUT.write("{}\n".format(read["header"]))
        OUT.write("{}\n".format(read["sequence"]))
        OUT.write("+\n")
        OUT.write("{}\n".format(read["qscore"]))
    global total_reads
    global div
    total_reads += 1
    if not args.quiet:
        if total_reads % div == 0:
            print("processed reads: %d" % total_reads)
            sys.stdout.flush()
        # don't make div larger than 500K
        if total_reads >= div*10 and div <= 50000:
            div = div*10




def submit_read(args, iq, rq, address, config, params, N):
    """
    submit a read to the basecall server
    """
    if args.profile:
        pr = cProfile.Profile()
        pr.enable()
    skipped = []
    read_counter = 0
    SPLIT_PASS = False
    if args.qscore:
        SPLIT_PASS = True
        qs_cutoff = float(args.qscore)
    done = 0
    bcalled_list = []

    client_sub = PyGuppyClient(address=address, config=config)
    client_sub.set_params(params)
    # submit a batch of reads to be basecalled
    with client_sub as client:
        while True:
            read_store = {}
            batch = iq.get()
            if batch is None:
                break
            for read in batch:
                read_id = read['read_id']
                if args.seq_sum:
                    read_store[read_id] = read
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
                            skipped.append(read_id)
                            break
                if result:
                    read_counter += 1
        
            # now collect the basecalled reads
            while done < read_counter:
                bcalled = client.get_completed_reads()
                if not bcalled:
                    time.sleep(client.throttle)
                    continue
                else:
                    for calls in bcalled:
                        done += 1
                        split_reads = False
                        if len(calls) > 1:
                            split_reads = True
                        for call in calls:
                            # print(call)
                            try:
                                bcalled_read = {}
                                bcalled_read["sam_record"] = ""
                                read_id = call['metadata']['read_id']
                                bcalled_read["parent_read_id"] = read_id
                                if split_reads:
                                    bcalled_read["read_id"] = call['metadata']['strand_id']
                                else:
                                    bcalled_read["read_id"] = read_id
                                bcalled_read["read_qscore"] = call['metadata']['mean_qscore']
                                bcalled_read["int_read_qscore"] = int(call['metadata']['mean_qscore'])
                                bcalled_read["header"] = "@{} parent_read_id={} model_version_id={} mean_qscore={}".format(bcalled_read["read_id"], bcalled_read["parent_read_id"], call['metadata']['model_version_id'], bcalled_read["int_read_qscore"])
                                bcalled_read["sequence"] = call['datasets']['sequence']
                                bcalled_read["qscore"] = call['datasets']['qstring']
                                if args.moves_out:
                                    bcalled_read["move_table"] = call['datasets']['movement']
                                    bcalled_read["model_stride"] = call['metadata']['model_stride']
                                if args.call_mods:
                                    try:
                                        bcalled_read["sam_record"] = call['metadata']['alignment_sam_record']
                                    except:
                                        # TODO: add warning that mods model not being used, and exit
                                        bcalled_read["sam_record"] = ""
                                if args.do_read_splitting:
                                    bcalled_read["num_samples"] = None
                                    bcalled_read["trimmed_samples"] = None
                                else:
                                    raw_num_samples = len(call['datasets']['raw_data'])
                                    bcalled_read["trimmed_samples"] = call['metadata']['trimmed_samples']
                                    trimmed_duration = call['metadata']['trimmed_duration']
                                    bcalled_read["num_samples"] = trimmed_duration + bcalled_read["trimmed_samples"]
                                    if bcalled_read["num_samples"] != raw_num_samples:
                                        print("WARNING: {} ns:i:{} != raw_num_samples:{}".format(bcalled_read["read_id"], bcalled_read["num_samples"], raw_num_samples))
                                if SPLIT_PASS:
                                    if bcalled_read["read_qscore"] >= qs_cutoff:
                                        # pass
                                        bcalled_read["out"] = "pass"
                                        passes_filtering = "TRUE"
                                    else:
                                        # fail
                                        bcalled_read["out"] = "fail"
                                        passes_filtering = "FALSE"
                                else:
                                    bcalled_read["out"] = "single"
                                    passes_filtering = "."
                                
                                # do barcoding
                                if args.barcode_kits:
                                    bcalled_read["barcode_arrangement"] = call['metadata']["barcode_arrangement"]
                                    
                                
                                 # create summary data
                                if args.seq_sum:
                                    minknow_events = call['metadata']['num_minknow_events']
                                    duration = call['metadata']['duration']
                                    num_events = call['metadata']['num_events']
                                    median = round(call['metadata']['median'], 6)
                                    med_abs_dev = round(call['metadata']['med_abs_dev'], 6)
                                    pore_type = read_store[read_id]["header_array"]['pore_type']
                                    experiment_id = read_store[read_id]["header_array"]['protocol_group_id']
                                    run_id = read_store[read_id]["header_array"]["run_id"]
                                    sample_id = read_store[read_id]["header_array"]["sample_id"]
                                    strand_score_template = round(call['metadata']['call_score'], 6)
                                    sequence_length = call['metadata']['sequence_length']
                                    channel = read_store[read_id]["aux_data"]['channel_number']
                                    mux = read_store[read_id]["aux_data"]['start_mux']
                                    start_time = read_store[read_id]["aux_data"]['start_time']
                                    end_reason_val = read_store[read_id]["aux_data"]['end_reason']
                                    end_reason = read_store[read_id]["aux_data"]['end_reason_labels'][end_reason_val]
                                    output_name = ""
                                    sum_out = "\t".join([str(i) for i in [read_store[read_id]["slow5_filename"], bcalled_read["parent_read_id"], read_id, run_id, channel, mux, minknow_events,
                                            start_time, duration, passes_filtering, ".", num_events, ".",
                                            sequence_length, round(bcalled_read["read_qscore"], 6), strand_score_template, median, med_abs_dev, pore_type,
                                            experiment_id, sample_id, end_reason]])
                                    bcalled_read["sum_out"] = sum_out

                                # create barcode summary data
                                if args.barcode_kits:
                                    bc_keys = ["barcode_arrangement", "barcode_full_arrangement", "barcode_kit", "barcode_variant", "barcode_score",
                                                "barcode_front_id", "barcode_front_score", "barcode_front_refseq", "barcode_front_foundseq", "barcode_front_foundseq_length",
                                                "barcode_front_begin_index", "barcode_rear_id", "barcode_rear_score", "barcode_rear_refseq", "barcode_rear_foundseq", "barcode_rear_foundseq_length",
                                                "barcode_rear_end_index"]
                                    bc_sum_out = "\t".join([bcalled_read["parent_read_id"]]+[bcalled_read["read_id"]]+[str(call['metadata'][i]) for i in bc_keys])
                                    bcalled_read["bc_sum_out"] = bc_sum_out


                                bcalled_list.append(bcalled_read)
                            except Exception as error:
                                # handle the exception
                                print("An exception occurred:", type(error).__name__, "-", error)
                                sys.exit(1)

                    
            read_counter = 0
            done = 0

            # TODO: make a skipped queue to handle skipped reads
            rq.put(bcalled_list)
            iq.task_done()
            bcalled_list = []
    
    if args.profile:
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        with open("submit_worker_{}.txt".format(N), 'w') as f:
            print(s.getvalue(), file=f)


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

def get_slow5_batch(args, slow5_obj, reads, size=4096, slow5_filename=None, header_array=None):
    """
    re-batchify slow5 output
    """
    batch = []
    for read in reads:
        if args.seq_sum:
            # get header once for each read group
            read_group = read["read_group"]
            # if read_group not in header_array:
            #     header_array[read_group] = slow5_obj.get_all_headers(read_group=read_group)
            # get aux data for ead read
            

            aux_data = {"channel_number": read["channel_number"], 
                                "start_mux": read["start_mux"],
                                "start_time": read["start_time"],
                                "read_number": read["read_number"],
                                "end_reason": read["end_reason"],
                                "median_before": read["median_before"],
                                "end_reason_labels": slow5_obj.get_aux_enum_labels('end_reason')
                                }
            read["aux_data"] = aux_data
            read["header_array"] = header_array[read_group]
            read["slow5_filename"] = slow5_filename
        
        batch.append(read)
        if len(batch) >= size:
            yield batch
            batch = []
    if len(batch) > 0:
        yield batch


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
    epilog="Citation: Hiruna Samarakoon, James M Ferguson, Hasindu Gamaarachchi, Ira W Deveson, Accelerated nanopore basecalling with SLOW5 data format, Bioinformatics, Volume 39, Issue 6, June 2023, btad352, https://doi.org/10.1093/bioinformatics/btad352",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    # Args for the wrapper, and then probably best to just have free form args for guppy
    parser.add_argument("-i", "--input", required=True,
                        help="input blow5 file or directory for basecalling")
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
    parser.add_argument("--procs", type=int, default=4,
                        help="Number of worker processes to use processing reads")
    parser.add_argument("--slow5_batchsize", type=int, default=4000,
                        help="Number of reads to process at a time reading slow5")
    parser.add_argument("--quiet", action="store_true",
                        help="Don't print progress")
    parser.add_argument("--max_read_queue_size", type=int, default=20000,
                        help="Number of reads to process at a time reading slow5")
    parser.add_argument("--log", default="buttery_guppy_logs",
                        help="guppy log folder path")
    parser.add_argument("--moves_out", action="store_true",
                        help="output move table (sam format only)")

    # read splitting
    parser.add_argument("--do_read_splitting", action="store_true",
                        help="Perform read splitting based on mid-strand adapter detection")
    parser.add_argument("--min_score_read_splitting", type=float, default=50.0,
                        help="Minimum mid-strand adapter score for reads to be split")
    
    # Adapter trimming
    parser.add_argument("--detect_adapter", action="store_true",
                        help="Enable detection of adapters at the front and rear of the sequence")
    parser.add_argument("--min_score_adapter", type=float, default=60.0,
                        help="Minimum score for a front or rear adapter to be classified. Default is 60.")
    parser.add_argument("--trim_adapters", action="store_true",
                        help="Flag indicating that adapters should be trimmed. Default is False.")
    parser.add_argument("--detect_mid_strand_adapter", action="store_true",
                        help="Flag indicating that read will be marked as unclassified if the adapter sequence appears within the strand itself. Default is False.")

    # Sequencing Summary file
    parser.add_argument("--seq_sum", action="store_true",
                        help="[Experimental] - Write out sequencing_summary.txt file")
    
    # barcode demultiplexing/trimming
    parser.add_argument("--barcode_kits", action="append",
                        help="Strings naming each barcode kit to use. Default is to not do barcoding.")
    parser.add_argument("--enable_trim_barcodes", action="store_true",
                        help="Flag indicating that barcodes should be trimmed.")
    parser.add_argument("--require_barcodes_both_ends", action="store_true",
                        help="Flag indicating that barcodes must be at both ends.")
    parser.add_argument("--detect_mid_strand_barcodes", action="store_true",
                        help="Flag indicating that read will be marked as unclassified if barcodes appear within the strand itself.")
    parser.add_argument("--min_score_barcode_front", type=float, default=60.0,
                        help="Minimum score for a front barcode to be classified")
    parser.add_argument("--min_score_barcode_rear", type=float, default=60.0,
                        help="Minimum score for a rear barcode to be classified")
    parser.add_argument("--min_score_barcode_mid", type=float, default=60.0,
                        help="Minimum score for mid barcodes to be detected")
    # parser.add_argument("--max_queued_reads", default="2000",
    #                     help="Number of reads to send to guppy server queue")
    # parser.add_argument("--chunk_size", default="2000",
    #                     help="signal chunk size, lower this for lower VRAM GPUs")
    parser.add_argument("--profile", action="store_true",
                        help="run cProfile on all processes - for debugging benchmarking")
    parser.add_argument("-v", "--version", action='version', version="buttery-eel - wrapping guppy for SLOW5 basecalling version: {}".format(VERSION),
                        help="Prints version")
    # parser.add_argument("--debug", action="store_true",
    #                     help="Set logging to debug mode")

    # args = parser.parse_args()
    # This collects known and unknown args to parse the server config options
    args, other_server_args = parser.parse_known_args()


    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    if args.slow5_batchsize > args.max_read_queue_size:
        print("slow5_batchsize > max_read_queue_size, please alter args so max_read_queue_size is the larger value")
        parser.print_help(sys.stderr)
        sys.exit(1)
    
    print()
    print("               ~  buttery-eel - SLOW5 Guppy Basecalling  ~")
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
            print("ERROR: Please use guppy and ont-pyguppy-client-lib version 6.3.0 or higher for modification calling")
            print()
            sys.exit(1)

    # ==========================================================================
    # Start guppy_basecall_server
    # ==========================================================================
    print("\n")
    print("==========================================================================\n  Starting Guppy Basecalling Server\n==========================================================================")
    with start_guppy_server_and_client(args, other_server_args) as client_one:
        client, address, config, params = client_one
        print(client)
        print("guppy_basecall_server started...")
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
        # print(client.get_barcode_kits("127.0.0.1:{}".format(args.port), 10))
        # print(client.get_protocol_version())
        # print(client.get_server_information("127.0.0.1:{}".format(args.port), 10))
        # print(client.get_software_version())

        print("\n")

        # ==========================================================================
        # Read signal file
        # ==========================================================================
        print("==========================================================================\n  Files\n==========================================================================")
        print("Reading from: {}".format(args.input))
        if args.input.split(".")[-1] not in ["slow5", "blow5"]:
            print("input file is not a slow5 or blow5 file")
            parser.print_help(sys.stderr)
            sys.exit(1)

        print("Output: {}".format(args.output))
        if args.output.split(".")[-1] not in ["fastq", "sam"]:
            print("output file is not a fastq or sam file")
            parser.print_help(sys.stderr)
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
                # PASS = open(pass_file, 'w') 
                # FAIL = open(fail_file, 'w')
                OUT = {"pass": pass_file, "fail": fail_file}
                # sam_header(PASS)
                # sam_header(FAIL)
                print("Writing to: {}".format(pass_file))
                print("Writing to: {}".format(fail_file))
            else:
                # single = open(args.output, 'w')
                OUT = {"single": args.output}
                # sam_header(single)
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
                # PASS = open(pass_file, 'w') 
                # FAIL = open(fail_file, 'w')
                OUT = {"pass": pass_file, "fail": fail_file}
                print("Writing to: {}".format(pass_file))
                print("Writing to: {}".format(fail_file))
            else:
                # single = open(args.output, 'w')
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
        skipped = []
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
