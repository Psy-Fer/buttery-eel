import sys
from io import StringIO
import numpy as np
import time
from contextlib import contextmanager, redirect_stdout


try:
    from pybasecall_client_lib.pyclient import PyBasecallClient as pclient
    from pybasecall_client_lib import helper_functions

except ImportError:
    # maybe i can do a version check or something? hard to do from this side of the software
    print("Could not load pybasecall, trying for version earlier versions <=7.2.15 pyguppy lib")
    try:
        from pyguppy_client_lib.pyclient import PyGuppyClient as pclient
        from pyguppy_client_lib import helper_functions
    except ImportError:
        print("Can't import pybasecall_client_lib or pyguppy_client_lib, please check environment and try again.")
        sys.exit(1)

import cProfile, pstats, io


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
    # the high priority queue uses a different batch size which alters the basecalls when called with dorado
    # leaving this on default should set it to medium and give 'correct' results
    # funny enough, it will call R9.4.1 data at higher accuracy, and the opposite impact on R10.4.1
    # params = {"priority": PyGuppyClient.high_priority}
    params = {}

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
    client = pclient(address=address, config=args.config)


    print("Setting params...")
    client.set_params(params)
    print("Connecting...")
    try:
        with client:
            yield [client, address, args.config, params]
    finally:
        server.terminate()


def calibration(digitisation, range):
    """
    input:
        digitisation: float
        range: float
    output:
        scale: float
    """
    return range / digitisation



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

    client_sub = pclient(address=address, config=config)
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
                # if args.seq_sum:
                read_store[read_id] = read
                # calculate scale
                scale = calibration(read['digitisation'], read['range'])
                result = False
                tries = 0
                while not result:
                    if args.above_7310:
                        result = client.pass_read(helper_functions.package_read(
                                    raw_data=np.frombuffer(read['signal'], np.int16),
                                    read_id=read_id,
                                    start_time=read['start_time'],
                                    daq_offset=read['offset'],
                                    daq_scaling=scale,
                                    sampling_rate=read['sampling_rate'],
                                    mux=read['start_mux'],
                                    channel=int(read["channel_number"]),
                                    run_id=read_store[read_id]["header_array"]["run_id"],
                                    duration=read['len_raw_signal'],
                                ))
                    else:
                        result = client.pass_read(helper_functions.package_read(
                                    raw_data=np.frombuffer(read['signal'], np.int16),
                                    read_id=read_id,
                                    daq_offset=read['offset'],
                                    daq_scaling=scale,
                                ))
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
                    model_id = ".".join(args.config.split(".")[:-1])
                    for calls in bcalled:
                        done += 1
                        if not isinstance(calls, list):
                            calls = [calls]
                        split_reads = False
                        if len(calls) > 1:
                            split_reads = True
                        for call in calls:
                            # for i in call:
                            #     if isinstance(call[i], dict):
                            #         for j in call[i]:
                            #             print("{}: {}".format(j, call[i][j]))
                            #     else:
                            #         print("{}: {}".format(i, call[i]))
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
                                bcalled_read["header"] = "@{} parent_read_id={} model_version_id={} mean_qscore={}".format(bcalled_read["read_id"], bcalled_read["parent_read_id"], call['metadata'].get('model_version_id', model_id), bcalled_read["int_read_qscore"])
                                bcalled_read["sequence"] = call['datasets']['sequence']
                            except Exception as error:
                                # handle the exception
                                print("An exception occurred in stage 1:", type(error).__name__, "-", error)
                                sys.exit(1)
                            try:
                                if len(bcalled_read["sequence"]) == 0:
                                    print("read_id: {} has a sequence length of zero, skipping".format(read_id))
                                    continue
                                bcalled_read["qscore"] = call['datasets']['qstring']
                                if args.moves_out:
                                    bcalled_read["move_table"] = call['datasets']['movement']
                                    bcalled_read["model_stride"] = call['metadata']['model_stride']
                                if args.call_mods:
                                    try:
                                        bcalled_read["sam_record"] = call['metadata']['alignment_sam_record']
                                    except Exception as error:
                                        # handle the exception
                                        print("An exception occurred getting sam_record/alignment_sam_record for {}:", read_id, type(error).__name__, "-", error)
                                        bcalled_read["sam_record"] = ""
                                        continue
                                if args.do_read_splitting:
                                    bcalled_read["num_samples"] = None
                                    bcalled_read["trimmed_samples"] = None
                                else:
                                    raw_num_samples = len(call['datasets']['raw_data'])
                                    bcalled_read["trimmed_samples"] = call['metadata']['trimmed_samples']
                                    trimmed_duration = call['metadata']['trimmed_duration']
                                    bcalled_read["num_samples"] = trimmed_duration + bcalled_read["trimmed_samples"]
                                    # turning this warning off, as it was put here to alert us to ns tag changing.
                                    # ONT has changed the definition of this field, and this warning told us about it. 
                                    # if bcalled_read["num_samples"] != raw_num_samples:
                                        # print("WARNING: {} ns:i:{} != raw_num_samples:{}".format(bcalled_read["read_id"], bcalled_read["num_samples"], raw_num_samples))
                            except Exception as error:
                                # handle the exception
                                print("An exception occurred in stage 2:", type(error).__name__, "-", error)
                                sys.exit(1)
                            try:
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
                                    minknow_events = call['metadata'].get('num_minknow_events', ".")
                                    sample_rate = float(read_store[read_id]["sampling_rate"])
                                    duration = round(float(call['metadata']['duration'] / sample_rate), 6)
                                    num_events = call['metadata']['num_events']
                                    median = round(call['metadata']['median'], 6)
                                    med_abs_dev = round(call['metadata']['med_abs_dev'], 6)
                                    # pore_type = read_store[read_id]["header_array"].get('pore_type', 'not_set')
                                    experiment_id = read_store[read_id]["header_array"]['protocol_group_id']
                                    run_id = read_store[read_id]["header_array"]["run_id"]
                                    sample_id = read_store[read_id]["header_array"]["sample_id"]
                                    strand_score_template = round(call['metadata']['call_score'], 6)
                                    sequence_length = call['metadata']['sequence_length']
                                    channel = read_store[read_id]["aux_data"]['channel_number']
                                    mux = read_store[read_id]["aux_data"]['start_mux']
                                    start_time = round(float(read_store[read_id]["aux_data"]['start_time']) / sample_rate, 6)
                                    end_reason_val = read_store[read_id]["aux_data"]['end_reason']
                                    end_reason = read_store[read_id]["aux_data"]['end_reason_labels'][end_reason_val]
                                    output_name = ""
                                    sum_out = "\t".join([str(i) for i in [read_store[read_id]["slow5_filename"], bcalled_read["parent_read_id"], read_id, run_id, channel, mux, minknow_events,
                                            start_time, duration, passes_filtering, ".", num_events, ".",
                                            sequence_length, round(bcalled_read["read_qscore"], 6), strand_score_template, median, med_abs_dev,
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
                                print("An exception occurred in stage 3:", type(error).__name__, "-", error)
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
