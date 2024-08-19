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

# region start basecaller
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
        if args.above_7412:
            params["move_enabled"] = True
        else:
            params["move_and_trace_enabled"] = True
    
    if args.call_mods:
        if args.above_7412:
            params["move_enabled"] = True
        else:
            params["move_and_trace_enabled"] = True
    
    if args.do_read_splitting and not args.above_7310:
        params["do_read_splitting"] = True
        params["min_score_read_splitting"] = args.min_score_read_splitting
    
    if args.detect_adapter and not args.above_7412:
        params["detect_adapter"] = True
        params["min_score_adapter"] = args.min_score_adapter
        
    if args.detect_mid_strand_adapter and not args.above_7412:
        params["detect_mid_strand_adapter"] = True

    if args.trim_adapters:
        params["trim_adapters"] = True
    
    if args.barcode_kits:
        params["barcode_kits"] = args.barcode_kits
        params["enable_trim_barcodes"] = args.enable_trim_barcodes
        params["require_barcodes_both_ends"] = args.require_barcodes_both_ends
        if not args.above_7412:
            params["min_score_barcode_front"] = args.min_score_barcode_front
            params["min_score_barcode_rear"] = args.min_score_barcode_rear
            params["min_score_barcode_mid"] = args.min_score_barcode_mid
            # docs are a bit wonky on this, enable_trim_barcodes vs barcode_trimming_enabled
            params["detect_mid_strand_barcodes"] = args.detect_mid_strand_barcodes
    
    if args.duplex:
        if args.above_7412:
            params["pair_by_channel"] = True
        else:
            raise RuntimeError("Duplex calling not avilable for versions lower than dorado-server 7.4.12")

    # This function has it's own prints that may want to be suppressed
    with redirect_stdout(StringIO()) as fh:
        server, port = helper_functions.run_server(server_args, bin_path=args.basecaller_bin)

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

# region submit reads
def submit_reads(args, client, sk, batch):
    '''
    Submit batch of reads to basecaller
    '''
    skipped = []
    read_counter = 0
    read_store = {}
    for read in batch:
        read_id = read['read_id']
        # TODO: remove the signal from the read_store to reduce memory usage
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
                    skipped.append([read_id, "stage-0", "timed out trying to submit read to client"])
                    break
        if result:
            read_counter += 1
    if args.duplex:
        # send fake reads with different channels to trick server to flushing cache
        for chhh in range(90000, 90011):
            result = client.pass_read(helper_functions.package_read(
                                raw_data=np.frombuffer(read['signal'], np.int16),
                                read_id=str(chhh),
                                start_time=read['start_time'],
                                daq_offset=read['offset'],
                                daq_scaling=scale,
                                sampling_rate=read['sampling_rate'],
                                mux=read['start_mux'],
                                channel=int(chhh),
                                run_id=read_store[read_id]["header_array"]["run_id"],
                                duration=read['len_raw_signal'],
                            ))
            if result:
                read_counter += 1
    if len(skipped) > 0:
        for i in skipped:
            sk.put(i)
    
    return read_counter, read_store


# region get reads
def get_reads(args, client, read_counter, sk, read_store):
    '''
    Get reads from the basecaller and process them
    '''
    SPLIT_PASS = False
    if args.qscore:
        SPLIT_PASS = True
        qs_cutoff = float(args.qscore)
    done = 0
    bcalled_list = []
    skipped_list = []

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
                    # if int(call['metadata']['channel']) > 20000:
                    #     print("Fake read:", call['metadata']['channel'])
                    #     continue
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
                        skipped_list.append([read_id, "stage-1", "Failed to get initial sequence data from read record"])
                        continue
                    try:
                        if len(bcalled_read["sequence"]) == 0:
                            print("read_id: {} has a sequence length of zero, skipping".format(read_id))
                            skipped_list.append([read_id, "stage-1", "Sequence length of zero"])
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
                                skipped_list.append([read_id, "stage-1", "Failed to get sam_record/alignment_sam_record"])
                                continue
                        if args.do_read_splitting or args.above_7310:
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
                        skipped_list.append([read_id, "stage-2", "Error getting data related to sam outpout"])
                        continue
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
                            sum_out = "\t".join([str(i) for i in [read_store[read_id]["slow5_filename"], bcalled_read["parent_read_id"], bcalled_read["read_id"], run_id, channel, mux, minknow_events,
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
                        skipped_list.append([read_id, "stage-3", "Error splitting barcodes or sequencing summary"])
                        continue
    
    if len(skipped_list) > 0:
        for i in skipped_list:
            sk.put(i)
    
    read_counter = 0
    done = 0
    return bcalled_list

# region entry point
def basecaller_proc(args, iq, rq, sk, address, config, params, N):
    """
    submit a read to the basecall server
    """
    if args.profile:
        pr = cProfile.Profile()
        pr.enable()
    
    client_sub = pclient(address=address, config=config)
    client_sub.set_params(params)
    # submit a batch of reads to be basecalled
    with client_sub as client:
        while True:
            batch = iq.get()
            if batch is None:
                break
            # print("[BASECALLER] - submitting channel: {}".format(batch[0]["channel_number"]))
            # Submit to be basecalled
            read_counter, read_store = submit_reads(args, client, sk, batch)
            # now collect the basecalled reads
            # print("[BASECALLER] - getting basecalled channel: {}".format(batch[0]["channel_number"]))
            bcalled_list = get_reads(args, client, read_counter, sk, read_store)
            # TODO: make a skipped queue to handle skipped reads
            # print("[BASECALLER] - writing channel: {}".format(batch[0]["channel_number"]))
            rq.put(bcalled_list)
            iq.task_done()
    
    if args.profile:
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        with open("submit_worker_{}.txt".format(N), 'w') as f:
            print(s.getvalue(), file=f)
