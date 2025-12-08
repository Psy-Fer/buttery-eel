import sys
from io import StringIO
import numpy as np
import time
from contextlib import contextmanager, redirect_stdout
import re

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
    basecaller_bin = args.basecaller_bin
    found = False
    gotten = False
    model_path = False
    mod_path = False
    tmp_args = []
    for arg in server_args:
        if not gotten:
            if arg == "--guppy_bin":
                found = True
                continue
            if found:
                basecaller_bin = arg
                gotten = True
                continue
        if arg in ["--do_read_splitting", "--detect_adapter", "--detect_mid_strand_adapter", "--detect_mid_strand_barcodes"]:
            print("==========================================================================\n  Depricated Arguments Detected\n==========================================================================")
            print("{} no longer used, please check docs to learn more".format(arg))
            print("\n")
            continue
        tmp_args.append(arg)
    
    if basecaller_bin is None:
        print("-g/--basecaller_bin/--guppy_bin is a required argument")
        sys.exit(1)

    server_args = tmp_args

    server_args.extend(["--log_path", args.log,
                        # "--port", args.port,
                        # "--max_queued_reads", args.max_queued_reads,
                        # "--chunk_size", args.chunk_size,
                        ])
    # if no model path given, add --config back in
    if args.model:
        # TODO: add version based backward compatible stuff here, because it was hinky for a bit before 7.11.2
        server_args.extend(["--model", args.model])
        model_path = True
        if args.modbase_models:
            server_args.extend(["--modbase_models", args.modbase_models])
            mod_path = True
    elif args.config:
        server_args.extend(["--config", args.config])
    else:
        print("ERROR: No model or config detected. Exiting")
        sys.exit(1)
    
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
    
    # if args.above_798:
    #         params["move_enabled"] = True
    
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
    if not args.above_7412:
        params["detect_adapter"] = True
        params["min_score_adapter"] = args.min_score_adapter
    
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
    if args.above_768:
        if args.estimate_poly_a:
            params["estimate_poly_a"] = True
            if args.poly_a_config is not None:
                params["poly_a_config"] = args.poly_a_config
    
    if args.duplex:
        if args.above_7412:
            params["pair_by_channel"] = True
        else:
            raise RuntimeError("Duplex calling not avilable for versions lower than dorado-server 7.4.12")

    # This function has it's own prints that may want to be suppressed
    with redirect_stdout(StringIO()) as fh:
        server, port = helper_functions.run_server(server_args, bin_path=basecaller_bin)

    if port == "ERROR":
        raise RuntimeError("Server couldn't be started")

    if port.startswith("ipc"):
        address = "{}".format(port)
    else:
        address = "localhost:{}".format(port)
    if model_path:
        # create the model set <simplex_model>|<mod_models>|<duplex_model> (must include the ||)
        # takes a single argument, but can be a comma sep list
        # it takes the simplex model, and appends the mod models to it
        # so if simplex is dna_r10.4.1_e8.2_400bps_hac@v4.3.0
        # and mod_path is 5mC_5hmC@v1,6mA@v2
        # dna_r10.4.1_e8.2_400bps_hac@v4.3.0_5mC_5hmC@v1,dna_r10.4.1_e8.2_400bps_hac@v4.3.0_6mA@v2
        if args.modbase_models:
            mod_models = args.modbase_models
        else:
            mod_models = ""

        # excluding this given duplex is dead
        duplex_model = ""
        model_set = "{}|{}|{}".format(args.model, mod_models, duplex_model)
        client = pclient(address=address, config=model_set)
    else:
        client = pclient(address=address, config=args.config)



    print("Setting params...")
    client.set_params(params)
    print("Connecting...")
    try:
        with client:
            if model_path:
                yield [client, address, model_set, params]
            else:
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
                            end_reason=read["aux_data"]["end_reason_labels"][read['end_reason']],
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
            if tries >= 2000:
                if not result:
                    print("Skipped a read: {}".format(read_id))
                    skipped.append([read_id, "stage-0", "timed out trying to submit read to client"])
                    break
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
            model_id = client.get_basecalling_config()[0]["model_version_id"]
            for calls in bcalled:
                done += 1
                if not isinstance(calls, list):
                    calls = [calls]
                split_reads = False
                if len(calls) > 1:
                    split_reads = True
                for call in calls:
                    # if split_reads and int(call['metadata']['split_point']) > 0:
                    #     print(call)
                    #     sys.exit(1)
                    # for i in call:
                    #     if isinstance(call[i], dict):
                    #         for j in call[i]:
                    #             print("{}: {}".format(j, call[i][j]))
                    #     else:
                    #         print("{}: {}".format(i, call[i]))
                    #     time.sleep(5)
                        # sys.exit(1)
                    if args.duplex:
                        if int(call['metadata']['channel']) > 20000:
                            # print("Fake read:", call['metadata']['channel'])
                            # print("Channel: {} - fake read {}/{}".format(call['metadata']['channel'], done, read_counter))
                            done -= 1
                            # print("ending set")
                            # ending = True
                            continue
                    # print("Channel: {} - read {}/{}".format(call['metadata']['channel'], done, read_counter))
                    # if call['metadata']['read_id'] == "6f57d0a2-cfe4-4678-b6bc-299527eec9cc":
                    try:
                        bcalled_read = {}
                        bcalled_read["split_read"] = False
                        bcalled_read["sam_record"] = ""
                        read_id = call['metadata']['read_id']
                        bcalled_read["parent_read_id"] = read_id
                        if split_reads:
                            bcalled_read["split_read"] = True
                            bcalled_read["read_id"] = call['metadata']['strand_id']
                            bcalled_read["split_point"] = int(call['metadata']['split_point'])
                        else:
                            bcalled_read["read_id"] = read_id
                        bcalled_read["read_qscore"] = call['metadata']['mean_qscore']
                        bcalled_read["float_read_qscore"] = round(float(call['metadata']['mean_qscore']), 3)
                        if args.call_mods:
                            bcalled_read["header"] = "@{} parent_read_id={} model_version_id={} modbase_model_version_id={} mean_qscore={}".format(bcalled_read["read_id"], bcalled_read["parent_read_id"], call['metadata'].get('model_version_id', model_id), call['metadata'].get('modbase_model_version_id', model_id), bcalled_read["float_read_qscore"])
                        else:
                            bcalled_read["header"] = "@{} parent_read_id={} model_version_id={} mean_qscore={}".format(bcalled_read["read_id"], bcalled_read["parent_read_id"], call['metadata'].get('model_version_id', model_id), bcalled_read["float_read_qscore"])
                        bcalled_read["sequence"] = call['datasets']['sequence']
                        if args.U2T:
                            seq = []
                            bcalled_read["sequence"] = re.sub("U", "T", bcalled_read["sequence"])
                        if args.above_768:
                            if args.estimate_poly_a:
                                bcalled_read["poly_tail_length"] = call['metadata'].get('poly_tail_length', 0)
                        if args.duplex:
                            bcalled_read["duplex_parent"] = call['metadata']['is_duplex_parent']
                            bcalled_read["duplex_strand_1"] = call['metadata'].get('duplex_strand_1', None)
                            bcalled_read["duplex_strand_2"] = call['metadata'].get('duplex_strand_2', None)

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
                        # if args.moves_out or args.above_798:
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
                            if len(bcalled_read["sam_record"]) > 0 and args.U2T:
                                splitrec = bcalled_read["sam_record"].split("\t")
                                splitrec[9] = re.sub("U", "T", splitrec[9])
                                bcalled_read["sam_record"] = "\t".join(splitrec)
                        if args.do_read_splitting and not args.above_7310:
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
                        # do it for every read now, but only write it if the flag is on
                        # if args.seq_sum:
                        minknow_events = call['metadata'].get('num_minknow_events', ".")
                        sample_rate = float(read_store[read_id]["sampling_rate"])
                        duration = round(float(call['metadata']['duration'] / sample_rate), 6)
                        bcalled_read["duration"] = duration
                        num_events = call['metadata']['num_events']
                        median = round(call['metadata']['median'], 6)
                        med_abs_dev = round(call['metadata']['med_abs_dev'], 6)
                        bcalled_read["scaling_median"] = round(float(call['metadata']['scaling_median']), 3)
                        bcalled_read["scaling_med_abs_dev"] = round(float(call['metadata']['scaling_med_abs_dev']), 8)
                        bcalled_read["scaling_version"] = call['metadata']['scaling_version']
                        # pore_type = read_store[read_id]["header_array"].get('pore_type', 'not_set')
                        experiment_id = read_store[read_id]["header_array"].get('protocol_group_id', ".")
                        run_id = read_store[read_id]["header_array"]["run_id"]
                        sample_id = read_store[read_id]["header_array"].get("sample_id", ".")
                        strand_score_template = round(call['metadata'].get('call_score', 0.0), 6)
                        sequence_length = call['metadata']['sequence_length']
                        bcalled_read["sequence_length"] = sequence_length
                        channel = read_store[read_id]["aux_data"]['channel_number']
                        bcalled_read["channel"] = channel
                        mux = read_store[read_id]["aux_data"]['start_mux']
                        bcalled_read["mux"] = int(mux)
                        start_time = round(float(read_store[read_id]["aux_data"]['start_time']) / sample_rate, 6)
                        end_reason_val = read_store[read_id]["aux_data"].get('end_reason', 0)
                        end_reason = read_store[read_id]["aux_data"].get('end_reason_labels', ["unknown"])[end_reason_val]
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
                            bc_sum_out = "\t".join([bcalled_read["parent_read_id"]]+[bcalled_read["read_id"]]+[str(call['metadata'].get(i, ".")) for i in bc_keys])
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


# region get reads
def get_reads2(args, client, bcalled, sk, read_store):
    '''
    Get reads from the basecaller and process them
    '''
    SPLIT_PASS = False
    if args.qscore:
        SPLIT_PASS = True
        qs_cutoff = float(args.qscore)
    bcalled_list = []
    skipped_list = []
    read_id_set = set()
    model_id = client.get_basecalling_config()[0]["model_version_id"]


    for calls in bcalled:
        if not isinstance(calls, list):
            calls = [calls]
        split_reads = False
        if len(calls) > 1:
            split_reads = True
        for call in calls:
            # print(call)
            # print("\n-------------------------------------------------------------------\n")
            try:
                bcalled_read = {}
                bcalled_read["split_read"] = False
                bcalled_read["sam_record"] = ""
                read_id = call['metadata']['read_id']
                bcalled_read["parent_read_id"] = read_id
                read_id_set.add(read_id)
                if split_reads:
                    bcalled_read["split_read"] = True
                    bcalled_read["read_id"] = call['metadata']['strand_id']
                    bcalled_read["split_point"] = int(call['metadata']['split_point'])
                else:
                    bcalled_read["read_id"] = read_id
                bcalled_read["read_qscore"] = call['metadata']['mean_qscore']
                bcalled_read["float_read_qscore"] = round(float(call['metadata']['mean_qscore']), 3)
                if args.call_mods:
                    bcalled_read["header"] = "@{} parent_read_id={} model_version_id={} modbase_model_version_id={} mean_qscore={}".format(bcalled_read["read_id"], bcalled_read["parent_read_id"], call['metadata'].get('model_version_id', model_id), call['metadata'].get('modbase_model_version_id', model_id), bcalled_read["float_read_qscore"])
                else:
                    bcalled_read["header"] = "@{} parent_read_id={} model_version_id={} mean_qscore={}".format(bcalled_read["read_id"], bcalled_read["parent_read_id"], call['metadata'].get('model_version_id', model_id), bcalled_read["float_read_qscore"])
                bcalled_read["sequence"] = call['datasets']['sequence']
                if args.U2T:
                    seq = []
                    bcalled_read["sequence"] = re.sub("U", "T", bcalled_read["sequence"])
                if args.above_768:
                    if args.estimate_poly_a:
                        bcalled_read["poly_tail_length"] = call['metadata'].get('poly_tail_length', 0)
                        bcalled_read["poly_tail_info"] = "{}|{}|{}|{}|{}".format(call['metadata'].get('poly_tail_anchor', -1),
                                                                                 call['metadata'].get('poly_tail_start_1', -1),
                                                                                 call['metadata'].get('poly_tail_end_1', -1),
                                                                                 call['metadata'].get('poly_tail_start_2', -1),
                                                                                 call['metadata'].get('poly_tail_end_2', -1),
                                                                                 )

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
                # if args.moves_out or args.above_798:
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
                    if len(bcalled_read["sam_record"]) > 0 and args.U2T:
                        splitrec = bcalled_read["sam_record"].split("\t")
                        splitrec[9] = re.sub("U", "T", splitrec[9])
                        bcalled_read["sam_record"] = "\t".join(splitrec)
                if args.do_read_splitting and not args.above_7310:
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
                # do it for every read now, but only write it if the flag is on
                # if args.seq_sum:
                minknow_events = call['metadata'].get('num_minknow_events', ".")
                sample_rate = float(read_store[read_id]["sampling_rate"])
                duration = round(float(call['metadata']['duration'] / sample_rate), 6)
                bcalled_read["duration"] = duration
                num_events = call['metadata']['num_events']
                median = round(call['metadata']['median'], 6)
                med_abs_dev = round(call['metadata']['med_abs_dev'], 6)
                bcalled_read["scaling_median"] = round(float(call['metadata']['scaling_median']), 3)
                bcalled_read["scaling_med_abs_dev"] = round(float(call['metadata']['scaling_med_abs_dev']), 8)
                bcalled_read["scaling_version"] = call['metadata']['scaling_version']
                # pore_type = read_store[read_id]["header_array"].get('pore_type', 'not_set')
                experiment_id = read_store[read_id]["header_array"].get('protocol_group_id', ".")
                run_id = read_store[read_id]["header_array"]["run_id"]
                sample_id = read_store[read_id]["header_array"].get("sample_id", ".")
                strand_score_template = round(call['metadata'].get('call_score', 0.0), 6)
                sequence_length = call['metadata']['sequence_length']
                bcalled_read["sequence_length"] = sequence_length
                channel = read_store[read_id]["aux_data"]['channel_number']
                bcalled_read["channel"] = channel
                mux = read_store[read_id]["aux_data"]['start_mux']
                bcalled_read["mux"] = int(mux)
                start_time = round(float(read_store[read_id]["aux_data"]['start_time']) / sample_rate, 6)
                end_reason_val = read_store[read_id]["aux_data"].get('end_reason', 0)
                end_reason = read_store[read_id]["aux_data"].get('end_reason_labels', ["unknown"])[end_reason_val]
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
                    bc_sum_out = "\t".join([bcalled_read["parent_read_id"]]+[bcalled_read["read_id"]]+[str(call['metadata'].get(i, ".")) for i in bc_keys])
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
    
    return bcalled_list, read_id_set

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
        if args.duplex:
            fake_channel_start = 90000  # we are going to increment this so it doesn't try storing the same channels
            read_counter = 0
            read_store = {}
            while True:
                batch = iq.get()
                if batch is None:
                    break
                if batch == "end":
                    # send 10 fake reads with different channels to trick server to flushing cache
                    # server flushes after 10
                    read_id = list(read_store.keys())[0]
                    read = read_store[read_id]
                    print("[BASECALLER] - Sending fake reads to trick basecaller for channel: {}".format(ch))
                    for chhh in range(fake_channel_start, fake_channel_start+10):
                        result = client.pass_read(helper_functions.package_read(
                                raw_data=np.frombuffer(read['signal'], np.int16),
                                read_id=str(chhh),
                                start_time=read['start_time'],
                                daq_offset=read['offset'],
                                daq_scaling=calibration(read['digitisation'], read['range']),
                                sampling_rate=read['sampling_rate'],
                                mux=read['start_mux'],
                                channel=int(chhh),
                                run_id=read_store[read_id]["header_array"]["run_id"],
                                duration=read['len_raw_signal'],
                            ))
                        # if result:
                        #     read_counter += 0
                        # else:
                        #     # TODO: put in some throttle/while loop for this
                        #     print("failed to stuff in fake reads")
                        #     sys.exit(1)
                    # increase counter by 1 to get the 1 fake read but not the other 10
                    bcalled_list = get_reads(args, client, read_counter, sk, read_store)
                    print("[BASECALLER] - writing channel: {}".format(ch))
                    rq.put(bcalled_list)
                    read_counter = 0
                    read_store = {}
                    fake_channel_start += 10
                else:
                    print("[BASECALLER] - submitting channel: {}".format(batch[0]["channel_number"]))
                    # Submit to be basecalled
                    rc, rs = submit_reads(args, client, sk, batch)
                    read_counter += rc
                    read_store.update(rs)
                    # now collect the basecalled reads
                    print("[BASECALLER] - getting basecalled channel: {}".format(batch[0]["channel_number"]))
                    ch = batch[0]["channel_number"]
                    iq.task_done()
        else:
            # get and submit first batch
            batch = iq.get()
            if batch is None:
                return
            
            bcalled_count = 0
            batch_left = 0
            none_batch = False # this detects when a None comes in from the queue to trigger shut down
            last_submited = False # this checks for the last batch being submitted for basecalling
            # Submit to be basecalled
            read_counter, read_store = submit_reads(args, client, sk, batch)
            while True:
                bcalled = client.get_completed_reads()
                if not bcalled:
                    time.sleep(client.throttle)
                    continue
                else:
                    bcalled_count = len(bcalled)
                    # process basecalled reads
                    bcalled_list, read_id_set = get_reads2(args, client, bcalled, sk, read_store)
                    # push to write queue
                    rq.put(bcalled_list)
                    if bcalled_count != len(read_id_set):
                        print("bcalled_count != len(read_id_set): {} vs {}".format(bcalled_count, len(read_id_set)))
                    read_counter -= bcalled_count
                    # remove read_store values already basecalled
                    for key in read_id_set:
                        del read_store[key]
                    # if number of reads basecalled > reads left in batch, get another batch
                    sub_batch = []
                    if batch_left < bcalled_count and not none_batch:
                        # store left over reads
                        if batch_left > 0:
                            sub_batch = [i for i in batch]
                        # mark old batch as done
                        iq.task_done()
                        # get new batch
                        batch = iq.get()
                        if batch is None:
                            none_batch = True
                    if not none_batch:
                        # pull same number of reads that were just basecalled
                        for _ in range(bcalled_count-len(sub_batch)):
                            if batch:  # Check if the list is not empty before popping
                                sub_batch.append(batch.pop())
                            # else:
                                # print("Batch is empty?!?!")
                                # print("batch:", batch)
                        if batch:
                            batch_left = len(batch)
                        else:
                            batch_left = 0
                    # handle last batch
                    elif none_batch and not last_submited:
                        for _ in range(batch_left):
                            if batch:  # Check if the list is not empty before popping
                                sub_batch.append(batch.pop())
                            # else:
                                # print("Batch is empty?!?!")
                        last_submited = True
                        if batch:
                            batch_left = len(batch)
                        else:
                            batch_left = 0
                    else:
                        # we are on the last batch, waiting for it to be completed
                        if read_counter == 0:
                            break
                        else:
                            continue
                    # get sub batch from batch and submit reads, update read_store and adjust counter
                    sub_read_counter, sub_read_store = submit_reads(args, client, sk, sub_batch)
                    read_store.update(sub_read_store)
                    read_counter += sub_read_counter
                    

        
    if args.profile:
        pr.disable()
        s = io.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        with open("submit_worker_{}.txt".format(N), 'w') as f:
            print(s.getvalue(), file=f)
