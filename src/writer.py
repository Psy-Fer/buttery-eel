import sys, os

from ._version import __version__

import cProfile, pstats, io


try:
    import pybasecall_client_lib
except ImportError:
    try:
        import pyguppy_client_lib
    except ImportError:
        print("Can't import pybasecall_client_lib or pyguppy_client_lib, please check environment and try again.")
        sys.exit(1)

# constants
total_reads = 0
div = 50
skipped = 0

def write_summary(summary, data):
    """
    write summary file output
    """
    summary.write("{}\n".format(data))

def sam_header(OUT, model_version_id, sep='\t'):
    """
    Format a string sam header.
    This is taken from Bonito by Chris Seymour at ONT.
    https://github.com/nanoporetech/bonito/blob/master/bonito/io.py#L103
    
    ReadGroups would require writing the header and reads separate then merging at the end
    Most importantly we want the basecalling/mod models
    So I will put them in the DS tag of the PG2 tag
    RG	ID	<runid>
        PU	<flow_cell_id>
        PM	<device_id>
        DT	<exp_start_time>
        PL	ONT
        DS	basecall_model=<basecall_model_name> modbase_models=<modbase_model_names> runid=<run_id>
        LB	<sample_id>
        SM	<sample_id>

    Sam tags:
    RG:Z:	<runid>_<basecalling_model>_<barcode_arrangement>
    qs:f:	mean basecall qscore
    ts:i:	the number of samples trimmed from the start of the signal
    ns:i:	the basecalled sequence corresponds to the interval signal[ts : ns]
            the move table maps to the same interval.
            note that ns reflects trimming (if any) from the rear
            of the signal.
    mx:i:	read mux
    ch:i:	read channel
    rn:i:	read number
    st:Z:	read start time (in UTC)
    du:f:	duration of the read (in seconds)
    fn:Z:	file name
    sm:f:	scaling midpoint/mean/median (pA to ~0-mean/1-sd)
    sd:f:	scaling dispersion (pA to ~0-mean/1-sd)
    sv:Z:	scaling version
    mv:B:c	sequence to signal move table (optional)
    dx:i:	bool to signify duplex read (only in duplex mode)
    pi:Z:	parent read id for a split read
    sp:i:	start coordinate of split read in parent read signal
    pt:i:	estimated poly(A/T) tail length in cDNA and dRNA reads
    bh:i:	number of detected bedfile hits (only if alignment was performed with a specified bed-file)
    MN:i:	Length of sequence at the time MM and ML were produced
    """
    try:
        basecaller_version = pybasecall_client_lib.__version__
    except:
        basecaller_version = pyguppy_client_lib.__version__

    HD = sep.join([
        '@HD',
        'VN:1.5',
        'SO:unknown',
    ])
    PG1 = sep.join([
        '@PG',
        'ID:basecaller',
        'PN:ont basecaller',
        'VN:%s' % basecaller_version,
    ])
    PG2 = sep.join([
        '@PG',
        'ID:wrapper',
        'PN:buttery-eel',
        'VN:%s' % __version__,
        'CL:buttery-eel %s' % ' '.join(sys.argv[1:]),
        'DS:ont basecaller wrapper basecall_model={}'.format(model_version_id),
    ])
    OUT.write("{}\n".format(HD))
    OUT.write("{}\n".format(PG1))
    OUT.write("{}\n".format(PG2))


def write_worker(args, q, files, SAM_OUT, model_version_id):
    '''
    single threaded worker to process results queue
    '''
    if args.profile:
        pr = cProfile.Profile()
        pr.enable()
    
    if args.seq_sum:
        try:
            if "/" in args.output:
                SUMMARY = open("{}/sequencing_summary.txt".format("/".join(args.output.split("/")[:-1])), "w")
                print("Writing summary file to: {}/sequencing_summary.txt".format("/".join(args.output.split("/")[:-1])))
            else:
                SUMMARY = open("./sequencing_summary.txt", "w")
                print("Writing summary file to: ./sequencing_summary.txt")
        except Exception as error:
            # handle the exception
            print("ERROR: An exception occurred file opening:", type(error).__name__, "-", error)
            sys.exit(1)

        SUMMARY_HEADER = "\t".join(["filename_out", "filename_slow5", "parent_read_id",
                                    "read_id", "run_id", "channel", "mux", "minknow_events", "start_time", "duration",
                                    "passes_filtering", "template_start", "num_events_template", "template_duration",
                                    "sequence_length_template", "mean_qscore_template", "strand_score_template",
                                    "median_template", "mad_template", "experiment_id", "sample_id", "end_reason"])
        write_summary(SUMMARY, SUMMARY_HEADER)
    else:
        SUMMARY = None
    
    if args.barcode_kits:
        bc_files = {}
        try:
            if "/" in args.output:
                BARCODE_SUMMARY = open("{}/barcoding_summary.txt".format("/".join(args.output.split("/")[:-1])), "w")
                print("Writing summary file to: {}/barcoding_summary.txt".format("/".join(args.output.split("/")[:-1])))
            else:
                BARCODE_SUMMARY = open("./barcoding_summary.txt", "w")
                print("Writing summary file to: ./barcoding_summary.txt")
        except Exception as error:
            # handle the exception
            print("ERROR: An exception occurred file opening:", type(error).__name__, "-", error)
            sys.exit(1)
        BARCODE_SUMMARY_HEADER = "\t".join(["parent_read_id", "read_id", "barcode_arrangement", "barcode_full_arrangement", "barcode_kit", "barcode_variant", "barcode_score",
                                            "barcode_front_id", "barcode_front_score", "barcode_front_refseq", "barcode_front_foundseq", "barcode_front_foundseq_length",
                                            "barcode_front_begin_index", "barcode_rear_id", "barcode_rear_score", "barcode_rear_refseq", "barcode_rear_foundseq", "barcode_rear_foundseq_length",
                                            "barcode_rear_end_index"])
        write_summary(BARCODE_SUMMARY, BARCODE_SUMMARY_HEADER)

    try:
        if SAM_OUT:
            if args.qscore:
                PASS = open(files["pass"], 'w') 
                FAIL = open(files["fail"], 'w')
                sam_header(PASS, model_version_id)
                sam_header(FAIL, model_version_id)
                OUT = {"pass": PASS, "fail": FAIL}
            else:
                single = open(files["single"], 'w')
                sam_header(single, model_version_id)
                OUT = {"single": single}
        else:
            if args.qscore:
                PASS = open(files["pass"], 'w') 
                FAIL = open(files["fail"], 'w')
                OUT = {"pass": PASS, "fail": FAIL}
            else:
                single = open(files["single"], 'w')
                OUT = {"single": single}
    except Exception as error:
        # handle the exception
        print("ERROR: An exception occurred file opening:", type(error).__name__, "-", error)
        sys.exit(1)

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
                    try:
                        bc_files[barcode_name] = open(bcod_file, 'w')
                    except Exception as error:
                        # handle the exception
                        print("ERROR: An exception occurred file opening:", type(error).__name__, "-", error)
                        sys.exit(1)
                    if SAM_OUT:
                        bc_writer = bc_files[barcode_name]
                        sam_header(bc_writer, model_version_id)

                # write the barcode split sam/fastq
                bc_writer = bc_files[barcode_name]
                if SAM_OUT:
                    # TODO: Add duplex calling to the barcoded output
                    if args.above_7412:
                        sam_tags = "MN:i:{}\tqs:f:{}\tmx:i:{}\tch:i:{}\tns:i:{}\tts:i:{}\tsm:f:{}\tsd:f:{}\tsv:Z:{}\tdu:f:{}".format(read["sequence_length"],
                                                                                                                        read["float_read_qscore"],
                                                                                                                        read["mux"],
                                                                                                                        read["channel"],
                                                                                                                        read["num_samples"],
                                                                                                                        read["trimmed_samples"],
                                                                                                                        read["scaling_median"],
                                                                                                                        read["scaling_med_abs_dev"],
                                                                                                                        read["scaling_version"],
                                                                                                                        read["duration"])
                        if args.call_mods:
                            bc_writer.write("{}\tpi:Z:{}\t{}\tBC:Z:{}\n".format(read["sam_record"], read["parent_read_id"], sam_tags, barcode))
                        elif args.moves_out:
                            m = read["move_table"].tolist()
                            move_str = ','.join(map(str, m))
                            bc_writer.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tmv:B:c,{},{}\t{}\tpi:Z:{}\tBC:Z:{}\n".format(read["read_id"], read["sequence"], read["qscore"], read["model_stride"], move_str, sam_tags, read["parent_read_id"], barcode))
                        else:
                            bc_writer.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\t{}\tpi:Z:{}\tBC:Z:{}\n".format(read["read_id"], read["sequence"], read["qscore"], sam_tags, read["parent_read_id"], barcode))

                    else:
                        if args.call_mods:
                            if args.do_read_splitting:
                                bc_writer.write("{}\tpi:Z:{}\tBC:Z:{}\n".format(read["sam_record"], read["parent_read_id"], barcode))
                            else:
                                bc_writer.write("{}\tBC:Z:{}\n".format(read["sam_record"], barcode))
                        elif args.moves_out:
                            m = read["move_table"].tolist()
                            move_str = ','.join(map(str, m))
                            if args.do_read_splitting:
                                bc_writer.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tmv:B:c,{},{}\tqs:f:{}\tpi:Z:{}\tBC:Z:{}\n".format(read["read_id"], read["sequence"], read["qscore"], read["model_stride"], move_str, read["float_read_qscore"], read["parent_read_id"], barcode))
                            else:
                                # do ns and ts tags
                                bc_writer.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tmv:B:c,{},{}\tqs:f:{}\tns:i:{}\tts:i:{}\tBC:Z:{}\n".format(read["read_id"], read["sequence"], read["qscore"], read["model_stride"], move_str, read["float_read_qscore"], read["num_samples"], read["trimmed_samples"], barcode))
                        else:
                            if args.do_read_splitting:
                                bc_writer.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tqs:f:{}\tpi:Z:{}\tBC:Z:{}\n".format(read["read_id"], read["sequence"], read["qscore"], read["float_read_qscore"], read["parent_read_id"], barcode))
                            else:
                                # do ns and ts tags
                                bc_writer.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tqs:f:{}\tns:i:{}\tts:i:{}\tBC:Z:{}\n".format(read["read_id"], read["sequence"], read["qscore"], read["float_read_qscore"], read["num_samples"], read["trimmed_samples"], barcode))

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
    
    print("Total reads: {}".format(total_reads))
    
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
        if args.above_7412:
            sam_tags = "MN:i:{}\tqs:f:{}\tmx:i:{}\tch:i:{}\tns:i:{}\tts:i:{}\tsm:f:{}\tsd:f:{}\tsv:Z:{}\tdu:f:{}".format(read["sequence_length"],
                                                                                                                        read["float_read_qscore"],
                                                                                                                        read["mux"],
                                                                                                                        read["channel"],
                                                                                                                        read["num_samples"],
                                                                                                                        read["trimmed_samples"],
                                                                                                                        read["scaling_median"],
                                                                                                                        read["scaling_med_abs_dev"],
                                                                                                                        read["scaling_version"],
                                                                                                                        read["duration"])
            if args.duplex:
                duplex_tag = "0"
                if read["duplex_strand_1"] is not None:
                    duplex_tag = "1"
                if read["duplex_parent"]:
                    duplex_tag = "-1"
                if args.call_mods:
                    OUT.write("{}\tpi:Z:{}\tdx:i:{}\n".format(read["sam_record"], read["parent_read_id"], duplex_tag))
                elif args.moves_out:
                    m = read["move_table"].tolist()
                    move_str = ','.join(map(str, m))
                    OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tmv:B:c,{},{}\tpi:Z:{}\t{}\tdx:i:{}\n".format(read_id, read["sequence"], read["qscore"], read["model_stride"], move_str, read["parent_read_id"], sam_tags, duplex_tag))
                else:
                    OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tpi:Z:{}\tqs:f:{}\tdx:i:{}\n".format(read_id, read["sequence"], read["qscore"], read["parent_read_id"], read["float_read_qscore"], duplex_tag))
            else:
                if args.call_mods:
                    OUT.write("{}\tpi:Z:{}\n".format(read["sam_record"], read["parent_read_id"]))
                elif args.moves_out:
                    m = read["move_table"].tolist()
                    move_str = ','.join(map(str, m))
                    OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tmv:B:c,{},{}\tpi:Z:{}\t{}\n".format(read_id, read["sequence"], read["qscore"], read["model_stride"], move_str, read["parent_read_id"], sam_tags))
                else:
                    OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tpi:Z:{}\t{}\n".format(read_id, read["sequence"], read["qscore"], read["parent_read_id"], sam_tags))
        else:
            if args.call_mods:
                if args.do_read_splitting:
                    OUT.write("{}\tpi:Z:{}\n".format(read["sam_record"], read["parent_read_id"]))
                else:
                    OUT.write("{}\n".format(read["sam_record"]))
            elif args.moves_out:
                m = read["move_table"].tolist()
                move_str = ','.join(map(str, m))
                if args.do_read_splitting:
                    OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tmv:B:c,{},{}\tpi:Z:{}\tqs:f:{}\n".format(read_id, read["sequence"], read["qscore"], read["model_stride"], move_str, read["parent_read_id"], read["float_read_qscore"]))
                else:
                    # do ns and ts tags
                    OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tmv:B:c,{},{}\tqs:f:{}\tns:i:{}\tts:i:{}\n".format(read_id, read["sequence"], read["qscore"], read["model_stride"], move_str, read["float_read_qscore"], read["num_samples"], read["trimmed_samples"]))
            else:
                if args.do_read_splitting:
                    OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tpi:Z:{}\tqs:f:{}\n".format(read_id, read["sequence"], read["qscore"], read["parent_read_id"], read["float_read_qscore"]))
                else:
                    # do ns and ts tags
                    OUT.write("{}\t4\t*\t0\t0\t*\t*\t0\t0\t{}\t{}\tqs:f:{}\tns:i:{}\tts:i:{}\n".format(read_id, read["sequence"], read["qscore"], read["float_read_qscore"], read["num_samples"], read["trimmed_samples"]))
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
            # TODO: Add duplex read count here as well with (%)
            sys.stdout.flush()
        # don't make div larger than 500K
        if total_reads >= div*10 and div <= 50000:
            div = div*10
