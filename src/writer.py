import sys

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

def sam_header(OUT, sep='\t'):
    """
    Format a string sam header.
    This is taken from Bonito by Chris Seymour at ONT.
    https://github.com/nanoporetech/bonito/blob/master/bonito/io.py#L103
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
        'DS:ont basecaller wrapper',
    ])
    OUT.write("{}\n".format(HD))
    OUT.write("{}\n".format(PG1))
    OUT.write("{}\n".format(PG2))


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
                        if args.do_read_splitting:
                            bc_writer.write("{}\tpi:Z:{}\tBC:Z:{}\n".format(read["sam_record"], read["parent_read_id"], barcode))
                        else:
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
            if args.do_read_splitting:
                OUT.write("{}\tpi:Z:{}\n".format(read["sam_record"], read["parent_read_id"]))
            else:
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
