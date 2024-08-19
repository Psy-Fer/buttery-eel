import argparse
import sys
from pathlib import Path
from ._version import __version__



class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def get_args(above_7310_flag, above_7412_flag):

    VERSION = __version__
    parser = MyParser(description="buttery-eel - wrapping ONT basecallers (guppy/dorado) for SLOW5 basecalling",
    epilog="Citation: Hiruna Samarakoon, James M Ferguson, Hasindu Gamaarachchi, Ira W Deveson, Accelerated nanopore basecalling with SLOW5 data format, Bioinformatics, Volume 39, Issue 6, June 2023, btad352, https://doi.org/10.1093/bioinformatics/btad352",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    if above_7310_flag:
        run_options = parser.add_argument_group("Run Options")
        seq_sum = parser.add_argument_group("Sequencing summary Options")
        read_splitting = parser.add_argument_group("Read splitting Options")
        adapter_trimming = parser.add_argument_group("Adapter trimming Options")
        barcode_dmux = parser.add_argument_group("Barcode demultiplexing Options")
        duplex = parser.add_argument_group("Duplex Options")

        # Args for the wrapper, and then probably best to just have free form args for basecaller
        run_options.add_argument("-i", "--input", required=True,
                            help="input blow5 file or directory for basecalling")
        run_options.add_argument("-o", "--output", required=True,
                            help="output .fastq or unaligned .sam file to write")
        run_options.add_argument("-g", "--basecaller_bin", type=Path, required=True,
                            help="path to basecaller bin folder, eg: ont-dorado-server/bin")
        run_options.add_argument("--config", default="dna_r10.4.1_e8.2_400bps_5khz_hac.cfg", required=True,
                            help="basecalling model config")
        run_options.add_argument("--call_mods", action="store_true",
                            help="output MM/ML tags for methylation - will output sam - use with appropriate mod config")
        run_options.add_argument("-q", "--qscore", type=int,
                            help="A mean q-score to split fastq/sam files into pass/fail output")
        run_options.add_argument("--slow5_threads", type=int, default=4,
                            help="Number of threads to use reading slow5 file")
        run_options.add_argument("--procs", type=int, default=4,
                            help="Number of worker processes to use processing reads")
        run_options.add_argument("--slow5_batchsize", type=int, default=4000,
                            help="Number of reads to process at a time reading slow5")
        run_options.add_argument("--quiet", action="store_true",
                            help="Don't print progress")
        run_options.add_argument("--max_read_queue_size", type=int, default=20000,
                            help="Number of reads to process at a time reading slow5")
        run_options.add_argument("--log", default="buttery_basecaller_logs",
                            help="basecaller log folder path")
        run_options.add_argument("--moves_out", action="store_true",
                            help="output move table (sam format only)")

        # read splitting
        # read_splitting.add_argument("--do_read_splitting", action="store_true",
        #                     help="Perform read splitting based on mid-strand adapter detection - On by default dorado-server >= v7.3.10")
        # read_splitting.add_argument("--min_score_read_splitting", type=float, default=50.0,
        #                     help="Minimum mid-strand adapter score for reads to be split")
        
        # Adapter trimming
        # adapter_trimming.add_argument("--detect_adapter", action="store_true",
        #                     help="Enable detection of adapters at the front and rear of the sequence")
        # adapter_trimming.add_argument("--min_score_adapter", type=float, default=60.0,
        #                     help="Minimum score for a front or rear adapter to be classified. Default is 60.")
        adapter_trimming.add_argument("--trim_adapters", action="store_true",
                            help="Flag indicating that adapters should be trimmed. Default is False.")
        # adapter_trimming.add_argument("--detect_mid_strand_adapter", action="store_true",
        #                     help="Flag indicating that read will be marked as unclassified if the adapter sequence appears within the strand itself. Default is False.")

        # Sequencing Summary file
        seq_sum.add_argument("--seq_sum", action="store_true",
                            help="Write out sequencing_summary.txt file")
        
        # barcode demultiplexing/trimming
        barcode_dmux.add_argument("--barcode_kits", action="append",
                            help="Strings naming each barcode kit to use. Default is to not do barcoding.")
        barcode_dmux.add_argument("--enable_trim_barcodes", action="store_true",
                            help="Flag indicating that barcodes should be trimmed.")
        barcode_dmux.add_argument("--require_barcodes_both_ends", action="store_true",
                            help="Flag indicating that barcodes must be at both ends.")
        # barcode_dmux.add_argument("--detect_mid_strand_barcodes", action="store_true",
        #                     help="Flag indicating that read will be marked as unclassified if barcodes appear within the strand itself.")
        # barcode_dmux.add_argument("--min_score_barcode_front", type=float, default=60.0,
        #                     help="Minimum score for a front barcode to be classified")
        # barcode_dmux.add_argument("--min_score_barcode_rear", type=float, default=60.0,
        #                     help="Minimum score for a rear barcode to be classified")
        # barcode_dmux.add_argument("--min_score_barcode_mid", type=float, default=60.0,
        #                     help="Minimum score for mid barcodes to be detected")
        
        # Duplex
        duplex.add_argument("--duplex", action="store_true",
                            help="Turn on duplex calling - channel based - NOT WORKING JUST YET")
        duplex.add_argument("--single", action="store_true",
                            help="use only a single proc for testing - DUPLEX TESTING")
        old_args = argparse.Namespace(
            do_read_splitting=False,
            detect_adapter=False,
            detect_mid_strand_adapter=False,
            detect_mid_strand_barcodes=False,
            )
        
    else:
        run_options = parser.add_argument_group("Run Options")
        seq_sum = parser.add_argument_group("Sequencing summary Options")
        read_splitting = parser.add_argument_group("Read splitting Options")
        adapter_trimming = parser.add_argument_group("Adapter trimming Options")
        barcode_dmux = parser.add_argument_group("Barcode demultiplexing Options")
        duplex = parser.add_argument_group("Duplex Options")

        # Args for the wrapper, and then probably best to just have free form args for basecaller
        run_options.add_argument("-i", "--input", required=True,
                            help="input blow5 file or directory for basecalling")
        run_options.add_argument("-o", "--output", required=True,
                            help="output .fastq or unaligned .sam file to write")
        run_options.add_argument("-g", "--basecaller_bin", type=Path, required=True,
                            help="path to basecaller bin folder, eg: ont-dorado-server/bin")
        run_options.add_argument("--config", default="dna_r9.4.1_450bps_fast.cfg", required=True,
                            help="basecalling model config")
        run_options.add_argument("--call_mods", action="store_true",
                            help="output MM/ML tags for methylation - will output sam - use with appropriate mod config")
        run_options.add_argument("-q", "--qscore", type=int,
                            help="A mean q-score to split fastq/sam files into pass/fail output")
        run_options.add_argument("--slow5_threads", type=int, default=4,
                            help="Number of threads to use reading slow5 file")
        run_options.add_argument("--procs", type=int, default=4,
                            help="Number of worker processes to use processing reads")
        run_options.add_argument("--slow5_batchsize", type=int, default=4000,
                            help="Number of reads to process at a time reading slow5")
        run_options.add_argument("--quiet", action="store_true",
                            help="Don't print progress")
        run_options.add_argument("--max_read_queue_size", type=int, default=20000,
                            help="Number of reads to process at a time reading slow5")
        run_options.add_argument("--log", default="buttery_basecaller_logs",
                            help="basecaller log folder path")
        run_options.add_argument("--moves_out", action="store_true",
                            help="output move table (sam format only)")

        # read splitting
        read_splitting.add_argument("--do_read_splitting", action="store_true",
                            help="Perform read splitting based on mid-strand adapter detection - On by default dorado-server >= v7.3.10")
        read_splitting.add_argument("--min_score_read_splitting", type=float, default=50.0,
                            help="Minimum mid-strand adapter score for reads to be split")
        
        # Adapter trimming
        adapter_trimming.add_argument("--detect_adapter", action="store_true",
                            help="Enable detection of adapters at the front and rear of the sequence")
        adapter_trimming.add_argument("--min_score_adapter", type=float, default=60.0,
                            help="Minimum score for a front or rear adapter to be classified. Default is 60.")
        adapter_trimming.add_argument("--trim_adapters", action="store_true",
                            help="Flag indicating that adapters should be trimmed. Default is False.")
        adapter_trimming.add_argument("--detect_mid_strand_adapter", action="store_true",
                            help="Flag indicating that read will be marked as unclassified if the adapter sequence appears within the strand itself. Default is False.")

        # Sequencing Summary file
        seq_sum.add_argument("--seq_sum", action="store_true",
                            help="Write out sequencing_summary.txt file")
        
        # barcode demultiplexing/trimming
        barcode_dmux.add_argument("--barcode_kits", action="append",
                            help="Strings naming each barcode kit to use. Default is to not do barcoding.")
        barcode_dmux.add_argument("--enable_trim_barcodes", action="store_true",
                            help="Flag indicating that barcodes should be trimmed.")
        barcode_dmux.add_argument("--require_barcodes_both_ends", action="store_true",
                            help="Flag indicating that barcodes must be at both ends.")
        barcode_dmux.add_argument("--detect_mid_strand_barcodes", action="store_true",
                            help="Flag indicating that read will be marked as unclassified if barcodes appear within the strand itself.")
        barcode_dmux.add_argument("--min_score_barcode_front", type=float, default=60.0,
                            help="Minimum score for a front barcode to be classified")
        barcode_dmux.add_argument("--min_score_barcode_rear", type=float, default=60.0,
                            help="Minimum score for a rear barcode to be classified")
        barcode_dmux.add_argument("--min_score_barcode_mid", type=float, default=60.0,
                            help="Minimum score for mid barcodes to be detected")
        
        # Duplex
        # duplex.add_argument("--duplex", action="store_true",
        #                     help="Turn on duplex calling - channel based - NOT WORKING JUST YET")
        # duplex.add_argument("--single", action="store_true",
        #                     help="use only a single proc for testing - DUPLEX TESTING")
        old_args = argparse.Namespace(
            duplex=False,
            single=False,
            )


    # parser.add_argument("--max_queued_reads", default="2000",
    #                     help="Number of reads to send to guppy server queue")
    # parser.add_argument("--chunk_size", default="2000",
    #                     help="signal chunk size, lower this for lower VRAM GPUs")
    parser.add_argument("--profile", action="store_true",
                        help="run cProfile on all processes - for debugging benchmarking")
    parser.add_argument("-v", "--version", action='version', version="buttery-eel - wrapping ONT basecallers (guppy/dorado) for SLOW5 basecalling version: {}".format(VERSION),
                        help="Prints version")
    # parser.add_argument("--debug", action="store_true",
    #                     help="Set logging to debug mode")

    # args = parser.parse_args()
    # This collects known and unknown args to parse the server config options
    args, other_server_args = parser.parse_known_args()

    # now merge them. This will all get printed into the arg print below which also helps with troubleshooting
    args = argparse.Namespace(**vars(args), **vars(old_args))

    # add super sneaky hidden flags the user can't interact with but makes global sharing easier
    extra_args = argparse.Namespace(
        above_7310=above_7310_flag, # is the version >= 7.3.* where the name and inputs change?
        above_7412=above_7412_flag,
    )

    # now merge them. This will all get printed into the arg print below which also helps with troubleshooting
    args = argparse.Namespace(**vars(args), **vars(extra_args))

    return args, other_server_args, parser.print_help


