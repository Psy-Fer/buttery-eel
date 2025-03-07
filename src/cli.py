import argparse
import sys
from pathlib import Path
from ._version import __version__



class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('ERROR: %s\n' % message)
        self.print_help()
        sys.exit(2)


def get_args():

    VERSION = __version__
    parser = MyParser(description="buttery-eel - An ONT dorado-server basecalling wrapper for slow5 files",
    epilog="Citation: Hiruna Samarakoon, James M Ferguson, Hasindu Gamaarachchi, Ira W Deveson, Accelerated nanopore basecalling with SLOW5 data format, Bioinformatics, Volume 39, Issue 6, June 2023, btad352, https://doi.org/10.1093/bioinformatics/btad352",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    run_options = parser.add_argument_group("Run Options")
    seq_sum = parser.add_argument_group("Sequencing summary Options")
    adapter_trimming = parser.add_argument_group("Adapter trimming Options")
    barcode_dmux = parser.add_argument_group("Barcode demultiplexing Options")
    rna = parser.add_argument_group("RNA options")

    # Args for the wrapper, and then probably best to just have free form args for basecaller
    run_options.add_argument("-i", "--input", required=True,
                        help="input blow5 file or directory for basecalling")
    run_options.add_argument("-o", "--output", required=True,
                        help="output .fastq or unaligned .sam file to write")
    run_options.add_argument("-g", "--basecaller_bin", type=Path,
                        help="path to basecaller bin folder, eg: ont-dorado-server/bin")
    run_options.add_argument("--config", default="dna_r10.4.1_e8.2_400bps_5khz_hac.cfg",
                        help="basecalling model config")
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
    run_options.add_argument("--max_batch_time", type=int, default=1200,
                        help="Maximum seconds to wait for batch to be basecalled before killing basecalling. Used to detect locked states/hung servers. Default=1200 (20min)")
    run_options.add_argument("--resume", default=None,
                        help="Resume a sequencing run. fastq or sam input.")
    run_options.add_argument("--moves_out", action="store_true",
                        help="output move table (sam format only)")
    

    # Methylation calling
    run_options.add_argument("--call_mods", action="store_true",
                        help="output MM/ML tags for methylation - will output sam - use with appropriate mod config")
    
    # Adapter trimming
    adapter_trimming.add_argument("--trim_adapters", action="store_true",
                        help="Flag indicating that adapters should be trimmed. Default is False.")

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
    
    # RNA
    rna.add_argument("--U2T", action="store_true",
                        help="Convert Uracil (U) to Thymine (T) in direct RNA output")
    rna.add_argument("--estimate_poly_a", action="store_true",
                        help="Perform polyA/T tail length estimation")
    rna.add_argument("--poly_a_config",
                        help="Filename of a custom polyA/T configuration to use for estimation.")

    # General
    parser.add_argument("--profile", action="store_true",
                        help="run cProfile on all processes - for debugging benchmarking")
    parser.add_argument("-v", "--version", action='version', version="buttery-eel - wrapping ONT basecallers (guppy/dorado) for SLOW5 basecalling version: {}".format(VERSION),
                        help="Prints version")

    old_args = argparse.Namespace(
        do_read_splitting=True,
        detect_adapter=False,
        detect_mid_strand_adapter=False,
        detect_mid_strand_barcodes=False,
        )
    # args = parser.parse_args()
    # This collects known and unknown args to parse the server config options
    args, other_server_args = parser.parse_known_args()

    # now merge them. This will all get printed into the arg print below which also helps with troubleshooting
    args = argparse.Namespace(**vars(args), **vars(old_args))

    if "--dorado_model_path" in other_server_args:
        dorado_model_path_flag = True
    else:
        dorado_model_path_flag = False
    # add super sneaky hidden flags the user can't interact with but makes global sharing easier
    extra_args = argparse.Namespace(
        resume_run=False,
        dorado_model_path_flag=dorado_model_path_flag,
    )

    # now merge them. This will all get printed into the arg print below which also helps with troubleshooting
    args = argparse.Namespace(**vars(args), **vars(extra_args))

    return args, other_server_args, parser.print_help


