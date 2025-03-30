# buttery-eel
<div style="width: 80%; height: 80%">
  
  ![](/docs/buttery-eel-mascot-banner.svg)
  
</div>




## The buttery eel - A slow5 guppy/dorado basecaller wrapper

`buttery-eel` is a wrapper for `guppy` and `dorado`. It allows us to read [`SLOW5` files](https://github.com/hasindu2008/slow5tools), and send that data to [`guppy`] or `dorado` server (https://community.nanoporetech.com/downloads) to basecall. It requires matching versions of [`guppy/dorado`](https://community.nanoporetech.com/downloads) and [`ont-pyguppy-client-lib`](https://pypi.org/project/ont-pyguppy-client-lib/)/[`ont-pybasecall-client-lib`](https://pypi.org/project/ont-pybasecall-client-lib/) to work.

You can download guppy or dorado server here: https://community.nanoporetech.com/downloads. An ONT login is required to access that page, sorry no easy way around that one without legal headaches.

- Currently, the main branch is the multi-process version (parallel processes to communicate to/from Guppy client) that enables performance scaling for multi-GPU setups, especially for FAST basecalling or shorter reads. A simple single-process version (one process to communicate to/from the Guppy client) that works well for HAC and SUP models is available in the `singleproc` branch for learning purposes. 
- Before v0.3.3, the main branch was the single-process version (`singleproc` branch) and the multi=process version was under the `multiproc` branch.


# Quickstart

Using python3, preferably python3.8 to 3.12. Older versions of guppy will work with 3.7, but not with 3.10 or higher.

Install a version of `guppy` (something higher than 4) where `GUPPY_VERSION` is the version, for example, `6.3.8`. Alternatively, you can install a version of `dorado server` too.
Download: https://community.nanoporetech.com/downloads

The `guppy/dorado` and `ont-pyguppy-client-lib/ont-pybasecall-client-lib` versions need to match
```
git clone https://github.com/Psy-Fer/buttery-eel.git
cd buttery-eel
python3 -m venv venv3
source ./venv3/bin/activate
pip install --upgrade pip
pip install --upgrade setuptools wheel

# if your slow5 file uses zstd compression and you have zstd installed
# see slow5lib for more info
# set this first to ensure pyslow5 installs with zstd:
# export PYSLOW5_ZSTD=1
# newer versions of pyslow5 have zstd libs included by default

# if GUPPY_VERSION=6.3.8
# modify requirements.txt to have:
#   ont-pyguppy-client-lib==6.3.8
# if using DORADO_SERVER_VERSION=7.4.12
#   ont-pybasecall-client-lib==7.4.12

# Install using pip
pip install .

buttery-eel --help

```

Suppose the name of the virtual environment you created is venv3 and resides directly in the root of the cloned buttery-eel git repository. In that case, you can use the wrapper script available under `/path/to/repository/scripts/eel` for conveniently executing buttery-eel. This script will automatically source the virtual environment, find a free port, execute the buttery-eel with the parameters you specified and finally deactivate the virtual environment. If you add the path of `/path/to/repository/scripts/` to your PATH environment variable, you can simply use buttery-eel as:
```
eel -g /path/to/ont-guppy/bin/ --config dna_r10.4.1_e8.2_400bps_hac.cfg --device cuda:all -i reads.blow5 -o reads.reads # and any other parameters
```

Alternatively, you can manually execute buttery-eel if you have sourced the virtual environment. You must provide `--port PORT --use_tcp` parameters manually in this case. Example:
```
buttery-eel -g /path/to/ont-guppy/bin/ --config dna_r10.4.1_e8.2_400bps_hac.cfg --device cuda:all -i reads.blow5 -o reads.reads.fastq --port 5000 --use_tcp  # and any other parameters
```

furthermore, if you are using the latest version, and are using dorado-server backend, then simply set the port argumnet to `--port auto` and it will automatically find a free port for you.

## Breaking changes

all:
- `-g/--guppy_bin` is now `-g/--basecaller_bin` in buttery-eel. However `--guppy_bin` will still work.

v7.3.10
- ONT have removed the `--do_read_splitting` option, and it is now on by default
- `--detect_adapter` and a number of other barcode/adapter options were removed


## Duplex calling

#### Duplex looks to be depricated - leaving this for legacy sake

The duplex calling does work, so long as you provide a duplex model for `--config` and the `--duplex` flag.

However there are some things to note:
- Use a single blow5 file rather than many smaller ones.
- The basecalling server stores all the reads for 10 channels, then on the 11th, it releases the first. Buttery-eel sends 1 channel per client connection, controlled by `--procs`, and in order to force the basecaller to release the data, it sends 10 "fake" reads to the basecaller with channel numbers >9000. This is mostly due to the poor implementation of duplex in the ONT library, so I can't really do much about that.
- You should write duplex data out using `.sam`. This will mean you get the duplex tags, dx:i:N where N=0 is simplex, N=-1 is a parent of a duplex read, and N=1 is a duplex read.
- There is a bug in the ONT library, where if a read is split, and two reads from that split read are parents of a duplex read, one of those parent reads won't be flagged with dx:i:-1, but dx:i:0 instead. I have told ONT and they said they will fix it (Bug present in `ont-pybasecall-client-lib v7.4.12`, now fixed in `ont-pybasecall-client-lib v7.6.8`)
- When duplex first starts, it sequentially reads the whole blow5 file to create the channel groups. This can take a while, so please be patient.

I wouldn't recommend using duplex just yet because of the issues and poor performance.


# Usage

The `--help` shown will be different for different versions of the ont library installed.

```
usage: buttery-eel [-h] -i INPUT -o OUTPUT [-g BASECALLER_BIN] [--config CONFIG] [--call_mods] [-q QSCORE] [--slow5_threads SLOW5_THREADS] [--procs PROCS] [--slow5_batchsize SLOW5_BATCHSIZE] [--quiet]
                   [--max_read_queue_size MAX_READ_QUEUE_SIZE] [--log LOG] [--moves_out] [--max_batch_time MAX_BATCH_TIME] [--resume RESUME] [--trim_adapters] [--seq_sum] [--barcode_kits BARCODE_KITS]
                   [--enable_trim_barcodes] [--require_barcodes_both_ends] [--U2T] [--estimate_poly_a] [--poly_a_config POLY_A_CONFIG] [--duplex] [--single] [--profile] [-v]

buttery-eel - wrapping ONT basecallers (guppy/dorado) for SLOW5 basecalling

options:
  -h, --help            show this help message and exit
  --profile             run cProfile on all processes - for debugging benchmarking (default: False)
  -v, --version         Prints version

Run Options:
  -i INPUT, --input INPUT
                        input blow5 file or directory for basecalling (default: None)
  -o OUTPUT, --output OUTPUT
                        output .fastq or unaligned .sam file to write (default: None)
  -g BASECALLER_BIN, --basecaller_bin BASECALLER_BIN
                        path to basecaller bin folder, eg: ont-dorado-server/bin (default: None)
  --config CONFIG       basecalling model config (default: dna_r10.4.1_e8.2_400bps_5khz_hac.cfg)
  --call_mods           output MM/ML tags for methylation - will output sam - use with appropriate mod config (default: False)
  -q QSCORE, --qscore QSCORE
                        A mean q-score to split fastq/sam files into pass/fail output (default: None)
  --slow5_threads SLOW5_THREADS
                        Number of threads to use reading slow5 file (default: 4)
  --procs PROCS         Number of worker processes to use processing reads (default: 4)
  --slow5_batchsize SLOW5_BATCHSIZE
                        Number of reads to process at a time reading slow5 (default: 4000)
  --quiet               Don't print progress (default: False)
  --max_read_queue_size MAX_READ_QUEUE_SIZE
                        Number of reads to process at a time reading slow5 (default: 20000)
  --log LOG             basecaller log folder path (default: buttery_basecaller_logs)
  --moves_out           output move table (sam format only) (default: False)
  --max_batch_time MAX_BATCH_TIME
                        Maximum seconds to wait for batch to be basecalled before killing basecalling. Used to detect locked states/hung servers. Default=5000 (1.3h) (default: 5000)
  --resume RESUME       Resume a sequencing run. fastq or sam input. (default: None)

Sequencing summary Options:
  --seq_sum             Write out sequencing_summary.txt file (default: False)

Adapter trimming Options:
  --trim_adapters       Flag indicating that adapters should be trimmed. Default is False. (default: False)

Barcode demultiplexing Options:
  --barcode_kits BARCODE_KITS
                        Strings naming each barcode kit to use. Default is to not do barcoding. (default: None)
  --enable_trim_barcodes
                        Flag indicating that barcodes should be trimmed. (default: False)
  --require_barcodes_both_ends
                        Flag indicating that barcodes must be at both ends. (default: False)

RNA options:
  --U2T                 Convert Uracil (U) to Thymine (T) in direct RNA output (default: False)
  --estimate_poly_a     Perform polyA/T tail length estimation (default: False)
  --poly_a_config POLY_A_CONFIG
                        Filename of a custom polyA/T configuration to use for estimation. (default: None)

Duplex Options:
  --duplex              Turn on duplex calling - channel based - See README for information (default: False)
  --single              use only a single proc for testing - DUPLEX TESTING (default: False)

```

Set up flags needed and run (`--use_tcp` is needed but not forced in these early versions):

    buttery-eel -g ont-guppy-6.3.8/bin --use_tcp -x "cuda:all" --config dna_r9.4.1_450bps_fast.cfg --port 5558 -i PAF25452_pass_bfdfd1d8_11.blow5 -o test.fastq

To call modifications, provide a `modbases` model and the `--call_mods` flag. Output will now be unaligned-sam containing the `MM/ML` tags. It will also provide the move table

You must use guppy 6.3.0 or higher for mod calling

    buttery-eel -g ont-guppy-6.3.8/bin --use_tcp -x "cuda:all" --config dna_r9.4.1_450bps_modbases_5hmc_5mc_cg_fast.cfg --call_mods --port 5558 -i PAF25452_pass_bfdfd1d8_11.blow5 -o test.mod.sam

the `--config` file can be found using this command with guppy `guppy_basecaller --print_workflows` and looking up the appropriate kit and flowcell type. Specify the format like this `--config dna_r10.4.1_e8.2_400bps_5khz_hac.cfg` ending in `.cfg`

You can parse in any guppy/dordao flag, and buttery-eel will parse it through to the basecaller.

Some common args are:

- `-x/--device`: Specify CPU or GPU device: 'cpu', 'cuda:all', 'auto' or 'cuda:<device_id>[,<device_id> ...]', the default is 'cpu'. Specifying 'auto' will choose either 'cpu' or 'cuda:all' depending on the presence of a cuda device.
- `--use_tcp`: Make connections on a tcp port instead of a Unix socket file. This flag has no effect on Windows as connections are always via tcp.
- `-p/--port`:Port for hosting service. Specify "auto" to make server automatically search for a free port.
- `--dorado_model_path`/`--dorado_modbase_models`: Can be provided instead of using `--config` the same way you can define a model path in dorado. The modbase models can be added like this `--dorado_modbase_models 5mC_5hmC@v1,6mA@v2` where 2 models are given separated by a comma with no spaces.

Fore more information on model paths and modbase models, please refer to the models table in the dorado documentation [here](https://github.com/nanoporetech/dorado/tree/release-v0.9?tab=readme-ov-file#dna-models)

## Ultra long reads

If you have a run with ultra long reads, there is a case where you can run out of memory and get an OOM error, which will throw some errors that look like client timeout errors. If you are running on a cluster, this is difficult to catch. (it took us a while to figure out)

We found decreasing the `--slow5_batch_size` and the number of `--procs` helpped keep the ram usage low when there were lots of long reads.
4000 reads with a length of 5k bases each is MUCH smaller than 4000 reads with a length of 100k+ with a few reads over 1M thrown in.

We have found on HPC, instead of running with a batch size of 4000 and 20 procs, reducing to a batch size of 2000 and 10 procs uses a quarter of the RAM, and so can help get through the run, with minor performance impact. We were using a system with 4x Tesla V100-SXM2-32GB cards and 385GB of RAM
Depending on your system ram and GPUs you may need to alter this further to find the right balance.

## Resume a run

If a run did not complete, crashed or was interupted for some reason, you can resume the run with the `--resume <failed_run.fastq/sam>` flag and providing the fastq/sam file of the failed run. If you have multiple files created in the previous run due to barcoding and quality splitting, you can use the pattern `--resume file1.fastq,file2.fastq` separated by a comma and no space.

If the last record in a file is malformed, it will be skipped and a warning will be displayed. If there are more than 5 malformed records in a file, buttery-eel will exit with an error. This will most likely be caused by reads not having the `parent_read_id` field in fastq files or the `pi:z:` field in sam files.

Make sure the new run `-o/--output` filename is different to the file used in `--resume`, otherwise it will overwrite it and you will lose the previous run of data.

Once the resumed run is completed, you can merge the output files from incomplete run with their corresponding files in the resumed run.


### Estimate polyT/A tails

The arguments `--estimate_poly_a` and  `--poly_a_config` work the same way they do in dorado-server. Please see their documentation [here](https://github.com/nanoporetech/dorado/blob/release-v0.9/documentation/PolyTailConfig.md) for more information.


## Usage for older versions < 7.3.10
```
usage: buttery-eel [-h] -i INPUT -o OUTPUT [-g BASECALLER_BIN] --config CONFIG [--call_mods] [-q QSCORE] [--slow5_threads SLOW5_THREADS] [--procs PROCS] [--slow5_batchsize SLOW5_BATCHSIZE] [--quiet]
                   [--max_read_queue_size MAX_READ_QUEUE_SIZE] [--log LOG] [--moves_out] [--max_batch_time MAX_BATCH_TIME] [--resume RESUME] [--do_read_splitting] [--min_score_read_splitting MIN_SCORE_READ_SPLITTING]
                   [--detect_adapter] [--min_score_adapter MIN_SCORE_ADAPTER] [--trim_adapters] [--detect_mid_strand_adapter] [--seq_sum] [--barcode_kits BARCODE_KITS] [--enable_trim_barcodes]
                   [--require_barcodes_both_ends] [--detect_mid_strand_barcodes] [--min_score_barcode_front MIN_SCORE_BARCODE_FRONT] [--min_score_barcode_rear MIN_SCORE_BARCODE_REAR]
                   [--min_score_barcode_mid MIN_SCORE_BARCODE_MID] [--profile] [-v]

buttery-eel - wrapping ONT basecallers (guppy/dorado) for SLOW5 basecalling

options:
  -h, --help            show this help message and exit
  --profile             run cProfile on all processes - for debugging benchmarking (default: False)
  -v, --version         Prints version

Run Options:
  -i INPUT, --input INPUT
                        input blow5 file or directory for basecalling (default: None)
  -o OUTPUT, --output OUTPUT
                        output .fastq or unaligned .sam file to write (default: None)
  -g BASECALLER_BIN, --basecaller_bin BASECALLER_BIN
                        path to basecaller bin folder, eg: ont-dorado-server/bin (default: None)
  --config CONFIG       basecalling model config (default: dna_r9.4.1_450bps_fast.cfg)
  --call_mods           output MM/ML tags for methylation - will output sam - use with appropriate mod config (default: False)
  -q QSCORE, --qscore QSCORE
                        A mean q-score to split fastq/sam files into pass/fail output (default: None)
  --slow5_threads SLOW5_THREADS
                        Number of threads to use reading slow5 file (default: 4)
  --procs PROCS         Number of worker processes to use processing reads (default: 4)
  --slow5_batchsize SLOW5_BATCHSIZE
                        Number of reads to process at a time reading slow5 (default: 4000)
  --quiet               Don't print progress (default: False)
  --max_read_queue_size MAX_READ_QUEUE_SIZE
                        Number of reads to process at a time reading slow5 (default: 20000)
  --log LOG             basecaller log folder path (default: buttery_basecaller_logs)
  --moves_out           output move table (sam format only) (default: False)
  --max_batch_time MAX_BATCH_TIME
                        Maximum seconds to wait for batch to be basecalled before killing basecalling. Used to detect locked states/hung servers. Default=5000 (1.3h) (default: 5000)
  --resume RESUME       Resume a sequencing run. fastq or sam input. (default: None)

Sequencing summary Options:
  --seq_sum             Write out sequencing_summary.txt file (default: False)

Read splitting Options:
  --do_read_splitting   Perform read splitting based on mid-strand adapter detection - On by default dorado-server >= v7.3.10 (default: False)
  --min_score_read_splitting MIN_SCORE_READ_SPLITTING
                        Minimum mid-strand adapter score for reads to be split (default: 50.0)

Adapter trimming Options:
  --detect_adapter      Enable detection of adapters at the front and rear of the sequence (default: False)
  --min_score_adapter MIN_SCORE_ADAPTER
                        Minimum score for a front or rear adapter to be classified. Default is 60. (default: 60.0)
  --trim_adapters       Flag indicating that adapters should be trimmed. Default is False. (default: False)
  --detect_mid_strand_adapter
                        Flag indicating that read will be marked as unclassified if the adapter sequence appears within the strand itself. Default is False. (default: False)

Barcode demultiplexing Options:
  --barcode_kits BARCODE_KITS
                        Strings naming each barcode kit to use. Default is to not do barcoding. (default: None)
  --enable_trim_barcodes
                        Flag indicating that barcodes should be trimmed. (default: False)
  --require_barcodes_both_ends
                        Flag indicating that barcodes must be at both ends. (default: False)
  --detect_mid_strand_barcodes
                        Flag indicating that read will be marked as unclassified if barcodes appear within the strand itself. (default: False)
  --min_score_barcode_front MIN_SCORE_BARCODE_FRONT
                        Minimum score for a front barcode to be classified (default: 60.0)
  --min_score_barcode_rear MIN_SCORE_BARCODE_REAR
                        Minimum score for a rear barcode to be classified (default: 60.0)
  --min_score_barcode_mid MIN_SCORE_BARCODE_MID
                        Minimum score for mid barcodes to be detected (default: 60.0)

```

## Aligning uSAM output and getting sorted bam using -y in minimap2

    samtools fastq -TMM,ML test.mod.sam | minimap2 -ax map-ont -y -Y ref.fa - | samtools sort - > test.aln.mod.bam

If you also wish to keep the quality scores in the unofficial qs tags or if mapping a regular unmapped sam the -T argument can be used in conjunction with minimap2 -y for example: `-TMM,ML,qs` or `-Tqs`. You can also get all sam tags with `-T'*'` but you need samtools of v1.16 or higher.

#### samtools v1.16 and higher

    samtools fastq -T '*' test.mod.sam | minimap2 -ax map-ont -y -Y ref.fa - | samtools sort - > test.aln.mod.bam


# Shutting down server

If everything goes right, the server will be terminated at the end of the basecalling.

However, sometimes things go wrong, and the wrapper will temrinate before it terminates the server.

I have mostly fixed this but sometimes it still happens. Here is how you check for the server and then kill it.

```
# check for guppy instanaces
ps -ef | grep guppy

# That might give you a result like this

# hasindu  27946 27905 99 19:31 pts/22   01:25:29 /install/ont-guppy-6.3.8/bin/guppy_basecall_server --log_path buttery_guppy_logs --config dna_r9.4.1_450bps_hac_prom.cfg --port 5558 --use_tcp -x cuda:all --max_queued_reads 2000 --chunk_size 2000

# using the --port to see that it is indeed the one you started.
# you can then kill the process with, where in this case, `PID=27946`

kill <PID>

# then you can try again
```


# Info

ONT have 4 basecallers.
- `Albacore` (archived)
- `Guppy`    (archived)
- `Bonito`   (research)
- `Dorado`   (current - production)

`Albacore` (not used anymore) and `Guppy` are closed source software, and researchers are required to sign a developer agreement to gain access to the source code. Any changes made by the researchers can't share that modified software, and there are a number of restrictions in place that are too boring to talk about.

`Bonito` and `Dorado` on the other hand, are MOSTLY open source. I say MOSTLY, because they have closed source dependencies, such as `koi`, and their license isn't really "true" open source. Either way, it's close enough. This is great, because any modifications a researcher makes to those software, can be shared with the wider community.

We have had an interest in the basecallers as we are the developers of `SLOW5`, a file format that is smaller and faster than `FAST5`. While ONT is developing `POD5`, there is a gap, where anyone using `SLOW5` for downstream analysis and storage, have to convert their `SLOW5` files back to `FAST5` to re-basecall with newer versions of `guppy`, bonito, or dorado. While we have built versions of `guppy` that can read `SLOW5`, we cannot distribute these due to developer agreement. We have built forks of both bonito and dorado and submitted pull requests on their official repositories, they were rejected by the ONT devs.

So there is a major road block for ease of use with `SLOW5` for re-basecalling data, especially with `guppy`, as that is the current production basecaller. I think this is something that stands in the way of more widespread adoption. While we can maintain our own forks of bonito and dorado as the changes to read `SLOW5` are minimal, there is no work around to share a `guppy` version that reads `SLOW5`....UNTIL NOW!

Adaptive sampling has really pushed the need for interfaces with basecalling in real-time, and a python interface has been built to help control `guppy_basecaller_client` and `guppy_basecall_server`. There is an open library from ONT called `pyguppyclient`, however I had some issues getting it to work, and support seems patchy. However it uses something called `ont-pyguppy-client-lib` which is built and released with every `guppy` release, and so is up to date with the latest versions. So, if we match the python library version, with the `guppy` release version, they should be compatible.

So using this library and some helpful hints from the readfish/ReadUntil basecalling scripts by Matt Loose and Alexander Payne, `buttery-eel` starts a `guppy_basecall_server`, connects to it as a client, reads a `SLOW5` file to get reads, converts them to a form to be accepted by `guppy`, and sends them to the `guppy_basecall_server`. It then collects the basecalled reads, and writes the data to a fastq. Once all the reads are processed, it sends a termination command to the `guppy_basecall_server` and that's it.

The best thing about this, is all of the libraries and code is open, and so we can share a method that doesn't cause any legal headaches, and also doesn't require the ONT devs to accept any pull requests or code changes. This just wraps around any `guppy` release like a "buttery eel", and you can use `SLOW5`.

There is now a new python library/API for the dorado-server builds called `ont-pybasecall-client-lib`. It has replaced `ont-pyguppy-client-lib`. There are a number of breaking changes with this change, and I have done my best to handle them. Please let me know if something does work for you.

# Acknowledgments

- Firstly, whoever maintains and develops the `ont-pyguppy-client-lib` at ONT, thank you. This library gives the minimum required to make this relatively easy to implement.
- Hasindu Gamaarachchi for having the idea to do this, and for the issue posted on the slow5tools repo by SziKayLeung
- My partner Hilary for coming up with the name.
- Matt Loose and Alexander Payne for having the basics of this all along in your readfish code and being awesome in general
- ONT and their open source code of bonito and dorado for handling uSAM writing.
- Mark Bicknell at ONT who has been very helpful and patient with me while figuring out how to use new released features.
- Lastly, I'd like to say i'm a little surprised this wasn't suggested to us by the devs at ONT when they were rejecting our pull requests on Guppy, Bonito, and Dorado. Oh well.

# Software used
- [slow5lib/pyslow5](https://github.com/hasindu2008/slow5lib)
- [ONT guppy](https://community.nanoporetech.com/downloads)
- [ONT dorado-server](https://community.nanoporetech.com/downloads)
- [ONT ont-pyguppy-client-lib](https://pypi.org/project/ont-pyguppy-client-lib/)
- [ONT ont-pyguppy-client-lib](https://pypi.org/project/ont-pybasecall-client-lib/)
- basecaller code and flow mostly follows the methods used in [readfish](https://github.com/LooseLab/readfish/blob/23dd37117bce576b99caf097e7711dc87d30fa0a/ru/basecall.py) by Matt Loose and Alexander Payne (thought not so much anymore from ~0.7.0)


# Citation
 Hiruna Samarakoon, James M Ferguson, Hasindu Gamaarachchi, Ira W Deveson, Accelerated nanopore basecalling with SLOW5 data format, Bioinformatics, Volume 39, Issue 6, June 2023, btad352, https://doi.org/10.1093/bioinformatics/btad352",