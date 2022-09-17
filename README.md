# buttery-eel

The buttery eel - A slow5 guppy basecaller wrapper

`buttery-eel` is a wrapper for `guppy`. It allows us to read [`SLOW5` files](https://github.com/hasindu2008/slow5tools), and send that data to [`guppy`](https://community.nanoporetech.com/downloads) to basecall. It requires matching versions of [`guppy`](https://community.nanoporetech.com/downloads) and [`ont-pyguppy-client-lib`](https://pypi.org/project/ont-pyguppy-client-lib/) to work.

You can download guppy here: https://community.nanoporetech.com/downloads. An ONT login is required to access that page, sorry no easy way around that one without legal headaches.


# Quick start

Using python3, preferably python3.7 to 3.9. Python 3.10 and higher does not yet have any pip wheel builds available for v6.3.7 and lower of guppy

Install a version of `guppy` (something higher than 4) where `GUPPY_VERSION` is the version, for example, `6.3.7`

Download: https://community.nanoporetech.com/downloads

The `guppy` and `ont-pyguppy-client-lib` versions need to match

    # if GUPPY_VERSION=6.3.7
    # modify requirements.txt to have:
    #   ont-pyguppy-client-lib==6.3.7


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

    python setup.py install

    buttery-eel --help


Usage:

    usage: buttery-eel [-h] -i INPUT -o OUTPUT -g GUPPY_BIN --config CONFIG [--call_mods] [--log LOG] [-v]

    buttery-eel - wrapping guppy for file agnostic basecalling

    optional arguments:
    -h, --help            show this help message and exit
    -i INPUT, --input INPUT
                            input blow5 file for basecalling (default: None)
    -o OUTPUT, --output OUTPUT
                            output .fastq or unaligned .sam file to write (default: None)
    -g GUPPY_BIN, --guppy_bin GUPPY_BIN
                            path to ont_guppy/bin folder (default: None)
    --config CONFIG       basecalling model config (default: dna_r9.4.1_450bps_fast.cfg)
    --call_mods           output MM tag for methylation - will output sam - use with appropriate mod config (default: False)
    --log LOG             guppy log folder path (default: buttery_guppy_logs)
    -v, --version         Prints version


Set up flags needed and run (`--use_tcp` is needed but not forced in these early versions):

    buttery-eel -g ont-guppy-6.3.7/bin --use_tcp -x "cuda:all" --config dna_r9.4.1_450bps_fast.cfg --port 5558 -i PAF25452_pass_bfdfd1d8_11.blow5 -o test.fastq

To call modifications, provide a `modbases` model and the `--call_mods` flag. Output will now be unaligned-sam containing the `MM/ML` tags. It will also provide the move table

You must use guppy 6.3.0 or higher for mod calling

    buttery-eel -g ont-guppy-6.3.7/bin --use_tcp -x "cuda:all" --config dna_r9.4.1_450bps_modbases_5hmc_5mc_cg_fast.cfg --call_mods --port 5558 -i PAF25452_pass_bfdfd1d8_11.blow5 -o test.mod.sam


the `--config` file can be found using this command with guppy `guppy_basecaller --print_workflows` and looking up the appropriate kit and flowcell type. Specify the format like this `--config dna_r9.4.1_450bps_fast.cfg` ending in `.cfg`

## Aligning uSAM output and getting sorted bam using -y in minimap2

    samtools fastq -TMM,ML test.mod.sam | minimap2 -ax map-ont -y ref.fa - | samtools view -Sb - | samtools sort - > test.aln.mod.bam


# Shutting down server

If everything goes right, the server will be terminated at the end of the basecalling.

However, sometimes things go wrong, and the wrapper will temrinate before it terminates the server.

I have mostly fixed this but sometimes it still happens. Here is how you check for the server and then kill it.

    # check for guppy instanaces
    ps -ef | grep guppy

    # That might give you a result like this

    # hasindu  27946 27905 99 19:31 pts/22   01:25:29 /install/ont-guppy-6.3.7/bin/guppy_basecall_server --log_path buttery_guppy_logs --config dna_r9.4.1_450bps_hac_prom.cfg --port 5558 --use_tcp -x cuda:all --max_queued_reads 2000 --chunk_size 2000

    # using the --port to see that it is indeed the one you started.
    # you can then kill the process with, where in this case, `PID=27946`

    kill <PID>

    # then you can try again


# Info

ONT have 4 basecallers.
- `Albacore` (archived)
- `Guppy`    (current - production)
- `Bonito`   (research)
- `Dorado`   (preview - future production)

`Albacore` (not used anymore) and `Guppy` are closed source software, and researchers are required to sign a developer agreement to gain access to the source code. Any changes made by the researchers can't share that modified software, and there are a number of restrictions in place that are too boring to talk about.

`Bonito` and `Dorado` on the other hand, are MOSTLY open source. I say MOSTLY, because they have closed source dependencies, such as `koi`, and their license isn't really "true" open source. Either way, it's close enough. This is great, because any modifications a researcher makes to those software, can be shared with the wider community.

We have had an interest in the basecallers as we are the developers of `SLOW5`, a file format that is smaller and faster than `FAST5`. While ONT is developing `POD5`, there is a gap, where anyone using `SLOW5` for downstream analysis and storage, have to convert their `SLOW5` files back to `FAST5` to re-basecall with newer versions of `guppy`, bonito, or dorado. While we have built versions of `guppy` that can read `SLOW5`, we cannot distribute these due to developer agreement. We have built forks of both bonito and dorado and submitted pull requests on their official repositories, they were rejected by the ONT devs.

So there is a major road block for ease of use with `SLOW5` for re-basecalling data, especially with `guppy`, as that is the current production basecaller. I think this is something that stands in the way of more widespread adoption. While we can maintain our own forks of bonito and dorado as the changes to read `SLOW5` are minimal, there is no work around to share a `guppy` version that reads `SLOW5`....UNTIL NOW!

Adaptive sampling has really pushed the need for interfaces with basecalling in real-time, and a python interface has been built to help control `guppy_basecaller_client` and `guppy_basecall_server`. There is an open library from ONT called `pyguppyclient`, however I had some issues getting it to work, and support seems patchy. However it uses something called `ont-pyguppy-client-lib` which is built and released with every `guppy` release, and so is up to date with the latest versions. So, if we match the python library version, with the `guppy` release version, they should be compatible.

So using this library and some helpful hints from the readfish/ReadUntil basecalling scripts by Matt Loose and Alexander Payne, `buttery-eel` starts a `guppy_basecall_server`, connects to it as a client, reads a `SLOW5` file to get reads, converts them to a form to be accepted by `guppy`, and sends them to the `guppy_basecall_server`. It then collects the basecalled reads, and writes the data to a fastq. Once all the reads are processed, it sends a termination command to the `guppy_basecall_server` and that's it.

The best thing about this, is all of the libraries and code is open, and so we can share a method that doesn't cause any legal headaches, and also doesn't require the ONT devs to accept any pull requests or code changes. This just wraps around any `guppy` release like a "buttery eel", and you can use `SLOW5`.


# Acknowledgments

- Firstly, whoever maintains and develops the `ont-pyguppy-client-lib` at ONT, thank you. This library gives the minimum required to make this relatively easy to implement.
- Hasindu Gamaarachchi for having the idea to do this, and for the issue posted on the slow5tools repo by SziKayLeung
- My partner Hilary for coming up with the name.
- Matt Loose and Alexander Payne for having the basics of this all along in your readfish code and being awesome in general
- ONT and their open source code of bonito and dorado for handling uSAM writing.
- Lastly, I'd like to say i'm a little surprised this wasn't suggested to us by the devs at ONT when they were rejecting our pull requests on Guppy, Bonito, and Dorado. Oh well.

# Software used
- [slow5lib/pyslow5](https://github.com/hasindu2008/slow5lib)
- [ONT guppy]()
- [ONT ont-pyguppy-client-lib](https://pypi.org/project/ont-pyguppy-client-lib/6.2.1/)
- basecaller code and flow mostly follows the methods used in [readfish](https://github.com/LooseLab/readfish/blob/23dd37117bce576b99caf097e7711dc87d30fa0a/ru/basecall.py) by Matt Loose and Alexander Payne
