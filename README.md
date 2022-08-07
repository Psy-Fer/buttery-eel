# buttery-eel
The buttery eel - Wrapping guppy for your file agnostic basecalling needs

`buttery-eel` is a wrapper for `guppy`. It allows us to read any data format, and send that data to `guppy` to basecall.

It requires matching versions of `guppy` and `ont-pyguppy-client-lib` to work.

You can download guppy here: https://community.nanoporetech.com/downloads

An ONT login is required to access that page, sorry no easy way around that one without legal headaches.


# Quick start

Using python3, preferably python3.7 or higher

Install a version of `guppy` (something higher than 4) where `GUPPY_VERSION` is the version, for example, `6.1.3`

Download: https://community.nanoporetech.com/downloads

The `guppy` and `ont-pyguppy-client-lib` versions need to match

    GUPPY_VERSION=6.1.3
    git clone https://github.com/Psy-Fer/buttery-eel.git
    cd butter-eel
    python3 -m venv venv3
    source ./venv3/bin/activate
    pip install --upgrade pip
    pip install pyslow5==0.5.0 ont-pyguppy-client-lib==${GUPPY_VERSION}


Set up flags needed and run

    ./buttery-eel.py --guppy_bin ont-guppy-6.1.3/bin --port 5558 -i ~/Data/bench/1_slow5/PAF25452_pass_bfdfd1d8_11.blow5 -o ~/Data/bench/buttery_test/test.fastq

# Info

ONT have 4 basecallers.
    - `Albacore` (archived)
    - `Guppy`    (current - production)
    - `Bonito`   (research)
    - `Dorado`   (preview - future production)

`Albacore` (not used anymore) and `Guppy` are closed source software, and researchers are required to sign a developer agreement to gain access to the source code. Any changes made by the researchers can't share that modified software, and there are a number of restrictions in place that are too boring to talk about.

`Bonito` and `Dorado` on the other hand, are MOSTLY open source. I say MOSTLY, because they have closed source dependencies, such as `koi`, and their license isn't really "true" open source. Either way, it's close enough. This is great, because any modifications a researcher makes to those software, can be shared with the wider community.

We have had an interest in the basecallers as we are the developers of `SLOW5`, a file format that is smaller and faster than `FAST5`. While ONT is developing `POD5`, there is a gap, where anyone using `SLOW5` for downstream analysis and storage, have to convert their `SLOW5` files back to `FAST5` to re-basecall with newer versions of `guppy`, bonito, or dorado. While we have built versions of `guppy` that can read `SLOW5`, we cannot distribute these due to developer agreement. We have built forks of both bonito and dorado and submitted pull requests on their official repositories, they were rejected by the ONT devs.

So there is a major road block for easy of use with `SLOW5` for re-basecalling data, especially with `guppy`, as that is the current production basecaller, and I think this is something that stands in the way of more widespread adoption. While we can maintain our own forks of bonito and dorado as the changes to read `SLOW5` are very minimal, there is no work around to share a `guppy` version that reads `SLOW5`....UNTIL NOW!

Adaptive sampling has really pushed the need for interfaces with basecalling in realtime, and a python interface has been built to help control `guppy_basecaller_client` and `guppy_basecall_server`. There is an open library from ONT called `pyguppyclient`, however I had some issues getting it to work, and support seems patchy. However it uses something called `ont-pyguppy-client-lib` which is is built and released with every `guppy` release, and so is up to date with the latest versions. So, if we match the python library version, with the `guppy` release version, they should be compatible.

So using this library and some helpful hints from the readfish/ReadUntil basecalling scripts by Matt Loose, `buttery-eel` starts a `guppy_basecall_server`, connects to it as a client, reads a `SLOW5` (or any format..more on that later) to get reads, converts them to a form to be accepted by `guppy`, and sends them to the `guppy_basecall_server`. It then collects the basecalled reads, and writes the data to a fastq. Once all the reads are processed, it sends a termination command to the `guppy_basecall_server` and that's it.

The best thing about this, is all of the libraries and code is open, and so we can share a method that doesn't cause any legal headaches, and also doesn't require the ONT devs to accept any pull requests/code changes. This just wraps around any `guppy` release like a "buttery eel", and you can use `SLOW5`.

# Other file formats

So currently, this only works with `SLOW5`. However, with only a few lines of code changed, it can easily work with `FAST5` or the new `POD5`. If there is sufficient interest in having access to more than just `SLOW5`, I can build these. (or feel free to submit a pull request)


# Acknowledgments

Firstly, whoever maintains and develops the `ont-pyguppy-client-lib` at ONT, thank you. This library gives the minimum required to make this relatively easy to implement.

Hasindu Gamaarachchi for having the idea to do this, and for the issue posted on the slow5tools repo by SziKayLeung

Matt Loose for having the basics of this all along in your readfish code and being awesome in general

Lastly, I'd like to say i'm a little surprised this wasn't suggested to us by the devs at ONT when they were rejecting our pull requests on Guppy, Bonito, and Dorado. Oh well.

# Software used
- slow5lib
- ONT guppy
- ONT ont-pyguppy-client-lib
- basecaller code and flow mostly follows the methods used here: https://github.com/LooseLab/readfish/blob/23dd37117bce576b99caf097e7711dc87d30fa0a/ru/basecall.py by Matt Loose
