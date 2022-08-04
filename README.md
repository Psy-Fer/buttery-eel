# buttery-eel
The buttery eel - Wrapping guppy for your file agnostic basecalling needs


# Quick start

clone
python3 -m venv venv3
source ./venv3/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

flags
run

    ./buttery-eel.py --guppy_bin ont-guppy-6.1.3/bin --port 5558 -i ~/Data/bench/1_slow5/PAF25452_pass_bfdfd1d8_11.blow5 -o ~/Data/bench/buttery_test/test.fastq

# Info

talk about need for tool


# Software used
- slow5lib
- ONT guppy
- ONT ont-pyguppy-client-lib
- basecaller code and flow mostly follows the methods used here: https://github.com/LooseLab/readfish/blob/23dd37117bce576b99caf097e7711dc87d30fa0a/ru/basecall.py by Matt Loose
