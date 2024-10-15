#!/bin/bash

# MIT License

# Copyright (c) 2023 Hasindu Gamaarachchi
# Copyright (c) 2023 James Ferguson

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

die(){
    echo "$@" >&2
    exit 1
}

GUPPY_VERSION=7.1.4
CURRENT_GUPPY=$(grep "ont-pyguppy-client-lib" requirements.txt | cut -d "=" -f 3)
GUPPY_LIB=$(grep "ont-pyguppy-client-lib" requirements.txt | cut -d "=" -f 1)
DORADO_LIB=$(grep "ont-pybasecall-client-lib" requirements.txt | cut -d "=" -f 1)
test -z ${CURRENT_GUPPY} && die "ont-pyguppy-client-lib version not found in requirements.txt"
test -z ${GUPPY_LIB} && die "ont-pyguppy-client-lib not found in requirements.txt"
sed -i "s/${GUPPY_LIB}/ont-pyguppy-client-lib/" requirements.txt || die "sed failed"
sed -i "s/${CURRENT_GUPPY}/${GUPPY_VERSION}/" requirements.txt || die "sed failed"
sed -i "s/${DORADO_LIB}/#${DORADO_LIB}/" requirements.txt || die "sed failed"

test -z $EEL_PYTHON3 && EEL_PYTHON3=python3
rm -rf venv3
${EEL_PYTHON3} -m venv venv3 || die "venv failed"
source venv3/bin/activate || die "venv activate failed"
pip3 install --upgrade pip || die "pip upgrade failed"
pip3 install . || die "pip install failed"
buttery-eel --version || die "buttery-eel --version failed"
deactivate || die "venv deactivate failed"

