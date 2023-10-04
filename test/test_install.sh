#!/bin/bash

die(){
    echo "$@" >&2
    exit 1
}

test -z $EEL_PYTHON3 && EEL_PYTHON3=python3
rm -rf venv3
${EEL_PYTHON3} -m venv venv3 || die "venv failed"
source venv3/bin/activate || die "venv activate failed"
pip3 install --upgrade pip || die "pip upgrade failed"
pip3 install . || die "pip install failed"
buttery-eel --version || die "buttery-eel --version failed"
deactivate || die "venv deactivate failed"

