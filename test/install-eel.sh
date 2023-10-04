#!/bin/bash

die(){
    echo "$@" >&2
    exit 1
}

GUPPY_VERSION=7.1.4 #>6 is dorado

test -z $EEL_PYTHON3 && EEL_PYTHON3=python3
git clone https://github.com/Psy-Fer/buttery-eel -b main || die "git clone failed"
mv buttery-eel buttery-eel-0.3.1+${GUPPY_VERSION} || die "mv failed"
cd buttery-eel-0.3.1+${GUPPY_VERSION} || die "cd failed"
#git checkout 855fff2 || die "git checkout failed"
CURRENT_GUPPY=$(grep "ont-pyguppy-client-lib" requirements.txt | cut -d "=" -f 3)
test -z ${CURRENT_GUPPY} && die "ont-pyguppy-client-lib not found in requirements.txt"
sed -i "s/${CURRENT_GUPPY}/${GUPPY_VERSION}/" requirements.txt || die "sed failed"

${EEL_PYTHON3} -m venv venv3 || die "venv failed"
source venv3/bin/activate || die "venv activate failed"
pip3 install --upgrade pip || die "pip upgrade failed"
pip3 install . || die "pip install failed"
buttery-eel --version || die "buttery-eel --version failed"
deactivate || die "venv deactivate failed"
cd .. || die "cd failed"
