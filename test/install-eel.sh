#!/bin/bash

die(){
    echo "$@" >&2
    exit 1
}

GUPPY_VERSION=7.1.4 #>6 is dorado

git clone https://github.com/Psy-Fer/buttery-eel -b main || die "git clone failed"
mv buttery-eel buttery-eel-0.3.1+${GUPPY_VERSION} || die "mv failed"
cd buttery-eel-0.3.1+${GUPPY_VERSION} || die "cd failed"
#git checkout 855fff2 || die "git checkout failed"
sed -i "s/6.3.8/${GUPPY_VERSION}/" requirements.txt || die "sed failed"

python3 -m venv venv3 || die "venv failed"
source venv3/bin/activate || die "venv activate failed"
pip3 install --upgrade pip || die "pip upgrade failed"
pip3 install . || die "pip install failed"
buttery-eel --version || die "buttery-eel --version failed"
deactivate || die "venv deactivate failed"
cd .. || die "cd failed"
