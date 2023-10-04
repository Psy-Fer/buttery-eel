#!/bin/bash

die() {
    echo "Error: $@" >&2
    exit 1
}

instal_venv(){
    test -d ${PATH_TO_EEL_VENV} && rm -r ${PATH_TO_EEL_VENV}
    python3.8 -m venv ${PATH_TO_EEL_VENV} || die "Failed to create venv"
    source ${PATH_TO_EEL_VENV} || die "Failed to source ${PATH_TO_EEL_VENV}"
    pip install --upgrade pip || die "Failed to upgrade pip"
    pip install --upgrade setuptools wheel || die "Failed to upgrade setuptools"
    python setup.py install || die "Failed to install buttery-eel"
    deactivate
}

CURRENT_GUPPY=$(grep "ont-pyguppy-client-lib" requirements.txt | cut -d "=" -f 3)
test -z ${CURRENT_GUPPY} && die "ont-pyguppy-client-lib not found in requirements.txt"
export PATH_TO_GUPPY=/install/ont-guppy-${CURRENT_GUPPY}/bin/
export GUPPY_OUT_TMP=ont-guppy-tmp
export EEL_OUT_TMP=buttery_eel_tmp.fastq


export PATH_TO_EEL_VENV=./venv3/bin/activate


export PATH_TO_IDENTITY=/install/biorand/bin/identitydna.sh
export REFIDX=/genome/hg38noAlt.idx

echo "R9.4.1 DNA - FAST model - 20k reads"
export PATH_TO_FAST5=/data/slow5-testdata/NA12878_prom_subsubsample/fast5/
export PATH_TO_BLOW5=/data/slow5-testdata/NA12878_prom_subsubsample/reads.blow5
export MODEL=dna_r9.4.1_450bps_fast_prom.cfg
test/test.sh || die "test failed"
echo ""
echo "********************************************************************"

echo "R10.4.1 DNA - HAC model - 20k reads - split qscore"
export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/fast5/
export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/reads.blow5
export MODEL=dna_r10.4.1_e8.2_400bps_hac_prom.cfg
test/test_qscore_split.sh || die "test failed"
echo ""
echo "********************************************************************"

echo "split reads"
export OPTS_GUPPY="--detect_mid_strand_adapter --trim_adapters --detect_adapter --do_read_splitting --trim_strategy dna"
export OPTS_EEL=$OPTS_GUPPY
test/test.sh || die "test failed"
echo ""
echo "********************************************************************"

echo "R10.4.1 DNA - FAST model - 500k reads"
export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/fast5/
export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/reads.blow5
export MODEL=dna_r10.4.1_e8.2_400bps_fast_prom.cfg
test/test.sh || die "test failed"
echo ""
echo "********************************************************************"

# remora?
# move table


#rna
export PATH_TO_IDENTITY=/install/biorand/bin/identityrna.sh
export REFIDX=/genome/gencode.v40.transcripts.fa
export PATH_TO_FAST5=/data/hasindu/hasindu2008.git/f5c/test/rna/
export PATH_TO_BLOW5=/data/hasindu/hasindu2008.git/f5c/test/rna/reads.blow5
export MODEL=rna_r9.4.1_70bps_fast_prom.cfg
test/test.sh || die "test failed"
