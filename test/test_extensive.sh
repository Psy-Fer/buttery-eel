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

die() {
    echo "Error: $@" >&2
    exit 1
}


CURRENT_GUPPY=$(grep "ont-pyguppy-client-lib" requirements.txt | cut -d "=" -f 3)
test -z ${CURRENT_GUPPY} && die "ont-pyguppy-client-lib not found in requirements.txt"

export PATH_TO_GUPPY=/install/ont-guppy-${CURRENT_GUPPY}/bin/
export GUPPY_OUT_TMP=ont-guppy-tmp
export EEL_OUT_TMP=buttery_eel_tmp.fastq


export PATH_TO_EEL_VENV=./venv3/bin/activate

export PATH_TO_IDENTITY=/install/biorand/bin/identitydna.sh
export REFIDX=/genome/hg38noAlt.idx

echo "Installation"
test/test_install.sh || die "test failed"
echo ""
echo "********************************************************************"

echo "R9.4.1 DNA - FAST model - 20k reads"
export PATH_TO_FAST5=/data/slow5-testdata/NA12878_prom_subsubsample/fast5/
export PATH_TO_BLOW5=/data/slow5-testdata/NA12878_prom_subsubsample/reads.blow5
export MODEL=dna_r9.4.1_450bps_fast_prom.cfg
test/test.sh || die "test failed"
echo ""
echo "********************************************************************"

echo "R10.4.1 DNA - HAC model - 20k reads - split qscore inbuilt"
export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/fast5/
export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/reads.blow5
export MODEL=dna_r10.4.1_e8.2_400bps_hac_prom.cfg
test/test_qscore_split.sh || die "test failed"
echo ""
echo "********************************************************************"

echo "R10.4.1 DNA - FAST model - 20k reads - split qscore script"
export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/fast5/
export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/reads.blow5
export MODEL=dna_r10.4.1_e8.2_400bps_fast_prom.cfg
test/test_qscore_split2.sh || die "test failed"
echo ""
echo "********************************************************************"

echo "read splitting"
export OPTS_GUPPY="--detect_mid_strand_adapter --trim_adapters --detect_adapter --do_read_splitting --trim_strategy dna"
export OPTS_EEL=$OPTS_GUPPY
test/test.sh || die "test failed"
echo ""
echo "********************************************************************"

#rna
export PATH_TO_IDENTITY=/install/biorand/bin/identityrna.sh
export REFIDX=/genome/gencode.v40.transcripts.fa
export PATH_TO_FAST5=/data/hasindu/hasindu2008.git/f5c/test/rna/
export PATH_TO_BLOW5=/data/hasindu/hasindu2008.git/f5c/test/rna/reads.blow5
export MODEL=rna_r9.4.1_70bps_fast_prom.cfg
test/test.sh || die "test failed"


echo "R10.4.1 DNA - FAST model - 500k reads"
export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/fast5/
export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/reads.blow5
export MODEL=dna_r10.4.1_e8.2_400bps_fast_prom.cfg
test/test.sh || die "test failed"
echo ""
echo "********************************************************************"

# remora?
# move table


