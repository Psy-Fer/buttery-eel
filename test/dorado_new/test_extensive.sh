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
    #exit 1
}

rm -f *.log

echo "Installation"
test/dorado_new/test_install.sh &> install.log || die "test failed. see install.log for details"
echo ""
echo "********************************************************************"

CURRENT_GUPPY=$(grep "ont-pybasecall-client-lib" requirements.txt | cut -d "=" -f 3)
test -z ${CURRENT_GUPPY} && die "ont-pybasecall-client-lib not found in requirements.txt"

export PATH_TO_GUPPY=/install/ont-dorado-server-${CURRENT_GUPPY}/bin/
export GUPPY_OUT_TMP=ont-guppy-tmp
export EEL_OUT_TMP=buttery_eel_tmp

export PATH_TO_EEL_VENV=./venv3/bin/activate

export PATH_TO_IDENTITY=/install/biorand/bin/identitydna.sh
export REFIDX=/genome/hg38noAlt.idx


echo "R9.4.1 DNA - FAST model - 20k reads"
export PATH_TO_FAST5=/data/slow5-testdata/NA12878_prom_subsubsample/fast5/
export PATH_TO_BLOW5=/data/slow5-testdata/NA12878_prom_subsubsample/reads.blow5
export MODEL=dna_r9.4.1_450bps_fast.cfg
test/dorado_new/test.sh &> r9_dna_fast.log || die "test failed. see r9_dna_fast.log for details"
echo ""
echo "********************************************************************"

echo "R10.4.1 DNA - FAST model - 20k reads - resume"
export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/pod5/
export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/reads.blow5
export MODEL=dna_r10.4.1_e8.2_400bps_5khz_fast.cfg
test/dorado_new/test_resume.sh &> r10_resume.log || die "test failed. see r10_resume.log for details"
echo ""
echo "********************************************************************"

echo "R10.4.1 DNA - HAC model - 20k reads - split qscore inbuilt"
export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/pod5/
export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/reads.blow5
export MODEL=dna_r10.4.1_e8.2_400bps_5khz_fast.cfg
test/dorado_new/test_qscore_split.sh &> r10_split1.log || die "test failed. see r10_split1.log for details"
echo ""
echo "********************************************************************"

echo "R10.4.1 DNA - FAST model - 20k reads - split qscore script"
export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/pod5/
export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/reads.blow5
export MODEL=dna_r10.4.1_e8.2_400bps_5khz_fast.cfg
test/dorado_new/test_qscore_split2.sh &> r10_split2.log || die "test failed. See r10_split2.log for details"
echo ""
echo "********************************************************************"

echo "SAM format qscore split script"
echo "Not yet implemented :("
echo ""
echo "********************************************************************"

echo "adapater trimming"
export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/pod5/
export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/reads.blow5
export OPTS_GUPPY="--trim_adapters"
export OPTS_EEL=$OPTS_GUPPY
test/dorado_new/test.sh &> r10_adaptertrim.log  || echo "test failed. See r10_adaptertrim.log for details"
unset OPTS_GUPPY
unset OPTS_EEL
echo ""
echo "********************************************************************"

# echo "read splitting"
# export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/fast5/
# export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/reads.blow5
# export OPTS_GUPPY="--do_read_splitting --min_score_read_splitting 50"
# export OPTS_EEL=$OPTS_GUPPY
# test/dorado_new/test.sh &> r10_readsplit.log  || echo "test failed. See r10_readsplit.log for details"
# unset OPTS_GUPPY
# unset OPTS_EEL
# echo ""
# echo "********************************************************************"

# echo "adapter trimming with read splitting"
# export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/fast5/
# export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/reads.blow5
# export OPTS_GUPPY="--detect_mid_strand_adapter --trim_adapters --detect_adapter --min_score_adapter 60 --do_read_splitting --min_score_read_splitting 50"
# export OPTS_EEL=$OPTS_GUPPY
# test/dorado_new/test.sh &> r10_readsplittrim.log  || echo "test failed. See r10_readsplittrim.log for details"
# unset OPTS_GUPPY
# unset OPTS_EEL
# echo ""
# echo "********************************************************************"

echo "seqsum"
export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/pod5/
export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/reads.blow5
test/dorado_new/test_seqsum.sh &> seqsum.log || die "test failed. See seqsum.log for details"
echo ""
echo "********************************************************************"

echo "seqsum - multiple BLOW5"
export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/pod5/
export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/reads.blow5
test/dorado_new/test_seqsum.sh &> seqsum_multiblow.log
echo ""
echo "********************************************************************"

echo "demux - FASTQ and SAM"
export PATH_TO_FAST5=/data/slow5-testdata/barcode_test/fast5/
export PATH_TO_BLOW5=/data/slow5-testdata/barcode_test/merged_rand.blow5
export MODEL=dna_r10.4.1_e8.2_400bps_fast.cfg
test/dorado_new/test_demux.sh &> demux.log || die "test failed. See demux.log for details"
echo ""
echo "********************************************************************"

echo "demux - qscore - FASTQ and SAM"
export PATH_TO_FAST5=/data/slow5-testdata/barcode_test/fast5/
export PATH_TO_BLOW5=/data/slow5-testdata/barcode_test/merged_rand.blow5
test/dorado_new/test_demux_qscore_split.sh &> demux_qscore.log  || die "test failed. See demux_qscore.log for details"
echo ""
echo "********************************************************************"

echo "demux - qscore - FASTQ and SAM - BARCODE+adapter trimming"
export PATH_TO_FAST5=/data/slow5-testdata/barcode_test/fast5/
export PATH_TO_BLOW5=/data/slow5-testdata/barcode_test/merged_rand.blow5
export OPTS_GUPPY="--trim_adapters "
export OPTS_BARCODER="--enable_trim_barcodes"
export OPTS_EEL=$OPTS_GUPPY" "$OPTS_BARCODER
test/dorado_new/test_demux_qscore_split.sh &> demux_qscore_trim.log  || die "test failed. See demux_qscore_trim.log for details"
unset OPTS_GUPPY
unset OPTS_EEL
echo ""
echo "********************************************************************"

echo "move table"
echo "Not yet implemented :("
echo ""
echo "********************************************************************"

echo "move table when adaptor/barcode trimming"
echo "Not yet implemented :("
echo ""
echo "********************************************************************"

echo "remora"
export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/pod5/
export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsubsample/reads.blow5
export MODEL=dna_r10.4.1_e8.2_400bps_5khz_modbases_5hmc_5mc_hac.cfg
test/dorado_new/test_remora.sh &> remora.log || die "test failed. See remora.log for details"
echo ""
echo "********************************************************************"

echo "remora with qscore split and dumux"
echo "Not yet implemented :("
echo ""
echo "********************************************************************"

echo "remora with adaptor/barcode trimming"
echo "Not yet implemented :("
echo ""
echo "********************************************************************"

echo "R10.4.1 DNA - FAST model - 500k reads"
export PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsample/pod5/
export PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_5khz_subsample/reads.blow5
export MODEL=dna_r10.4.1_e8.2_400bps_5khz_fast.cfg
test/dorado_new/test.sh &> dna_500k.log || die "test failed. See dna_500k.log for details"
echo ""
echo "********************************************************************"

echo "R9.4.1 RNA - FAST model"
export PATH_TO_IDENTITY=/install/biorand/bin/identityrna.sh
export REFIDX=/genome/gencode.v40.transcripts.fa
export PATH_TO_FAST5=/data/slow5-testdata/uhr_prom_rna002_subsubsample/fast5/
export PATH_TO_BLOW5=/data/slow5-testdata/uhr_prom_rna002_subsubsample/PRPN119035_reads_20k.blow5
export MODEL=rna_r9.4.1_70bps_fast.cfg
test/dorado_new/test.sh &> rna.log || die "test failed. See rna.log for details"

echo "RNA004 RNA - rna_rp4_130bps_fast"
export PATH_TO_IDENTITY=/install/biorand/bin/identityrna.sh
export REFIDX=/genome/gencode.v40.transcripts.fa
export PATH_TO_FAST5=/data/slow5-testdata/uhr_prom_rna004_subsubsample/pod5/
export PATH_TO_BLOW5=/data/slow5-testdata/uhr_prom_rna004_subsubsample/PNXRXX240011_reads_20k.blow5
export MODEL=rna_rp4_130bps_hac.cfg
test/dorado_new/test.sh &> rna.log || die "test failed. See rna.log for details"

echo "RNA004 RNA - rna_rp4_130bps_modbases_m6a_drach_sup.cfg"
echo "Not yet implemented :("
echo ""
echo "********************************************************************"
