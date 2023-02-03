#!/bin/bash

die() {
    echo "Error: $@" >&2
    exit 1
}

PATH_TO_GUPPY=/install/ont-guppy-6.4.2/bin/
PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/fast5/
PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/reads.blow5
PATH_TO_IDENTITY=/install/biorand/bin/identitydna.sh
PATH_TO_EEL_VENV=./venv3-multi-guppy-6.4.2/bin/activate
MODEL=dna_r10.4.1_e8.2_400bps_hac_prom.cfg
REFIDX=/genome/hg38noAlt.idx
GUPPY_OUT_TMP=ont-guppy-tmp
EEL_OUT_TMP=buttery_eel_tmp.fastq

test -e ${PATH_TO_GUPPY}/guppy_basecaller || die  "${PATH_TO_GUPPY}/guppy_basecaller not found"
test -d ${PATH_TO_FAST5} || die  "${PATH_TO_FAST5} not found"
test -e ${PATH_TO_BLOW5} || die  "${PATH_TO_BLOW5} not found"
test -e ${PATH_TO_IDENTITY} || die  "${PATH_TO_IDENTITY} not found"

test -d ${GUPPY_OUT_TMP} && rm -r ${GUPPY_OUT_TMP}
test -e ${EEL_OUT_TMP%.*}.pass.fastq && rm ${EEL_OUT_TMP%.*}.pass.fastq
test -e ${EEL_OUT_TMP%.*}.fail.fastq && rm ${EEL_OUT_TMP%.*}.fail.fastq

source ${PATH_TO_EEL_VENV} || die "Failed to source ${PATH_TO_EEL_VENV}"

echo "Running guppy"
${PATH_TO_GUPPY}/guppy_basecaller -c ${MODEL}  -i ${PATH_TO_FAST5} -s ${GUPPY_OUT_TMP}  -x cuda:all --recursive  --min_qscore 7
cat ${GUPPY_OUT_TMP}/pass/* > ${GUPPY_OUT_TMP}/reads_tmp.fastq
${PATH_TO_IDENTITY} ${REFIDX} ${GUPPY_OUT_TMP}/reads_tmp.fastq  | cut -f 2-> ${GUPPY_OUT_TMP}/reads_tmp.identity

echo "Running buttery-eel"
buttery-eel  -g ${PATH_TO_GUPPY}  --config ${MODEL} --device 'cuda:all' -i  ${PATH_TO_BLOW5} -o  ${EEL_OUT_TMP} --port 5555  --use_tcp --qscore 7
${PATH_TO_IDENTITY} ${REFIDX} ${EEL_OUT_TMP%.*}.pass.fastq | cut -f 2-> ${EEL_OUT_TMP}.identity

echo "Comparing results"
diff ${GUPPY_OUT_TMP}/reads_tmp.identity ${EEL_OUT_TMP}.identity || die "Results differ"

