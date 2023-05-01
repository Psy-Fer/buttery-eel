#!/bin/bash

die() {
    echo "Error: $@" >&2
    exit 1
}

get_port(){
#from https://unix.stackexchange.com/questions/55913/whats-the-easiest-way-to-find-an-unused-local-port
PORT=$(netstat -aln | awk '
  $6 == "LISTEN" {
    if ($4 ~ "[.:][0-9]+$") {
      split($4, a, /[:.]/);
      port = a[length(a)];
      p[port] = 1
    }
  }
  END {
    for (i = 5000; i < 65000 && p[i]; i++){};
    if (i == 65000) {exit 1};
    print i
  }
  ')
echo $PORT
}


#defaults if not set
test -z $PATH_TO_GUPPY && PATH_TO_GUPPY=/install/ont-guppy-6.4.2/bin/
test -z $PATH_TO_FAST5 && PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/fast5/
test -z $PATH_TO_BLOW5 && PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_subsubsample/reads.blow5
test -z $PATH_TO_IDENTITY && PATH_TO_IDENTITY=/install/biorand/bin/identitydna.sh
test -z $PATH_TO_EEL_VENV && PATH_TO_EEL_VENV=./venv3-multi-guppy-6.4.2/bin/activate
test -z $MODEL && MODEL=dna_r10.4.1_e8.2_400bps_hac_prom.cfg
test -z $REFIDX && REFIDX=/genome/hg38noAlt.idx
test -z $GUPPY_OUT_TMP && GUPPY_OUT_TMP=ont-guppy-tmp
test -z $EEL_OUT_TMP && EEL_OUT_TMP=buttery_eel_tmp.fastq

#checks
test -e ${PATH_TO_GUPPY}/guppy_basecaller || die  "${PATH_TO_GUPPY}/guppy_basecaller not found"
test -d ${PATH_TO_FAST5} || die  "${PATH_TO_FAST5} not found"
test -e ${PATH_TO_BLOW5} || die  "${PATH_TO_BLOW5} not found"
test -e ${PATH_TO_IDENTITY} || die  "${PATH_TO_IDENTITY} not found"

#cleanups
test -d ${GUPPY_OUT_TMP} && rm -r ${GUPPY_OUT_TMP}
test -e ${EEL_OUT_TMP%.*}.pass.fastq && rm ${EEL_OUT_TMP%.*}.pass.fastq
test -e ${EEL_OUT_TMP%.*}.fail.fastq && rm ${EEL_OUT_TMP%.*}.fail.fastq

#sourcing venv
source ${PATH_TO_EEL_VENV} || die "Failed to source ${PATH_TO_EEL_VENV}"

echo "Running guppy"
${PATH_TO_GUPPY}/guppy_basecaller -c ${MODEL}  -i ${PATH_TO_FAST5} -s ${GUPPY_OUT_TMP}  -x cuda:all --recursive  --min_qscore 7 ${OPTS_GUPPY}
cat ${GUPPY_OUT_TMP}/pass/* > ${GUPPY_OUT_TMP}/reads_tmp.fastq
${PATH_TO_IDENTITY} ${REFIDX} ${GUPPY_OUT_TMP}/reads_tmp.fastq  | cut -f 2-> ${GUPPY_OUT_TMP}/reads_tmp.identity

echo "Running buttery-eel"
PORT=$(get_port)
/usr/bin/time -v buttery-eel  -g ${PATH_TO_GUPPY}  --config ${MODEL} --device 'cuda:all' -i  ${PATH_TO_BLOW5} -o  ${EEL_OUT_TMP} --port ${PORT}  --use_tcp --qscore 7 ${OPTS_EEL} &> eel.log
MEM=$(grep "Maximum resident set size" eel.log | cut -d " " -f 6)
if [ $MEM -gt 8000000 ]; then
    die "Memory usage is too high: $MEM"
else
    echo "Memory usage is OK: $MEM"
fi
cat eel.log
${PATH_TO_IDENTITY} ${REFIDX} ${EEL_OUT_TMP%.*}.pass.fastq | cut -f 2-> ${EEL_OUT_TMP}.identity

echo "Comparing results"
diff ${GUPPY_OUT_TMP}/reads_tmp.identity ${EEL_OUT_TMP}.identity || die "Results differ"

echo "Test passed"
cat ${GUPPY_OUT_TMP}/reads_tmp.identity
cat ${EEL_OUT_TMP}.identity
