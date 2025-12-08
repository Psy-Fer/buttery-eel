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

# script to execute buttery-eel and guppy on a test dataset and compare the results

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
  }')
echo $PORT
}

CURRENT_GUPPY=$(grep "ont-pybasecall-client-lib" requirements.txt | cut -d "=" -f 3)
test -z ${CURRENT_GUPPY} && die "ont-pybasecall-client-lib not found in requirements.txt"

#defaults if not set
test -z $PATH_TO_GUPPY && PATH_TO_GUPPY=/install/ont-dorado-server-${CURRENT_GUPPY}/bin/
test -z $PATH_TO_FAST5 && PATH_TO_FAST5=/data/slow5-testdata/NA12878_prom_subsubsample/fast5/
test -z $PATH_TO_BLOW5 && PATH_TO_BLOW5=/data/slow5-testdata/NA12878_prom_subsubsample/reads.blow5
test -z $PATH_TO_IDENTITY && PATH_TO_IDENTITY=/install/biorand/bin/identitydna.sh
test -z $PATH_TO_EEL_VENV && PATH_TO_EEL_VENV=./venv3/bin/activate
test -z $MODEL && MODEL=dna_r10.4.1_e8.2_400bps_fast@v5.2.0
test -z $REFIDX && REFIDX=/genome/hg38noAlt.idx
test -z $GUPPY_OUT_TMP && GUPPY_OUT_TMP=ont-guppy-tmp
test -z $EEL_OUT_TMP && EEL_OUT_TMP=buttery_eel_tmp

#check if files exist
test -e ${PATH_TO_GUPPY}/dorado_basecall_server || die  "${PATH_TO_GUPPY}/dorado_basecall_server not foundd"
test -e ${PATH_TO_GUPPY}/ont_basecall_client || die  "${PATH_TO_GUPPY}/ont_basecall_client not found"
test -d ${PATH_TO_FAST5} || die  "${PATH_TO_FAST5} not found"
test -e ${PATH_TO_BLOW5} || die  "${PATH_TO_BLOW5} not found"
test -e ${PATH_TO_IDENTITY} || die  "${PATH_TO_IDENTITY} not found"

#clean up
test -d ${GUPPY_OUT_TMP} && rm -r ${GUPPY_OUT_TMP}
test -d ${EEL_OUT_TMP} && rm -r ${EEL_OUT_TMP}
mkdir ${EEL_OUT_TMP} || die "Failed to create ${EEL_OUT_TMP}"

#sourcing venv
source ${PATH_TO_EEL_VENV} || die "Failed to source ${PATH_TO_EEL_VENV}"

echo "Running server"
LOGPATH=$(mktemp -d)
${PATH_TO_GUPPY}/dorado_basecall_server --model ${MODEL} --port 5000 --use_tcp -x cuda:all --log_path ${LOGPATH} &
pid=$!
echo "Running client"
${PATH_TO_GUPPY}/ont_basecall_client --model ${MODEL} -i ${PATH_TO_FAST5} -s ${GUPPY_OUT_TMP} --recursive ${OPTS_GUPPY} --port 5000 --use_tcp
kill $pid

cat ${GUPPY_OUT_TMP}/pass/* ${GUPPY_OUT_TMP}/fail/* > ${GUPPY_OUT_TMP}/reads_tmp.fastq
${PATH_TO_IDENTITY} ${REFIDX} ${GUPPY_OUT_TMP}/reads_tmp.fastq | cut -f 2- >  ${GUPPY_OUT_TMP}/reads_tmp.identity

echo "Running buttery-eel"
PORT=$(get_port)
/usr/bin/time -v buttery-eel -g ${PATH_TO_GUPPY} --model ${MODEL} --device 'cuda:all' -i ${PATH_TO_BLOW5} -o ${EEL_OUT_TMP}/reads.fastq --port ${PORT} --use_tcp ${OPTS_EEL} &> eel.log
cat eel.log
MEM=$(grep "Maximum resident set size" eel.log | cut -d " " -f 6)
if [ $MEM -gt 8000000 ]; then
    die "Memory usage is too high: $MEM"
else
    echo "Memory usage is OK: $MEM"
fi
${PATH_TO_IDENTITY} ${REFIDX} ${EEL_OUT_TMP}/reads.fastq | cut -f 2-> ${EEL_OUT_TMP}/reads.identity
DUPLI=$(awk '{if(NR%4==1) {print $1}}' ${EEL_OUT_TMP}/reads.fastq  | tr -d '@' | sort | uniq -c | sort -nr -k1,1 | head -1 | awk '{print $1}')
test -z $DUPLI && die "Error in extracting reads ids"
test $DUPLI -gt 1 && die "Duplicate reads found"

echo "Comparing results"
diff ${GUPPY_OUT_TMP}/reads_tmp.identity ${EEL_OUT_TMP}/reads.identity || die "Results differ"

echo "Test passed"
cat ${GUPPY_OUT_TMP}/reads_tmp.identity
cat ${EEL_OUT_TMP}/reads.identity