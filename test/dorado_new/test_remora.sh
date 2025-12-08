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

# script to execute buttery-eel with remora and guppy on a test dataset and compare the results

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
test -z $PATH_TO_MODKIT && PATH_TO_MODKIT=/install/modbam2bed/modbam2bed
test -z $PATH_TO_FAST5 && PATH_TO_FAST5=/data/slow5-testdata/hg2_prom_lsk114_chr22/pod5/
test -z $PATH_TO_BLOW5 && PATH_TO_BLOW5=/data/slow5-testdata/hg2_prom_lsk114_chr22/PGXX22394_reads_chr22.blow5
test -z $BIS && BIS=/data/slow5-testdata/hg2_prom_lsk114_chr22/chr22_bi.tsv
test -z $PATH_TO_IDENTITY && PATH_TO_IDENTITY=/install/biorand/bin/identitydna.sh
test -z $PATH_TO_EEL_VENV && PATH_TO_EEL_VENV=./venv3/bin/activate
test -z $MODEL && MODEL=dna_r10.4.1_e8.2_400bps_hac@v5.2.0
test -z $MOD_MODEL && MOD_MODEL=dna_r10.4.1_e8.2_400bps_hac@v5.2.0_5mCG_5hmCG@v2
test -z $REFIDX && REFIDX=/genome/hg38noAlt.idx
test -z $REF && REF=/genome/hg38noAlt.fa
test -z $GUPPY_OUT_TMP && GUPPY_OUT_TMP=ont-guppy-tmp
test -z $EEL_OUT_TMP && EEL_OUT_TMP=buttery_eel_tmp
test -z $SAMTOOLS && SAMTOOLS=samtools
test -z $MINIMAP2 && MINIMAP2=minimap2
test -z $COMPARE_METH && COMPARE_METH=/install/scripts/compare_methylation.py

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
${PATH_TO_GUPPY}/dorado_basecall_server  --model ${MODEL} --modbase_models ${MOD_MODEL} --port 5000 --use_tcp -x cuda:all --log_path ${LOGPATH} &
pid=$!
echo "Running client"
${PATH_TO_GUPPY}/ont_basecall_client --model ${MODEL} --modbase_models ${MOD_MODEL} -i ${PATH_TO_FAST5} -s ${GUPPY_OUT_TMP}  --recursive --port 5000 --use_tcp ${OPTS_GUPPY}
kill $pid
${SAMTOOLS} cat ${GUPPY_OUT_TMP}/pass/*.bam ${GUPPY_OUT_TMP}/fail/*.bam -o ${GUPPY_OUT_TMP}/reads_tmp.bam
${PATH_TO_IDENTITY} ${REFIDX} ${GUPPY_OUT_TMP}/reads_tmp.bam | cut -f 2- >  ${GUPPY_OUT_TMP}/reads_tmp.identity

#meth align
${SAMTOOLS} fastq -TMM,ML ${GUPPY_OUT_TMP}/reads_tmp.bam  |  ${MINIMAP2} -x map-ont -a -t 32 -y --secondary=no ${REFIDX} - | samtools sort  - > ${GUPPY_OUT_TMP}/reads_tmp_sorted.bam || die "remora mapping failed"
${SAMTOOLS} index ${GUPPY_OUT_TMP}/reads_tmp_sorted.bam || die "samtools index failed"
${PATH_TO_MODKIT} --cpg -m 5mC -t 32 ${REF} ${GUPPY_OUT_TMP}/reads_tmp_sorted.bam -r chr22 |  grep -v nan > ${GUPPY_OUT_TMP}/remora.bedmethyl || die "modbam2bed failed"
python3 ${COMPARE_METH} ${BIS} ${GUPPY_OUT_TMP}/remora.bedmethyl > ${GUPPY_OUT_TMP}/remora.tsv || die "compare failed"
cat ${GUPPY_OUT_TMP}/remora.tsv | tail -n+2 | cut -f3,5 | datamash ppearson 1:2 > ${GUPPY_OUT_TMP}/remora.compare || die "compare failed"

echo "Running buttery-eel"
PORT=$(get_port)
/usr/bin/time -v buttery-eel  -g ${PATH_TO_GUPPY}  --model ${MODEL} --modbase_models ${MOD_MODEL} --device 'cuda:all' -i  ${PATH_TO_BLOW5} -o  ${EEL_OUT_TMP}/reads.sam --port ${PORT} --call_mods --use_tcp ${OPTS_EEL} &> eel.log
cat eel.log
MEM=$(grep "Maximum resident set size" eel.log | cut -d " " -f 6)
if [ $MEM -gt 8000000 ]; then
    die "Memory usage is too high: $MEM"
else
    echo "Memory usage is OK: $MEM"
fi
${PATH_TO_IDENTITY} ${REFIDX} ${EEL_OUT_TMP}/reads.sam  | cut -f 2-> ${EEL_OUT_TMP}/reads.identity

#dupli check
DUPLI=$(${SAMTOOLS} view ${EEL_OUT_TMP}/reads.sam | sort | uniq -c | sort -nr -k1,1 | head -1 | awk '{print $1}')
test -z $DUPLI && die "Error in extracting reads ids "
test $DUPLI -gt 1 && die "Duplicate reads found"

#meth align
${SAMTOOLS} fastq -TMM,ML ${EEL_OUT_TMP}/reads.sam   |  ${MINIMAP2} -x map-ont -a -t 32 -y --secondary=no ${REFIDX} - | samtools sort  - > ${EEL_OUT_TMP}/reads_tmp_sorted.bam || die "remora mapping failed"
${SAMTOOLS} index ${EEL_OUT_TMP}/reads_tmp_sorted.bam || die "samtools index failed"
${PATH_TO_MODKIT} --cpg -m 5mC -t 32 ${REF} ${EEL_OUT_TMP}/reads_tmp_sorted.bam -r chr22 |  grep -v nan > ${EEL_OUT_TMP}/remora.bedmethyl || die "modbam2bed failed"
python3 ${COMPARE_METH} ${BIS} ${EEL_OUT_TMP}/remora.bedmethyl > ${EEL_OUT_TMP}/remora.tsv || die "compare failed"
cat ${EEL_OUT_TMP}/remora.tsv | tail -n+2 | cut -f3,5 | datamash ppearson 1:2 > ${EEL_OUT_TMP}/remora.compare || die "compare failed"


echo "Comparing Identity results"
diff ${GUPPY_OUT_TMP}/reads_tmp.identity ${EEL_OUT_TMP}/reads.identity || die "Results differ"

echo "Identity Test passed"
cat ${GUPPY_OUT_TMP}/reads_tmp.identity
cat ${EEL_OUT_TMP}/reads.identity

echo "Comparing methylation results"
diff ${GUPPY_OUT_TMP}/remora.compare ${EEL_OUT_TMP}/remora.compare || die "Results differ"
cat ${GUPPY_OUT_TMP}/remora.compare
cat ${EEL_OUT_TMP}/remora.compare

