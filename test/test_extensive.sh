#!/bin/bash

set -e

rm -f *.log
rm -rf logs_guppy logs_dorado logs_dorado_new
test/guppy/test_extensive.sh
mkdir logs_guppy logs_dorado logs_dorado_new
mv *.log logs_guppy
test/dorado/test_extensive.sh
mv *.log logs_dorado
test/dorado_new/test_extensive.sh
mv *.log logs_dorado_new





