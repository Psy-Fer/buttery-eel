#!/bin/bash

set -e

rm -f *.log
rm -rf logs_guppy logs_dorado
test/guppy/test_extensive.sh 
mkdir logs_guppy 
mv *.log logs_guppy
test/dorado/test_install.sh
mv *.log logs_dorado





