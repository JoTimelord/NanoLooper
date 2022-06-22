#!/bin/bash

PROCESSES="VBSOSWWH_C2V_3 \
VBSWWH_C2V_3 \
VBSWZH_C2V_3 \
VBSZZH_C2V_3 \
DYJETSbkg \
ttbar"

mkdir -p hadds/

rm -f .hadd.txt

for PROCESS in ${PROCESSES}; do
    echo "hadd -f hadds/${PROCESS}.root outputs/${PROCESS}*.root > hadds/${PROCESS}.log 2>&1" >> .hadd.txt
done

xargs.sh .hadd.txt
