#!/bin/bash

PROCESSES="DYJETSbkg \
WWinclusive \
VBSOSWWH_C2V_4 \
VBSOSWWH_C2V_3 \
VBSWZH_C2V_3 \
VBSWZH_C2V_4 \
WWdilep"

mkdir -p hadds/

rm -f .hadd.txt

for PROCESS in ${PROCESSES}; do
    echo "hadd -f hadds/${PROCESS}.root outputs/${PROCESS}*.root > hadds/${PROCESS}.log 2>&1" >> .hadd.txt
done

xargs.sh .hadd.txt
