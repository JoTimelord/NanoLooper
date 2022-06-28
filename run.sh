#!/bin/bash

SAMPLEDIR=/home/users/joytzphysics/

PROCESSES="VBSOSWWH_C2V_3 \
DYJETSbkg \
ttbar"

rm -f .jobs.txt
mkdir -p logfiles
mkdir -p outputs

for PROCESS in ${PROCESSES}; do
    if [[ ${PROCESS} == *"VBSOSWWH_C2V_3" ]]; then SCALE1FB=8.752681e-6; fi
    if [[ ${PROCESS} == *"VBSOSWWH_C2V_4" ]]; then SCALE1FB=9.940178296e-6; fi
    if [[ ${PROCESS} == *"VBSWWH_C2V_3" ]]; then SCALE1FB=5.580557e-6; fi
    if [[ ${PROCESS} == *"VBSWZH_C2V_3" ]]; then SCALE1FB=5.797057e-6; fi
    if [[ ${PROCESS} == *"VBSWZH_C2V_4" ]]; then SCALE1FB=1.206246359e-5; fi
    if [[ ${PROCESS} == *"VBSZZH_C2V_3" ]]; then SCALE1FB=4.663566e-6; fi
    if [[ ${PROCESS} == *"ttbar" ]]; then SCALE1FB=0.002466; fi
    if [[ ${PROCESS} == *"DYJETSbkg" ]]; then SCALE1FB=0.2036172; fi

    IFILE=0
    for INPUTFILE in $(ls ${SAMPLEDIR}/${PROCESS}); do
        IFILE=$((IFILE+1))
        echo "./runNanoLooper --input ${SAMPLEDIR}/${PROCESS}/${INPUTFILE}  --output outputs/${PROCESS}_${IFILE}.root --scale1fb ${SCALE1FB} > logfiles/${PROCESS}_${IFILE}.log 2>&1" >> .jobs.txt
    done
done

xargs.sh .jobs.txt
