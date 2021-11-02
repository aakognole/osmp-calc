#!/bin/bash

source ../../setenv
source ../../inputs

method=<method>

for r in $(seq 1 3); do
    if [ $method == "drude" ]; then
	${PYTHONDIR}/python openmm_drude.py <mol1> <mol2> <CONC> ${r} > openmm_drude.${r}.out
	if [ $DISKCLEAN == 1 ] && [ -e <mol1>_<mol2>_<CONC>.${r}.nowat.dcd ]; then
	    rm <mol1>_<mol2>_<CONC>.${r}.dcd
	fi
    elif [ $method == "c36" ]; then
	${PYTHONDIR}/python openmm_c36.py <mol1> <mol2> <CONC> ${r} > openmm_c36.${r}.out
        if [ $DISKCLEAN == 1 ] && [ -e <mol1>_<mol2>_<CONC>.${r}.nowat.dcd ]; then
            rm <mol1>_<mol2>_<CONC>.${r}.dcd
        fi
    fi
done

wait

exit
