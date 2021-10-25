#!/bin/bash

source ../../setenv
source ../../inputs

if [ $method == "drude" ] || [ $method == "both" ]; then
    ${PYTHONDIR}/python openmm_drude.py <mol1> <mol2> <CONC> > openmm_drude.out
elif [ $method == "c36" ] || [ $method == "both" ]; then
    ${PYTHONDIR}/python openmm_c36.py <mol1> <mol2> <CONC> > openmm_c36.out
fi

wait

exit
