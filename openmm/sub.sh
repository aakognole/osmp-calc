#!/bin/bash

source ../../setenv

${PYTHONDIR}/python openmm_drude.py <mol1> <mol2> <CONC> 
wait

exit
