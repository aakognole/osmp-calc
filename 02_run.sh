#!/bin/bash

cwd=`pwd`
if [ -e setenv ]; then source setenv; else ./setenv.sh; source setenv; fi
source inputs

printf "\nRunning simulations for ( $mol1 $mol2 ) systems\n"
printf "\nPress ENTER to continue OR ctrl+C to cancel\n\n>>> "
read -t 15 rep

cd run_${mol1}_${mol2}
if [ ${method} == "drude" ] || [ ${method} == "both" ]; then
    for c in ${concs}; do
	cd drude_at_$c
	printf ">>> drude: $mol1 $mol2 $c \n"
	./run.sh
	cp -rfp ${cwd}/openmm/calc_osmp.py ./
	${PYTHONDIR}/python calc_osmp.py $mol1 $mol2 $c
	cd ..; done
elif [ ${method} == "c36" ] || [ ${method} == "both" ]; then
    for c in ${concs}; do
	cd c36_at_$c
	printf ">>> c36: $mol1 $mol2 $c \n"
	./run.sh
        cp -rfp ${cwd}/openmm/calc_osmp.c36.py ./
	${PYTHONDIR}/python calc_osmp.c36.py $mol1 $mol2 $c
	cd ..; done
fi
