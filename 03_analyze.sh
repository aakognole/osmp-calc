#!/bin/bash

cwd=`pwd`
source setenv
source inputs

printf "\nAnalyzing simulations for ( $mol1 $mol2 ) systems\n"
printf "\nPress ENTER to continue OR ctrl+C to cancel\n\n>>> "
read rep

cd run_${mol1}_${mol2}

if [ ${method} == "drude" ] || [ ${method} == "both" ]; then
    for c in ${concs}; do
	cd drude_at_$c
	rsync -auz ${cwd}/openmm/calc_osmp.py ./
	if [ ! -e osmp.${mol1}_${mol2}_${c}.1.dat ]; then
	    printf "\n>>> drude: $mol1 $mol2 $c ... "
	    ${PYTHONDIR}/python calc_osmp.py $mol1 $mol2 $c
	    printf "done!\n"
	fi
	cd ..; done
    printf "\n"
fi
if [ ${method} == "c36" ] || [ ${method} == "both" ]; then
    for c in ${concs}; do
	cd c36_at_$c
        rsync -auz ${cwd}/openmm/calc_osmp.c36.py ./
        if [ ! -e osmp.${mol1}_${mol2}_${c}.1.dat ]; then
	    printf "\n>>> c36: $mol1 $mol2 $c ... "
	    ${PYTHONDIR}/python calc_osmp.c36.py $mol1 $mol2 $c
	    printf "done!\n"
	fi
	cd ..; done
    printf "\n"
fi

rsync -auz ${cwd}/openmm/plot_osmp.py ./
${PYTHONDIR}/python plot_osmp.py ${mol1} ${mol2}

exit
