#!/bin/bash

charmm=`which charmm`

if [ ! $1 ]; then
    echo "Usage: sh step1.sh <mol1> <mol2>"
    exit
elif [ ! $2 ]; then
    echo "Usage: sh step1.sh $1 <mol2>"
    exit
fi

ion=$1
mol=$2

declare -a conc=( '0.125' '0.25' '0.5' '1' '1.5' '2' )

for ion in $1
do
    for mol in $2
    do
	for i in 0 1 2 3 4 5
	do
	    echo $ion $mol ${conc[$i]}
	    $mpirun $charmm -i charmm_mini.inp ion=${ion} mol=${mol} conc=${conc[$i]} > charmm_mini.${ion}_${mol}_${conc[$i]}.out
	done
    done
done
