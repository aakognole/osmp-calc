#!/bin/bash

cwd=`pwd`
if [ -e setenv ]; then source setenv; else ./setenv.sh; fi

update_inputs () {
    printf "\n# Run calculations using: drude OR c36 OR both \n" > inputs
    #printf "\nRun calculations using: drude OR c36 OR both \n>>> "
    #read method; echo "method=${method}" >> inputs
    method=drude; echo "method=${method}" >> inputs
    printf "\n# Enter residue name for first molecule \n" >> inputs
    printf "\nEnter residue name for first molecule \n>>> "
    read mol1; echo "mol1=${mol1}" >> inputs
    printf "\n# Enter residue name for second molecule \n" >> inputs
    printf "\nEnter residue name for second molecule \n>>> "
    read mol2; echo "mol2=${mol2}" >> inputs
    printf "\n# Enter ratio e.g. for MgCl2 enter >>> 1 2 \n" >> inputs
    printf "\nEnter ratio e.g. for MgCl2 enter >>> 1 2 \n>>> "
    read ratio; echo "ratio=\"${ratio}\"" >> inputs
    printf "\n# Enter concentrations (in M) you want separated by space e.g. 0.125 0.25 0.5 1 \n" >> inputs
    printf "\nEnter concentrations (in M) you want separated by space e.g. 0.125 0.25 0.5 1 \n>>> "
    read concs; echo "concs=\"${concs}\"" >> inputs
}

if [ -e inputs ]; then
    printf "\nFound inputs as:\n-----------------------------------\n"
    cat inputs
    printf "\n-----------------------------------\n"
    printf "Press ENTER to continue OR press any other key to update\n>>> "
    read rep
    if [ ! $rep ]; then
	source inputs
    else
	update_inputs
    fi
else
    update_inputs
fi

dir="${cwd}/run_${mol1}_${mol2}"; mkdir -p ${dir}
declare -a nmols; x=0; for r in $ratio;do nmols[$x]=$r; x=$((x+1)); done
declare -a conc
x=0; for c in $concs;do
         conc[$x]=$c
         if [ $method == "drude" ] || [ $method == "both" ]; then mkdir -p ${dir}/drude_at_${c}; fi
         if [ $method == "c36" ] || [ $method == "both" ]; then mkdir -p ${dir}/c36_at_${c}; fi
         x=$((x+1)); done
nconc=$((x-1))
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

#cd charmm
## calculate molecular volumes
#${CHARMMDIR}/charmm -i charmm_volume.inp resi=$mol1 > charmm.out
#v_mol1=`cat fort.99`
#${CHARMMDIR}/charmm -i charmm_volume.inp resi=$mol2 >> charmm.out
#v_mol2=`cat fort.99`
#v_wat=25.2
#cd ..

echo "Calculating compositions ..."
waters=`echo "1" | awk '{printf "%d", 6*0.0001*(48*48*96)*55}'`
declare -a nmol1; declare -a nmol2; declare -a nwat
for i in $(seq 0 $nconc);do
    nmol=`echo ${conc[$i]} | awk '{printf "%d", 6*0.0001*(48*48*96)*$1}'`
    nmol1[$i]=`echo ${nmol} ${nmols[0]} | awk '{printf "%d", $1*$2}'`
    nmol2[$i]=`echo ${nmol} ${nmols[1]} | awk '{printf "%d", $1*$2}'`
    #nwat[$i]=$waters
    nwat[$i]=`echo ${waters} ${nmol1[$i]} ${nmol2[$i]} | awk '{printf "%d",$1-$2-$3}'`
    #nwat[$i]=`echo ${nmol1[$i]} ${nmol2[$i]} $waters ${v_mol1} ${v_mol2} ${v_wat} | awk '{printf "%d", $3-((($4*$1)+($5*$2))/($6))}'`
    echo ">>> For ${conc[$i]} M solution: ${nmol1[$i]} ${mol1}, ${nmol2[$i]} ${mol2} and ${nwat[$i]} waters."
done
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

echo "Building systems for ( $mol1 $mol2 ) using Packmol ..."
cd packmol
for i in $(seq 0 $nconc);do
    printf ">>> ${conc[$i]} M ..."
    if [ $method == "drude" ] || [ $method == "both" ]; then
	printf " drude..."
	outpdb="${dir}/${mol1}_${mol2}_${conc[$i]}.pdb"
	outcrd="${dir}/${mol1}_${mol2}_${conc[$i]}.crd"
	sed -e "s~OPDB~${outpdb}~g" -e "s~OCRD~${outcrd}~g" -e "s~III~${mol1}~g" \
	    -e "s~MMM~${mol2}~g" -e "s~N_I~${nmol1[$i]}~g" -e "s~N_M~${nmol2[$i]}~g" \
	    -e "s~N_W~${nwat[$i]}~g" packmol.tmpl > ${dir}/packmol.${mol1}.${mol2}.${conc[$i]}.inp
	${packmol} < ${dir}/packmol.${mol1}.${mol2}.${conc[$i]}.inp > ${dir}/packmol.${mol1}.${mol2}.${conc[$i]}.out
	${convpdb} -renumber 1 -crdext ${outpdb} > ${outcrd}
    elif [ $method == "drude" ] || [ $method == "both" ]; then
	printf " c36..."
	outpdb="${dir}/${mol1}_${mol2}_${conc[$i]}.c36.pdb"
	outcrd="${dir}/${mol1}_${mol2}_${conc[$i]}.c36.crd"
	sed -e "s~OPDB~${outpdb}~g" -e "s~OCRD~${outcrd}~g" -e "s~III~${mol1}~g" \
	    -e "s~MMM~${mol2}~g" -e "s~N_I~${nmol1[$i]}~g" -e "s~N_M~${nmol2[$i]}~g" \
	    -e "s~N_W~${nwat[$i]}~g" packmol.c36.tmpl > ${dir}/packmol.${mol1}.${mol2}.${conc[$i]}.c36.inp
#	${packmol} < ${dir}/packmol.${mol1}.${mol2}.${conc[$i]}.c36.inp > ${dir}/packmol.${mol1}.${mol2}.${conc[$i]}.c36.out 
#	${convpdb} -renumber 1 -crdext ${outpdb} > ${outcrd}
	#sed -i -e "s~${mol1^^} A~${mol1^^}  ~g" -e "s~${mol1^^}A~${mol1^^} ~g" -e "s~${mol2^^} B~${mol2^^}  ~g" -e "s~${mol2^^}B~${mol2^^}~g" -e "s~SWM4C~SWM4 ~g" -e "s~BULK~SOLV~g" ${outpdb}
    fi
    printf " done!\n"
done
cd ..
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

echo "Running minimization and equilibration for ( $mol1 $mol2 ) using charmm ..."
echo "(Estimated time: 2 hrs per system on charmm-serial)"
cd charmm
for i in $(seq 0 $nconc);do
    printf ">>> ${conc[$i]} M ..."
    if [ $method == "drude" ] || [ $method == "both" ]; then
	printf " drude..."
	${CHARMMDIR}/charmm -i charmm_psf.inp mol1=${mol1} mol2=${mol2} nmol1=${nmol1[$i]} nmol2=${nmol2[$i]} nwater=${nwat[$i]} conc=${conc[$i]} >> ${dir}/charmm.out
	${CHARMMDIR}/charmm -i charmm_mini.inp mol1=${mol1} mol2=${mol2} conc=${conc[$i]} > ${dir}/charmm_mini.${mol1}_${mol2}_${conc[$i]}.out
	${CHARMMDIR}/charmm -i write.omm.inp mol1=${mol1} mol2=${mol2} conc=${conc[$i]} > ${dir}/drude_at_${conc[$i]}/write.omm.out
    elif [ $method == "c36" ] || [ $method == "both" ]; then
	printf " c36..."
#	${CHARMMDIR}/charmm -i charmm_psf.c36.inp mol1=${mol1} mol2=${mol2} nmol1=${nmol1[$i]} nmol2=${nmol2[$i]} nwater=${nwat[$i]} conc=${conc[$i]} >> ${dir}/charmm.c36.out
#	${CHARMMDIR}/charmm -i charmm_mini.c36.inp mol1=${mol1} mol2=${mol2} conc=${conc[$i]} > ${dir}/charmm_mini.${mol1}_${mol2}_${conc[$i]}.c36.out
#       ${CHARMMDIR}/charmm -i write.omm.c36.inp mol1=${mol1} mol2=${mol2} conc=${conc[$i]} > ${dir}/c36_at_${conc[$i]}/write.omm.out
    fi
    printf " done!\n"
done
cd ..
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

echo "Setting up openmm simulations ..."
for i in $(seq 0 $nconc);do
    printf ">>> ${conc[$i]} M ..."
    if [ $method == "drude" ] || [ $method == "both" ]; then
        printf " drude..."
	dir2="${cwd}/run_${mol1}_${mol2}/drude_at_${conc[$i]}"; mkdir -p ${dir2}
	cp -rfp openmm/* ${dir2}/
	cd ${dir2}
	${PYTHONDIR}/python atoms_for_plumed.py ${mol1} ${mol2} ${conc[$i]} 2>> ${cwd}/error.out
	sed -i -e "s~<mol1>~${mol1}~g" -e "s~<mol2>~${mol2}~g" -e "s~<CONC>~${conc[$i]}~g" run.sh 2>> ${cwd}/error.out
	cd ${cwd}; printf " done!\n"
    elif [ $method == "c36" ] || [ $method == "both" ]; then
        printf " c36..."
        dir2="${cwd}/run_${mol1}_${mol2}/c36_at_${conc[$i]}"; mkdir -p ${dir2}
        cp -rfp openmm/* ${dir2}/
        cd ${dir2}
        ${PYTHONDIR}/python atoms_for_plumed.py ${mol1} ${mol2} ${conc[$i]} 2>> ${cwd}/error.out
        sed -i -e "s~<mol1>~${mol1}~g" -e "s~<mol2>~${mol2}~g" -e "s~<CONC>~${conc[$i]}~g" run.sh 2>> ${cwd}/error.out
        cd ${cwd}; printf " done!\n"
    fi
done

exit
