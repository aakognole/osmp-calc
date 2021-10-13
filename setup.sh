#!/bin/bash

cwd=`pwd`
cd toppar
tar -xzf toppar_drude.tgz
tar -xzf toppar_c36.tgz
cd ..

charmm=`which charmm`
printf "Enter path to charmm binary i.e. \${CHARMMDIR}/charmm (Found: ${charmm:0:-7})\n"
printf "Press ENTER to continue or specify path \n>>> "
read rep; if [ $rep ]; then charmm=${rep}; else charmm=${charmm:0:-7}; fi
export CHARMMDIR=${charmm}; echo -e "\nCHARMMDIR set to ${charmm}"
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

export convpdb=${cwd}/bin/toolset/perl/convpdb.pl
if [ ! -e $convpdb ]; then
    printf "Installing MMTSB Toolset...\n---------------------------\n"
    cd bin
    ./install_mmtsb.sh
    cd ..
    if [ ! -e $convpdb ]; then
	echo "Something went wrong! MMTSB Toolset could not be installed. Exiting now..."
	echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
	exit
    else
	printf "\n---------------------------\nsuccess"
	echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
    fi
fi

export packmol=${cwd}/bin/packmol/packmol
if [ ! -e $packmol ]; then
    printf "Installing Packmol...\n---------------------\n"
    cd bin
    ./install_packmol.sh
    cd ..
    if [ ! -e $packmol ]; then
	echo "Something went wrong! Packmol could not be installed. Exiting now..."
	echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
        exit
    else
        printf "\n---------------------\nsuccess"
	echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"
    fi
fi

printf "Run calculations using: drude OR c36 OR both \n>>> "
read method
printf "\nEnter residue name for first molecule \n>>> "
read mol1
printf "\nEnter residue name for second molecule \n>>> "
read mol2
dir="${cwd}/run_${mol1}_${mol2}"
mkdir -p ${dir}

printf "\nEnter ratio e.g. for MgCl2 enter >>> 1 2 \n>>> "
read ratio
declare -a nmols; x=0; for r in $ratio;do nmols[$x]=$r; x=$((x+1)); done
printf "\nEnter concentrations (in M) you want separated by space e.g. 0.125 0.25 0.5 1 \n>>> "
read concs
echo $concs > ${dir}/concs.dat
declare -a conc; x=0; for c in $concs;do conc[$x]=$c; x=$((x+1)); done; nconc=$((x-1))
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
    nwat[$i]=$waters
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
	${packmol} < ${dir}/packmol.${mol1}.${mol2}.${conc[$i]}.c36.inp > ${dir}/packmol.${mol1}.${mol2}.${conc[$i]}.c36.out
	${convpdb} -renumber 1 -crdext ${outpdb} > ${outcrd}
	#sed -i -e "s~${mol1^^} A~${mol1^^}  ~g" -e "s~${mol1^^}A~${mol1^^} ~g" -e "s~${mol2^^} B~${mol2^^}  ~g" -e "s~${mol2^^}B~${mol2^^}~g" -e "s~SWM4C~SWM4 ~g" -e "s~BULK~SOLV~g" ${outpdb}
    fi
    printf " done!\n"
done
cd ..
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

echo "Running minimization and equilibration for ( $mol1 $mol2 ) using charmm ..."
cd charmm
for i in $(seq 0 $nconc);do
    printf ">>> ${conc[$i]} M ..."
    if [ $method == "drude" ] || [ $method == "both" ]; then
	printf " drude..."
	${CHARMMDIR}/charmm -i charmm_psf.inp ion=${mol1} mol=${mol2} nion=${nmol1[$i]} nmol=${nmol2[$i]} nwater=${nwat[$i]} conc=${conc[$i]} >> ${dir}/charmm.out
	${CHARMMDIR}/charmm -i charmm_mini.inp ion=${mol1} mol=${mol2} conc=${conc[$i]} > ${dir}/charmm_mini.${mol1}_${mol2}_${conc[$i]}.out
    elif [ $method == "c36" ] || [ $method == "both" ]; then
	printf " c36..."
	${CHARMMDIR}/charmm -i charmm_psf.c36.inp ion=${mol1} mol=${mol2} nion=${nmol1[$i]} nmol=${nmol2[$i]} nwater=${nwat[$i]} conc=${conc[$i]} >> ${dir}/charmm.c36.out
	${CHARMMDIR}/charmm -i charmm_mini.c36.inp ion=${mol1} mol=${mol2} conc=${conc[$i]} > ${dir}/charmm_mini.${mol1}_${mol2}_${conc[$i]}.c36.out
    fi
    printf " done!\n"
done
cd ..
echo -e "\n-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-\n"

exit
