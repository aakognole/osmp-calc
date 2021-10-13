if [ ! $1 ]; then
    echo "Usage: sh step1.sh <mol1> <mol2>"
elif [ ! $2 ]; then
    echo "Usage: sh step1.sh $1 <mol2>"
fi

charmm=`which charmm`
python=`which python`
#charmm=/usr/local/bin/charmm
#charmm=/home/akognole/bin/c44b2
#python=/home/akognole/modules/miniconda2/bin/python

declare -a conc=('0.125' '0.25' '0.5' '1' '1.5' '2')

for ion in $1
do
    for mol in $2
    do
	for i in 0 1 2 3 4 5
	do
	    echo $ion $mol ${conc[$i]}
	    dir="../$1_$2_files/at_${conc[$i]}"
	    mkdir -p ${dir}
	    cp -rfp write.omm.inp ${dir}/
	    cp -rfp atoms_for_plumed.py ${dir}/
	    cp -rfp template ${dir}/
	    cd ${dir}
	    $charmm -i write.omm.inp ion=${ion} mol=${mol} conc=${conc[$i]} > write.omm.out
	    $python atoms_for_plumed.py ${ion} ${mol} ${conc[$i]}
	    sed -i -e "s~<ION>~${ion}~g" -e "s~<MOL>~${mol}~g" -e "s~<CONC>~${conc[$i]}~g" template/sub.sh
	    cd ../../03_openmm/
	done
    done
done
