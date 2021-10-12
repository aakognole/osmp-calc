#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -o jobout
#$ -j y
#$ -N <ION>_<MOL>_<CONC>
#$ -l h_data=1000M,h_rt=12:00:00,gpu=1
#$ -l h=(n036|n042|n043|n044|n045|n046|n047|n048|n049|n050|n051|n052|n053|n054|n055|n056|n057|n058|n059|n060|n061|n062|n063|n064|n065|n066|n067|n068|n069|n070|n071|n072|n073|n074|n075|n076|n077|n078|n079|n080)
#$ -pe smp 2
#$ -R y

echo $HOSTNAME $JOB_ID `date`

export PYTHONPATH="/opt/openmm/latest/pythonapi:$PYTHONPATH"
export LD_LIBRARY_PATH=/opt/cuda/latest/lib64:/opt/openmm/7.4.0/plugins/plumed/plumed_in_openmm:/opt/openmm/7.4.0/plugins/plumed/lib

ion=<ION>
mol=<MOL>
conc=<CONC>

cwd=`pwd`
cd /tmp
mkdir -p akognole
cd akognole
mkdir job_${JOB_ID}
cd job_${JOB_ID}
rsync -auz ${cwd}/* ./

ln -s ${cwd}/../../../toppar_drude
ln -s ${cwd}/../../${ion}_${mol}_${conc}.xplor.psf

/opt/openmm/latest/bin/python run.py ${ion} ${mol} ${conc}

rsync -auz ./* ${cwd}/
cd ..
rm -rf job_${JOB_ID}
cd ${cwd}
