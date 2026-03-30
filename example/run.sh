#!/bin/bash

#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --mem=5000M
#SBATCH --array=1-2048
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH -L sps
#SBATCH -A juno

source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh

gibuu_base=/sps/juno/yqiyu/GiBUU2025
GIBUU=$gibuu_base/release/objects/GiBUU.x

mkdir ${SLURM_ARRAY_TASK_ID} -p
cd ${SLURM_ARRAY_TASK_ID}

if [ -f success ]; then
    echo "This task has already been completed successfully."
    exit 0
fi


sed -e "s|@RANDON_SEED@|$RANDOM|g" -e "s|@GIBUU_BASE@|$gibuu_base|g" ../run1.job.template >> run1.job
ln -sf ../run1.inp run1.inp
sed -e "s|@RANDON_SEED@|$RANDOM|g" -e "s|@GIBUU_BASE@|$gibuu_base|g" ../run2.job.template >> run2.job
ln -sf ../epi.json epi.json

mkdir init_state -p
cd init_state
ln -s ../run1.inp run1.inp
$GIBUU < ../run1.job | gzip -c > init_state.log.gz
cd ..
/sps/juno/yqiyu/pdk/gibuu-pdk/build/gibuu-pdk epi.json
mkdir final_state -p
cd final_state
ln -sf ../epi.inp run2.inp
$GIBUU < ../run2.job | gzip -c > final_state.log.gz
cd ..
/sps/juno/yqiyu/pdk/Pert_to_root/build/pert-to-root epi.detail.txt final_state/PertParticles_Final_mom.dat result.root \
&& touch success \
&& rm -rf final_state init_state \
|| touch fail


