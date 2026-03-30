#!/bin/bash

# Local execution script (non-SLURM version of run.sh)
# Usage: ./run_local.sh <run_index>
# Example: ./run_local.sh 1

if [ $# -lt 1 ]; then
    echo "Usage: $0 <run_index>"
    echo "Example: $0 1"
    exit 1
fi

RUN_INDEX=$1

# Set up paths for local execution
gibuu_base=/var/home/yan/code/GiBUU
GIBUU=$gibuu_base/release/objects/GiBUU.x

# Use local build executables (adjust paths as needed)
# Can be overridden by setting PDK_BUILD_DIR environment variable
PDK_BUILD_DIR=$(realpath ${PDK_BUILD_DIR:-../build})
GIBUU_PDK=$PDK_BUILD_DIR/gibuu-pdk
PERT_TO_ROOT=$PDK_BUILD_DIR/pert_to_root

mkdir -p ${RUN_INDEX}
cd ${RUN_INDEX}

if [ -f success ]; then
    echo "This task has already been completed successfully."
    exit 0
fi

# Generate job files from templates (replace placeholders)
sed -e "s|@RANDON_SEED@|$RANDOM|g" -e "s|@GIBUU_BASE@|$gibuu_base|g" ../run1.job.template >> run1.job
ln -sf ../run1.inp run1.inp
sed -e "s|@RANDON_SEED@|$RANDOM|g" -e "s|@GIBUU_BASE@|$gibuu_base|g" ../run2.job.template >> run2.job
ln -sf ../epi.json epi.json

mkdir -p init_state
cd init_state
ln -sf ../run1.inp run1.inp
$GIBUU < ../run1.job | gzip -c > init_state.log.gz
cd ..

$GIBUU_PDK epi.json

mkdir -p final_state
cd final_state
ln -sf ../epi.inp run2.inp
$GIBUU < ../run2.job | gzip -c > final_state.log.gz
cd ..

$PERT_TO_ROOT epi.detail.txt final_state/PertParticles_Final_mom.dat result.root \
&& touch success \
&& rm -rf final_state init_state \
|| touch fail
