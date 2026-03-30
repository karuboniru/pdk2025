#!/bin/bash
# Local (non-SLURM) version of run.sh.
# Runs a single GiBUU CC job for one flavor × target and converts the output
# to ROOT with gibuu2root.
#
# Usage: ./run_local.sh [target] [flavor_pdg]
#   target     : O or H          (default: O)
#   flavor_pdg : 14 12 -14 -12   (default: 14)

target=${1:-O}
flavor=${2:-14}

gibuu_base=${GIBUU_BASE:-/var/home/yan/code/GiBUU}
gibuu=$gibuu_base/release/objects/GiBUU.x

PDK_BUILD_DIR=$(realpath ${PDK_BUILD_DIR:-../build})
gibuu2root=$PDK_BUILD_DIR/gibuu2root

flv2name() {
    case $1 in
      14)  echo "numu"    ;;
      12)  echo "nue"     ;;
      -14) echo "numubar" ;;
      -12) echo "nuebar"  ;;
    esac
}

flv2process_ID() {
    case $1 in
      14|12)   echo "2"  ;;
      -14|-12) echo "-2" ;;
    esac
}

flv2flavor_ID() {
    case $1 in
      14|-14) echo "2" ;;
      12|-12) echo "1" ;;
    esac
}

flv2T() {
    case $1 in
      14|12)   echo "0" ;;
      -14|-12) echo "1" ;;
    esac
}

name=$(flv2name $flavor)
process_ID=$(flv2process_ID $flavor)
flavor_ID=$(flv2flavor_ID $flavor)
T=$(flv2T $flavor)

run_dir=${name}_${target}_local
mkdir -p $run_dir
cd $run_dir

sed \
  -e "s/@PROCESS@/$process_ID/g" \
  -e "s/@FLAVOR@/$flavor_ID/g" \
  -e "s/@T@/$T/g" \
  -e "s|@FLUXFILE@|$(realpath ../$name)|g" \
  -e "s|path_to_input='[^']*'|path_to_input='$gibuu_base/buuinput'|g" \
  ../run_${target}.job.template > run.job

$gibuu < run.job > gibuu.log 2>&1
$gibuu2root FinalEvents.dat record.root $flavor
