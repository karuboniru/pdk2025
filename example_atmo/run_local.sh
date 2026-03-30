#!/bin/bash
# run_local.sh -- Local (non-SLURM) version of run.sh.
# Runs a single GiBUU CC job for one flavor × target combination and converts
# the output to ROOT with gibuu2root.
#
# Usage: ./run_local.sh [target] [flavor]
#   target : O (oxygen-16) or H (hydrogen-1)  [default: O]
#   flavor : neutrino PDG ID: 14, 12, -14, -12  [default: 14]
#
# Environment variables (override defaults):
#   GIBUU          path to GiBUU.x binary
#   GIBUU_INPUT    path to GiBUU buuinput data directory
#   PDK_BUILD_DIR  path to this project's build directory (for gibuu2root)

target=${1:-O}
flavor=${2:-14}

# ---------------------------------------------------------------------------
# Path configuration -- adjust or export these before running
# ---------------------------------------------------------------------------
GIBUU=${GIBUU:-$(command -v GiBUU.x)}
GIBUU_INPUT=${GIBUU_INPUT:-/usr/share/GiBUU/buuinput}
PDK_BUILD_DIR=$(realpath ${PDK_BUILD_DIR:-../build})
GIBUU2ROOT=$PDK_BUILD_DIR/gibuu2root

if [ -z "$GIBUU" ] || [ ! -x "$GIBUU" ]; then
    echo "Error: GiBUU.x not found. Set GIBUU=/path/to/GiBUU.x"
    exit 1
fi
if [ ! -x "$GIBUU2ROOT" ]; then
    echo "Error: gibuu2root not found at $GIBUU2ROOT. Set PDK_BUILD_DIR."
    exit 1
fi

# ---------------------------------------------------------------------------
# Flavor → GiBUU parameter mappings
# ---------------------------------------------------------------------------
flv2name() {
    case $1 in
      14) echo "numu" ;;
      12) echo "nue" ;;
      -14) echo "numubar" ;;
      -12) echo "nuebar" ;;
      *) echo "unknown" ;;
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

if [ "$name" = "unknown" ]; then
    echo "Error: unknown flavor '$flavor'. Use 14, 12, -14, or -12."
    exit 1
fi

flux=$(readlink -f $name)
if [ ! -f "$flux" ]; then
    echo "Error: flux file '$name' not found in $(pwd)"
    exit 1
fi

template=$(readlink -f run_${target}.job.template)
if [ ! -f "$template" ]; then
    echo "Error: template run_${target}.job.template not found."
    exit 1
fi

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
run_dir=${name}_${target}_local
mkdir -p $run_dir
cd $run_dir

sed \
  -e "s/@PROCESS@/$process_ID/g" \
  -e "s/@FLAVOR@/$flavor_ID/g" \
  -e "s/@T@/$T/g" \
  -e "s|@FLUXFILE@|$flux|g" \
  -e "s|/sps/juno/yqiyu/GiBUU2025/buuinput|$GIBUU_INPUT|g" \
  $template > run.job

echo "Running GiBUU (target=$target, flavor=$name, process_ID=$process_ID)..."
$GIBUU < run.job > gibuu.log 2>&1 && echo "GiBUU done." || { echo "GiBUU failed; see gibuu.log"; exit 1; }

echo "Converting FinalEvents.dat → record.root ..."
$GIBUU2ROOT FinalEvents.dat record.root $flavor && echo "Done: $run_dir/record.root" || { echo "gibuu2root failed"; exit 1; }
