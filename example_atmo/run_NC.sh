#!/bin/bash
# run_NC.sh -- Launch SLURM array jobs to generate atmospheric neutrino NC
# events with GiBUU for numu and numubar on both O16 and H1 targets.
#
# NC jobs use process_ID=±3 and share the combined "neutrino"/"antineutrino"
# flux files (numu+nue summed), hence only two flavors are looped here.
# Array indices 1-512 are used (complementing the CC run which uses 513-1024).
#
# Prerequisites: same as run.sh.  Flux files "neutrino" and "antineutrino"
# (the combined-flavor files) must exist in the working directory.

tools_base=/sps/juno/yqiyu/TKI_analysis/build
gibuu_base=/sps/juno/yqiyu/GiBUU2025

build_jobcard=$tools_base/build_jobcard
gibuu2root=$tools_base/gibuu2root
gibuu=$gibuu_base/release/objects/GiBUU.x
gibuu_input=$gibuu_base/buuinput


flv2name() {
    case $1 in
      14) echo "neutrino" ;;
      12) echo "neutrino" ;;
      -14) echo "antineutrino" ;;
      -12) echo "antineutrino" ;;
    esac
}

flv2process_ID() {
    case $1 in
      14) echo "3" ;;
      12) echo "3" ;;
      -14) echo "-3" ;;
      -12) echo "-3" ;;
    esac
}

flv2flavor_ID() {
    case $1 in
      14) echo "2" ;;
      12) echo "1" ;;
      -14) echo "2" ;;
      -12) echo "1" ;;
    esac
}

flv2T() {
    case $1 in
      14) echo "0" ;;
      12) echo "0" ;;
      -14) echo "1" ;;
      -12) echo "1" ;;
    esac
}
for target in O H; do
for flavor in 14 -14; do
    name=$(flv2name $flavor)
    process_ID=$(flv2process_ID $flavor)
    flavor_ID=$(flv2flavor_ID $flavor)
    T=$(flv2T $flavor)
    flux=$(readlink -f $name)

    echo "Processing flavor: $name"
    echo "Process ID: $process_ID"
    echo "Flavor ID: $flavor_ID"

    sedcommand="s/@PROCESS@/$process_ID/g;s/@FLAVOR@/$flavor_ID/g;s/@T@/$T/g;s|@FLUXFILE@|$flux|g"

    template=$(readlink -f run_${target}_nc.job.template)

    (
      mkdir -p ${name}_${target}_run
      cd ${name}_${target}_run
      sed "$sedcommand" $template > run.job
      jobcard=$(readlink -f run.job)
      mkdir result
      result_dir=$(readlink -f result)
      cat > run.sh <<EOF
#!/bin/bash
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH --mem=800M
#SBATCH --array=1-512
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH -L sps
#SBATCH -A juno

source /sps/juno/yqiyu/GiBUUGEN/env.sh
cd \$TMPDIR
( $gibuu < $jobcard |& gzip > gibuu.log.gz ) && touch done || touch fail
( $gibuu2root FinalEvents.dat record.root |& gzip > record.log.gz )  && touch record-done || touch record-fail
cp record.root $result_dir/\$SLURM_ARRAY_TASK_ID.root
EOF
      chmod +x run.sh
      sbatch run.sh
    )

done
done
