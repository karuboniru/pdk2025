#!/bin/bash
# run.sh -- Launch SLURM array jobs to generate atmospheric neutrino CC events
# with GiBUU for all four flavors (numu, nue, numubar, nuebar) on both O16 and
# H1 targets.
#
# Each flavor×target combination gets its own sub-directory (e.g. numu_O_run/)
# containing the instantiated GiBUU jobcard and a SLURM batch script.  GiBUU
# runs in $TMPDIR, converts FinalEvents.dat to ROOT via gibuu2root, and copies
# the result to a local result/ sub-directory.
#
# Prerequisites:
#   - Adjust tools_base and gibuu_base to point to your build/install paths.
#   - Flux files numu, nue, numubar, nuebar must exist in the working directory
#     (symlinks or copies of the columnar [E, flux] files shipped in example_atmo/).
#   - The SLURM directives (-A, -L) may need adapting for your cluster.

tools_base=/sps/juno/yqiyu/TKI_analysis/build
gibuu_base=/sps/juno/yqiyu/GiBUU2025

build_jobcard=$tools_base/build_jobcard
gibuu2root=$tools_base/gibuu2root
gibuu=$gibuu_base/release/objects/GiBUU.x
gibuu_input=$gibuu_base/buuinput


flv2name() {
    case $1 in
      14) echo "numu" ;;
      12) echo "nue" ;;
      -14) echo "numubar" ;;
      -12) echo "nuebar" ;;
    esac
}

flv2process_ID() {
    case $1 in
      14) echo "2" ;;
      12) echo "2" ;;
      -14) echo "-2" ;;
      -12) echo "-2" ;;
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
for flavor in 14 12 -14 -12; do
    name=$(flv2name $flavor)
    process_ID=$(flv2process_ID $flavor)
    flavor_ID=$(flv2flavor_ID $flavor)
    T=$(flv2T $flavor)
    flux=$(readlink -f $name)

    echo "Processing flavor: $name"
    echo "Process ID: $process_ID"
    echo "Flavor ID: $flavor_ID"

    sedcommand="s/@PROCESS@/$process_ID/g;s/@FLAVOR@/$flavor_ID/g;s/@T@/$T/g;s|@FLUXFILE@|$flux|g"

    template=$(readlink -f run_${target}.job.template)

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
#SBATCH --array=513-1024
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
