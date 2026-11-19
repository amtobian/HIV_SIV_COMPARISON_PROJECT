#!/bin/bash
#SBATCH --job-name=ALIGN_REGIONS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1000
#SBATCH --output=output/ALIGN_REGIONS_%j.out
#SBATCH --error=output/ALIGN_REGIONS_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=amtobian@wm.edu

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$PROJECT_ROOT"

module load mafft/gcc-11.4.1/7.525

if ! command -v mafft >/dev/null 2>&1; then
    echo "ERROR: MAFFT is unavailable after loading its module." >&2
    exit 1
fi

threads="${SLURM_CPUS_PER_TASK:-1}"

for region in env pol gag; do
    input="results/${region}.fasta"
    output="results/${region}.aln.fasta"

    if [[ ! -s "$input" ]]; then
        echo "ERROR: Expected nonempty input file $input." >&2
        exit 1
    fi

    if [[ -s "$output" ]]; then
        echo "Skipping existing alignment $output."
        continue
    fi
    if [[ -e "$output" ]]; then
        echo "ERROR: $output exists but is empty; inspect or remove it before rerunning." >&2
        exit 1
    fi

    temporary="${output}.tmp.${SLURM_JOB_ID:-$$}"
    echo "Aligning $region with $threads thread(s)..."
    if ! mafft --auto --thread "$threads" --inputorder "$input" > "$temporary"; then
        rm -f "$temporary"
        echo "ERROR: MAFFT failed for $region." >&2
        exit 1
    fi
    if [[ ! -s "$temporary" ]]; then
        rm -f "$temporary"
        echo "ERROR: MAFFT produced an empty alignment for $region." >&2
        exit 1
    fi
    if [[ -e "$output" ]]; then
        rm -f "$temporary"
        echo "ERROR: $output appeared while this job was running; refusing to overwrite it." >&2
        exit 1
    fi
    mv "$temporary" "$output"
done

echo "Alignment files are in $PROJECT_ROOT/results."
