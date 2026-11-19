#!/bin/bash
#SBATCH --job-name=HIV_DIVERSITY_STATS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=output/DIVERSITY_STATS_%j.out
#SBATCH --error=output/DIVERSITY_STATS_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=amtobian@wm.edu

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
WORKFLOW_SCRIPT="$PROJECT_ROOT/scripts/10_statistical_diversity.sh"

if [[ ! -f "$WORKFLOW_SCRIPT" ]]; then
    echo "ERROR: Could not find $WORKFLOW_SCRIPT." >&2
    echo "Submit this job from the project root directory." >&2
    exit 1
fi

cd "$PROJECT_ROOT"
bash "$WORKFLOW_SCRIPT"
