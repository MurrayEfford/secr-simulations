#!/bin/bash -e

#SBATCH --job-name      simtasks
#SBATCH --account	    uoo03368
#SBATCH --time          01:00:00
#SBATCH --cpus-per-task	20
#SBATCH --mem-per-cpu   512MB
#SBATCH --output        simtasks.run_%a.out # Include the task ID %a in the 
#SBATCH --error         simtasks.run_%a.err # names of the output and error files
#SBATCH --array         5

# 2024-06-06

module load R/4.3.2-foss-2023a
module load pandoc

echo "Executing R ..."
srun  Rscript onesimtask.R
echo "R finished."
