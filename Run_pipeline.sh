#!/bin/bash
#SBATCH --job-name=Run_pipeline_py
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --ntasks 1
#SBATCH --mem 110000
#SBATCH --time 7-0:0:0
#SBATCH --cpus-per-task 20

module load python
python main.py