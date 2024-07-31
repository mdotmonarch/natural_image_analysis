#!/bin/bash
#SBATCH --job-name=natural_images_analysis
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-16
#SBATCH --output=job-natural_images_analysis_%A_%a.out
#SBATCH --error=job-natural_images_analysis_%A_%a.err
#SBATCH --mem-per-cpu=8G

julia src/main.jl $SLURM_ARRAY_TASK_ID