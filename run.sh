#!/bin/bash
#SBATCH --job-name=methylation_analysis
#SBATCH --output=methylation_analysis_%j.out
#SBATCH --error=methylation_analysis_%j.err
#SBATCH --time=12:00:00  
#SBATCH --mem=50GB  
#SBATCH --cpus-per-task=10  

source activate environment_shannon 
python /home/eharpu/methylation_analysis/run.py
