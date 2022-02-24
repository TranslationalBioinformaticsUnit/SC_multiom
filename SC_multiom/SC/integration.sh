#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH -o /home/mlagasag/Scripts/jobs/%x_%j.out
#SBATCH -e /home/mlagasag/Scripts/jobs/%x_%j.err
#SBATCH --job-name=integration
#SBATCH -p short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=90G

module load R/4.0.5-foss-2020b


#module load R/3.5.1-foss-2018b
#Code_Dir=$1


#Work Directory
#Work_Dir=$2

Rscript /home/mlagasag/Scripts/SC/harmony.R

#R CMD BATCH /home/mlagasag/Scripts/SC/integration.R


