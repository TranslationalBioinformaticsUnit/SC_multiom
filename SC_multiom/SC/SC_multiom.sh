#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH -o /home/mlagasag/Scripts/jobs/%x_%j.out
#SBATCH -e /home/mlagasag/Scripts/jobs/%x_%j.err
#SBATCH --job-name=multiom
#SBATCH -p short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=90G

module load R/4.0.5-foss-2020b

#Code_Dir=$1


#Work Directory
#Work_Dir=$2

Rscript /home/mlagasag/Scripts/SC/SC_multiom.R