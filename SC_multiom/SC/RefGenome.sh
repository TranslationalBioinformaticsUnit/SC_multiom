#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH -o /home/mlagasag/Scripts/jobs/%x_%j.out
#SBATCH -e /home/mlagasag/Scripts/jobs/%x_%j.err
#SBATCH --job-name=RefGenome
#SBATCH -p short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G


module load CellRanger/5.0.1

abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

set -e
trap 'abort' 0



cellranger mkref \
--genome=hg38 \
--fasta=/home/mlagasag/RefGenomes/Homo_Sapiens38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
--genes=/home/mlagasag/RefGenomes/Homo_Sapiens38/Homo_sapiens.GRCh38.93.filtered.gtf

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'