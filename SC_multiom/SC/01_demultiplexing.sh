#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH -o /home/mlagasag/Scripts/jobs/%x_%j.out
#SBATCH -e /home/mlagasag/Scripts/jobs/%x_%j.err
#SBATCH --job-name=demultiplexing
#SBATCH -p short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G


module load CellRanger/5.0.1
module load CellRanger-ATAC/1.0.1
module load bcl2fastq2/2.20.0-foss-2018b



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

# Get genome name
#FastqID=$1

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2

/home/mlagasag/cellranger-arc-2.0.0/bin/cellranger-arc  mkfastq --id=MM_HR_907052_RE\
                     --run=/datos/intercambio/mlagasag/211104_VH00461_83_AAAH2MKHV \
                     --csv=/home/mlagasag/MM_HR_907052/MM_HR_907052.csv \
                     --filter-dual-index \

#cellranger mkfastq --id=Fanconi_2006\
#                     --run=/datos/intercambio/Miren-Amaia/210816_VH00461_61_AAAG3FCHV \
#                     --csv=/home/mlagasag/fanconi_2006.csv \
#                     --filter-single-index \
                 
                     
printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'
