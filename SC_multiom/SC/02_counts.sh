#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -o /home/mlagasag/Scripts/jobs/%x_%j.out
#SBATCH -e /home/mlagasag/Scripts/jobs/%x_%j.err
#SBATCH --job-name=counts
#SBATCH -p short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=90G


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

# Get genome name

#Work_dir=$1

FastqID=$1

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2


#/home/mlagasag/cellranger-arc-2.0.0/bin/cellranger-arc count --id=${FastqID}_counts \
#                                       --reference=/home/mlagasag/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
#                                       --libraries=/home/mlagasag/libraries_MM_HR_907052_RE.csv \
                                       
cellranger count --id=${FastqID}_job \
                 --transcriptome=/home/mlagasag/hg38/ \
                 --fastqs=/home/mlagasag/MM_HR_907052_RE/outs/fastq_path/AAAH2MKHV/ \
                 --sample=${FastqID}

printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'