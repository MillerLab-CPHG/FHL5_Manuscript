#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -n 4
#SBATCH -p standard
#SBATCH -J sra_download
#SBATCH -A cphg-millerlab
#SBATCH --output=sra_download.out
#SBATCH --mem 28000

module load sratoolkit/2.10.5

file=$1

fasterq-dump ${file}
