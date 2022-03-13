#!/bin/bash
#SBATCH -t 10:00:00
#SBATCH -n 4
#SBATCH -p standard
#SBATCH -J quant_salmon
#SBATCH -A cphg-millerlab
#SBATCH --output=salmon_quant.out
#SBATCH --mem 28000

module load  gcc/9.2.0 salmon/1.2.1
FASTQ_DIR=/scratch/dw2vr/HCASMC_fastq
INDEX_DIR=/scratch/dw2vr/Cosmo_Validation/GRCh38_salmon_index
OUT_DIR=/scratch/dw2vr/Cosmo_Validation/ENCODE_RNAseq/salmon_out
n=$1

#SRR748303_1.fastq

salmon quant \
 --libType A \
 -i ${INDEX_DIR} \
 -1 ${FASTQ_DIR}/${n}_1.fastq \
 -2 ${FASTQ_DIR}/${n}_2.fastq \
 --output ${OUT_DIR}/${n}

