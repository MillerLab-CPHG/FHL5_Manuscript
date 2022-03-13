#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -n 2 
#SBATCH -p standard
#SBATCH -J quant_salmon
#SBATCH -A cphg-millerlab 
#SBATCH --output=salmon_quant.out
#SBATCH --mem 8000

module load gcc star fastqc samtools bedtools

module load salmon/0.11.2
 

FASTQ_DIR=/scratch/dw2vr/RNAseq_v2/fastq

INDEX_DIR=/scratch/dw2vr/RNAseq_v2/fastq/salmon_outHomo_sapiens.GRCh38.cdna.salmon_idx

FASTQC_REPORTS_DIR=/scratch/dw2vr/RNAseq_v2/fastq_reports

BAM_DIR=/scratch/dw2vr/RNAseq_v2/bam

OUT_DIR=/scratch/dw2vr/RNAseq_v2/salmon_outdir_2

FILE=$1

salmon quant --threads 8 \
 --libType A \
 -i ${INDEX_DIR} \
 -1 ${FASTQ_DIR}/${FILE}_R1_001.fastq \
 -2 ${FASTQ_DIR}/${FILE}_R2_001.fastq \
 --output ${OUT_DIR}/${FILE} 
