#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 50G
#SBATCH --time 12:00:00



#splitbam.sh

FILESTEM=${1}

DATADIR='INPUT/BAM_FILES_RENAMED'
OUTDIR='INPUT/REMADE_FASTQ_FROM_BAM'




bam="${DATADIR}/${FILESTEM}.final.bam"
r1="${OUTDIR}/${FILESTEM}_R1_001.fq"
r2="${OUTDIR}/${FILESTEM}_R2_001.fq"

echo "input bam is $bam"

# modules
module purge
module load GATK/4.4
#module load samtools/1.17



# Convert from bam to FASTQ
echo "INFO: Converting bam to fastq with SamToFastq"
echo "      gatk SamToFastq I=${bam} FASTQ=${r1} SECOND_END_FASTQ=${r2}"
gatk SamToFastq I=${bam} FASTQ=${r1} SECOND_END_FASTQ=${r2}




# Split FASTQ into 4

module load seqkit
seqkit split2 \
-1 ${r1} \
-2 ${r2} \
-p 4 \
-j 8 \
-e .gz \
-O ${OUTDIR}/${FILESTEM}_split 