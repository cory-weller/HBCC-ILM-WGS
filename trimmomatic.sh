#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks 5
#SBATCH --mem 20G
#SBATCH --time 2:59:00
#SBATCH --partition quick

module load trimmomatic

cd /vf/users/CARDPB/users/wellerca/HBCC-ILM-WGS/INPUT/FASTQ_FILES_RENAMED


java -jar ${TRIMMOJAR} PE -threads 5 \
    82053_L003_R1_001.fastq.gz \
    82053_L003_R2_001.fastq.gz \
    82053_L003_R1_paired.fastq.gz \
    82053_L003_R1_unpaired.fastq.gz \
    82053_L003_R2_paired.fastq.gz \
    82053_L003_R2_unpaired.fastq.gz \
    LEADING:3 TRAILING:3 MINLEN:36

# rename files to use only paired R1 and R2
rm 82053_L003_R1_001.fastq.gz
rm 82053_L003_R2_001.fastq.gz
rm 82053_L003_R1_unpaired.fastq.gz
rm 82053_L003_R2_unpaired.fastq.gz
mv 82053_L003_R1_paired.fastq.gz 82053_L003_R1_001.fastq.gz
mv 82053_L003_R2_paired.fastq.gz 82053_L003_R2_001.fastq.gz
