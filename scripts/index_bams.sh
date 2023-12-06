#!/usr/bin/env bash
#SBATCH --mem 20G
#SBATCH --time 1:00:00
#SBATCH --ntasks-per-node 5
#SBATCH --nodes 1
#SBATCH --partition quick,norm
module purge

module load samtools

bam=$(sed -n ${SLURM_ARRAY_TASK_ID}p need_index.txt)

samtools index -@ 5 -c ${bam}
