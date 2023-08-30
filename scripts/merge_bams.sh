#!/usr/bin/env bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 20G
#SBATCH --time 24:00:00

IDS=($@)
echo "iterating over ids ${IDS[@]}"
export OUTDIR='OUTPUT/MERGED_BAMS'
mkdir -p ${OUTDIR}
module purge
module load samtools/1.17

mergebam() {
    local id=${1}
    local in_bams=$(ls OUTPUT/UNMERGED_BAMS/${id}_L00*.bam | tr '\n' ' ')
    echo "id is ${id}"
    echo "CMD: samtools merge --write-index -@ 8  ${OUTDIR}/${id}.merged.bam ${in_bams}"
    samtools merge --write-index -@ 8  ${OUTDIR}/${id}.merged.bam ${in_bams}
    echo "done with ${id}"
}

export -f mergebam


parallel -j 1 mergebam ::: ${IDS[@]}
