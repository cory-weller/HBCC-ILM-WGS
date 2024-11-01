#!/usr/bin/env bash
#SBATCH --gres lscratch:500
#SBATCH --mem 16G
#SBATCH --time 0:30:00
#SBATCH --partition quick
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8


PROJECT_DIR='/vf/users/CARDPB/data/HBCC/ILM_WGS_data'
BAM_DIR="${PROJECT_DIR}/BAM_FILES"
OUTPUT_DIR="${PROJECT_DIR}/REMADE_FASTQ_FROM_BAM"


# Define current job's file based on job array number
N=${SLURM_ARRAY_TASK_ID}
BAM_ARRAY=($(ls ${BAM_DIR}))
INFILE_BASENAME=${BAM_ARRAY["${N}"]}
FILESTEM=$(echo ${INFILE_BASENAME%.final.bam})
INFILE_PATH="${BAM_DIR}/${INFILE_BASENAME}"


# Work in lscratch
TMPDIR="/lscratch/${SLURM_JOB_ID}"
cd ${TMPDIR}

# Set up bam conversion table
cat << EOF > convert.txt
CMC-HBCC-ACC-DNA-5591 82017
CMC-HBCC-ACC-DNA-5785 82029
CMC-HBCC-DNA-ACC-4025 81953
CMC-HBCC-DNA-ACC-4042 81926-BAM
CMC-HBCC-DNA-ACC-4051 82054
CMC-HBCC-DNA-ACC-4057 82037
CMC-HBCC-DNA-ACC-4080 82036
CMC-HBCC-DNA-ACC-4127 82013
CMC-HBCC-DNA-ACC-4131 82042
CMC-HBCC-DNA-ACC-4135 82057
CMC-HBCC-DNA-ACC-4176 82035
CMC-HBCC-DNA-ACC-4178 82048
CMC-HBCC-DNA-ACC-4190 81944
CMC-HBCC-DNA-ACC-4218 81992
CMC-HBCC-DNA-ACC-4220 82015
CMC-HBCC-DNA-ACC-4231 82033
CMC-HBCC-DNA-ACC-4243 81941
CMC-HBCC-DNA-ACC-4249 82004
CMC-HBCC-DNA-ACC-4255 82074
CMC-HBCC-DNA-ACC-4257 81927
CMC-HBCC-DNA-ACC-4261 81991
CMC-HBCC-DNA-ACC-4275 82018
CMC-HBCC-DNA-ACC-4284 81943
CMC-HBCC-DNA-ACC-4290 82030
CMC-HBCC-DNA-ACC-4312 82063
CMC-HBCC-DNA-ACC-4314 81952
CMC-HBCC-DNA-ACC-5563 82012
CMC-HBCC-DNA-ACC-5569 81985
CMC-HBCC-DNA-ACC-5573 81925
CMC-HBCC-DNA-ACC-5578 82011
CMC-HBCC-DNA-ACC-5600 82016
CMC-HBCC-DNA-ACC-5602 81951
CMC-HBCC-DNA-ACC-5611 82038
CMC-HBCC-DNA-ACC-5643 81940
CMC-HBCC-DNA-ACC-5649 81990
CMC-HBCC-DNA-ACC-5655 82003
CMC-HBCC-DNA-ACC-5657 81983
CMC-HBCC-DNA-ACC-5669 82008
CMC-HBCC-DNA-ACC-5677 81947
CMC-HBCC-DNA-ACC-5682 81938
CMC-HBCC-DNA-ACC-5689 81959
CMC-HBCC-DNA-ACC-5754 82032
CMC-HBCC-DNA-ACC-5761 81977
CMC-HBCC-DNA-ACC-5765 81929
CMC-HBCC-DNA-ACC-5769 81930
CMC-HBCC-DNA-ACC-5771 82025
CMC-HBCC-DNA-ACC-5789 82026
CMC-HBCC-DNA-ACC-5791 81928
CMC-HBCC-DNA-ACC-5793 82027
CMC-HBCC-DNA-ACC-5797 81936
CMC-HBCC-DNA-ACC-5799 81982
CMC-HBCC-DNA-ACC-5805 82023
CMC-HBCC-DNA-ACC-5809 82082
CMC-HBCC-DNA-ACC-6001 82051
CMC-HBCC-DNA-ACC-6007 82049
CMC-HBCC-DNA-ACC-6009 82060
CMC-HBCC-DNA-ACC-6011 81997
CMC-HBCC-DNA-ACC-6013 82014
CMC-HBCC-DNA-ACC-6017 82061
CMC-HBCC-DNA-ACC-6031 82052
CMC-HBCC-DNA-ACC-6033 82039
CMC-HBCC-DNA-ACC-6035 82069
CMC-HBCC-DNA-ACC-6041 82059
CMC-HBCC-DNA-ACC-6058 82066
CMC-HBCC-DNA-ACC-6756 81974
CMC-HBCC-DNA-ACC-6758 82020
CMC-HBCC-DNA-ACC-6760 82050
CMC-HBCC-DNA-ACC-6767 81980
CMC-HBCC-DNA-ACC-6768 81973
CMC-HBCC-DNA-ACC-6769 82041
CMC-HBCC-DNA-ACC-6774 81968
CMC-HBCC-DNA-ACC-6776 81963
CMC-HBCC-DNA-ACC-6779 82040
CMC-HBCC-DNA-ACC-6781 82047
CMC-HBCC-DNA-ACC-6783 81984
CMC-HBCC-DNA-ACC-6785 82009
CMC-HBCC-DNA-ACC-6802 81948
CMC-HBCC-DNA-ACC-6809 81958
CMC-HBCC-DNA-ACC-6811 81969
CMC-HBCC-DNA-ACC-6812 81965
CMC-HBCC-DNA-ACC-6814 81949
CMC-HBCC-DNA-ACC-6818 81964
CMC-HBCC-DNA-ACC-6872 82021
HBCC-DNA-CER-13122 81932
HBCC-DNA-CER-13126 82075
HBCC-DNA-CER-13147 82081
HBCC-DNA-CER-13162 82071
HBCC-DNA-CER-13171 82010
HBCC-DNA-CER-13177 82080
HBCC-DNA-CER-13191 82070
HBCC-DNA-CER-13197 82079
HBCC-DNA-CER-13206 82131
HBCC-DNA-CER-13212 81935
HBCC-DNA-CER-13218 82005
HBCC-DNA-CER-13224 81934
HBCC-DNA-CER-13251 81987
HBCC-DNA-CER-13269 82078
HBCC-DNA-CER-13272 82086
HBCC-DNA-CER-13287 81962
HBCC-DNA-CER-13296 81961
HBCC-DNA-PFC-13349 81954
HBCC-DNA-PFC-13351 81931
HBCC-DNA-PFC-13353 81945
HBCC-DNA-PFC-13357 82045
EOF

NEW_ID=$(awk -v s=${FILESTEM} '$1 == s {print $2}' convert.txt)

# Convert to uncompressed BAM
module load samtools/1.17

echo "Converting to uncompressed bam"
samtools view --threads 8 -u ${INFILE_PATH} > ${NEW_ID}.bam

echo "Extracting reads from bam -> fastq.gz"
samtools fastq \
    --threads 8 \
    -1 ${NEW_ID}_R1_001.fastq.gz \
    -2 ${NEW_ID}_R2_001.fastq.gz \
    ${NEW_ID}.bam

# Copy to output dir
mkdir -p ${OUTPUT_DIR}
cp ${NEW_ID}_R1_001.fastq.gz ${OUTPUT_DIR}
cp ${NEW_ID}_R2_001.fastq.gz ${OUTPUT_DIR}

cd 