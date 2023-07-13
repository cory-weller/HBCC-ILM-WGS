#!/usr/bin/env bash

nrow='100000'
testdir='INPUT/TEST'
mkdir -p ${testdir}

# Tiny version of (original) fastq
infilestem='INPUT/FASTQ_FILES_RENAMED/82084_L001'
outfilestem="${testdir}/tinyfastq_L001"

zcat "${infilestem}_R1_001.fastq.gz" | \
    head -n ${nrow} | \
    gzip -c > "${outfilestem}_R1_001.fastq.gz"

zcat "${infilestem}_R2_001.fastq.gz" | \
    head -n ${nrow} | \
    gzip -c > "${outfilestem}_R2_001.fastq.gz"

# Tiny version of regenerated (from bam) fastq
infilestem='INPUT/REMADE_FASTQ_FROM_BAM/81925'
outfilestem="${testdir}/tinyBAMfastq"

zcat "${infilestem}_R1_001.fastq.gz" | \
    head -n ${nrow} | \
    gzip -c > "${outfilestem}_R1_001.fastq.gz"

zcat "${infilestem}_R2_001.fastq.gz" | \
    head -n ${nrow} | \
    gzip -c > "${outfilestem}_R2_001.fastq.gz"
