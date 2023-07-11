#!/usr/bin/env python3

import os
import glob
import pandas


def rename_files(df, olddir, newdir, suffix):
    os.makedirs(newdir, exist_ok=True)
    for oldpattern, newpattern in zip(df.OLD_ID.tolist(), df.NEW_ID.tolist()):
        files = glob.glob(f"{olddir}/{oldpattern}{suffix}")
        for oldname in files:
            newname = oldname.replace(oldpattern, newpattern).replace(olddir, newdir)
            #print(oldname, newname)
            os.symlink(oldname, newname)

FASTQS = pandas.read_csv('fastq_samples.tsv', sep='\t', dtype=str)
BAMS = pandas.read_csv('bam_samples.tsv', sep='\t', dtype=str)


rename_files(df=FASTQS, 
             olddir='/vf/users/CARDPB/data/HBCC/ILM_WGS_data/FASTQ_FILES',
             newdir='/vf/users/CARDPB/users/wellerca/HBCC-ILM-WGS/INPUT/FASTQ_FILES_RENAMED',
             suffix='_L*.fastq.gz')

rename_files(df=BAMS, 
             olddir='/vf/users/CARDPB/data/HBCC/ILM_WGS_data/BAM_FILES',
             newdir='/vf/users/CARDPB/users/wellerca/HBCC-ILM-WGS/INPUT/BAM_FILES_RENAMED',
             suffix='.final.bam')
