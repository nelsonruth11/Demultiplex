#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=part3_unzipped.py
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --time=20:00:00

conda activate bgmp_py39

/usr/bin/time -v ./part_3.py -f1 "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz" \
-f2 "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz" \
-f3 "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz" \
-f4 "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz" \
-i "/projects/bgmp/shared/2017_sequencing/indexes.txt"