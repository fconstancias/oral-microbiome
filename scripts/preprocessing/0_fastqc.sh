#!/bin/bash

#$ -q normal.q
#$ -N fastqc
#$ -M florentin.constancias@cirad.fr
#$ -pe parallel_smp 6
#$ -l mem_free=6G
#$ -V
#$ -cwd


# JOB BEGIN


module purge
module load bioinfo/FastQC/0.11.3 


IN=raw/
mkdir $IN'qc'

fastqc \
-o $IN'qc' \
-t $NSLOTS \
$IN*gz

# JOB END
date

exit 0

