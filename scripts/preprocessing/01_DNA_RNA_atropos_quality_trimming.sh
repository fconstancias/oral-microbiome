#!/bin/bash

#$ -q long.q
#$ -N trim_pipe
#$ -M florentin.constancias@cirad.fr
#$ -pe parallel_smp 4
#$ -l mem_free=8G
#$ -V
#$ -cwd


# JOB BEGIN


SAMPLE=samples_names.txt

echo "#### 01 Cutadapt/atropos #####"

module load system/python/3.4.3

IN=raw/
OUT=01_trimmed/
mkdir $OUT

for NAME in `awk '{print $1}' DNA_samples_names.txt`
do


echo "#### Analyzing " $NAME

/usr/local/bioinfo/python/3.4.3_build2/bin/atropos \
            -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
	    --report-file $OUT$NAME"atropos_log.txt" \
	    -m 150 \
	    -e 0.1 -O 1 \
            --quality-base 33 \
            -q 20,20 --max-n 0 \
	    -T $NSLOTS \
            -o $OUT$NAME"_trimmed_R1_.fastq.gz" -p $OUT$NAME"_trimmed_R2_.fastq.gz" \
            -pe1 $IN$NAME*_R1*fastq* -pe2 $IN$NAME*_R2*fastq*

done

for NAME in `awk '{print $1}' RNA_samples_names.txt`
do


echo "#### Analyzing " $NAME

/usr/local/bioinfo/python/3.4.3_build2/bin/atropos \
            -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
	    --report-file $OUT$NAME"atropos_log.txt" \
	    -m 70 \
	    -e 0.1 -O 1 \
            --quality-base 33 \
            -q 20,20 --max-n 0 \
	    -T $NSLOTS \
            -o $OUT$NAME"_trimmed_R1_.fastq.gz" -p $OUT$NAME"_trimmed_R2_.fastq.gz" \
            -pe1 $IN$NAME*_R1*fastq* -pe2 $IN$NAME*_R2*fastq*

done


module purge
module load bioinfo/FastQC/0.11.3 

mkdir $OUT'qc'

fastqc \
-o $OUT'qc' \
-t $NSLOTS \
$OUT*gz

echo "####  Cutadapt/atropos DONE #####"
echo "#### 02 bowtie PhiX screening #####"

module purge
module load bioinfo/bowtie2/2.3.4.1

IN=01_trimmed/
OUT=02_phix/
mkdir $OUT

REF=~/db/refgenomes_bowtie2/PhiX/genome


for NAME in `awk '{print $1}' cat all_samples_names.txt`
do

bowtie2 --un-conc-gz $OUT$NAME"_phix_tr.fastq.gz" \
--sensitive --dovetail -p $NSLOTS \
-x $REF -1 $IN$NAME*R1*.fastq.gz -2 $IN$NAME*R2*.fastq.gz \
-S $OUT$NAME"phix_mapping.sam" > $OUT$NAME"_bowtie_report.txt" 2>&1

rm -f $OUT$NAME"phix_mapping.sam"
rm 01_trimmed/$NAME"_trimmed_R1_.fastq.gz" 01_trimmed/$NAME"_trimmed_R2_.fastq.gz"

done

module purge
module load bioinfo/FastQC/0.11.3 

mkdir $OUT'qc'

fastqc \
-o $OUT'qc' \
-t $NSLOTS \
$OUT*gz

echo "#### 02 bowtie PhiX screening DONE #####"

# JOB END
date

exit 0


