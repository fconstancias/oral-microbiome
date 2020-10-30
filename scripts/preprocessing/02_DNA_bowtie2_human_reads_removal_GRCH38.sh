#!/bin/bash

#$ -q long.q
#$ -N bowtie2_DNA_GRCH38
#$ -M florentin.constancias@cirad.fr
#$ -pe parallel_smp 20
#$ -l mem_free=8G
#$ -V
#$ -cwd
#$ -V

# JOB BEGIN

module load bioinfo/bowtie2/2.3.4.1 bioinfo/samtools/1.3 bioinfo/bedtools/2.24.0 system/parallel/20150822 

IN=02_phix/
REF=/homedir/constancias/db/refgenomes_bowtie2/human/300718/GRCh38_latest_genomic_bt22341
OUT='03_DNA_GRCH38/'
mkdir $OUT


for NAME in `awk '{print $1}' DNA_samples_names.txt`
do

echo "#### Analyzing " $NAME


    ls $IN$NAME*QUALITY*_R1*fastq* $IN$NAME*QUALITY*_R2*fastq*
    echo "###Mapping Sample" $NAME" Start###"

bowtie2 -x $REF \
        -1 $IN$NAME*QUALITY*_R1*fastq* -2 $IN$NAME*QUALITY*_R2*fastq* \
        --dovetail -p $NSLOTS -S $OUT$NAME".sam" > $OUT$NAME"_bowtie_report.txt" 2>&1

samtools view -bS $OUT$NAME.sam -@ $NSLOTS > $OUT$NAME'.bam'
rm $OUT$NAME.sam

samtools view -b -f 12 -F 256 $OUT$NAME.bam -@ $NSLOTS > $OUT$NAME'_bothEndsUnmapped.bam'

samtools sort -@ $NSLOTS -n $OUT$NAME'_bothEndsUnmapped.bam' -o $OUT$NAME'_bothEndsUnmapped_sorted.bam'
rm $OUT$NAME'_bothEndsUnmapped.bam'
bamToFastq -i $OUT$NAME'_bothEndsUnmapped_sorted.bam' \
-fq $OUT$NAME'_tr_hostout_R1_.fastq' -fq2 $OUT$NAME'_tr_hostout_R2_.fastq'

ls  $OUT$NAME*fastq | parallel -t gzip

rm $OUT$NAME'_bothEndsUnmapped_sorted.bam'; date
echo "#### Analyzing done for " $NAME


done

# JOB END
date

exit 0
