#!/bin/bash

#$ -q normal.q
#$ -N hisat_RNA_gen
#$ -M florentin.constancias@cirad.fr
#$ -pe parallel_smp 10
#$ -l mem_free=6G
#$ -V
#$ -cwd
#$ -V

# JOB BEGIN
module purge
module load bioinfo/hisat2/2.1.0 system/parallel/20150822 bioinfo/samtools/1.3 bioinfo/bedtools/2.24.0  # latest

##https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5032908/pdf/nihms816842.pdf
##issue with whist and human rRNA ? 
##https://github.com/infphilo/hisat2/issues/101

IN=03_RNA_smrna/
##REF=/homedir/constancias/db/refgenomes_bowtie2/human/300718/GRCh38_latest_tran
REF=/homedir/constancias/db/refgenomes_bowtie2/human/gencodegenes_08082018/GRCh38.primary_assembly.genome_tran
OUT='04_RNA_hisat2_gencodegenes/'
mkdir $OUT


for NAME in `awk '{print $1}' RNA_samples_names.txt`
do

echo "#### Analyzing " $NAME 

    ls $IN$NAME*_R1_*mRNA*fastq* $IN$NAME*_R2_*mRNA*fastq*
    echo "###Mapping Sample" $NAME" Start###"

hisat2 -p $NSLOTS --dta -x $REF -1 $IN$NAME*_R1_*mRNA*fastq* -2 $IN$NAME*_R2_*mRNA*fastq* \
-S $OUT"/"$NAME".sam" --un-conc-gz $OUT"/"$NAME"_OUT.fastq.gz" --al-conc-gz $OUT"/"$NAME"_IN.fastq.gz" \
--summary-file $OUT$NAME"_hisat2_report.txt"

samtools view -bS $OUT$NAME.sam -@ $NSLOTS > $OUT$NAME'.bam'
rm $OUT$NAME.sam

samtools view -b -f 12 -F 256 $OUT$NAME.bam -@ $NSLOTS > $OUT$NAME'_bothEndsUnmapped.bam'

samtools sort -@ $NSLOTS -n $OUT$NAME'_bothEndsUnmapped.bam' -o $OUT$NAME'_bothEndsUnmapped_sorted.bam'
rm $OUT$NAME'_bothEndsUnmapped.bam'
bamToFastq -i $OUT$NAME'_bothEndsUnmapped_sorted.bam' \
-fq $OUT$NAME'_tr_hostout_R1_.fastq' -fq2 $OUT$NAME'_tr_hostout_R2_.fastq'

ls $OUT$NAME*fastq | parallel -t gzip

rm $OUT$NAME'_bothEndsUnmapped_sorted.bam'; date
echo "#### Analyzing done for " $NAME

##http://www.metagenomics.wiki/tools/samtools/number-of-reads-in-bam-file
##counting only mapped (primary aligned) reads
samtools view  -c -F 260 $OUT$NAME'.bam' > $OUT$NAME'_hisat_primlary_alignment_counts.txt'

done

# JOB END
date

exit 0

