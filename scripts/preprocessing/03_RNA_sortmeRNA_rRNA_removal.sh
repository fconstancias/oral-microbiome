#!/bin/bash

#$ -q long.q
#$ -N sortmeRNA_para_re
#$ -M florentin.constancias@cirad.fr
#$ -pe parallel_smp 20
#$ -l mem_free=6G
#$ -V
#$ -cwd
#$ -V

# JOB BEGIN
module load system/git/2.8.3 system/parallel/20150822 bioinfo/samtools/1.3

IN=02_phix/
OUT=03_RNA_smrna_rev/
mkdir $OUT

for NAME in `awk '{print $1}' RNA_samples_names_reverse.txt`
do

echo "### sortmeRNA Sample" $NAME "Start###"

R1input=$IN$NAME'_phix_tr.fastq.1'
R2input=$IN$NAME'_phix_tr.fastq.2'

inR1R2=$OUT$NAME"_R1R2_tr.fastq"
outmRNAR1R2=$OUT$NAME"_R1R2_mRNA_tr"
outrRNAR1R2=$OUT$NAME"_R1R2_rRNA_tr"

echo "##  unzipping inputs : start ##"


ls $R1input".gz" $R2input".gz" | parallel -t gunzip

echo "##  merging R1 and R2 : start ##"


/gs7k1/home/constancias/tools/sortmernav2.1/scripts/merge-paired-reads.sh $R1input $R2input $inR1R2

echo "##  smRNA : start ##"


/gs7k1/home/constancias/tools/sortmernav2.1/bin/sortmerna  \
--ref /homedir/constancias/tools/sortmerna-2.1/rRNA_databases/silva-bac-16s-id90.fasta,\
/homedir/constancias/tools/sortmerna-2.1/rRNA_databases/silva-bac-16s-db:\
/homedir/constancias/tools/sortmerna-2.1/rRNA_databases/silva-bac-23s-id98.fasta,\
/homedir/constancias/tools/sortmerna-2.1/rRNA_databases/silva-bac-23s-db:\
/homedir/constancias/tools/sortmerna-2.1/rRNA_databases/silva-euk-18s-id95.fasta,\
/homedir/constancias/tools/sortmerna-2.1/rRNA_databases/silva-euk-18s-db:\
/homedir/constancias/tools/sortmerna-2.1/rRNA_databases/silva-euk-28s-id98.fasta,\
/homedir/constancias/tools/sortmerna-2.1/rRNA_databases/silva-euk-28s:\
/homedir/constancias/tools/sortmerna-2.1/rRNA_databases/rfam-5s-database-id98.fasta,\
/homedir/constancias/tools/sortmerna-2.1/rRNA_databases/rfam-5s-db \
--num_alignments 1 \
--sam \
--reads $inR1R2 \
--paired_in \
--other $outmRNAR1R2 \
--aligned $outrRNAR1R2 \
--log \
-a $NSLOTS \
-m 12000 \
-v \
--fastx

echo "##  unmerging R1 and R2 mRNA ##"
date

/gs7k1/home/constancias/tools/sortmernav2.1/scripts/unmerge-paired-reads.sh \
$outmRNAR1R2".fastq" $OUT$NAME"_R1_tr_mRNA.fastq" $OUT$NAME"_R2_tr_mRNA.fastq"

echo "##  unmerging R1 and R2 rRNA ##"
date

/gs7k1/home/constancias/tools/sortmernav2.1/scripts/unmerge-paired-reads.sh \
$outrRNAR1R2".fastq" $OUT$NAME"_R1_tr_rRNA.fastq" $OUT$NAME"_R2_tr_rRNA.fastq"

echo "##  removing merged input, m/r RNA and sam outputs ##"

rm $inR1R2
rm $outmRNAR1R2".fastq"
rm $outrRNAR1R2".fastq"

samtools view -bS $outrRNAR1R2".sam" -@ $NSLOTS > $outrRNAR1R2".bam"

rm $outrRNAR1R2".sam"


echo "##  gzipping rRNA and input ##"

ls $OUT$NAME*.fastq $R1input $R2input | parallel -t gzip

done

# JOB END
date

exit 0

