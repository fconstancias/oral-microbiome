#source ~/.bash_profile
conda activate humann2

cd /media/DataDrive05/Flo/Saliva/Saliva2

NSLOTS=26

INPUT=03_DNA_GRCH38/
OUT=humann2_tuned5_DNA/
mkdir $OUT

#for NAME in K7_SUB K8_SALIVA 
for NAME in `awk '{print $1}' samples_DNA_RNA_S2.txt`

do

echo "#### humann2 DNA " $NAME

cat $INPUT$NAME*_R1_*fastq* $INPUT$NAME*_R2_*fastq*  > $INPUT$NAME"_R1R2_cat_fastq.gz" #| gunzip | sed "s/ /:/g" 

humann2 \
--input $INPUT$NAME"_R1R2_cat_fastq.gz" \
--output $OUT \
--memory-use maximum \
--threads $NSLOTS \
--protein-database /media/DataDrive05/Flo/db/humann2/uniref/ \
--prescreen-threshold 0.01 \
--output-format tsv \
--taxonomic-profile metaphlan2/tuned5/${NAME}*_profile.txt

# not sure prescreen-threshold 0.0001 make sense ust add mapping noise bequce more tqxq

##humann2_tuned5_DNA/K10_SALIVA_R1R2_cat_mod_fastq_humann2_temp/K10_SALIVA_R1R2_cat_mod_fastq_bowtie2_aligned.sam
#*bt2

samtools view -bS ${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/${NAME}_R1R2_cat_fastq_bowtie2_aligned.sam -@ $NSLOTS \
 -o ${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/${NAME}_R1R2_cat_fastq_bowtie2_aligned.bam

rm ${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/${NAME}_R1R2_cat_fastq_bowtie2_aligned.sam

#humann2_tuned5_DNA/K10_SALIVA_R1R2_cat_mod_fastq_humann2_temp

ls ${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/*tsv \
${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/*unaligned.fa \
${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/*ffn \
${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/*bt2 | parallel -t gzip

##ls ${OUT}/${NAME}*_humann2_temp/*fa | parallel -t gzip

rm $INPUT$NAME"_R1R2_cat_fastq.gz"
done

echo "#### summarize #####"
date

##https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-standard-operating-procedure-v2
humann2_join_tables -s --input ${OUT} --file_name pathabundance --output ${OUT}/DNA_humann2_pathabundance.tsv
humann2_join_tables -s --input ${OUT} --file_name pathcoverage --output ${OUT}/DNA_humann2_pathcoverage.tsv
humann2_join_tables -s --input ${OUT} --file_name genefamilies --output ${OUT}/DNA_humann2_genefamilies.tsv

humann2_renorm_table --input ${OUT}DNA_humann2_pathabundance.tsv --update-snames --units relab --output ${OUT}DNA_humann2_pathabundance_relab.tsv
humann2_renorm_table --input ${OUT}DNA_humann2_genefamilies.tsv --update-snames --units relab --output ${OUT}DNA_humann2_genefamilies_relab.tsv

humann2_split_stratified_table --input ${OUT}DNA_humann2_pathabundance_relab.tsv --output ${OUT}
humann2_split_stratified_table --input ${OUT}DNA_humann2_genefamilies_relab.tsv --output ${OUT}
humann2_split_stratified_table --input ${OUT}DNA_humann2_pathcoverage.tsv --output ${OUT}

##humann2_join_tables --input $OUT/ --file_name pathabundance --output $OUT/humann2_DNA_pathabundance.tsv
##humann2_renorm_table --input $OUT/humann2_pathabundance.tsv --units relab --output $OUT/humann2_DNA_pathabundance_relab.tsv
##humann2_split_stratified_table --input $OUT/humann2_pathabundance_relab.tsv --output $OUT/


INPUT=04_RNA_hisat2_gencodegenes/

OUT=humann2_tuned5_RNA/
mkdir $OUT

#for NAME in K7_SUB K8_SALIVA 
for NAME in `awk '{print $1}' samples_DNA_RNA_S2.txt`

do

echo "#### humann2 RNA " $NAME

cat $INPUT$NAME*_R1_*fastq* $INPUT$NAME*_R2_*fastq*  > $INPUT$NAME"_R1R2_cat_fastq.gz" #| gunzip | sed "s/ /:/g" 

humann2 \
--input $INPUT$NAME"_R1R2_cat_fastq.gz" \
--output $OUT \
--memory-use maximum \
--threads $NSLOTS \
--protein-database /media/DataDrive05/Flo/db/humann2/uniref/  \
--prescreen-threshold 0.01 \
--output-format tsv \
--taxonomic-profile metaphlan2/tuned5/${NAME}*_profile.txt

# not sure prescreen-threshold 0.0001 make sense ust add mapping noise bequce more tqxq

samtools view -bS ${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/${NAME}_R1R2_cat_fastq_bowtie2_aligned.sam -@ $NSLOTS \
 -o ${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/${NAME}_R1R2_cat_fastq_bowtie2_aligned.bam

 rm ${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/${NAME}_R1R2_cat_fastq_bowtie2_aligned.sam

#humann2_tuned5_DNA/K10_SALIVA_R1R2_cat_mod_fastq_humann2_temp

ls ${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/*tsv \
${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/*unaligned.fa \
${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/*ffn \
${OUT}${NAME}_R1R2_cat_fastq_humann2_temp/*bt2 | parallel -t gzip

##ls ${OUT}/${NAME}*_humann2_temp/*fa | parallel -t gzip

rm $INPUT$NAME"_R1R2_cat_fastq.gz"
done

echo "#### summarize #####"
date
##https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-standard-operating-procedure-v2
humann2_join_tables -s --input ${OUT} --file_name pathabundance --output ${OUT}RNA_humann2_pathabundance.tsv
humann2_join_tables -s --input ${OUT} --file_name pathcoverage --output ${OUT}RNA_humann2_pathcoverage.tsv
humann2_join_tables -s --input ${OUT} --file_name genefamilies --output ${OUT}RNA_humann2_genefamilies.tsv

humann2_renorm_table --input ${OUT}RNA_humann2_pathabundance.tsv --update-snames --units relab --output ${OUT}RNA_humann2_pathabundance_relab.tsv
humann2_renorm_table --input ${OUT}RNA_humann2_genefamilies.tsv --update-snames --units relab --output ${OUT}RNA_humann2_genefamilies_relab.tsv

humann2_split_stratified_table --input ${OUT}RNA_humann2_pathabundance_relab.tsv --output ${OUT}
humann2_split_stratified_table --input ${OUT}RNA_humann2_genefamilies_relab.tsv --output ${OUT}
humann2_split_stratified_table --input ${OUT}RNA_humann2_pathcoverage.tsv --output ${OUT}

##humann2_join_tables --input $OUT/ --file_name pathabundance --output $OUT/humann2_pathabundance.tsv
##humann2_renorm_table --input $OUT/humann2_pathabundance.tsv --units relab --output $OUT/humann2_pathabundance_relab.tsv
##humann2_split_stratified_table --input $OUT/humann2_pathabundance_relab.tsv --output $OUT/

# laplace bullshit et cie

# ? On all .tsv outputs? on combined allready?

