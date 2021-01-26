
#!/bin/bash

#export PATH=~/miniconda3/bin/:$PATH
#source ~/.bashrc 
conda activate biobakery3

cd /media/DataDrive05/Flo/Saliva/Saliva2

INPUT=/datadrive05/Flo/Saliva/Saliva2/03_DNA_GRCH38/
OUT=/datadrive05/Flo/Saliva/Saliva2/humann3_DNA/
mkdir $OUT

NSLOTS=20

# for NAME in `awk '{print $1}' samples_DNA_RNA_S2.txt`

for NAME in K5_TONGUE K6_SALIVA K6_SUB K6_TONGUE K7_SALIVA K7_SUB K7_TONGUE K8_SALIVA K8_SUB K8_TONGUE K9_SALIVA K9_SUB K9_TONGUE MP10_SALIVA MP10_SUB MP10_TONGUE MP11_SALIVA MP11_SUB MP11_TONGUE MP1_SALIVA MP1_TONGUE MP2_SALIVA MP2_TONGUE MP3_SALIVA MP3_SUB MP3_TONGUE MP4_SALIVA MP4_SUB MP4_TONGUE MP5_SALIVA MP5_SUB MP5_TONGUE MP6_SALIVA MP6_SUB MP6_TONGUE MP7_SALIVA MP7_SUB MP7_TONGUE MP8_SALIVA MP8_SUB MP8_TONGUE MP9_SALIVA MP9_SUB

do

    echo "#### humann2 DNA " $NAME

    humann \
    --input ${INPUT}${NAME}_DNA_cat.fastq.gz \
    --output ${OUT} \
    --memory-use maximum \
    --threads ${NSLOTS} \
    --protein-database /datadrive05/Flo/db/humann3/uniref \
    --search-mode uniref90 
# \
# --taxonomic-profile /datadrive05/Flo/Saliva/Saliva2/metaphlan31/default/${NAME}_profiled_metagenome.txt

samtools view -bS ${OUT}${NAME}_DNA_cat_humann_temp/${NAME}_DNA_cat_bowtie2_aligned.sam -@ $NSLOTS \
-o ${OUT}${NAME}_DNA_cat_humann_temp/${NAME}_DNA_cat_bowtie2_aligned.bam

rm ${OUT}${NAME}_DNA_cat_humann_temp/${NAME}_DNA_cat_bowtie2_aligned.sam

ls ${OUT}${NAME}_DNA_cat_humann_temp/*tsv \
${OUT}${NAME}_DNA_cat_humann_temp/*unaligned.fa \
${OUT}${NAME}_DNA_cat_humann_temp/*ffn \
${OUT}${NAME}_DNA_cat_humann_temp/*bt2 | parallel -t gzip

mkdir ${OUT}/genefamilies ${OUT}/genefamilies_relab ${OUT}/genefamilies_cpm ${OUT}/pathabundance ${OUT}/pathabundance_relab ${OUT}/pathcoverage ${OUT}/pathabundance_cpm

mv ${OUT}${NAME}_DNA_cat_genefamilies.tsv ${OUT}/genefamilies/${NAME}_genefamilies.tsv
mv ${OUT}${NAME}_DNA_cat_pathabundance.tsv ${OUT}/pathabundance/${NAME}_pathabundance.tsv
mv ${OUT}${NAME}_DNA_cat_pathcoverage.tsv ${OUT}/pathcoverage/${NAME}_pathcoverage.tsv

humann_renorm_table \
--input ${OUT}/genefamilies/${NAME}_genefamilies.tsv \
--output ${OUT}/genefamilies_relab/${NAME}_genefamilies.tsv \
-p \
--units relab

humann_renorm_table \
--input ${OUT}/genefamilies/${NAME}_genefamilies.tsv \
--output ${OUT}/genefamilies_cpm/${NAME}_genefamilies.tsv \
-p \
--units cpm

humann_renorm_table \
--input ${OUT}/pathabundance/${NAME}_pathabundance.tsv \
--output ${OUT}/pathabundance_relab/${NAME}_pathabundance.tsv \
-p \
--units relab

humann_renorm_table \
--input ${OUT}/pathabundance/${NAME}_pathabundance.tsv \
--output ${OUT}/pathabundance_cpm/${NAME}_pathabundance.tsv \
-p \
--units cpm

done



echo "#### summarize #####"
date

##https://github.com/LangilleLab/microbiome_helper/wiki/Metagenomics-standard-operating-procedure-v2
# humann_join_tables -s --input ${OUT} --file_name pathabundance --output ${OUT}/DNA_humann2_pathabundance.tsv
# humann_join_tables -s --input ${OUT} --file_name pathcoverage --output ${OUT}/DNA_humann2_pathcoverage.tsv
# humann_join_tables -s --input ${OUT} --file_name genefamilies --output ${OUT}/DNA_humann2_genefamilies.tsv

#humann2_renorm_table --input ${OUT}DNA_humann2_pathabundance.tsv --update-snames --units relab --output ${OUT}DNA_humann2_pathabundance_relab.tsv
#humann2_renorm_table --input ${OUT}DNA_humann2_genefamilies.tsv --update-snames --units relab --output ${OUT}DNA_humann2_genefamilies_relab.tsv

#humann2_split_stratified_table --input ${OUT}DNA_humann2_pathabundance_relab.tsv --output ${OUT}
#humann2_split_stratified_table --input ${OUT}DNA_humann2_genefamilies_relab.tsv --output ${OUT}
#humann2_split_stratified_table --input ${OUT}DNA_humann2_pathcoverage.tsv --output ${OUT}

##humann2_join_tables --input $OUT/ --file_name pathabundance --output $OUT/humann2_DNA_pathabundance.tsv
##humann2_renorm_table --input $OUT/humann2_pathabundance.tsv --units relab --output $OUT/humann2_DNA_pathabundance_relab.tsv
##humann2_split_stratified_table --input $OUT/humann2_pathabundance_relab.tsv --output $OUT/


#!/bin/bash

#export PATH=~/miniconda3/bin/:$PATH
#source ~/.bashrc 
conda activate biobakery3

cd /media/DataDrive05/Flo/Saliva/Saliva2
NSLOTS=20

INPUT=04_RNA_hisat2_gencodegenes/

OUT=humann2_RNA/
mkdir $OUT

#for NAME in K7_SUB K8_SALIVA 
for NAME in `awk '{print $1}' samples_DNA_RNA_S2.txt`

do

    echo "#### humann2 RNA " $NAME

# cat $INPUT$NAME*_R1_*fastq* $INPUT$NAME*_R2_*fastq*  > $INPUT$NAME"_R1R2_cat_fastq.gz" #| gunzip | sed "s/ /:/g" 

humann \
--input ${INPUT}${NAME}_cat.fastq.gz \
--output ${OUT} \
--memory-use maximum \
--threads ${NSLOTS} \
--protein-database /datadrive05/Flo/db/humann3/uniref \
--search-mode uniref90 \
--taxonomic-profile /datadrive05/Flo/Saliva/Saliva2/humann3_DNA/${NAME}_DNA_cat_humann_temp/${NAME}_DNA_cat_metaphlan_bugs_list.tsv

# not sure prescreen-threshold 0.0001 make sense ust add mapping noise bequce more tqxq

samtools view -bS ${OUT}${NAME}_cat_humann_temp/${NAME}_cat_bowtie2_aligned.sam -@ $NSLOTS \
-o ${OUT}${NAME}_cat_humann_temp/${NAME}_cat_bowtie2_aligned.bam

rm ${OUT}${NAME}_cat_humann_temp/${NAME}_cat_bowtie2_aligned.sam

ls ${OUT}${NAME}_cat_humann_temp/*tsv \
${OUT}${NAME}_cat_humann_temp/*unaligned.fa \
${OUT}${NAME}_cat_humann_temp/*ffn \
${OUT}${NAME}_cat_humann_temp/*bt2 | parallel -t gzip

mkdir ${OUT}/genefamilies ${OUT}/genefamilies_relab ${OUT}/genefamilies_cpm ${OUT}/pathabundance ${OUT}/pathabundance_relab ${OUT}/pathcoverage ${OUT}/pathabundance_cpm

mv ${OUT}${NAME}_cat_genefamilies.tsv ${OUT}/genefamilies/${NAME}_genefamilies.tsv
mv ${OUT}${NAME}_cat_pathabundance.tsv ${OUT}/pathabundance/${NAME}_pathabundance.tsv
mv ${OUT}${NAME}_cat_pathcoverage.tsv ${OUT}/pathcoverage/${NAME}_pathcoverage.tsv

humann_renorm_table \
--input ${OUT}/genefamilies/${NAME}_genefamilies.tsv \
--output ${OUT}/genefamilies_relab/${NAME}_genefamilies.tsv \
-p \
--units relab

humann_renorm_table \
--input ${OUT}/genefamilies/${NAME}_genefamilies.tsv \
--output ${OUT}/genefamilies_cpm/${NAME}_genefamilies.tsv \
-p \
--units cpm

humann_renorm_table \
--input ${OUT}/pathabundance/${NAME}_pathabundance.tsv \
--output ${OUT}/pathabundance_relab/${NAME}_pathabundance.tsv \
-p \
--units relab

humann_renorm_table \
--input ${OUT}/pathabundance/${NAME}_pathabundance.tsv \
--output ${OUT}/pathabundance_cpm/${NAME}_pathabundance.tsv \
-p \
--units cpm

done

###################################
# Then transfert locally.

# Then: join, regroup / rename PWAYs
merge_metaphlan_tables.py humann3_DNA/*/*bug* > humann3_DNA/merged_metaphlan3_humann3.tsv

for NAME in humann3_DNA humann3_RNA
do
    for toto in genefamilies genefamilies_relab genefamilies_cpm pathabundance pathabundance_cpm pathabundance_relab pathcoverage
    do 
humann_join_tables --input ${NAME}/${toto}  --output ${NAME}/${toto}/${toto}_joined_tables.tsv #--file_name ${toto}
done
done


#It is recommended to renorm before regrouping.
#https://groups.google.com/forum/#!topic/humann-users/tlkjsp22iPg
#If, after regrouping, you want all samples to sum to the same value again, you would run renorm_table again.
#To improve statistical power, we usually limit the number of features we test. One easy way to do this is with a variance filter: e.g. selecting the top 10% most variable features across your samples. 
#This tends to isolate features that are fairly abundant and fairly prevalent, while discounting super-conserved features (e.g. ribosomal proteins) that are unlikely to differ in abundance across samples.
#It's almost always going to be a good idea to normalize the RPK outputs using the renorm_table script (to adjust for differences in sequencing depth between samples). 
#In some cases it's useful to have the RPK numbers, which is why they are left unnormalized by default. For example, if you want to talk about the coverage of a genome in your sample, you need unnormalized RPK values.

for NAME in humann3_DNA humann3_RNA
do
    for toto in genefamilies genefamilies_relab genefamilies_cpm
    do
        for cate in uniref90_rxn uniref90_go uniref90_ko uniref90_level4ec uniref90_pfam uniref90_eggnog
        do
            humann_regroup_table -i ${NAME}/${toto}/${toto}_joined_tables.tsv -g ${cate} -f sum -e 100 -u N -o ${NAME}/${toto}/${toto}_joined_tables_${cate}.tsv
        done
    done
done

for NAME in humann3_DNA humann3_RNA
do
    for toto in genefamilies genefamilies_relab genefamilies_cpm
    do
        for cate in uniref90_rxn uniref90_go uniref90_ko uniref90_level4ec uniref90_pfam uniref90_eggnog
        do
            for catenames in kegg-orthology kegg-pathway kegg-module ec metacyc-rxn metacyc-pwy pfam eggnog go infogo1000
            do
                humann_rename_table -i ${NAME}/${toto}/${toto}_joined_tables_${cate}.tsv -n ${catenames} -s -o ${NAME}/${toto}/${toto}_joined_tables_${cate}_renamed_${catenames}.tsv
            done
        done
    done
done

# then for pway: replace biosynthesis by bs. ...

