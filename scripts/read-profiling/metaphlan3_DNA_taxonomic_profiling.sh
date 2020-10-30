#!/bin/sh
NSLOTS=10
eval "$(conda shell.bash hook)" #https://github.com/conda/conda/issues/7980

conda --version

conda activate metaphlan3

INPUT=03_DNA_GRCH38/


for NAME in `awk '{print $1}' samples_DNA_RNA_S1-S2.txt`
do

	echo "#### Analyzing " $NAME


	ls ${INPUT}${NAME}*_R1_*fastq* ${INPUT}${NAME}*_R2_*fastq*

	echo "#### Starting Default ####"

	OUT=metaphlan3/default/
	mkdir -p ${OUT}

	R1R2s=`ls -m  ${INPUT}${NAME}*.*fastq* | sed 's/ //g' | tr -d '\n'`
	echo ${R1R2s}

	metaphlan ${R1R2s} \
		--bowtie2out ${OUT}/${NAME}_profile.bowtie2.bz2 --nproc ${NSLOTS} \
		--input_type fastq \
		--force --bt2_ps very-sensitive-local \
		--bowtie2out ${OUT}/${NAME}_profile.bowtie2.bz2 \
		--min_ab 0.000001 \
		--legacy-output \
		--min_alignment_len 100 \
		--samout ${OUT}/${NAME}_profile.bowtie2.sam.bz2 \
		--add_viruses --tax_lev a --min_cu_len 2000 --sample_id ${NAME} --nproc ${NSLOTS} \
		--read_min_len 100 -o ${OUT}/${NAME}_profiled_metagenome.txt

	echo "#### Starting tuned ####"

	OUT=metaphlan3/tuned/
	mkdir -p ${OUT}

	metaphlan ${R1R2s} \
		--bowtie2out ${OUT}/${NAME}_profile.bowtie2.bz2 --nproc ${NSLOTS} \
		--input_type fastq --force --bt2_ps very-sensitive \
		--min_ab 0.0001 --samout ${OUT}/${NAME}_profile.bowtie2.sam.bz2 \
		--legacy-output \
		--add_viruses --tax_lev a --min_cu_len 2000 --sample_id ${NAME} --nproc ${NSLOTS} \
		--read_min_len 100 -o ${OUT}/${NAME}_profiled_metagenome.txt

	OUT=metaphlan3/tuned2/
	mkdir -p ${OUT}

	metaphlan ${R1R2s} \
		--bowtie2out ${OUT}/${NAME}_profile.bowtie2.bz2 --nproc ${NSLOTS} \
		--input_type fastq --force --bt2_ps very-sensitive-local \
		--min_ab 0.000001 --samout ${OUT}/${NAME}_profile.bowtie2.sam.bz2 \
		--legacy-output --avoid_disqm \
		--min_alignment_len 100 \
		--add_viruses --tax_lev a --min_cu_len 2000  --sample_id ${NAME} --nproc ${NSLOTS} \
		--read_min_len 70 -o ${OUT}/${NAME}_profiled_metagenome.txt
done

date


conda activate humann2.8.1  
for tag in default tuned tuned2
do
humann2_join_tables -i ${tag} --file_name profiled_metagenome.txt -o ${tag}_merged_metagenome.txt
# merge_metaphlan_tables.py ${tag}/*profiled_metagenome.txt > ${tag}_merged_metagenome.txt
grep -E "(s__)|(^#SampleID)" ${tag}_merged_metagenome.txt| grep -v "t__" > species_${tag}.txt
done
