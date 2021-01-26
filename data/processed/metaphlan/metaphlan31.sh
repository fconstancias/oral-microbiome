#!/bin/sh
NSLOTS=12
eval "$(conda shell.bash hook)" #https://github.com/conda/conda/issues/7980

conda --version

conda activate metaphlan31
# conda install -c bioconda bowtie2=2.3.5.1
INPUT=03_DNA_GRCH38/


for NAME in `awk '{print $1}' samples_DNA_RNA_S1-S2.txt`
do

	echo "#### Analyzing " $NAME


	ls ${INPUT}${NAME}*_R1_*fastq* ${INPUT}${NAME}*_R2_*fastq*

	echo "#### Starting Default ####"

	OUT=metaphlan31/default/
	mkdir -p ${OUT}

	R1R2s=`ls -m  ${INPUT}${NAME}*_R*_*fastq* | sed 's/ //g' | tr -d '\n'`
	echo ${R1R2s}

	metaphlan ${R1R2s} \
	--bowtie2out ${OUT}/${NAME}_profile.bowtie2.bz2 --nproc ${NSLOTS} \
	--input_type fastq \
	--force \
	--bowtie2out ${OUT}/${NAME}_profile.bowtie2.bz2 \
	--min_ab 0.000001 \
	--legacy-output \
	--samout ${OUT}/${NAME}_profile.bowtie2.sam.bz2 \
	--add_viruses --tax_lev a --min_cu_len 2000 \
	--sample_id ${NAME} --nproc ${NSLOTS} \
	--read_min_len 100 -o ${OUT}/${NAME}_profiled_metagenome.txt

	echo "#### Starting tuned ####"

	OUT=metaphlan31/tuned/
	mkdir -p ${OUT}

	metaphlan ${R1R2s} \
	--bowtie2out ${OUT}/${NAME}_profile.bowtie2.bz2 --nproc ${NSLOTS} \
	--input_type fastq \
	--force --bt2_ps sensitive-local \
	--bowtie2out ${OUT}/${NAME}_profile.bowtie2.bz2 \
	--min_ab 0.000001 \
	--legacy-output \
	--min_alignment_len 75 \
	--samout ${OUT}/${NAME}_profile.bowtie2.sam.bz2 \
	--add_viruses --tax_lev a --min_cu_len 2000 \
	--sample_id ${NAME} --nproc ${NSLOTS} \
	--read_min_len 100 -o ${OUT}/${NAME}_profiled_metagenome.txt

	OUT=metaphlan31/tuned2/
	mkdir -p ${OUT}

	metaphlan ${R1R2s} \
	--bowtie2out ${OUT}/${NAME}_profile.bowtie2.bz2 --nproc ${NSLOTS} \
	--input_type fastq --force --bt2_ps very-sensitive-local \
	--min_ab 0.000001 --samout ${OUT}/${NAME}_profile.bowtie2.sam.bz2 \
	--legacy-output \
	--min_alignment_len 75 \
	--add_viruses --tax_lev a --min_cu_len 2000 \
	--sample_id_key ${NAME} --sample_id ${NAME} --nproc ${NSLOTS} \
	--read_min_len 100 -o ${OUT}/${NAME}_profiled_metagenome.txt

	OUT=metaphlan31/tuned3/
	mkdir -p ${OUT}

	metaphlan ${R1R2s} \
	--bowtie2out ${OUT}/${NAME}_profile.bowtie2.bz2 --nproc ${NSLOTS} \
	--input_type fastq --force --bt2_ps sensitive-local \
	--min_ab 0.000001 --samout ${OUT}/${NAME}_profile.bowtie2.sam.bz2 \
	--legacy-output \
	--min_alignment_len 100 \
	--add_viruses --tax_lev a --min_cu_len 2000 \
	--sample_id_key ${NAME} --sample_id ${NAME} --nproc ${NSLOTS} \
	--read_min_len 100 -o ${OUT}/${NAME}_profiled_metagenome.txt
done

date

																		# JOB END


																		exit 0


# for toto in default tuned1 tuned2 tuned3 tuned4 tuned5; do grep -E "(s__)|(^ID)" ${toto}_merged.txt | grep -v "t__" > species_${toto}.txt; done
