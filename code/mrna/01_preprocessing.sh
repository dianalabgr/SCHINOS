#!/bin/bash
min_quality=20
min_length=80
THREADS=8

IN="/mnt/raid1b/philip/schinos_seq/data/mrna/batch3/merged_fastq_files"
OUT="/mnt/raid1b/philip/schinos_seq/data/mrna/batch3/processed_fastq_files"
sample_names=/mnt/raid1b/philip/schinos_seq/data/mrna/batch3/samples_list.txt #

while [ "$1" != "" ]; do
    case $1 in
		--cores | --threads )
			shift
			THREADS=$1
			shift
			;;
		--in | --input )
			shift
			IN=$1
			shift
			;;
		--out | --output )
			shift
			OUT=$1
			shift
			;;
		--sample_list )
			shift
			sample_names=$1
			shift
			;;
    esac
done

echo $sample_names
samples=`cat $sample_names|cut -f1`
! type "foo" > /dev/null 2>&1;
#starting_path=$(pwd)
cd $IN
mkdir $OUT
for sample in $samples; do
	echo $sample
	echo "$sample , running trim_galore"
	trim_galore -j 2 --cores $THREADS --paired --length $min_length \
	-q $min_quality $sample*_R1.fastq.gz $sample*_R2.fastq.gz --output_dir $OUT
	if [ $? -ne 0 ]
	then
	  echo "$sample $sample trim_galore did not work" >> $starting_path"/error_log.txt"
	fi
	mv $OUT/$sample*"R1_val_1.fq.gz" $OUT/$sample"_R1_processed.fastq.gz"
	mv $OUT/$sample*"R2_val_2.fq.gz" $OUT/$sample"_R2_processed.fastq.gz"
done
