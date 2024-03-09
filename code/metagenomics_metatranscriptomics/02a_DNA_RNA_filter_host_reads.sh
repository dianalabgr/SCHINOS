#!/usr/bin/env bash
#conda activate schinos
THREADS=8
USER_DIR="/mnt/raid1b/philip/"
IN=$USER_DIR"schinos_seq/data/metagenomics/batch2/processed_fastq_files/"
OUT=$USER_DIR"schinos_seq/data/metagenomics/batch2/filtered_1_fastq_files/"
sample_names=$USER_DIR"schinos_seq/data/metagenomics/batch2/samples_list.txt" #
index=$USER_DIR"Homo_sapiens/grch38_hisat2_index/genome"
NAME="G13GF0DNA_S49"
created_list=""


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
		--out )
			shift
			OUT=$1
			shift
			;;
		--sample )
			shift
			echo $1 > ./$1.txt
			sample_names=$(pwd)"/$1.txt"
			created_list=$(pwd)"/$1.txt"
			shift
			;;
		--sample_list )
			shift
			sample_names=$1
			shift
			;;
		--index )
			shift
			index=$1
			shift
			;;
    esac
done

mkdir -p $OUT

echo "########### running hisat2 ##########"

while read NAME
do
	echo  "running hisat2 for $NAME"
hisat2 -t -p $THREADS -x $index -1 $IN$NAME"_processed_R1.fq.gz" -2 $IN$NAME"_processed_R2.fq.gz" -S $OUT$NAME"_to_human.sam" --summary-file $OUT$NAME"_summary_to_human.txt" --un-conc-gz $OUT$NAME"_unmapped_%.fq.gz"
    rm -f $OUT$NAME"_to_human.sam"
    rm -f $created_list

done < $sample_names
#

