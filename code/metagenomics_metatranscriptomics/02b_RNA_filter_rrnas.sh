#!/usr/bin/env bash
#conda activate schinos
THREADS=8
USER_DIR="/mnt/raid1b/philip/"
IN=$USER_DIR"schinos_seq/data/metagenomics/batch2/filtered_1_fastq_files/"
OUT=$USER_DIR"schinos_seq/data/metagenomics/batch2/filtered_2_fastq_files/"
sample_names=$USER_DIR"schinos_seq/data/metagenomics/batch2/samples_list.txt" #
#index=$USER_DIR"Homo_sapiens/grch38_hisat2_index/genome"
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
echo "###################################### SortmeRNA removal of rRNA reads from RNA data ##########"
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda config --set channel_priority strict
#conda create --name sortmerna_env
#conda activate sortmerna_env
#conda install sortmerna
conda activate sortmerna_env
for NAME in `awk '{print $1}' $sample_names`
do
  echo "################### run sortmeRNA for Sample: " $NAME " ###"

  R1input=$IN$NAME"_unmapped_1.fq.gz"
  R2input=$IN$NAME"_unmapped_2.fq.gz"
  outDNA=$OUT$NAME"_filtered"

  rm -f -r $OUT$NAME"/kvdb"
  rm -f -r $OUT$NAME"/readb"
  mkdir -p $OUT"aligned_with_sortmerna"
  sortmerna \
    --ref $USER_DIR"software/sortmerna-4.3.6/data/rRNA_databases/silva-bac-16s-id90.fasta" \
    --ref $USER_DIR"software/sortmerna-4.3.6/data/rRNA_databases/silva-bac-23s-id98.fasta" \
    --ref $USER_DIR"software/sortmerna-4.3.6/data/rRNA_databases/silva-euk-18s-id95.fasta" \
    --ref $USER_DIR"software/sortmerna-4.3.6/data/rRNA_databases/silva-euk-28s-id98.fasta" \
    --ref $USER_DIR"software/sortmerna-4.3.6/data/rRNA_databases/rfam-5s-database-id98.fasta" \
    -workdir $OUT$NAME \
	--idx-dir $USER_DIR"software/sortmerna-4.3.6/indexes" \
	-reads $R1input -reads $R2input \
	-paired_in -other $outDNA \
	-aligned $OUT"aligned_with_sortmerna/"$NAME \
	--no-best --num_alignments 1 \
	--zip-out 1 \
    --out2 --threads $THREADS -v -fastx
  rm -f $OUT"aligned_with_sortmerna/"$NAME*
  rm -f -r $OUT$NAME
  rm -f $created_list

done
