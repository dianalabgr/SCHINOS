#!/bin/bash
#conda activate py3.11
THREADS=8
rsem_reference="/mnt/raid1b/philip/Homo_sapiens/rsem/Homo_sapiens"
genomeDir="/mnt/raid1b/philip/Homo_sapiens/STAR_index"
IN="/mnt/raid1b/philip/schinos_seq/data/mrna/batch3/processed_fastq_files"
sample_names="/mnt/raid1b/philip/schinos_seq/data/mrna/batch3/samples_list.txt"
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
        --sample )
			shift
			echo $1 > ./$1.txt
			sample_names=$(pwd)"/$1.txt"
			created_list=$(pwd)"/$1.txt"
			shift
			;;
    esac
done
echo $sample_names
samples=`cat $sample_names|cut -f1`
! type "foo" > /dev/null 2>&1;
cd $IN
mkdir ../mapped
mkdir ../counts
for sample in $samples; do
	echo $sample
	echo "$sample , running STAR"
#    STAR --genomeDir $genomeDir --readFilesIn $sample"_R1_processed.fastq.gz" $sample"_R2_processed.fastq.gz" \
#		--outFileNamePrefix ../mapped/$sample ulimit -n 10000 --runThreadN $THREADS \
#		--readFilesCommand zcat --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate
    STAR --genomeDir $genomeDir --readFilesIn $sample"_R1_processed.fastq.gz" $sample"_R2_processed.fastq.gz" \
		--outFileNamePrefix ../mapped_to_genes/$sample ulimit -n 10000 --runThreadN $THREADS \
		--readFilesCommand zcat --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate

#	rsem-calculate-expression -p $THREADS --paired-end --bam ../mapped/$sample"Aligned.toTranscriptome.out.bam" $rsem_reference ../counts/$sample
    rm -f $created_list
#     echo "$sample , running hisat2"
#     hisat2 -x "/mnt/raid1/philip/Homo_sapiens/grch38_tran/genome_tran" -1 $sample"_R1_filt.fastq.gz" -2 $sample"_R2_filt.fastq.gz" \
#     -S "../mapped/hisat2/"$sample"_hisat2.sam" -p $THREADS \
#     --summary-file  "../mapped/hisat2/"$sample"_hisat2_summary.txt"
done
