#!/usr/bin/env bash
THREADS=20
#conda activate schinos
IN=/mnt/raid1b/philip/schinos_seq/data/metagenomics/batch2/merged_fastq_files
OUT=/mnt/raid1b/philip/schinos_seq/data/metagenomics/batch2/processed_fastq_files
sample_names=/mnt/raid1b/philip/schinos_seq/data/metagenomics/batch2/samples_list.txt #
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
    esac
done

 echo "####################################"
#
mkdir -p $OUT
#
for NAME in `awk '{print $1}' $sample_names`
do
    echo "#### Trimming with ATROPOS the sample: " $NAME
    if [ -f $IN/$NAME*"_R1"*".fastq.gz" ]; then
        atropos -a CTGTCTCTTATA -A CTGTCTCTTATA \
            --report-file $OUT/$NAME"_atropos3_log.txt" \
            --minimum-length 60 -e 0.1 \
            --quality-base 33 \
            --quality-cutoff 20,20 \
            --max-n 0 \
            --threads $THREADS \
            -o $OUT/$NAME"_processed_R1.fq.gz" \
            --paired-output $OUT/$NAME"_processed_R2.fq.gz" \
            -pe1 $IN/$NAME*"_R1"*".fastq.gz" -pe2 $IN/$NAME*"_R2"*".fastq.gz"
    else
        echo 'File does not exist in input folder !'
    fi
done
rm $created_list


# echo "##################################################################                   Quality analysis of TRIMMED DATA by FastQC         ######################"
# mkdir $OUT"fastqc"
# echo "#### FastQC for sample :" $NAME
# fastqc -o $OUT"fastqc" -t $THREADS $OUT*fastq.gz #both pairs will be analysed
#
#
# #checking the quality score results together with multiqc
# cd $OUT"fastqc"
# multiqc --outdir multiqc .
#
# echo "#### Quality analysis done"


