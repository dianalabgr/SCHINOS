#!/usr/bin/env bash

#eval "$(conda shell.bash hook)" #https://github.com/conda/conda/issues/7980
#conda activate biobakery3 #The env containing MetaPhlan4 and human 3.6 compatible with metaphlan 4

# # #Download the protein(uniref_90_diamond) and chocophlan databases
# # #humann_databases --download uniref uniref90_diamond /mnt/raid1/armen/Schinos/Beltsrom_2021/humann/databases/
# # #humann_databases --download chocophlan full /mnt/raid1/armen/Schinos/Beltsrom_2021/humann/databases

NSLOTS=30

IN=/mnt/raid1b/philip/schinos_seq/data/metatranscriptomics/batch_8-2-2024/filtered_2_fastq_files/
OUT=/mnt/raid1b/philip/schinos_seq/data/metatranscriptomics/batch_8-2-2024/humann_output/
#NAME=$1
mkdir $OUT

echo "**************************************************************** RUNNING HUMANN********************************"
echo "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ analyzing RNA data @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"

for NAME in `awk '{print $1}' /mnt/raid1b/philip/schinos_seq/data/metatranscriptomics/batch_8-2-2024/samples_list_2.txt`
do
    echo "########################################################################### Humann sample: " $NAME

    cat $IN$NAME'_filtered_1.fq.gz' $IN$NAME'_filtered_2.fq.gz'  > $IN$NAME"_R1R2_fastq.gz"

    mkdir $OUT$NAME

    # before running this we should index the whole chocophlan database

    humann \
    --input $IN$NAME"_R1R2_fastq.gz" \
    -v\
    --output $OUT$NAME \
    --memory-use minimum \
    --threads $NSLOTS \
    --bypass-nucleotide-index\
    --nucleotide-database  /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/chocophlan/full_chocophlan_index/bowtie2_index\
    --output-format tsv \
    --protein-database /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/humann/databases/uniref/

    samtools view -bS $OUT$NAME"/"$NAME'_R1R2_fastq_humann_temp/'$NAME'_R1R2_fastq_bowtie2_aligned.sam' -@ $NSLOTS \
    -o $OUT$NAME"/"$NAME'_R1R2_fastq_humann_temp/'$NAME'_R1R2_fastq_bowtie2_aligned.bam'

    rm $OUT$NAME"/"$NAME'_R1R2_fastq_humann_temp/'$NAME'_R1R2_fastq_bowtie2_aligned.sam'

    ls $OUT$NAME"/"$NAME'_R1R2_fastq_humann_temp/'*tsv \
    $OUT$NAME"/"$NAME'_R1R2_fastq_humann_temp/'*unaligned.fa \
    $OUT$NAME"/"$NAME'_R1R2_fastq_humann_temp/'*ffn \
    $OUT$NAME"/"$NAME'_R1R2_fastq_humann_temp/'*bt2 | parallel -t gzip

    rm $IN$NAME"_R1R2_fastq.gz"

    echo "************* done with "$NAME "*****************"
    date
done
