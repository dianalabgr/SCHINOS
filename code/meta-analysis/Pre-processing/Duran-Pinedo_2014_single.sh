#!/usr/bin/env bash

IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/RawData/
OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/RawData/QC/
RNA_sample_names=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/RNA_single.txt
THREADS=40

echo "######################################################################################              ATROPOS TRIMMING                   #################################"

OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/atropos_trimmed/
mkdir $OUT
for NAME in `awk '{print $1}' $RNA_sample_names`
do
echo "#### Trimming with ATROPOS the RNA sample: " $NAME
mkdir $OUT$NAME

atropos -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a CTGTCTCTTATACACATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  -A CTGTCTCTTATACACATCT \
        --report-file $OUT$NAME"/"$NAME"_atropos_log.txt" \
	    --minimum-length 30 \
	    --threads $THREADS \
        -o $OUT$NAME"/"$NAME"_trimmed.fasta.gz" \
        -se $IN$NAME*fasta*
done
date


echo "######################################################    Filtering out the human host reads   #####################################################"


echo "################################################################################################## SortmeRNA removal of rRNA reads from RNA data ##########"

#sortmeRNA user Manual:
#https://bioinfo.lifl.fr/RNA/sortmerna/code/SortMeRNA-user-manual-v2.1.pdf
# #2 commands: 'sortmeRNA' (main command) and 'indexdb_rna' (to index databases)

# the RNA database fasta files should be indexed first. We use command 'indexdb_rna'.
echo "************##  Indexing of the databases for smrna ##"
indexdb_rna  \
--ref /mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-bac-16s-id90.fasta,\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-bac-16s-db:\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-bac-23s-id98.fasta,\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-bac-23s-db:\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-euk-18s-id95.fasta,\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-euk-18s-db:\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-euk-28s-id98.fasta,\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-euk-28s:\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/rfam-5s-database-id98.fasta,\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/rfam-5s-db \
-m 12000 \
-v

echo "##  indexing done ##"

IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/atropos_trimmed/
OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/sortmeRNA/
RNA_DBs=/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/

mkdir $OUT

for NAME in `awk '{print $1}' $RNA_sample_names`
do
mkdir $OUT$NAME
echo "### sortmeRNA Sample" $NAME "Start###"
# first we unzip as the scripts accepts only fastq files
echo "##  unzipping input : start ##"
gunzip $IN$NAME"/"$NAME'_trimmed.fasta.gz'

echo "##  smRNA : start ##"

sortmerna  \
--ref /mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-bac-16s-id90.fasta,\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-bac-16s-db:\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-bac-23s-id98.fasta,\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-bac-23s-db:\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-euk-18s-id95.fasta,\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-euk-18s-db:\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-euk-28s-id98.fasta,\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/silva-euk-28s:\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/rfam-5s-database-id98.fasta,\
/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/rfam-5s-db \
--num_alignments 1 \
--sam \
--SQ \
--reads $IN$NAME"/"$NAME'_trimmed.fasta' \
--other $OUT$NAME"/"$NAME"_mRNA_tr" \
--aligned $OUT$NAME"/"$NAME"_rRNA_tr" \
--log \
-a $THREADS \
-m 12000 \
-v \
--fastx


samtools view -bS $OUT$NAME"/"$NAME"_rRNA_tr.sam" -@ $THREADS > $OUT$NAME"/"$NAME"_rRNA_tr.bam"

rm $OUT$NAME"/"$NAME"_rRNA_tr.sam"


echo "##  gzipping ##"

#ls $IN$NAME"/"$NAME'_trimmed.fasta' $OUT$NAME"/"$NAME"_mRNA_tr.fasta" $OUT$NAME"/"$NAME"_rRNA_tr.fasta" | parallel -t gzip
done

date


echo "############################### Align SortmeRNA output reads to Human genome by Hisat2 to remove host related reads  ##############################"

IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/sortmeRNA/
OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/host_filtered/
REF=/mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/GRCh38.primary_assembly.genome.fa

#echo "*******Building index files for human genome GRCh38 with hisat2********"
#mkdir /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/index_hisat2
#mkdir /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/index_hisat2/GRCh38
#hisat2-build $REF /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/index_hisat2/GRCh38

for NAME in `awk '{print $1}' $RNA_sample_names`
do

mkdir $OUT$NAME
echo "#### Analyzing " $NAME

ls $IN$NAME"/"$NAME*_R1_*mRNA*fastq.gz $IN$NAME"/"$NAME*_R2_*mRNA*fastq.gz
echo "###Mapping Sample" $NAME" Start###"

hisat2 -p $THREADS --dta -x /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/index_hisat2/GRCh38 -f -U $IN$NAME"/"$NAME'_mRNA_tr.fasta.gz' -S $OUT$NAME"/"$NAME".sam" --un-gz $OUT$NAME"/"$NAME"_OUT.fasta.gz" --al-gz $OUT$NAME"/"$NAME"_IN.fasta.gz" --summary-file $OUT$NAME"/"$NAME"_hisat2_report.txt"

samtools view -bS $OUT$NAME"/"$NAME'.sam' -@ $THREADS > $OUT$NAME"/"$NAME'.bam'
rm $OUT$NAME"/"$NAME'.sam'

done

date

echo "################################################################################# Running Agamemnon for abundances ########################################"


IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/host_filtered/
OUT_RNA=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/For_Agamemnon/RNA/
#move all the pretreated files to one file and change names
mkdir /mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/For_Agamemnon/
mkdir $OUT_RNA

for NAME in `awk '{print $1}' $RNA_sample_names`
do
mv $IN$NAME"/"$NAME"_OUT.fasta.gz"  $OUT_RNA$NAME'.fasta.gz'

done




