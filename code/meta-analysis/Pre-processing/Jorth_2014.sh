#!/usr/bin/env bash

echo  "######################################################## Quality analysis of RAW DATA ###################################"
IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Jorth_2014/RawData/RNA/
OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Jorth_2014/RawData/RNA/QC/
RNA_sample_names=/mnt/raid1/armen/Schinos/META-ANALYSIS/Jorth_2014/SRR_Acc_List.txt
THREADS=32

echo "******############# Quality analysis by FastQC ##########******"
mkdir $OUT
for NAME in `awk '{print $1}' $RNA_sample_names`
do

echo "#### FastQC of RAW DATA for sample :" $NAME
mkdir $OUT$NAME
fastqc -o $OUT$NAME -t $THREADS $IN$NAME"/"*fastq*
done

# checking the quality score results together
cd $OUT
multiqc .

echo "#### Quality analysis done"

date



echo "######################################              ATROPOS TRIMMING                   #################################"

OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Jorth_2014/atropos_trimmed/
mkdir $OUT


# adapters: https://emea.support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html

for NAME in `awk '{print $1}' $RNA_sample_names`
do
echo "#### Trimming with ATROPOS the RNA sample: " $NAME
mkdir $OUT$NAME

atropos -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a CTGTCTCTTATACACATCT \
        --report-file $OUT$NAME"/"$NAME"_atropos_log.txt" \
	    --minimum-length 50 -e 0.1 \
	    --quality-base 33 \
        --quality-cutoff 20,20 \
        --max-n 0 \
	    --threads $THREADS \
        -o $OUT$NAME"/"$NAME"_trimmed.fastq.gz" \
        -se $IN$NAME"/"$NAME*fastq*

done

echo "#### Trimming of the RNA samples is done"

date

echo "#########################################          Quality analysis of TRIMMED DATA by FastQC         ######################"

OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Jorth_2014/atropos_trimmed/
mkdir $OUT"QC"

for NAME in `awk '{print $1}' $RNA_sample_names`
do

echo "#### FastQC for sample :" $NAME
mkdir $OUT"QC/"$NAME
fastqc -o $OUT"QC/"$NAME -t $THREADS $OUT$NAME"/"*fastq*

done

#checking the quality score results together with multiqc
cd $OUT
multiqc .

echo "#### Quality analysis done"




echo "##############################################        Filtering out the human host reads     ######################################################"

echo "################################################################################################## SortmeRNA removal of rRNA reads from RNA data ##########"

#sortmeRNA user Manual:
#https://bioinfo.lifl.fr/RNA/sortmerna/code/SortMeRNA-user-manual-v2.1.pdf

#2 commands: 'sortmeRNA' (main command) and 'indexdb_rna' (to index databases)

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

IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Jorth_2014/atropos_trimmed/
OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Jorth_2014/sortmeRNA/
RNA_DBs=/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/

mkdir $OUT

for NAME in `awk '{print $1}' $RNA_sample_names`
do

mkdir $OUT$NAME

echo "### sortmeRNA Sample" $NAME "Start###"


# first we unzip as the scripts accepts only fastq files
echo "##  unzipping input : start ##"
gunzip $IN$NAME"/"$NAME'_trimmed.fastq.gz'

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
--reads $IN$NAME"/"$NAME'_trimmed.fastq' \
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

ls $IN$NAME"/"$NAME'_trimmed.fastq' $OUT$NAME"/"$NAME"_mRNA_tr.fastq" $OUT$NAME"/"$NAME"_rRNA_tr.fastq" | parallel -t gzip


done

date

echo "################################################################################# Align SortmeRNA output reads to Human genome by Hisat2 to remove host related reads########################################"

IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Jorth_2014/sortmeRNA/
OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Jorth_2014/host_filtered/
REF=/mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/GRCh38.primary_assembly.genome.fa

echo "*******Building index files for human genome GRCh38 with hisat2********"
mkdir /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/index_hisat2
mkdir /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/index_hisat2/GRCh38
hisat2-build $REF /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/index_hisat2/GRCh38
mkdir $OUT

for NAME in `awk '{print $1}' $RNA_sample_names`
do

mkdir $OUT$NAME
echo "#### Analyzing " $NAME

echo "###Mapping with HISAT2 the Sample" $NAME" Start###"

hisat2 -p $THREADS --dta -x /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/index_hisat2/GRCh38 -U $IN$NAME"/"$NAME"_mRNA_tr.fastq.gz" -S $OUT$NAME"/"$NAME".sam" --un-gz $OUT$NAME"/"$NAME"_OUT.fastq.gz" --al-gz $OUT$NAME"/"$NAME"_IN.fastq.gz" --summary-file $OUT$NAME"/"$NAME"_hisat2_report.txt"

samtools view -bS $OUT$NAME"/"$NAME'.sam' -@ $THREADS > $OUT$NAME"/"$NAME'.bam'
rm $OUT$NAME"/"$NAME.sam

gzip $OUT$NAME"/"$NAME*fastq


done

date


echo "################################################################################# Running Agamemnon for abundances ########################################"


IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Jorth_2014/host_filtered/
OUT_RNA=/mnt/raid1/armen/Schinos/META-ANALYSIS/Jorth_2014/For_Agamemnon/RNA/
#move all the pretreated files to one file and change names
mkdir /mnt/raid1/armen/Schinos/META-ANALYSIS/Jorth_2014/For_Agamemnon/
mkdir $OUT_RNA

for NAME in `awk '{print $1}' $RNA_sample_names`
do

mv $IN$NAME"/"$NAME'_OUT.fastq.gz'  $OUT_RNA$NAME'.fastq.gz'

done
