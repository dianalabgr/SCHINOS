#!/usr/bin/env bash

echo  "###############################################################################################################################  Quality analysis of RAW DATA ###################################"
IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/RawData/
OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/RawData/QC/
RNA_sample_names=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/RNA.txt
THREADS=32

echo "######################################################################################              ATROPOS TRIMMING                   #################################"

OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/atropos_trimmed/
mkdir $OUT

echo "#### Trimming with ATROPOS the RNA sample: " $NAME

mkdir $OUT$NAME

atropos -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a CTGTCTCTTATACACATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT  -A CTGTCTCTTATACACATCT \
        --report-file $OUT$NAME"/"$NAME"_atropos_log.txt" \
	    --minimum-length 30 \
        --threads $THREADS \
        -o $OUT$NAME"/"$NAME"_trimmed_R1_.fasta.gz" \
        --paired-output $OUT$NAME"/"$NAME"_trimmed_R2_.fasta.gz" \
        -pe1  $IN$NAME*_R1.fasta* -pe2  $IN$NAME*_R2.fasta* \

done

date


echo "################################################################################################## SortmeRNA removal of rRNA reads from RNA data ##########"

# #sortmeRNA user Manual:
# #https://bioinfo.lifl.fr/RNA/sortmerna/code/SortMeRNA-user-manual-v2.1.pdf
#
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

IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/atropos_trimmed/
OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/sortmeRNA/
RNA_DBs=/mnt/raid1/armen/sortmerna-4.3.6/data/rRNA_databases/

mkdir $OUT

for NAME in `awk '{print $1}' $RNA_sample_names`
do
mkdir $OUT$NAME
echo "### sortmeRNA Sample" $NAME "Start###"

R1input=$IN$NAME"/"$NAME'_trimmed_R1_.fasta'
R2input=$IN$NAME"/"$NAME'_trimmed_R2_.fasta'
inR1R2=$OUT$NAME"/"$NAME"_R1R2_tr.fasta"
outmRNAR1R2=$OUT$NAME"/"$NAME"_R1R2_mRNA_tr"
outrRNAR1R2=$OUT$NAME"/"$NAME"_R1R2_rRNA_tr"

# sortmeRNA accepts only paired reads in one file, thus we have to merge pairs together.
#The merge script was not available in the latest version of sortmeRNA, it was taken from older version:https://github.com/biocore/sortmerna/releases/tag/2.0
#first we unzip as the scripts accepts only fastq files
echo "##  unzipping inputs : start ##"
ls $R1input'.gz' $R2input'.gz' | parallel -t gunzip
echo "##  merging R1 and R2 : start ##"
/mnt/raid1/armen/sortmerna-4.3.6/scripts/merge-paired-reads.sh $R1input $R2input $inR1R2
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
--reads $inR1R2 \
--paired_in \
--other $outmRNAR1R2 \
--aligned $outrRNAR1R2 \
--log \
-a $THREADS \
-m 12000 \
-v \
--fastx

echo "##  unmerging R1 and R2 mRNA ##"

/mnt/raid1/armen/sortmerna-4.3.6/scripts/unmerge-paired-reads.sh \
$outmRNAR1R2".fasta" $OUT$NAME"/"$NAME"_R1_tr_mRNA.fasta" $OUT$NAME"/"$NAME"_R2_tr_mRNA.fasta"

echo "##  unmerging R1 and R2 rRNA ##"

/mnt/raid1/armen/sortmerna-4.3.6/scripts/unmerge-paired-reads.sh \
$outrRNAR1R2".fasta" $OUT$NAME"/"$NAME"_R1_tr_rRNA.fasta" $OUT$NAME"/"$NAME"_R2_tr_rRNA.fasta"

echo "##  removing merged input, m/r RNA and sam outputs ##"

rm $inR1R2
rm $outmRNAR1R2".fasta"
rm $outrRNAR1R2".fasta"

samtools view -bS $outrRNAR1R2".sam" -@ $THREADS > $outrRNAR1R2".bam"

rm $outrRNAR1R2".sam"

echo "##  gzipping rRNA and input ##"
ls $OUT$NAME"/"$NAME*.fasta $R1input $R2input | parallel -t gzip
done

date

echo "################################################################################# Align SortmeRNA output reads to Human genome by Hisat2 to remove host related reads########################################"

IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/sortmeRNA/
OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/host_filtered/
REF=/mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/GRCh38.primary_assembly.genome.fa

echo "*******Building index files for human genome GRCh38 with hisat2********"
mkdir /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/index_hisat2
mkdir /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/index_hisat2/GRCh38
hisat2-build $REF /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/index_hisat2/GRCh38

for NAME in `awk '{print $1}' $RNA_sample_names`
do

mkdir $OUT$NAME
echo "#### Analyzing " $NAME

ls $IN$NAME"/"$NAME*_R1_*mRNA*fastq.gz $IN$NAME"/"$NAME*_R2_*mRNA*fastq.gz

echo "###Mapping Sample" $NAME" Start###"

hisat2 -p $THREADS --dta -x /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/refGenomes/index_hisat2/GRCh38 -1 $IN$NAME"/"$NAME'_R1_tr_mRNA.fasta.gz' -f -2 $IN$NAME"/"$NAME'_R2_tr_mRNA.fasta.gz' -S $OUT$NAME"/"$NAME".sam" --un-conc-gz $OUT$NAME"/"$NAME"_OUT.fasta.gz" --al-conc-gz $OUT$NAME"/"$NAME"_IN.fasta.gz" --summary-file $OUT$NAME"/"$NAME"_hisat2_report.txt"

samtools view -bS $OUT$NAME"/"$NAME'.sam' -@ $THREADS > $OUT$NAME"/"$NAME'.bam'
rm $OUT$NAME"/"$NAME'.sam'

samtools view -b -f 12 -F 256 $OUT$NAME"/"$NAME'.bam' -@ $THREADS > $OUT$NAME"/"$NAME'_unmapped_pair.bam'

samtools sort -@ $THREADS -n $OUT$NAME"/"$NAME'_unmapped_pair.bam' -o $OUT$NAME"/"$NAME'_unmapped_pair_sorted.bam'

rm $OUT$NAME"/"$NAME'_unmapped_pair.bam'

samtools fasta -1 $OUT$NAME"/"$NAME'_tr_hostout_R1_.fasta'    -2  $OUT$NAME"/"$NAME'_tr_hostout_R2_.fasta' $OUT$NAME"/"$NAME'_unmapped_pair_sorted.bam'

bamToFastq -i $OUT$NAME"/"$NAME'_unmapped_pair_sorted.bam' -fq $OUT$NAME"/"$NAME'_tr_hostout_R1_.fastq' -fq2 $OUT$NAME"/"$NAME'_tr_hostout_R2_.fastq'

ls $OUT$NAME"/"$NAME*fasta | parallel -t gzip
rm $OUT$NAME"/"$NAME'_unmapped_pair_sorted.bam'; date


done

date

echo "################################################################################# Running Agamemnon for abundances ########################################"


IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/host_filtered/
OUT_DNA=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/For_Agamemnon/DNA/
OUT_RNA=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/For_Agamemnon/RNA/
#move all the pretreated files to one file and change names
mkdir /mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/For_Agamemnon/
mkdir $OUT_RNA



for NAME in `awk '{print $1}' $RNA_sample_names`
do

mv $IN$NAME"/"$NAME'_tr_hostout_R1_.fasta.gz'  $OUT_RNA$NAME'_1.fasta.gz'
mv $IN$NAME"/"$NAME'_tr_hostout_R2_.fasta.gz'  $OUT_RNA$NAME'_2.fasta.gz'

rm $IN$NAME"/"$NAME'_tr_hostout_R1_.fastq'
rm $IN$NAME"/"$NAME'_tr_hostout_R2_.fastq'

done




