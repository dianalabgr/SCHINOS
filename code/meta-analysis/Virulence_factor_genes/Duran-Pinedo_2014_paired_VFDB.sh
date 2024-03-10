#!/usr/bin/env bash

echo  "###############################################################################################################################  Mapping with bowtie2 to the VFDB ###################################"
IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/For_Agamemnon/
OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/VirGenes/
RNA_sample_names=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Paired/RNA.txt
THREADS=8

# #build index
# cd /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/VFDB
# bowtie2-build bowtie2-build VFDB_nt.fa VFDB_nt


mkdir $OUT

#Analysis for RNA samples
for NAME in `awk '{print $1}' $RNA_sample_names`

do

echo "#### mapping for RNA sample :" $NAME
mkdir $OUT$NAME
#THE INPUT FILES ARE FASTA (PUT -f)
bowtie2 -x /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/VFDB/index/VFDB_nt -f\
        -1 $IN"RNA/"$NAME"_1.fastq.gz" -2 $IN"RNA/"$NAME"_2.fastq.gz" \
        --al-conc-gz $OUT$NAME"/"$NAME"_VirPairs.fastq.gz"  -p $THREADS \
        -S $OUT$NAME"/"$NAME"_vir_mapping.sam" > $OUT$NAME"/"$NAME"_bowtie_report.txt" 2>&1

done
