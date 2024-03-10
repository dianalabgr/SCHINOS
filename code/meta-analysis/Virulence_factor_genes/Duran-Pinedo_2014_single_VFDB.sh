#!/usr/bin/env bash
echo  "###############################################################################################################################  Mapping with bowtie2 to the VFDB ###################################"

IN=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/For_Agamemnon/
OUT=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/VirGenes/
RNA_sample_names=/mnt/raid1/armen/Schinos/META-ANALYSIS/Duran-Pinedo_2014/Single/RNA_single.txt
THREADS=8

#this is done only once to index the database
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
bowtie2 -x /mnt/raid1/armen/Schinos/META-ANALYSIS/Databases/VFDB/index/VFDB_nt -f \
        -U $IN"RNA/"$NAME".fastq.gz"\
        --al-gz $OUT$NAME"/"$NAME"_VirReads.fastq.gz"  -p $THREADS \
        -S $OUT$NAME"/"$NAME"_vir_mapping.sam" > $OUT$NAME"/"$NAME"_bowtie_report.txt" 2>&1

done
