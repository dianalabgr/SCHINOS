# the working directory path
work_dir="/mnt/raid1b/philip/schinos_seq/data/metagenomics"

# need a text file that contains a list of all sample names, here called "sample_list.txt"
list=$work_dir"/sample_list.txt"

# step 1: preprocess fastq files
awk '{print $1}' $list | parallel -j 6 /mnt/raid1b/philip/schinos_seq/code/01_DNA_preprocess.sh \
    --input $work_dir"/merged_fastq_files" --threads 8 \
    --output $work_dir"/processed_fastq_files" --sample {}

# step 2: filter host reads from processed fastq files
awk '{print $1}' $list  | parallel -j 6 /mnt/raid1b/philip/schinos_seq/code/02a_DNA_RNA_filter_host_reads.sh \
    --threads 8 --sample {} \
    --in $work_dir"/processed_fastq_files/" \
    --out $work_dir"/filtered_1_fastq_files/"

# step 3: run agamemnon
cd /mnt/raid1b/philip/software/agamemnon/
snakemake --snakefile /mnt/raid1b/philip/software/agamemnon/AGAMEMNON \
    --cores 26 --rerun-incomplete --resources mem_mb=300000

