#
work_dir="/mnt/raid1b/philip/schinos_seq/data/mrna"

#sequencer_output is the directory where the raw fastq files are stored
cd $work_dir"sequencer_output"
#save all fastq file names into a text file
ls *_R1* > ../samples_list.txt
#convert them into sample ids
sed -i 's/_L001_R1_001\.fastq\.gz//g' ../samples_list.txt

#step 1, preprocessing
/mnt/raid1/philip/schinos_seq/code/mrna/01_preprocessing.sh --sample_list $work_dir"/samples_list.txt" \
    --input $work_dir"/sequencer_output" \
    --output $work_dir"/processed_fastq_files"

#step 2, mapping and counting with STAR
list=$work_dir"/samples_list.txt"
cd $work_dir
awk '{print $1}' $list | parallel -j 6 /mnt/raid1b/philip/schinos_seq/code/mrna/02_mapping_and_counts.sh \
    --sample {} \
    --input $work_dir"/processed_fastq_files" --threads 8



cd /mnt/raid1b/philip/schinos_seq/data/mrna
ls samples_list_* | parallel -j 4 /mnt/raid1b/philip/schinos_seq/code/mrna/01_preprocessing.sh \
    --sample_list {} \
    --input /mnt/raid1b/philip/schinos_seq/data/mrna/merged_fastq_files \
    --output /mnt/raid1b/philip/schinos_seq/data/mrna/processed_fastq_files --threads 8
ls samples_list_* | parallel -j 4 /mnt/raid1b/philip/schinos_seq/code/mrna/02_mapping_and_counts.sh \
    --sample_list {} \
    --input /mnt/raid1b/philip/schinos_seq/data/mrna/processed_fastq_files --threads 8

cd /mnt/raid1b/philip/schinos_seq/data/mrna/batch2
ls samples_list_* | parallel -j 3 /mnt/raid1b/philip/schinos_seq/code/mrna/02_mapping_and_counts.sh \
    --sample_list {} \
    --input /mnt/raid1b/philip/schinos_seq/data/mrna/batch2/processed_fastq_files --threads 8



####2-3-2024

work_dir="/mnt/raid1/philip/schinos_seq/data/mrna/"


#step 2
#/mnt/raid1b/philip/schinos_seq/code/mrna/02_mapping_and_counts.sh --sample_list /mnt/raid1b/philip/schinos_seq/data/mrna/batch3/samples_list.txt \
#    --input $work_dir"processed_fastq_files" --threads 8
#or
list=$work_dir"samples_list_all.txt"
awk '{print $1}' $list | parallel -j 6 /mnt/raid1/philip/schinos_seq/code/mrna/02_mapping_and_counts.sh \
    --sample {} \
    --input $work_dir"processed_fastq_files" --threads 9


