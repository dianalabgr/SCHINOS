# schinos-mirna-pipeline
This pipeline allows the user to perform miRNA NGS data analysis from raw sRNA-Seq libraries to quantification and Differential Expression Analysis.

step 1:
schinos_mirna_analysis_pipeline.r --user_config='schinos_mirna_analysis_config.r' <input_fastq_file> <output_directory>


step 2:
for each pair-wise comparison:
schinos_mirna_analysis_pipeline_DE.r --difexpression=<PATH_of_condition_table> <input_directory_with> <output_directory> <number_of_cores>

where <input_directory> contains sample*_Counts.txt files (it is te output) of step 1.