#!/usr/bin/R

## General paths
path_fastqc <- "/home/progs/FastQC/" # the FastQC executable folder
path_dnapi <- "/home/progs/DNApi/" # the DNApi executable folder
path_umi_tools <- "" # the UMI Tools executable folder
path_cutadapt <- "" # the Cutadapt executable folder
path_mirdeep <- "/home/progs/mirdeep2/bin/" # the Mirdeep2 executable folder
path_aligner_executable <- "/home/progs/mirdeep2/essentials/bowtie-1.1.1/bowtie" # the aligner executable folder, utilised the DNApi adapter inference through mapping.
# the path to mirdeep-2/essentials folder needs to be on the path variable

#######################################################################

## General info
exp_name = '' # Optional experiment name to be utilized in the naming of directories and files.
threads = 2 # The thread number for the parallelizable processes.

#######################################################################

## Analysis-specific variables and paths.

# Pre-processing (Quality Trimming and Adapter Removal)
preprocess = TRUE # Flag to perform the Pre-process step or not.
adapter = "RAW_INPUT" # The adapter to be removed. Empty quotes if adapter is unknown or "RAW_INPUT" if none. 
adapter_type = "a" # Provide "a" for 3' adapters or "g" for 5' adapters.
trim_quality = 10 # The minimum quality of a base allowed.
trim_max_err_rate = 0.1 # The allowed mismatch error rate between adapter and read sequence during a match.
trim_min_adapter_overlap = 3 # Minimum overlap between adapter and read sequence required to trim the sequence.
trim_min_len = 10 # Minimum read length allowed after trimming.
trim_max_len = 50 # Maximum read length allowed after trimming.
adapter_library = "/home/analysis/known_adapter_library.fa" # The path to the known adapter library fasta file.
adapter_identity_threshold = 90 # The min identity percetange of an adapter with a known adapter in the above library.
adapter_kmer_size = 10 # The kmer size an adapter to be split during the pre-processing loop process to identify and clean adapter fragments.

# Mapping
max_multimaps = 5 # Max number of allowed multimaps for reads.
path_ref_genome_prefix = "/home/analysis/Genomes/hsa_GRCh38/Index/hsa_GRCh38" # The reference genome directory.
large_index = TRUE # Bowtie large-index flag option, for more info refer to the bowtie manual online.

# Quantification
max_mismatch_to_precursors = 1 # Max allowed mismatched bases for the alignment of reads to known miRNA precursors.
species = "hsa" # The 3-letter abreviation of the species studied, used by miRBase.
path_to_hairpin = "/home/analysis/miRBase_miRNAs/hairpin_v22.fa" # The miRBase hairpin fasta file.
path_to_mature = "/home/analysis/miRBase_miRNAs/mature_v22.fa" # The miRBase mature fasta file.
