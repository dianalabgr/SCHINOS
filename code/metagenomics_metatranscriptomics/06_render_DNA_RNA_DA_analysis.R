library(this.path())
code_dir=this.path::here()
project_dir=this.path::here(..=1)

#make reports
for (type in c("metagenomics", "metatranscriptomics" )){
  agamemnon_output=paste0("/mnt/raid1b/philip/schinos_seq/results/",type,"/agamemnon_output/final/all")
  report_dir=file.path(this.path::here(..=1),"results",type,"reports")
  for (i in c("SA","SB","SP","GF")){
    rmarkdown::render(file.path(code_dir,'06_DNA_RNA_DA_analysis.Rmd'),
                      output_file = paste0(report_dir,'/SCHINOS_',type,'_',i,'_report.', Sys.Date(), sep=''),
                      params = list(
                        tissue=i,
                        counts_folder=agamemnon_output,
                        sequencing=type))
  }
}
# # make abundances tables
# #metagenomics
# code_dir=this.path::here()
# for (type in c("metagenomics" )){
#   output_folder=paste0("/mnt/raid1b/philip/schinos_seq/results/",type,"/abundances")
#   agamemnon_output=paste0("/mnt/raid1b/philip/schinos_seq/results/",type,"/agamemnon_output/final/all")
#     rmarkdown::render(file.path(code_dir,'04_create_abundance_tables.Rmd'),
#                       output_file = file.path(output_folder,"junk.pdf"),
#                       params = list(
#                         counts_folder=agamemnon_output,
#                         sequencing=type,
#                         output_folder=output_folder))
# }
# 
# #metatranscriptomics
# code_dir=this.path::here()
# for (type in c( "metatranscriptomics" )){
#   output_folder=paste0("/mnt/raid1b/philip/schinos_seq/results/",type,"/abundances")
#   print("")
#   agamemnon_output=paste0("/mnt/raid1b/philip/schinos_seq/results/",type,"/agamemnon_output/final/all")
#   for (i in c("SB","SP")){
#     print("")
#     rmarkdown::render(file.path(code_dir,'04_create_abundance_tables.Rmd'),
#                       output_file = file.path(output_folder,"junk.pdf"),
#                       params = list(
#                         tissue=i,
#                         counts_folder=agamemnon_output,
#                         sequencing=type,
#                         output_folder=output_folder))
#   }
# }
