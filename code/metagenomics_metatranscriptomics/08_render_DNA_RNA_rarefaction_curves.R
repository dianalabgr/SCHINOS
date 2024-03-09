library(this.path())
code_dir=this.path::here()
project_dir=this.path::here(..=1)

#make rarefaction curves
for (type in c("metatranscriptomics" )){
  agamemnon_output=paste0("/mnt/raid1b/philip/schinos_seq/results/",type,"/agamemnon_output/final/all")
  report_dir=file.path(this.path::here(..=1),"results",type,"rarefaction_curves")
  for (i in c("SB","SP")){
    rmarkdown::render(file.path(code_dir,'08_DNA_RNA_rarefaction_curves.Rmd'),
                      output_file = paste0(report_dir,'/rarefaction_curves_',i,'_',type,'_', Sys.Date(), sep=''),
                      params = list(
                        tissue=i,
                        counts_folder=agamemnon_output,
                        sequencing=type))
  }
}
for (type in c("metatranscriptomics" )){
  agamemnon_output=paste0("/mnt/raid1b/philip/schinos_seq/results/",type,"/agamemnon_output/final/all")
  report_dir=file.path(this.path::here(..=1),"results",type,"rarefaction_curves")
  for (i in c("SΑ")){
    rmarkdown::render(file.path(code_dir,'08_DNA_RNA_rarefaction_curves.Rmd'),
                      output_file = paste0(report_dir,'/rarefaction_curves_',i,'_',type,'_', Sys.Date(), sep=''),
                      params = list(
                        tissue=i,
                        counts_folder=agamemnon_output,
                        sequencing=type))
  }
}
for (type in c("metagenomics" )){
  agamemnon_output=paste0("/mnt/raid1b/philip/schinos_seq/results/",type,"/agamemnon_output/final/all")
  report_dir=file.path(this.path::here(..=1),"results",type,"rarefaction_curves")
  for (i in c("SB","SA","SP","GF")){
    rmarkdown::render(file.path(code_dir,'08_DNA_RNA_rarefaction_curves.Rmd'),
                      output_file = paste0(report_dir,'/rarefaction_curves_',i,'_',type,'_', Sys.Date(), sep=''),
                      params = list(
                        tissue=i,
                        counts_folder=agamemnon_output,
                        sequencing=type))
  }
}