###################################################################################################################
####### Perform secretome QC 
###################################################################################################################

library(artMS)
library('org.Dm.eg.db')
library(stringr)
library(tidyverse)
library(ggrepel)
library(ggthemes)

inputdir <-  file.path("/projectnb/mccall/sbandyadka/mmNCC/")
setwd(file.path(paste0(inputdir),'scripts','secretome'))
resultsdir <- file.path(paste0(inputdir),'secretome','analysis/')

artmsWriteConfigYamlFile(config_file_name = file.path(paste0(inputdir),'secretome','analysis', "artMS_config.yaml"), 
                           overwrite = TRUE)

evidence <- as.data.frame(read.table(file.path(paste0(inputdir),'secretome','data', 'evidence.txt'),
                          sep="\t",header=1,stringsAsFactors=FALSE))

phagoptosiskeys <- as.data.frame(read.table(file.path(paste0(inputdir),'scripts', 'secretome','artMS_phagoptosis_keys.txt')),
                               sep="\t",header=1,stringsAsFactors=FALSE)
phagoptosiskeys <- phagoptosiskeys[-1,]
colnames(phagoptosiskeys) <- c("Raw.file", "Condition", "IsotopeLabelType", "BioReplicate", "Run")


starvationkeys <- as.data.frame(read.table(file.path(paste0(inputdir),'scripts', 'secretome','artMS_starvation_keys.txt')),
                                            sep="\t",header=1,stringsAsFactors=FALSE)
starvationkeys <- starvationkeys[-1,]
colnames(starvationkeys) <- c("Raw.file", "Condition", "IsotopeLabelType", "BioReplicate", "Run")


runartMS <- function(contrastname,outdir,datadir){
  dir.create(file.path(outdir), showWarnings = FALSE)
  setwd(outdir)
  #expnames <- c("APMS","AB","UB","AC")
  expnames <- c("AB")
  for(exp in expnames){
    print(exp)
    artmsQualityControlEvidenceBasic( 
      evidence_file = evidence,
      keys_file = contrastname,
      prot_exp = exp, 
      plotPTMSTATS = TRUE,
      plotINTDIST = TRUE, plotREPRO = TRUE,
      plotCORMAT = TRUE, plotINTMISC = TRUE,
      printPDF = TRUE, verbose = TRUE, output_dir= outdir)
  }
  
  artmsQualityControlEvidenceExtended(
    evidence_file = evidence,
    keys_file = contrastname)
  
  artmsQualityControlSummaryExtended(summary_file = file.path(inputdir,'secretome','data','summary.txt'),
                                     keys_file = contrastname)
  

  
}


runartMS(phagoptosiskeys,paste0(resultsdir,"phagoptosis"),inputdir)
runartMS(starvationkeys,paste0(resultsdir,"starvation"),inputdir)




