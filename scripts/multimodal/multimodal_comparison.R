######################################################################
# This script takes in DE outputs from translatome, singlecell, secretome
# And outputs DE fold change and adjusted p-values across datasets 
# in phagoptosis and starvation comparisons separately 
# and plots  fold change comparisons of genes across datasets
######################################################################

library(tidyverse)
library(reshape2)
library("Hmisc")
library(corrplot)
library(GGally)
library("ggradar")
library(scales)
library("reactablefmtr")


project_dir <- file.path("/projectnb2/mccall/sbandyadka/mmNCC/")
setwd(file.path(project_dir,'scripts'))
source(file.path(project_dir,'scripts', 'translatome','diffExpression_deseq2', 'utils.R'))
output_directory <- file.path(paste0(project_dir),'comparison/')

genelists <- get_genelists(file.path(project_dir, 'translatome','analysis', 'glad_annotations'))

secretome_phagoptosis <- read.table(file.path(project_dir,'secretome','analysis','FedSFCs_vs_FedMBFCs.tsv'),sep="\t",header=TRUE)
translatome_phagoptosis <- read.table(file.path(project_dir,'translatome','analysis','deseq2','genotype_PG150_vs_GR1_lfc.tsv'),sep="\t",header=TRUE)
scrnaseq_phagoptosis <- readRDS(file.path(project_dir,'scrnaseq','analysis','mbfc_vs_sfc_trapequiv_markers.RDS'))

## Consolidate Log2FC and pvalues of SFC vs MBFC results from translatome, scRNAseq and secretome in one dataframe.
scrnaseq_phagoptosis$gene <- rownames(scrnaseq_phagoptosis)
secretome_phagoptosis$gene <- secretome_phagoptosis$genesymbol

trsc <- full_join(translatome_phagoptosis,scrnaseq_phagoptosis,
                                       by = 'gene', suffix = c(".translatome", ".scAtlas"))
phagoptosis_3way <- full_join(trsc,secretome_phagoptosis, 
                  by = 'gene', suffix = c(".trsc", ".sec"))

phagoptosis_3way <- phagoptosis_3way %>% select(gene,log2FoldChange,padj,avg_log2FC,p_val_adj,ratio,p.adj)
colnames(phagoptosis_3way) <- c("gene","translatomeL2FC","translatomePADJ","scrnaseqL2FC","scrnaseqPADJ","secretomeL2FC","secretomePADJ")
phagoptosis_3way$gene <- make.unique(as.character(phagoptosis_3way$gene), sep = "_")
phagoptosis_3way <- phagoptosis_3way[!is.na(phagoptosis_3way$gene),]

#rownames(phagoptosis_3way) <- phagoptosis_3way$gene
#phagoptosis_3way <- phagoptosis_3way %>% select(-gene)
head(phagoptosis_3way)

write.table(phagoptosis_3way,file = file.path(output_directory,"/phagoptosis_comparisonsummary.tsv"),sep="\t",row.names = FALSE)

## Consolidate starvation-induced L2FCs from translatome and secretome. 
secretome_starvation <- read.table(file.path(project_dir,'secretome','analysis','StarvedMBFCs_vs_FedMBFCs2.tsv'),sep="\t",header=TRUE)
translatome_starvation <- read.table(file.path(project_dir,'translatome','analysis','deseq2','treatment_starved_vs_fed_lfc.tsv'),sep="\t",header=TRUE)
secretome_starvation$gene <- secretome_starvation$genesymbol
starvation_combined <- full_join(translatome_starvation,secretome_starvation,
                  by = 'gene', suffix = c(".translatome", ".secretome"))
starvation_combined <- starvation_combined %>% select(gene,log2FoldChange,padj,ratio,p.adj)
colnames(starvation_combined) <- c("gene","translatomeL2FC","translatomePADJ","secretomeL2FC","secretomePADJ")
starvation_combined$gene <- make.unique(as.character(starvation_combined$gene), sep = "_")
starvation_combined <- starvation_combined[!is.na(starvation_combined$gene),]
write.table(starvation_combined,file = file.path(output_directory,"/starvation_comparisonsummary.tsv"),sep="\t",row.names = FALSE)




### How many genes are identified across all comparisons with padj < 0.05
phagoptosis_siginall <- phagoptosis_3way %>% filter(translatomePADJ < 0.05 & scrnaseqPADJ < 0.05 & secretomePADJ < 0.05) %>% 
  as_tibble() %>% select(-contains("PADJ")) 

starvation_siginall <- starvation_combined %>% filter(translatomePADJ < 0.05  & secretomePADJ < 0.05) %>% 
  as_tibble() %>% select(-contains("PADJ")) 



cytoskeletal <- phagoptosis_3way %>% filter(gene %in% genelists$cytoskeletal)



