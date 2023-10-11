###################################################################################################################
#######  Main script for running - 
########### 1. Differential expression analysis and sample QC on salmon quant files
########### 2. Compare translational profiles of starvation and phagoptosis 
########### 3. Perform GO/functional enrichment
####### Calls - utils.R
# Author: Shruthi Bandyadka, McCall Lab, Boston University (sbandya@bu.edu)
###################################################################################################################

library("tximport")
library("readr")
library("tximportData")
library("tximeta")
library("DESeq2")
library(org.Dm.eg.db)
library("tidyverse")
library("AnnotationDbi")
library("clusterProfiler")
library("ggrepel")
library("ggbeeswarm")
library("pheatmap")
library("vsn")
library("RColorBrewer")
library("gplots")
library("reshape2")
library(extrafont)
library(enrichplot)
library(ggnewscale)
library(ggthemes)
library(ggbreak) 
require(RColorBrewer)
require(ComplexHeatmap)
require(circlize)
require(digest)
require(cluster)
library(stringr)
library(MoMAColors)

## Set up directories
project_dir <- file.path("/projectnb/mccall/sbandyadka/mmNCC/")
setwd(file.path(paste0(project_dir),'scripts','translatome','diffExpression_deseq2'))
source("utils.R")

input_directory <- file.path(paste0(project_dir),'translatome','analysis','salmon_BDGP6_32_104')
output_directory <- file.path(paste0(project_dir),'translatome','analysis','deseq2')
dir.create(output_directory, showWarnings = FALSE)


## Prepare salmon counts data

samples <- read.table(file.path(input_directory,"design_matrix.txt"), sep="\t", header=TRUE)
rownames(samples) <- samples$sample

files <- file.path(input_directory,samples$sample, "quant.sf")
names(files) <- samples$sample
coldata <- samples
coldata$files <- files
coldata$names <- coldata$sample
se <- tximeta(coldata)
gse <- summarizeToGene(se)
gse <- addIds(gse, "REFSEQ", gene=TRUE)

dds <- DESeqDataSet(gse, design = ~ genotype + treatment)
colSums(counts(dds))
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
resultsNames(dds)


## Run differential expression analysis, GO-term/functional enrichment 
genotype_lfc_results <- get_contrast_result( "genotype_PG150_vs_GR1",dds,'lfc','ENSEMBL',output_directory,"#068E8C","#65BADA",5)
treatment_lfc_resutls <- get_contrast_result( "treatment_starved_vs_fed",dds,'lfc','ENSEMBL',output_directory,"#D86F27","#068E8C",2)

## Plot PCA of replicates after Rlog transformation
rld <- rlog(dds, blind=FALSE)

pcaData <- plotPCA(rld, intgroup=c("genotype", "sample"), returnData=TRUE)
pcaData <- pcaData %>% mutate(colors = ifelse(grepl('P', sample),"#65BADA",
                                       ifelse(grepl('GC', sample),"#068E8C","#D86F27")))
percentVar <- round(100 * attr(pcaData, "percentVar"))
rld_pca <- ggplot(pcaData, aes(PC1, PC2, color=colors)) +
  geom_point(size=6) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_colour_identity() + 
  coord_fixed() +
  theme_clean(base_size = 18) 
rld_pca
ggsave(rld_pca,filename = paste0(output_directory,"/pca.png"),width = 6,height=6)



## Plot counts for specific genes

majorpathwaygenes <- c("Dad","wgn","kek1","E(spl)mbeta-HLH","numb","PGRP-LC","PGRP-SC2","PGRP-SD","SPE","Ilp6")

majorpathwaygenes_counts <- c()
for(g in majorpathwaygenes){
  genesymbol <- mapIds(org.Dm.eg.db, keys=g, column="ENSEMBL",keytype="SYMBOL", multiVals="first")
  counts <- plotCounts(dds, gene=genesymbol, intgroup=c("genotype","treatment"), returnData=TRUE)
  counts <- counts %>% mutate(symbol=g) %>% mutate(sample=paste0(genotype,":",treatment)) %>%
    filter(sample=="PG150:fed" | sample=="GR1:fed") %>% select(symbol,count,sample)
  majorpathwaygenes_counts[[g]] <- counts
}
majorpathwaygenes_counts_df <- do.call(rbind, majorpathwaygenes_counts)

sfc_pathways_counts <- ggplot(majorpathwaygenes_counts_df, aes(x=sample, y=count, fill=sample)) + 
  geom_boxplot() +
  geom_line() +
  facet_grid(cols = vars(symbol)) +
  theme_clean(base_size = 18) + 
  ylab("aggregared gene-level read counts") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        strip.text.x = element_text(size = 11,face = "italic"),
        panel.border = element_rect(color = "#eeeeee", fill = NA, size = 1,linetype = "dotted")) + 
  scale_fill_manual(values=c("#068E8C","#65BADA"),guide="none") 

sfc_pathways_counts
ggsave(sfc_pathways_counts,filename = paste0(output_directory,"/sfc_pathwayscounts.png"),width = 15,height=7)

## Toll genes enriched in SFCs discovered using PANGEA
tollgenes_inSFCs <- c("PGRP-SD","PGRP-LC","SPE","Usp16-45","Pp2A-29B", "BomS1","BomS3","for",
                      "mnd","senju",
                      "Drsl4","Drs","Tsf1","CG16799","Tsf2",
                      "Listericin","Def","CG43236")
#tollgenes_inSFCs <- c("PGRP-SD","PGRP-LC","SPE","Drs")
tollgenes_inSFCs_counts <- c()
for(g in tollgenes_inSFCs){
  genesymbol <- mapIds(org.Dm.eg.db, keys=g, column="ENSEMBL",keytype="SYMBOL", multiVals="first")
  counts <- plotCounts(dds, gene=genesymbol, intgroup=c("genotype","treatment"), returnData=TRUE)
  counts <- counts %>% mutate(symbol=g) %>% mutate(sample=paste0(genotype,":",treatment)) %>%
    filter(sample=="PG150:fed" | sample=="GR1:fed") %>% select(symbol,count,sample)
  tollgenes_inSFCs_counts[[g]] <- counts
}
tollgenes_inSFCs_counts_df <- do.call(rbind, tollgenes_inSFCs_counts)
tollgenes_inSFCs_counts_plot <- ggplot(tollgenes_inSFCs_counts_df,
                                       aes(x=sample, y=count, fill=sample, color=sample)) + 
  geom_boxplot(width=0.3, alpha=0.5) +
  geom_line() +
  facet_wrap(~factor(symbol, levels=tollgenes_inSFCs),scales = "free_y",ncol=4)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
  theme_clean(base_size = 16) + 
  ylab("Read counts per gene") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y =element_text(size = 16),
        axis.title.y =element_text(size = 16),
        strip.text.x = element_text(size = 16,face = "italic"),
        panel.border = element_rect(color = "#eeeeee", fill = NA, size = 1,linetype = "dotted")) + 
  scale_fill_manual(values=c("#068E8C","#65BADA"),guide="none") +
  scale_color_manual(values=c("#068E8C","#65BADA"),guide="none")

tollgenes_inSFCs_counts_plot
ggsave(tollgenes_inSFCs_counts_plot,filename = paste0(output_directory,"/tollgenes_inSFCs_counts_plot.png"),width = 12,height=5)

## Plot summary of DE genes.

desummary_genotype <- genotype_lfc_results %>% as.data.frame() %>% mutate(enrichment=ifelse(log2FoldChange>0, "SFC","MBFC")) %>%
  mutate(Adj.p = ifelse(padj >= 0.05, ">= 0.05","< 0.05")) %>% filter(Adj.p=="< 0.05" | Adj.p==">= 0.05") %>%
  mutate(Adj.p = factor(Adj.p, levels=c(">= 0.05","< 0.05"))) %>% 
  ggplot(aes(x=enrichment,fill=Adj.p))+ 
  geom_bar(width=0.2)+  scale_y_break(c(1000, 4000)) +
  xlab("condition") +
  ylab("number of DE genes") +
  scale_fill_manual(values=c("#c7e4ef","#65BADA")) + 
  theme_classic(base_size = 16)
desummary_genotype
ggsave(desummary_genotype,filename = paste0(output_directory,"/genotype_desummary.png"),width = 5,height=5)

desummary_treatment <- treatment_lfc_resutls %>% as.data.frame() %>% mutate(enrichment=ifelse(log2FoldChange>0, "starved","fed")) %>%
  mutate(Adj.p = ifelse(padj >= 0.05, ">= 0.05","< 0.05")) %>% filter(Adj.p=="< 0.05" | Adj.p==">= 0.05") %>%
  mutate(Adj.p = factor(Adj.p, levels=c(">= 0.05","< 0.05"))) %>% 
  ggplot(aes(x=enrichment,fill=Adj.p))+ 
  geom_bar(width=0.2)+  scale_y_break(c(100, 5000)) +
  xlab("condition") +
  ylab("number of DE genes") +
  scale_fill_manual(values=c("#f2b083","#D86F27")) + 
  theme_classic(base_size = 16)
desummary_treatment
ggsave(desummary_treatment,filename = paste0(output_directory,"/treatment_desummary.png"),width = 5,height=5)


## Plot heatmap of specific gene groups
vsd <- vst(dds, blind=FALSE)
vsd_raw <- as.data.frame(assay(vsd))
vsd_scaled <- as.data.frame(t(scale(t(vsd_raw)))) 



vsd_scaled$gene_symbol <- mapIds(org.Dm.eg.db,
                                 keys=rownames(vsd_scaled),
                                 column = 'SYMBOL',
                                 keytype = 'FLYBASE',
                                 multiVals="first")

genelists <- get_genelists(file.path(paste0(project_dir),'translatome','analysis','glad_annotations'))
genotype_lfc_results_df <- genotype_lfc_results %>% as.data.frame() %>% filter(padj < 0.05) 
genotype_lfc_results_df$genesymbol <- mapIds(org.Dm.eg.db, keys=rownames(genotype_lfc_results_df), column="SYMBOL",keytype="ENSEMBL", multiVals="first")
sfc_genes <- genotype_lfc_results_df %>%  filter(log2FoldChange> 2) %>% pull(genesymbol)
fedmbfc_genes <- genotype_lfc_results_df %>%  filter(log2FoldChange < -2) %>% pull(genesymbol)

sfc_patterns <- sfc_genes[sfc_genes %in% genelists$transmembrane] %>% sort() 
fedmbfc_patterns <- fedmbfc_genes[fedmbfc_genes %in% genelists$transmembrane] %>% sort()
patterns <- append(sfc_patterns,fedmbfc_patterns)

#patterns <- patterns[!str_detect(patterns,pattern="^CG")]


genegroup_heatmap <- vsd_scaled %>% select(contains("P") | contains("GC") | contains("gene")) %>%
  filter(gene_symbol %in% patterns) %>% melt() %>% 
  ggplot(aes(gene_symbol,variable, fill= value)) + 
  geom_tile() + 
  coord_equal() + 
  ylab(" ") +
  xlab("Gene") +
  labs(fill='scaled \ncounts')  +
  scale_fill_moma_c("Picabia") +
  scale_x_discrete(limits = patterns) +
  theme_classic(base_size = 16) +
  theme(axis.text.y = element_blank(),
    axis.text.x = element_text(angle=45,vjust=0.9,
                               hjust=0.9))
genegroup_heatmap

ggsave(genegroup_heatmap,filename = paste0(figure_dir,"/genotype_genegroupheatmap_transmembrane.png"),width = 15,height=3)

tollgenes_heatmap <- vsd_scaled %>% select( contains("GC") | contains("P") | contains("gene")) %>% 
  filter(gene_symbol %in% tollgenes_inSFCs) %>% melt() %>%
  ggplot(aes(x=gene_symbol,y=variable, fill= value)) + 
  geom_tile() + 
  coord_equal() + 
  ylab(" ") +
  xlab("Gene") +
  labs(fill='scaled \ncounts')  +
  scale_fill_moma_c("Picabia") +
  scale_x_discrete(limits = tollgenes_inSFCs) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle=30,vjust=0.9, face = "bold",
                                   hjust=0.9), axis.text.y = element_blank()
        )
tollgenes_heatmap

ggsave(tollgenes_heatmap,filename = paste0(output_directory,"/genotype_genegroupheatmap_tollpathway.png"),width = 10,height=3)

treatment_lfc_resutls_df <- treatment_lfc_resutls %>% as.data.frame() %>% filter(padj < 0.05) 
treatment_lfc_resutls_df$genesymbol <- mapIds(org.Dm.eg.db, keys=rownames(treatment_lfc_resutls_df), column="SYMBOL",keytype="ENSEMBL", multiVals="first")
starved_genes <- treatment_lfc_resutls_df %>%  filter(log2FoldChange> 0) %>% pull(genesymbol)
fedmbfc_genes <- treatment_lfc_resutls_df %>%  filter(log2FoldChange < 0) %>% pull(genesymbol)

starved_patterns <- starved_genes[starved_genes %in% genelists$secreted] %>% sort() 
fedmbfc_patterns <- fedmbfc_genes[fedmbfc_genes %in% genelists$secreted] %>% sort()
patterns <- append(starved_patterns,fedmbfc_patterns)

genegroup_heatmap <- vsd_scaled %>% select(contains("GC") | contains("GS") | contains("gene")) %>%
  filter(gene_symbol %in% patterns) %>% melt() %>% 
  ggplot(aes(gene_symbol,variable, fill= value)) + 
  geom_tile() + 
  coord_equal() + 
  ylab(" ") +
  xlab("Gene") +
  labs(fill='scaled \ncounts')  +
  scale_fill_moma_c("Picabia") +
  scale_x_discrete(limits = patterns) +
  theme_classic(base_size = 18) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle=45,vjust=0.9,
                                   hjust=0.9))
genegroup_heatmap
ggsave(genegroup_heatmap,filename = paste0(output_directory,"/treatment_genegroupheatmap_secreted.png"),width = 15,height=4)


patterns <- c("^Hsp|^Hsf")
genegroup_heatmap <- vsd_scaled %>% 
  mutate(genegroupname = ifelse(grepl(paste(patterns, collapse="|"),gene_symbol),"genegroup","other")) %>% 
  select(contains("P") | contains("GC") | contains("GS") | contains("gene")) %>%
  filter(genegroupname == "genegroup") %>% melt() %>%
  ggplot(aes(gene_symbol,variable, fill= value)) + 
  geom_tile() + 
  coord_equal() + 
  xlab("Replicate") +
  ylab("Gene") +
  scale_fill_moma_c("Picabia") +
  theme_classic() +
  theme(#axis.text.y = element_blank(),
        axis.text.x = element_text(angle=45,vjust=0.7,hjust=0.7))
genegroup_heatmap


ht <- vsd_scaled %>% 
  mutate(genegroupname = ifelse(grepl(paste(patterns, collapse="|"),gene_symbol),"genegroup","other")) %>%
  select(contains("GC") | contains("GS") | contains("P") | contains("gene")) %>%
  filter(genegroupname == "genegroup") %>% select(-one_of("genegroupname"))
rownames(ht) <- ht$gene_symbol
ht <- ht %>% select(-one_of("gene_symbol"))

coul <- colorRampPalette(brewer.pal(8, "YlGnBu"))(25)
heatmap(as.matrix(ht), col = coul,Colv = NA)


# Compare DEGs in phagoptosis vs. starvation

## Plot LFC quadrants in phagoptosis vs starvation
compare_degs <- compare_degs_prep_dds(dds)
comparison_plot <- plot_comparison_degs(compare_degs)
comparison_plot
ggsave(comparison_plot, filename = paste0(output_directory,"/starvation_vs_phagoptosis_scatter.png"),width = 15,height=10)

## Plot counts of concordant/discordant categories
category_counts <- compare_degs %>% drop_na() %>% 
  filter(category !="Adj p >= 0.05 in both") %>%
  group_by(category) %>% summarise(n()) %>% 
  mutate(pie_labels=paste0(category,"-",`n()`)) 
  #mutate(perc=round(`n()`*100/sum(`n()`),digits=1)) %>%
  #mutate(pie_labels=paste0(category,"-",perc,"%"))
png(file=paste0(output_directory,"/starvation_vs_phagoptosis_categorycounts.png"),width = 7,height=7,units="in",res=300)
pie(labels=category_counts$pie_labels,x=category_counts$`n()`,
    col=c("#9270b0","#503174","#65BADA","#D86F27","#12759a","#873803"),border = NA)
dev.off()


compare_categories_df <- compare_degs %>% drop_na() %>% filter(category != "Adj p >= 0.05 in both")
terms <- compare_go_annotation(compare_categories_df,output_directory)
write.table(as.data.frame(terms),paste0(output_directory,"/starvation_vs_phagoptosis_goterms.tsv"),sep="\t",row.names = FALSE)
write.table(compare_degs,paste0(output_directory,"/starvation_vs_phagoptosis_categories.tsv"),sep="\t",row.names = FALSE)


compare_categories_df$geneid <- mapIds(org.Dm.eg.db, 
                                   keys=compare_categories_df$gene, 
                                   column="ENSEMBL", 
                                   keytype="SYMBOL",
                                   multiVals="first")

### Distribution of differences in LFC in each category

lfcdiff_dist <- compare_degs %>% drop_na() %>% filter(category != "Adj p >= 0.05 in both") %>%  
  mutate("fcdiff"= abs(log2FoldChange.phagoptosis - log2FoldChange.starvation)) %>%
  dplyr::select(fcdiff,category) %>%  
  ggplot(aes(x=abs(fcdiff),color=factor(category),fill=factor(category)))+
  geom_density(alpha=0.01,size=2)+
  scale_color_manual(values=c("#9270b0","#503174","#65BADA","#D86F27","#12759a","#873803"))+
  scale_fill_manual(values=c("#9270b0","#503174","#65BADA","#D86F27","#12759a","#873803"))+
  geom_vline(xintercept = 0,linetype="dotted")+
  xlab("|(PhagoptosisLFC - StarvationDeathLFC)|")+ 
  theme_classic(base_size=16)+
  theme(axis.text = element_text(color="black",size=20),
        legend.position="none")

lfcdiff_dist
ggsave(lfcdiff_dist, filename = paste0(output_directory,"/starvation_vs_phagoptosis_lfcdiff_dist.png"),width = 7,height=7)
write.table(compare_degs, file = paste0(output_directory,"/starvation_vs_phagoptosis_summary.tsv"),sep="\t",row.names = FALSE)


## Plot GLAD gene group enrichment
sfc_glad_genegroup_enrichment <- read.table(file.path(paste0(project_dir),'translatome','analysis','glad_annotations','sfc_genegroups_glad.txt'),sep="\t",header=TRUE)
mbfc_glad_genegroup_enrichment <- read.table(file.path(paste0(project_dir),'translatome','analysis','glad_annotations','mbfc_genegroups_glad.txt'),sep="\t",header=TRUE)

mbfc_sfc_genegroups <- mbfc_glad_genegroup_enrichment %>%
  full_join(sfc_glad_genegroup_enrichment, by = "Group", keep = TRUE,suffix=c(".MBFC",".SFC"))

sfc_glad_genegroup_enrichment <- sfc_glad_genegroup_enrichment %>% mutate(condition="SFC")%>%
  mutate(sfc_genelabels=genes )%>%
  mutate(mbfc_genelabels="" )
mbfc_glad_genegroup_enrichment <- mbfc_glad_genegroup_enrichment %>% mutate(condition="MBFC") %>%
  mutate(X..p.values.log10 = -X..p.values.log10 )%>%
  mutate(mbfc_genelabels = genes )%>%
  mutate(sfc_genelabels = "" )
mbfc_sfc_genegroups_row <- rbind(mbfc_glad_genegroup_enrichment,sfc_glad_genegroup_enrichment)

mbfc_sfc_genegroups_row_plot <- mbfc_sfc_genegroups_row %>% mutate("gene_fraction"= paste0(Count,"/",genes.in.list)) %>%
  ggplot(aes(x=X..p.values.log10,y=reorder(Group,X..p.values.log10),fill=factor(condition)))+ 
  geom_bar(stat="identity", width = 0.5, position="identity") +
  #geom_text(aes(label=gene_fraction),vjust=0.01,angle=0,hjust=0.1,size=3,color="black") +
  geom_text(aes(label=sfc_genelabels),vjust=0.5,angle=0,hjust=-0.1,size=3,fontface="italic",color="#65BADA") +
  geom_text(aes(label=mbfc_genelabels),vjust=0.5,angle=0,hjust=1.2,size=3,fontface="italic",color="#068E8C") +
  scale_fill_manual(values=c("#068E8C","#65BADA"))+
  xlab("-Log10(p-value)")+
  ylab("GLAD gene group")+
  scale_x_continuous(breaks = seq(-10, 10, by = 5),limits = c(-10, 10),labels=c(10,5,0,5,10)) +
  theme_classic(base_size = 15)+
  theme(legend.title=element_blank())

mbfc_sfc_genegroups_row_plot


sfc_glad_enrichment_plot <- sfc_glad_genegroup_enrichment %>% mutate("gene_fraction"= paste0(Count,"/",genes.in.list)) %>%
  ggplot(aes(x=X..p.values.log10,y=reorder(Group,X..p.values.log10))) + 
  geom_bar(stat="identity",fill="#65BADA",alpha=0.8) +
  geom_text(aes(label=gene_fraction),vjust=0.5,angle=0,hjust=1,size=3,color="black") +
  geom_text(aes(label=genes),vjust=0.5,angle=0,hjust=-0.1,size=4,fontface="italic",color="#65BADA") +
  xlab("-Log10(p-value)")+
  ylab("GLAD gene group")+
  scale_x_continuous(breaks = seq(0, 10, by = 5),limits = c(0, 10)) +
  theme_classic(base_size = 18)

sfc_glad_enrichment_plot
ggsave(mbfc_sfc_genegroups_row_plot,filename = paste0(output_directory,"/sfc_glad_enrichment.png"),width = 15,height=10)


################################################################
#### Heatmap of top DE genes - VST normalize counts and subset to significant DEGs #######
################################################################

vsd <- vst(dds, blind=FALSE)
vsd_raw <- as.data.frame(assay(vsd))
vsd_raw['gene_symbol'] <- mapIds(org.Dm.eg.db,
                                 keys=rownames(vsd_raw),
                                 column = 'SYMBOL',
                                 keytype = 'FLYBASE',
                                 multiVals="first")

vsd_raw$fbgn <- rownames(vsd_raw) 
no_names <- vsd_raw[is.na(vsd_raw$gene_symbol),] # 22 genes have no symbol mappings, but these don't pass when selecting for DE significance

genotype_lfc_siggenes <- as.data.frame(genotype_lfc_results) %>% filter(padj < 0.05 & abs(log2FoldChange) > 3)
treatment_lfc_siggenes <- as.data.frame(treatment_lfc_resutls) %>% filter(padj < 0.05 & abs(log2FoldChange) > 3)
sigGenes <- c(rownames(genotype_lfc_siggenes), rownames(treatment_lfc_siggenes))

vsd_df <- vsd_raw[rownames(vsd_raw) %in% sigGenes,] # Subset to signif DE
dim(vsd_df) ##164 x 11

sig_no_names <- vsd_df[is.na(vsd_df$gene_symbol),] #check: all sig DEG have gene name annotations
rownames(vsd_df) <- vsd_df$gene_symbol # Set gene symbols as row labels in final heatmap
vsd_df_columns <- vsd_df %>% select (-c(fbgn, gene_symbol)) #select only sample counts columns for heatmap

################################################################
#### Draw heatmap ##############################################
################################################################
heat <- t(scale(t(vsd_df_columns)))
myCol <- colorRampPalette(c('#FFE6B3','#FFC857','#235FA4', '#0A284B'))(100)
myBreaks <- seq(-3, 3, length.out = 100)
colsplit = rep(1:3, each = 3)


ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=c('#FFC857','#235FA4', '#0A284B'),
                                                  col=c('#FFC857','#235FA4', '#0A284B')), 
                                        labels_gp = gpar(col = "white", fontsize = 8),
                                        labels = c("GR1-Conditioned","GR1-Starved","PG150-Conditioned")))
pdf(paste0(output_directory,"heatmap_p005_lfc2.pdf"), width = 8, height = 8)
hmap_p005_lfc2 <- Heatmap(heat,
                          col = colorRamp2(myBreaks, myCol),
                          row_names_gp = gpar(fontsize = 2.5),
                          column_names_gp = gpar(fontsize = 0),
                          cluster_rows = TRUE,
                          cluster_columns = FALSE,
                          column_dend_reorder = FALSE,
                          row_gap = unit(1, "mm"),
                          row_title = NULL,
                          column_title = NULL,
                          row_dend_width = unit(15,'mm'),
                          heatmap_legend_param = list(
                            title =  "counts",
                            title_gp = gpar(fontsize = 8,fontface="bold"), 
                            color_bar = 'continuous',
                            labels_gp=gpar(fontsize = 6)),
                          column_split = colsplit, 
                          row_split = 3,
                          top_annotation = ha,
)


draw(hmap_p005_lfc2,
     heatmap_legend_side = 'right',
     annotation_legend_side = 'right',
     row_sub_title_side = 'left',
     padding = unit(c(10, 10, 10, 10), "mm")
)

dev.off()







