#############################################################################################################################

### Identify MBFC and SFC populations in integrated scRNAseq data and estimate differentially expressed genes  

############################################################################################################################


library(Seurat)
library(tidyverse)
library(grDevices)
library(RColorBrewer)
library("ggsci")
library("ggbeeswarm")
library("ggthemes")
library(viridis)
library(gridExtra)

library("clusterProfiler")
library("org.Dm.eg.db")
library("AnnotationHub")
library(enrichplot)
library(clusterProfiler)
library(phateR)
library(slingshot)
library(scater)
library(tradeSeq)
library(reshape2)

projectdir <- file.path("/projectnb/mccall/sbandyadka/mmNCC")
source(file.path(projectdir, "scripts", "scRNAseq_integration","integration_utils.R")) 
input_dir <- file.path(projectdir,  "scrnaseq","data/")
output_dir <- file.path(projectdir,  "scrnaseq","analysis/")
figures <- file.path(projectdir,  "scrnaseq","analysis","figures/")

dir.create(figures)

combined.clusters <- readRDS(file.path(output_dir,"integrated.RDS"))
allmarkers <- readRDS(file.path(output_dir,"allmarkers.RDS"))
genes.of.interest <- loadGeneLists()

integrated_clusters <- DimPlot(combined.clusters, label=FALSE,group.by="orig.ident",pt.size=0.5)+
  scale_color_manual(values=c("#B482C2","#48344d"))
ggsave(integrated_clusters, filename = file.path(figures,"integratedclusters.png"),width = 7,height=7)

DefaultAssay(combined.clusters) <- "SCT"
combined.clusters <- PrepSCTFindMarkers(object = combined.clusters)

trapseq_DEresults <- read.table(file.path(projectdir, "translatome", "analysis", "deseq2", "genotype_PG150_vs_GR1_lfc.tsv"),
                               sep="\t", header=TRUE)
trapseq_sfcgenes <- trapseq_DEresults %>% filter(padj < 0.05 & log2FoldChange > 4) %>% drop_na() %>% pull(gene) 
sfcgenes <- append(genes.of.interest$sfc_markers,trapseq_sfcgenes)
combined.clusters <- AddModuleScore(
  object = combined.clusters,
  features = sfcgenes,nbin=10,
  ctrl=5,
  name = 'SFCscore'
)

sfcscoreplot <- FeaturePlot(combined.clusters,features = "SFCscore1", pt.size = 0.5) + 
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu")) 
sfcscoreplot

ggsave(sfcscoreplot, filename = paste0(figures,"/sfcscore.png"),width = 7,height=7)


trapseq_mbfcgenes <- trapseq_DEresults %>% filter(padj < 0.05 & log2FoldChange < -5) %>% drop_na() %>% pull(gene) 
mbfcgenes <- c("Yp1","Sox14","br","CG3397","CG6403","Nach","CG11600","Cyp4p2","CG7080")
combined.clusters <- AddModuleScore(
  object = combined.clusters,
  features = mbfcgenes,nbin=10,
  ctrl=5,
  name = 'MBFCscore'
)

mbfcscoreplot <- FeaturePlot(combined.clusters,features = "MBFCscore1") +  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlGnBu"))
mbfcscoreplot
ggsave(mbfcscoreplot, filename = paste0(figures,"/mbfcscore.png"),width = 7,height=7)

#identify cells with high sfc and mbfc module scores and subset integrated dataset 
highsfcscore <- FetchData(combined.clusters,vars = "SFCscore1") %>% filter(SFCscore1>1) %>% rownames()
highmbfcscore <- FetchData(combined.clusters,vars = "MBFCscore1") %>% filter(MBFCscore1>2) %>% rownames()

combined.clusters$trap.equivalentcells<- ifelse(colnames(combined.clusters) %in% highsfcscore, "SFC",
                                    ifelse(colnames(combined.clusters) %in% highmbfcscore, "MBFC",
                                           ifelse(colnames(combined.clusters) %in% highmbfcscore &
                                                    colnames(combined.clusters) %in% highsfcscore, "both","other")))




mbfc_sfc_cells <- FetchData(combined.clusters,vars = "trap.equivalentcells") %>% filter(trap.equivalentcells=="MBFC"|trap.equivalentcells=="SFC") %>% rownames()
mbfc_sfc_subset <- subset(combined.clusters,cells=mbfc_sfc_cells)
DimPlot(mbfc_sfc_subset,label=TRUE,group.by = "trap.equivalentcells")  +
  scale_color_manual(values=c("#068E8C","#65BADA"))


DefaultAssay(mbfc_sfc_subset) <- "integrated"
trapequiv_integratedplot <- DimPlot(combined.clusters,group.by = "trap.equivalentcells") +
  scale_color_manual(values=c("#068E8C","grey","#65BADA"))
trapequiv_integratedplot
ggsave(trapequiv_integratedplot, filename = paste0(figures,"/trapequiv_integrated.png"),width = 7,height=7)

mbfc_sfc_subset <- RunPCA(mbfc_sfc_subset, verbose = TRUE)
mbfc_sfc_subset <- RunUMAP(mbfc_sfc_subset, reduction = "pca", dims = 1:15)
mbfc_sfc_subset <- FindNeighbors(mbfc_sfc_subset, reduction = "pca", dims = 1:15)
mbfc_sfc_subset <- FindClusters(mbfc_sfc_subset, resolution = 0.5,reduction = "umap",verbose=TRUE)

mbfc_sfc_subsetplot <- DimPlot(mbfc_sfc_subset,label=FALSE,group.by = "trap.equivalentcells",pt.size=0.5)+
  scale_color_manual(values=c("#068E8C","#65BADA")) + coord_flip()
mbfc_sfc_subsetplot
ggsave(mbfc_sfc_subsetplot, filename = paste0(figures,"/mbfc_sfc_subset.png"),width = 7,height=7)

DimPlot(mbfc_sfc_subset,label=TRUE,group.by = "known_sfc_cells")
DimPlot(mbfc_sfc_subset,label=TRUE,group.by = "rust2020.ident")
DimPlot(mbfc_sfc_subset,label=TRUE,group.by = "jevitt2020.ident")

# Check expression levels of markers in MBFCs and SFCs and identify differetially expressed genes in the two groups
DefaultAssay(mbfc_sfc_subset) <- "SCT"
Idents(mbfc_sfc_subset) <- mbfc_sfc_subset$trap.equivalentcells

modulescoregenes_dotplot <- DotPlot(mbfc_sfc_subset,features=append(mbfcgenes,sfcgenes))+RotatedAxis()+ coord_flip() +
  theme_classic(base_size=16) +
  theme(axis.text.y = element_text(size=18,color="black",face="italic"),
        axis.text.x = element_text(size=18,color="black"))  & 
  scale_colour_gradientn(colours = brewer.pal(n =8, name = "YlGnBu")) 
modulescoregenes_dotplot
ggsave(modulescoregenes_dotplot, filename = paste0(figures,"/mbfc_sfc_modulescoremarkers.png"),width = 7,height=9)

classicalmarkers <- DotPlot(mbfc_sfc_subset,features=c("br","Yp1","mirr","cv-2","dpp","eya"))+RotatedAxis()+ coord_flip()+
  theme_classic(base_size=16) +
  theme(axis.text.y = element_text(size=20,color="black"),
        axis.text.x = element_text(size=20,color="black"))  & 
  scale_colour_gradientn(colours = brewer.pal(n =3, name = "YlGnBu")) 
classicalmarkers 
ggsave(classicalmarkers, filename = paste0(figures,"/mbfc_sfc_classicalmarkers.png"),width = 7,height=7)


DefaultAssay(mbfc_sfc_subset) <- "SCT"
mbfc_sfc_subset <- PrepSCTFindMarkers(object = mbfc_sfc_subset)
mbfc_vs_sfc <- FindMarkers(mbfc_sfc_subset,ident.1="SFC", ident.2="MBFC",
                           only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "DESeq2")
saveRDS(mbfc_vs_sfc,file.path(paste0(output_dir,"/mbfc_vs_sfc_trapequiv_markers.RDS")))
#mbfc_vs_sfc <- readRDS(file.path(paste0(output_dir,"/mbfc_vs_sfc_trapequiv_markers.RDS")))



comparisondf <- compare_translatome_scrnaseq_foldchanges(mbfc_vs_sfc,trapseq_DEresults,scalexy = c(-5,5),output_dir = figures)

lfcdifferential_dist <- comparisondf %>% mutate("fcdiff"=log2FoldChange-avg_log2FC) %>%
  filter(signifcolors != "#FBDCC0") %>% 
  dplyr::select(fcdiff,signifcolors) %>% 
  ggplot(aes(x=fcdiff,fill=factor(signifcolors),color=factor(signifcolors)))+
  geom_density(alpha=0.7)+
  scale_fill_identity() +
  scale_color_identity() +
  xlab("Log2FC Differential (Translatome - scRNAseq)")+ 
  theme_classic(base_size=16)+
  theme(axis.text = element_text(color="black",size=20))

lfcdifferential_dist 
ggsave(lfcdifferential_dist, filename = paste0(figures,"/lfcdiff_dist.png"),width = 7,height=7)
