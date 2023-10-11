
library(ggrepel)
library(ggthemes)
library(clusterProfiler)

loadGeneLists <- function(){
  genelists <- list(
  `rust2020_markers` = c("Mhc", "vas", "upd1", "cas", "en", "Wnt4", "zfh1", "pnt", "Tl","Sox14", "cv-2","slbo","br", "eya"),
  `jevitt2020_markers` = c("Past1",  "peb",  "Ilp6",  "Glut4EF",  "Ilp8",  "Diap1",  "Mmp2",  "abd-A",  "Hml",  "ct",  "br",  "Mp20",  "Wnt4",  "slbo",  "mirr",  "psd",  "Vml",  "dec-1",  "osk",  "orb",  "upd1",  "bbg",  "Atf3",  "cas",  "Cad74A",  "yellow-g",  "Fcp3C",  "ttk"),
  `sfc_markers` = c("cv-2","dpp","eya","Vha16-1","Vha100-2","Past1", "peb","drpr","puc", "trol", "kay"),
  `vatpase_a_subunits` = c("Vha100-1","Vha100-2","Vha100-3","Vha100-4","Vha100-5"),
  `v_atpases` = c("VhaM8.9","Vha13","Vha14-1","Vha16-1","Vha16-5","Vha26","Vha36-1","Vha36-3","Vha44","Vha55","Vha68-1","Vha68-2","Vha68-3","Vha100-1","Vha100-2","Vha100-4",
                 "VhaAC39-1","VhaAC39-2","VhaAC45","VhaM9.7-a","VhaM9.7-b","VhaM9.7-c","VhaPPA1-1","VhaSFD"),
  `motor_proteins` = c("Dhc16F","Dhc36C","Dhc62B","Dhc98D","Dnah3","kl-2","CG3339","Dhc1","Dhc93AB","kl-3","kl-5","btv","Dhc64C","cana","CG10845","CG14535","cmet","cos","Khc","Khc-73","Kif3C","Kif19A","Klp3A","Klp10A","Klp31E","Klp54D","Klp59C","Klp59D","Klp61F","Klp64D","Klp67A","Klp68D","Klp98A","ncd","neb","nod","pav","sub","unc-104","CG10793","Fign","Kat60","kat-60L1","spas","ck","didum","jar","Mhc","Myo10A","Myo28B1","Myo31DF","Myo61F","Myo81F","Myo95E","ninaC","zip"),
  `integrins` = c("if","ItgaPS4","ItgaPS5","Itgbn","mew","mys","scb"),
  `adaptor_proteins` = c("g","rb","cm","or","AP-1gamma","AP-1mu","AP-1sigma","AP-1-2beta","AP-2sigma","AP-2alpha")
  )
  #genelists <- list(rust2020_markers,jevitt2020_markers,sfc_tfc_markers,vatpase_a_subunits,v_atpases,motor_proteins,integrins,adaptor_proteins)
  return(genelists)
}

### Function to tabulate gene expression in clusters for heatmap keeping cluster numbers - adapted from UpsetR function fromList
fromList_withclusters <- function(input){
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x){
    x <- as.vector(match(elements, x))
  }
  ))
  data[is.na(data)] <- as.integer(0); data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) !=0), ]
  names(data) <- names(input)
  rownames(data) <- elements
  return(data)
}  


plotUMAPs <- function(seurat.object,filepath){
  pdf(file = filepath,width = 10, height = 15,onefile = TRUE) 
  DimPlot(combined.clusters, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 1,label.size = 5)
  DimPlot(combined.clusters, reduction = "umap", label = TRUE, repel = TRUE,group.by = "orig.ident",pt.size = 1,label.size = 5)
  DimPlot(combined.clusters, reduction = "umap", label = TRUE, repel = TRUE,split.by = "orig.ident",pt.size = 1,label.size = 5)
  DimPlot(combined.clusters, reduction = "umap", label = TRUE, repel = TRUE,group.by = "dataset",pt.size = 1,label.size = 5)
  DimPlot(combined.clusters, reduction = "umap", label = TRUE, repel = TRUE,group.by = "rust2020.ident",pt.size = 1,label.size = 5)
  DimPlot(combined.clusters, reduction = "umap", label = TRUE, repel = TRUE,group.by = "rust2020.ident",split.by = "orig.ident",pt.size = 1,label.size = 5)
  DimPlot(combined.clusters, reduction = "umap", label = TRUE, repel = TRUE,group.by = "jevitt2020.ident",pt.size = 1,label.size = 5)
  DimPlot(combined.clusters, reduction = "umap", label = TRUE, repel = TRUE, group.by = "jevitt2020.ident",split.by = "orig.ident",pt.size = 1,label.size = 5)
  DimPlot(combined.clusters, reduction = "umap", label = TRUE, repel = TRUE, group.by = "known_sfc_cells",pt.size = 1,label.size = 6)
  DimPlot(combined.clusters, reduction = "umap", label = TRUE, repel = TRUE, group.by = "known_sfc_cells",split.by="orig.ident",pt.size = 1,label.size = 6)
  
  
  dev.off()
}

plotMarkers <- function(seurat.object,genelists,filepath){
  pdf(file = filepath ,width = 15, height = 10,onefile = TRUE) 
  for(genelist in genelists){
    p <- DotPlot(seurat.object,features =genelist)+RotatedAxis()
    print(p)
  }
  dev.off()
}

# Plot heatmap of clusters vs. genes for SFC markers
plotSFCMarkersHeatmap <- function(FindAllMarkers.ResultsTable,logFC.threshold){
  heatmap.data <- list(`cv-2` = FindAllMarkers.ResultsTable %>% filter(gene == 'cv-2') %>% filter(avg_log2FC > logFC.threshold & p_val_adj < 0.05) %>% pull(cluster),
                           `dpp` = FindAllMarkers.ResultsTable %>% filter(gene == 'dpp') %>% filter(avg_log2FC > logFC.threshold & p_val_adj < 0.05) %>% pull(cluster),
                           `Vha100-2` = FindAllMarkers.ResultsTable %>% filter(gene == 'Vha100-2') %>% filter(avg_log2FC > logFC.threshold & p_val_adj < 0.05) %>% pull(cluster),
                           `Vha16-1` = FindAllMarkers.ResultsTable %>% filter(gene == 'Vha16-1') %>% filter(avg_log2FC > logFC.threshold & p_val_adj < 0.05) %>% pull(cluster),
                           `drpr` = FindAllMarkers.ResultsTable %>% filter(gene == 'drpr') %>% filter(avg_log2FC > logFC.threshold & p_val_adj < 0.05) %>% pull(cluster),
                           `Past1` = FindAllMarkers.ResultsTable %>% filter(gene == 'Past1') %>% filter(avg_log2FC > logFC.threshold & p_val_adj < 0.05) %>% pull(cluster),
                           `peb` = FindAllMarkers.ResultsTable %>% filter(gene == 'peb') %>% filter(avg_log2FC > logFC.threshold & p_val_adj < 0.05) %>% pull(cluster),
                           `crq` = FindAllMarkers.ResultsTable %>% filter(gene == 'crq') %>% filter(avg_log2FC > logFC.threshold & p_val_adj < 0.05) %>% pull(cluster),
                           `kay` = FindAllMarkers.ResultsTable %>% filter(gene == 'kay') %>% filter(avg_log2FC > logFC.threshold & p_val_adj < 0.05) %>% pull(cluster),
                           `trol` = FindAllMarkers.ResultsTable %>% filter(gene == 'trol') %>% filter(avg_log2FC > logFC.threshold & p_val_adj < 0.05) %>% pull(cluster),
                           `eya` = FindAllMarkers.ResultsTable %>% filter(gene == 'eya') %>% filter(avg_log2FC > logFC.threshold & p_val_adj < 0.05) %>% pull(cluster)
                       )
  legend_low <- paste0("Avg.LFC < ",logFC.threshold," | Adj.p > 0.05" )
  legend_high <- paste0("Avg.LFC > ",logFC.threshold," & Adj.p < 0.05" )
  pheatmap::pheatmap(fromList_withclusters(heatmap.data),
                   color = c("wheat","orange"),
                   legend_breaks = 0:1, legend=T,
                   legend_labels = c(legend_low ,legend_high))

}

compare_translatome_scrnaseq_foldchanges <- function(scrnaseq_markers,deseq_results, genelist,scalexy = c(-3,3), output_dir){
  
  scrnaseq_markers$gene <- rownames(scrnaseq_markers)
  #print(head(scrnaseq_markers))
  
  deseq_results <- deseq_results %>%  mutate(gene_unique = make.unique(as.character(gene))) %>% drop_na()
  rownames(deseq_results) <- deseq_results$gene_unique
  
  #print(head(deseq_results))
  
  translatome_atlas_compare <- full_join(deseq_results,scrnaseq_markers,
                                         by = 'gene', suffix = c(".translatome", ".scAtlas"))

  
  translatome_atlas_compare <- translatome_atlas_compare %>% 
    mutate("diffavg_log2FC" = abs(log2FoldChange - avg_log2FC)) %>%
    mutate("signifcolors" = ifelse(padj < 0.05 & p_val_adj < 0.05, 'black',
                                   ifelse(padj < 0.05 & p_val_adj >= 0.05, "#DDAFB0",
                                          ifelse(padj >= 0.05 & p_val_adj < 0.05, "#8F6980","#FBDCC0")))) %>%
    mutate("genelabels" = ifelse((log2FoldChange > 1 & avg_log2FC > 0.02) | #quandrant 1
                                   (log2FoldChange > 1 & avg_log2FC > 1.5 & signifcolors == 'black') | #quandrant 1
                                 (log2FoldChange < -1.5 | avg_log2FC > 2) |#quandrant 2
                                 (log2FoldChange < -2 | avg_log2FC < -2),#quandrant 3
                                 gene,""))  %>%
    mutate("gstd" = ifelse(grepl("Obp",gene),gene,'')) %>%
    mutate("quadrant" = ifelse((log2FoldChange > 0 & avg_log2FC > 0), "plum3", #quandrant 1
                               ifelse((log2FoldChange < 0 & avg_log2FC > 0), "darkslategray2", #quandrant 2
                                        ifelse((log2FoldChange < 0 & avg_log2FC < 0), "lightpink2", #quandrant 3
                                               ifelse((log2FoldChange > 0 & avg_log2FC < 0),"aquamarine3" #quandrant 4
                                 ,"")))))  %>%
    mutate("quadrantcompare" = ifelse((log2FoldChange > 0 & avg_log2FC > 0), "concordant", #quandrant 1
                               ifelse((log2FoldChange < 0 & avg_log2FC > 0), "discordant", #quandrant 2
                                      ifelse((log2FoldChange < 0 & avg_log2FC < 0), "concordant", #quandrant 3
                                             ifelse((log2FoldChange > 0 & avg_log2FC < 0),"discordant" #quandrant 4
                                                    ,""))))) 
  translatome_atlas_compare <- translatome_atlas_compare %>% drop_na()
  
  #print((translatome_atlas_compare %>% filter(log2FoldChange > 2 & avg_log2FC > 2)))
  
  
  #### Fold change directions for SFC (+) and MBFC (-) are the same in both datasets.
  translatome_atlas_compare_plot <- ggplot(translatome_atlas_compare, aes(x=log2FoldChange, y=avg_log2FC)) +
    geom_point(size=1,color=translatome_atlas_compare$signifcolors) +
    ggtitle(" Translatome vs.scAtlas") +
    xlab("TRAP-seq log2(Fed SFC/Fed MBFC)") + 
    ylab(paste0("Integrated scRNAseq log2(SFC/MBFC)")) +
    geom_text_repel(label = translatome_atlas_compare$genelabels,
                    size=5,max.overlaps = 300, fontface = "italic",
                    nudge_x = 0.2, nudge_y = -0.1,force=1,
                    #fill=translatome_atlas_compare$signifcolors, color = "white",
                    color= translatome_atlas_compare$signifcolors
                    ) + 
    geom_hline(yintercept = 0,linetype=2,color="grey2") +
    geom_vline(xintercept = 0,linetype=2,color="grey2") +
    scale_x_continuous(limits=scalexy,breaks = seq(-5, 5, by = 1)) +
    scale_y_continuous(limits=c(-3,3.5),breaks = seq(-3, 3.5, by = 1)) +
    theme_clean(base_size = 16) 
  
  
  print(translatome_atlas_compare_plot)
  
  quantify_concordance <- translatome_atlas_compare  %>%
    ggplot(aes(x=quadrantcompare,fill=signifcolors))+
    geom_bar(width=0.2)+
    scale_fill_identity()+
    theme_classic(base_size = 16) +
    theme(axis.text=element_text(size=20,color="black"),
          axis.text.x=element_text(angle=30,vjust=0.8,hjust=0.7))
  
  print(quantify_concordance)

  ggsave(translatome_atlas_compare_plot,filename = paste0(output_dir,"/sc_trap_compare.png"),width=8,height=8)
  ggsave(quantify_concordance,filename = paste0(output_dir,"/sc_trap_concordance.png"),width=5,height=7)
  return(translatome_atlas_compare)
}

