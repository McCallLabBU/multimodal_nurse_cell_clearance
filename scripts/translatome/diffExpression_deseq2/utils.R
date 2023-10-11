###################################################################################################################
####### Helper functions for translatome Differential expression analysis and visualization
# Author: Shruthi Bandyadka, McCall Lab, Boston University (sbandya@bu.edu)
###################################################################################################################

get_contrast_result <- function(contrast_name, deseq_data, type,convertid_type,outdir,color1,color2,lfc_cutoff){
  if(type == 'lfc'){
    contrast_result <- lfcShrink(deseq_data, coef = contrast_name,type='apeglm')
    
  }else{
    contrast_result <- results(deseq_data, name = contrast_name)
  }
  
  deseq_result <- contrast_result
  contrast_result <- as.data.frame(contrast_result)
  
  ## Make volcano plots of each contrast
  contrast_result$transcript <- rownames(contrast_result)
  contrast_result['gene'] <- mapIds(org.Dm.eg.db,
                                    keys=contrast_result$transcript,
                                    column = 'SYMBOL',
                                    keytype = convertid_type,
                                    multiVals="first")
  
  
  contrast_result <- contrast_result %>% arrange(padj)
  
  contrast_result <- contrast_result %>% mutate(log10p = -log10(padj)) %>%
    mutate(gene1 = ifelse(transcript=="FBgn0039307","CR13656",  ##update labels for pseudogenes not automatically added
                          ifelse(transcript=="FBgn0053487","CR33487",gene))) %>% 
    mutate(genelabels = ifelse(padj < 0.05 & abs(log2FoldChange) > lfc_cutoff,  gene1,
                               ifelse(log10p > 18, gene1, "")))
  #print(head(contrast_result,40))

  #contrast_result$genelabels[1:30] <- as.character(contrast_result$gene[1:30])
  
  contrast_result <- contrast_result %>% mutate(transcriptlabels = "")
  contrast_result$transcriptlabels[1:30] <- as.character(contrast_result$transcript[1:30])
  
  plot_title <- paste0(contrast_name,"_",type)
  contrast_filepath <- file.path(outdir,paste0(plot_title,".tsv"))
  write.table(contrast_result,file=contrast_filepath,sep="\t",row.names = FALSE)
  
  volcplot_gene <- plotvolcano(contrast_result,plot_title,"gene",color1,color2) ## Saves plot with gene symbol / names
  volcplot_path_gene <- file.path(outdir,paste0(plot_title,"_gene.png"))
  ggsave(volcplot_gene,filename = volcplot_path_gene,width = 10, height = 10)
  
  
  compare_clusters_Enrichgo(deseq_result,contrast_name,outdir)
  
  return(deseq_result)
}


plotvolcano <- function(results_table,title,category, poscolor, negcolor){
  results_table <- results_table %>% 
    
    mutate("pointcolors" = ifelse(padj < 0.05 & log2FoldChange < 0, negcolor,
                                   ifelse(padj < 0.05 & log2FoldChange > 0, poscolor,
                                          ifelse(padj > 0.05 ,"grey","grey"))))  
    
  if (category == 'gene'){
    volcplot <- ggplot(results_table,aes(x=log2FoldChange,y=-log10(padj))) +
      geom_point(aes(colour = pointcolors), size=2)+ 
      scale_x_continuous(breaks = seq(-12, 12, by = 2))+
      scale_colour_manual(values = c(negcolor,poscolor,"grey"))+
      geom_text_repel(aes(label = genelabels),size = 6,
                      max.overlaps = 3000)+
      geom_vline(xintercept = 0,color="black",linetype="longdash") +
      ggtitle(title) +
      xlab("Log2 Fold Change") + 
      ylab("-Log10 Adjusted p-value") +
      theme_clean()+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            plot.title = element_text(size = 14),
            axis.title=element_text(size= 20),
            legend.position="none",
            axis.line.x.bottom=element_line(size=1),
            axis.line.y.left=element_line(size=1),
            axis.text = element_text(size= 20 )) 
  }else{
    volcplot <- ggplot(results_table,aes(x=log2FoldChange,y=-log10(padj))) +
      geom_point(aes(colour = pointcolors), size=0.5)+ 
      scale_colour_manual(values = c(negcolor,poscolor,"grey"))+
      geom_text_repel(aes(label = transcriptlabels),size=2)+
      ggtitle(title) +
      xlab("log2 fold change") + 
      ylab("-log10 adjusted p-value") +
      theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black"),
            plot.title = element_text(size = 10),
            axis.title=element_text(size=10)) 
  }
  return(volcplot)
}


compare_clusters_Enrichgo <- function(deseq_resultstable,contrast_name,resultsdir){

  deseq_resultstable$entrez <- mapIds(org.Dm.eg.db, 
                                    keys=rownames(deseq_resultstable), 
                                    column="ENTREZID", 
                                    keytype="ENSEMBL",
                                    multiVals="first")
  
  
  enrichment_df <- deseq_resultstable %>% as.data.frame() %>% filter(padj < 0.05) 
  mydf <- data.frame(Entrez=enrichment_df$entrez, FC=enrichment_df$log2FoldChange)
  mydf <- mydf[abs(mydf$FC) > 1,]
  mydf$group <- paste0(contrast_name," Numerator")
  mydf$group[mydf$FC < 0] <- paste0(contrast_name," Denominator")
  mydf$othergroup <- "< 2x LogFoldChange"
  mydf$othergroup[abs(mydf$FC) > 2] <- "> 2x LogFoldChange"
  
  mydf <- mydf %>% drop_na() 
  
  formula_res_bp <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichGO",OrgDb='org.Dm.eg.db',ont="BP")
  formula_res_mf <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichGO",OrgDb='org.Dm.eg.db',ont="MF")
  formula_res_cc <- compareCluster(Entrez~group+othergroup, data=mydf, fun="enrichGO",OrgDb='org.Dm.eg.db',ont="CC")
  
  theme_text_size = 16
  
  pdf(file.path(resultsdir,paste0(contrast_name,"enrichGO_BP",".pdf")), width = 15, height = 10)
  go_bp <- dotplot(formula_res_bp,  x="group", title = "GO: Biological Process") +
    facet_grid(~othergroup) +
    scale_color_continuous(high="#f2b083",low="#6d2f04") +
    scale_y_discrete(guide = guide_axis(check.overlap = TRUE)) +
    scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(size= 16 ),
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(colour = "#eeeeee",linetype = "dotted"),
          panel.ontop = FALSE)
  
  print(go_bp)
  ggsave(go_bp,filename=file.path(resultsdir,paste0(contrast_name,"enrichGO_BP",".png")),width = 15, height = 10)
  dev.off()
  
  pdf(file.path(resultsdir,paste0(contrast_name,"enrichGO_MF",".pdf")), width = 15, height = 10)
  go_mf <- dotplot(formula_res_mf,  x="group", title = "GO: Molecular Function") + 
    facet_grid(~othergroup) +theme_classic(base_size = theme_text_size) + 
    scale_color_continuous() +
    scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
    scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
    theme_clean(base_size = 16) +
    theme(axis.text = element_text(size= 16 ))
  print(go_mf)
  ggsave(go_mf,filename=file.path(resultsdir,paste0(contrast_name,"enrichGO_MF",".png")),width = 15, height = 10)
  dev.off()
  
  pdf(file.path(resultsdir,paste0(contrast_name,"enrichGO_CC",".pdf")), width = 15, height = 10)
  go_cc <- dotplot(formula_res_cc,  x="group", title = "GO: Cellular Component") + 
    facet_grid(~othergroup) +theme_classic(base_size = theme_text_size) + 
    scale_color_continuous()  +
    scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
    scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
    theme_clean(base_size = 16) +
    theme(axis.text = element_text(size= 16 ))
  print(go_cc)
  ggsave(go_cc,filename=file.path(resultsdir,paste0(contrast_name,"enrichGO_CC",".png")),width = 15, height = 10)
  dev.off()
  
  deseq_resultstable <- as.data.frame(deseq_resultstable)
  significant_genes <- deseq_resultstable %>% filter(padj < 0.05 & abs(log2FoldChange) >1 )
  ego <- enrichGO(gene          = significant_genes$entrez,
                  universe      = deseq_resultstable$entrez,
                  keyType       = 'ENTREZID',
                  OrgDb         = org.Dm.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  
  ego_table_filename <- file.path(resultsdir,paste0(contrast_name,"_EnrichGO.tsv"))
  write.table(ego,ego_table_filename,sep="\t",row.names = FALSE)
  geneList <- significant_genes$log2FoldChange
  names(geneList) <- significant_genes$entrez 
  geneList <- sort(geneList, decreasing = TRUE)
  
  enrich_heatplot <- heatplot(ego, foldChange=geneList) + 
    scale_fill_continuous() + 
    theme(axis.text = element_text(size=14)) + coord_equal() +
    scale_y_discrete(guide = guide_axis(check.overlap = TRUE))
  
  pdf(file.path(resultsdir,paste0(contrast_name,"enrichGO_heatplot",".pdf")), width = 20, height = 10)
  print(enrich_heatplot)
  ggsave(enrich_heatplot, filename = paste0(resultsdir,"/", contrast_name,"_enrichGOheatplot.svg"),width=20,height=10)
  
  dev.off()
}


plot_counts <- function(categorized_dds,groupname){

  fbgn_genelist <- categorized_dds %>% filter(category == groupname)  %>% drop_na() %>% pull(geneid) 
  count_tables <- list()
  
  for(i in fbgn_genelist){
    #print(i)
    try({
      genecounts <- plotCounts(dds, gene=i, intgroup="genotype", returnData=TRUE)
      genecounts <- mutate(genecounts, group=ifelse(grepl('^GC',rownames(genecounts)), "MBFC-fed",
                                                    ifelse(grepl('^GS',rownames(genecounts)), 'MBFC-starved', "SFC-fed")))
      meangc <- genecounts %>% group_by(group) %>%
        summarise_at(vars(count), list(name = mean))
      meangc <- mutate(meangc, colorpoints=ifelse(grepl('^MBFC-fed',group), "#068E8C",
                                                  ifelse(grepl('^MBFC-starved',group), '#D86F27', "#65BADA")))
      meangc <- mutate(meangc,gene=i)
      #print(meangc)
      genesymbol <- mapIds(org.Dm.eg.db, keys=i, column="SYMBOL",keytype="ENSEMBL", multiVals="first")
      meangc <- mutate(meangc,symbol = if_else(grepl('^SFC-fed',group), genesymbol,""))
      
      count_tables[[i]] <- meangc

    })
  }
  
  count_tables_df <- Reduce(full_join,count_tables)

  gp_withannot <- ggplot(count_tables_df, aes(x = group, y = name, color=colorpoints,group=gene)) +   scale_y_log10() + 
    ggtitle(groupname) + geom_beeswarm(size = 3) + 
    geom_line(size=0.5, color="black",linetype="dashed")+
    geom_label_repel(aes(label = symbol), size = 2) +
    scale_colour_identity() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  
  gp_withoutannot <- ggplot(count_tables_df, aes(x = group, y = name, color=colorpoints,group=gene)) +   scale_y_log10() + 
    ggtitle(groupname) + geom_beeswarm(size = 3) + 
    #geom_line(size=0.5, color="black",linetype="dashed")+
    scale_colour_identity() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
    
    return(list(gp_withannot,gp_withoutannot))

  }


visualize_data <-function(deseq_data,transformation,transformed_data,level, outdir){
  dipersion_plotname <- paste0(level,"_dispersion_plot.pdf")
  pdf(file.path(outdir,dipersion_plotname))
  plotDispEsts(deseq_data)
  dev.off()
  
  meansd_plotname <- paste0(level,"_",transformation,"_meansdplot.pdf")
  pdf(file.path(outdir,meansd_plotname))
  meanSdPlot(assay(transformed_data))
  dev.off()
  
  top20_counts_heatmap_plotname <- paste0(level,"_",transformation,"_top10counts_heatmap.pdf")
  select <- order(rowMeans(counts(deseq_data,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  top20heatmap <- ggplot(melt(assay(transformed_data)[select,]), aes(Var2,Var1, fill=value)) + 
    geom_tile() +
    theme(axis.text.x=element_text(angle=90,hjust=1)) + 
    scale_fill_distiller(palette = "Blues", direction = 1) +
    xlab("")+
    ylab("")
  ggsave(file.path(outdir,top20_counts_heatmap_plotname),top20heatmap)
  
  
  sampledist_plotname <- paste0(level,"_",transformation,"_sampledist.pdf")
  pdf(file.path(outdir,sampledist_plotname))
  sampleDists <- dist(t(assay(transformed_data)))
  
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(transformed_data$condition, sep="-")
  colnames(sampleDistMatrix) <- paste(transformed_data$condition, sep="-")
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  dev.off()
  
  pca_plotname <- paste0(level,"_",transformation,"_pca.pdf")
  pca_title <- paste0("PCA - ",transformation," - ",level)
  pcaData <- plotPCA(transformed_data, intgroup="genotype", returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  pca <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=1) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() + ggtitle(pca_title) + geom_text(aes(label = name),size=3,nudge_x = 2,nudge_y = 1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  ggsave(file.path(outdir,pca_plotname),pca)
  
}

contrast_heatmap <- function(deseq_data, deseq_lfc_result,contrast_name,outdir,foldchage_cutoff){
  
  #deseq_lfc_result <- as.data.frame(deseq_lfc_result)
  sigGenesdf <- deseq_lfc_result %>% drop_na() %>% filter(padj < 0.05 & abs(log2FoldChange) >= foldchage_cutoff) 
  sigGenes <- rownames(sigGenesdf) 
  
  ## VST Normalize counts
  vsd <- vst(deseq_data, blind=FALSE)
  vsd_raw <- as.data.frame(assay(vsd))
  vsd_raw['gene_symbol'] <- mapIds(org.Dm.eg.db,
                                   keys=rownames(vsd_raw),
                                   column = 'SYMBOL',
                                   keytype = 'FLYBASE',
                                   multiVals="first")
  
  ## Check if all FBgns have corresponding gene symbol mappings
  vsd_raw$fbgn <- rownames(vsd_raw) 
  no_names <- vsd_raw[is.na(vsd_raw$gene_symbol),] 
  
  vsd_df <- vsd_raw[rownames(vsd_raw) %in% sigGenes,] # Subset to signif DE
  #dim(vsd_df)
  
  sig_no_names <- vsd_df[is.na(vsd_df$gene_symbol),] #check: all sig DEG have gene name annotations
  
  rownames(vsd_df) <-  make.names(vsd_df$gene_symbol, unique=TRUE) # Set gene symbols as row labels in final heatmap
  
  keep_samples <- c()
  sample_colors <- c()
  sample_labels <- c()
  
  if(grepl("GR1_Conditioned_FC",contrast_name)){
    keep_samples <- append(keep_samples,c("GC218","GC220","GC226"))
    sample_colors <- append(sample_colors,'#FFC857')
    sample_labels <- append(sample_labels,"GR1-Conditioned")
  }
  if(grepl("GR1_starved_FC",contrast_name)){
    keep_samples <- append(keep_samples,c("GS209","GS211","GS217-225"))
    sample_colors <- append(sample_colors,'#235FA4')
    sample_labels <- append(sample_labels,"GR1-Starved")
  }
  if(grepl("PG150_conditioned_SFC",contrast_name)){
    keep_samples <- append(keep_samples,c("P203","P205","P207"))
    sample_colors <- append(sample_colors,'#0A284B')
    sample_labels <- append(sample_labels,"PG150-Conditioned")
  }
  
  vsd_df_remove_ids <- vsd_df %>% select (-c(fbgn, gene_symbol)) #select only sample counts columns for heatmap
  vsd_df_columns <- vsd_df_remove_ids %>% select (keep_samples) #subset to samples in contrast
  
  heat <- t(scale(t(vsd_df_columns)))
  myCol <- colorRampPalette(c('#FFE6B3','#FFC857','#235FA4', '#0A284B'))(100)
  myBreaks <- seq(-3, 3, length.out = 100)
  colsplit = rep(1:2, each = 3)
  
  ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill=sample_colors,
                                                    col=sample_colors), 
                                          labels_gp = gpar(col = "white", fontsize = 8),
                                          labels = sample_labels))
  
  heatmap_title <- paste0(contrast_name,"_heatmap_p005_lfc",foldchage_cutoff,".pdf")
  heatmap_filepath <- file.path(outdir,heatmap_title)
  
  pdf(heatmap_filepath, width = 8, height = 8)

  hmap <- Heatmap(heat,
                            col = colorRamp2(myBreaks, myCol),
                            row_names_gp = gpar(fontsize = 5),
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
                            top_annotation = ha,
  )
  
  
  draw(hmap,
       heatmap_legend_side = 'right',
       annotation_legend_side = 'right',
       row_sub_title_side = 'left',
       padding = unit(c(10, 10, 10, 10), "mm")
  )
  
  dev.off()
  
}
 
  
compare_degs_prep_dds <- function(deseq_data){
  ## Get VST Normalized read counts and their averages across replicates
  vsd <- vst(deseq_data, blind=FALSE)
  vsd_df <- as.data.frame(assay(vsd))
  vsd_df$gene <- rownames(vsd_df)
  vsd_df['gene_name'] <- mapIds(org.Dm.eg.db,
                                keys=vsd_df$gene,
                                column = 'SYMBOL',
                                keytype = 'ENSEMBL',
                                multiVals="first")
  
  phagoptosis_allgenes <- lfcShrink(dds, coef="genotype_PG150_vs_GR1", type="apeglm")
  starvation_allgenes <- lfcShrink(dds, coef="treatment_starved_vs_fed", type="apeglm")
  phagoptosis_allgenes <- as.data.frame(phagoptosis_allgenes)
  starvation_allgenes <- as.data.frame(starvation_allgenes)
  
  phagoptosis_allgenes['gene'] <- mapIds(org.Dm.eg.db,
                                   keys=rownames(phagoptosis_allgenes),
                                   column = 'SYMBOL',
                                   keytype = 'ENSEMBL',
                                   multiVals="first")
  phagoptosis_allgenes['entrez'] <- mapIds(org.Dm.eg.db,
                                     keys=rownames(phagoptosis_allgenes),
                                     column = 'ENTREZID',
                                     keytype = 'ENSEMBL',
                                     multiVals="first")
  
  
  starvation_allgenes['gene'] <- mapIds(org.Dm.eg.db,
                                 keys=rownames(starvation_allgenes),
                                 column = 'SYMBOL',
                                 keytype = 'ENSEMBL',
                                 multiVals="first")
  
  starvation_allgenes['entrez'] <- mapIds(org.Dm.eg.db,
                                   keys=rownames(starvation_allgenes),
                                   column = 'ENTREZID',
                                   keytype = 'ENSEMBL',
                                   multiVals="first")
  
  ## Replace NAs in gene names that weren't populated with FB ids
  
  #phagoptosis_allgenes <-phagoptosis_allgenes %>% mutate(gene = coalesce(gene,rownames(phagoptosis_allgenes)))
  #starvation_allgenes <-starvation_allgenes %>% mutate(gene = coalesce(gene,rownames(starvation_allgenes)))
  
  ## Check if there are any genes without identifiers
  phagoptosis_allgenes[is.na(phagoptosis_allgenes$gene),]
  starvation_allgenes[is.na(starvation_allgenes$gene),]
  
  ## Combine for further analysis
  starvation_vs_phagoptosis <- full_join(phagoptosis_allgenes,starvation_allgenes,
                                         by = 'gene', suffix = c(".phagoptosis", ".starvation"))
  
  ## Annotate genes and points for plotting 
  starvation_vs_phagoptosis <- starvation_vs_phagoptosis %>% 
    mutate("diffavg_log2FC" = abs(log2FoldChange.phagoptosis - log2FoldChange.starvation)) %>%
    mutate("signifcolors" = ifelse(padj.phagoptosis < 0.05 & padj.starvation < 0.05, "black",
                                   ifelse(padj.phagoptosis < 0.05 & (padj.starvation > 0.05 | is.na(padj.starvation)), "#65BADA",
                                          ifelse((padj.phagoptosis > 0.05 | is.na(padj.phagoptosis)) & padj.starvation < 0.05, "#D86F27","grey")))) %>%
    mutate("genelabels" = ifelse((log2FoldChange.phagoptosis > 0.5 & log2FoldChange.starvation < -0.5) | 
                                   (log2FoldChange.phagoptosis < -0.5 & log2FoldChange.starvation > 0.5) |
                                   (diffavg_log2FC > 5) |
                                   (diffavg_log2FC > 1.5 & signifcolors == 'black' ) |
                                   (log2FoldChange.phagoptosis > 0 & signifcolors == 'black' ) |
                                   (gene=='AttD') |
                                   (signifcolors=="#503174"), gene,"")) 
  
  starvation_vs_phagoptosis <- starvation_vs_phagoptosis %>% 
    mutate(category = 
             ifelse(padj.phagoptosis  >= 0.05 & padj.starvation >= 0.05, "Adj p >= 0.05 in both",
                    ifelse(padj.phagoptosis < 0.05 & log2FoldChange.phagoptosis > 0 & (log2FoldChange.starvation < 0 | padj.starvation >= 0.05), "up in phagoptosis only",
                           ifelse(padj.phagoptosis < 0.05 & log2FoldChange.phagoptosis < 0 & (log2FoldChange.starvation > 0 | padj.starvation >= 0.05), "down in phagoptosis only",
                                  ifelse( padj.starvation < 0.05 & log2FoldChange.starvation > 0 & (log2FoldChange.phagoptosis < 0 | padj.phagoptosis >= 0.05), "up in starvation-death only",
                                          ifelse( padj.starvation < 0.05 & log2FoldChange.starvation < 0 & (log2FoldChange.phagoptosis > 0 | padj.phagoptosis >= 0.05) , "down in starvation-death only",
                                                  ifelse( padj.starvation < 0.05 & padj.phagoptosis < 0.05 & log2FoldChange.starvation > 0 & log2FoldChange.phagoptosis > 0, "congruently up",
                                                          ifelse( padj.starvation < 0.05 & padj.phagoptosis < 0.05 & log2FoldChange.starvation < 0 & log2FoldChange.phagoptosis < 0, "congruently down","other"
                                                                  ))))))))
  return(starvation_vs_phagoptosis)
           
  
}

plot_comparison_degs <- function(compare_degs_table){
  
  x.axis.labels <- seq(-8,8,2) # positions of the subtle ticks
  y.axis.labels <- seq(-4,4,2) # positions of the subtle ticks
  starvation_vs_phagoptosis_plot <- ggplot(compare_degs_table, aes(x=log2FoldChange.phagoptosis, y=log2FoldChange.starvation,color=factor(category))) +
    geom_point(size=1) +
    ggtitle(" starvation_vs_phagoptosis") +
    xlab("Log2 Fold Change (Fed SFCs / Fed MBFCs)") + 
    ylab("Log2 Fold Change (Starved MBFC / Fed MBFCs)") +
    labs(fill='Category') +
    geom_text_repel(label = compare_degs_table$genelabels,fontface="italic",
                    size=6,nudge_x = .2,
                    box.padding = 0.5,
                    nudge_y = 0.2,
                    segment.curvature = -0.1,
                    segment.ncp = 3
                    #segment.angle = 20,
                   
                   ) + 
    scale_color_manual(values=c("grey","#9270b0","#503174","#65BADA","#D86F27","#12759a","#873803"))+
    geom_hline(yintercept = 0,linetype=2) +
    geom_vline(xintercept = 0,linetype=2) +
    scale_x_continuous(limits=c(-12,6)) +
    scale_y_continuous(limits=c(-6,3)) +
    theme(text = element_text(family="Arial"),
          axis.text.x=element_blank(), 
          axis.ticks.x=element_blank(),
          legend.position = "none") +
    theme_classic(base_size = 16) 
  
  return(starvation_vs_phagoptosis_plot)
}

compare_go_annotation <- function(congruency_lists, resultsdir){
  congruency_lists$Entrez = congruency_lists$entrez.starvation
  go_categories <- c("CC","MF","BP")
  pdf(file.path(resultsdir,"starvation_vs_phagoptosis_goterms.pdf"), width = 15, height = 10,onefile=TRUE)
  
  for(cat in go_categories){
    compare_go_enrich <- compareCluster(Entrez~category, data=congruency_lists, fun="enrichGO",OrgDb='org.Dm.eg.db', ont = cat)
    pairwise <- pairwise_termsim(compare_go_enrich)
    go_filtered <- simplify(compare_go_enrich, cutoff=0.7, by="p.adjust", select_fun=min)
    
    dp <- dotplot(go_filtered) +theme_classic(base_size = 16) + scale_color_continuous() +
      scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
      scale_x_discrete(guide = guide_axis(check.overlap = TRUE))
    print(dp)
    
    cnet_genes <- setReadable(go_filtered, 'org.Dm.eg.db', 'ENTREZID')
    cnp <- cnetplot(cnet_genes,circular=FALSE, colorEdge = TRUE)
    print(cnp)
    
    pairwise <- pairwise_termsim(go_filtered)
    tp <- treeplot(pairwise,geneClusterPanel="pie",label_format = 15,
                   offset=2,nWords=0,fontsize=2,nCluster=7, legend_n =5, 
                   group_color = rep("#333333",7))
    
    print(tp)
    emp <- emapplot(pairwise, layout="kk")
    print(emp)
    
    
  }
  dev.off()
   
  return(pairwise)
  
}


get_genelists <- function(dir){

  genelists <- list(
    `cytoskeletal` = read.table(paste0(dir,"/GLAD_cytoskeletal.csv"),sep=",",header=TRUE) %>% pull(Gene),
    `chaperoneheatshock` = read.table(paste0(dir,"/GLAD_chaperoneheatshock.csv"),sep=",",header=TRUE) %>% pull(Gene),
    `ionchannel` = read.table(paste0(dir,"/GLAD_ionchannel.csv"),sep=",",header=TRUE) %>% pull(Gene),
    `kinases` = read.table(paste0(dir,"/GLAD_kinases.csv"),sep=",",header=TRUE) %>% pull(Gene),
    `receptors` = read.table(paste0(dir,"/GLAD_receptors.csv"),sep=",",header=TRUE) %>% pull(Gene),
    `secreted` = read.table(paste0(dir,"/GLAD_secreted.csv"),sep=",",header=TRUE) %>% pull(Gene),
    `transmembrane` = read.table(paste0(dir,"/GLAD_transmembrane.csv"),sep=",",header=TRUE) %>% pull(Gene),
    `transporters` = read.table(paste0(dir,"/GLAD_transporters.csv"),sep=",",header=TRUE) %>% pull(Gene),
    `signalingpathways` = read.table(paste0(dir,"/GLAD_signalingpathways.csv"),sep=",",header=TRUE) %>% pull(Gene),
    `autophagy` = read.table(paste0(dir,"/GLAD_autophagy.csv"),sep=",",header=TRUE) %>% pull(Gene),
    `matrisome` = read.table(paste0(dir,"/GLAD_matrisome.csv"),sep=",",header=TRUE) %>% pull(Gene)
  )
  return(genelists)
  
  }


