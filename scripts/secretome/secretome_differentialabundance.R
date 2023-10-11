###################################################################################################################
#######  Main script for running secretome (LC-MS/MS) data analysis- 
########### 1. Differential abundance analysis and sample QC on Maxquant intensity 
########### 2. Compare secretome profiles of starvation and phagoptosis 
########### 3. Perform GO/functional enrichment
###################################################################################################################


library("DEP")
library("tidyverse")
library(httr)
library(ggrepel)
library(org.Dm.eg.db)
library(clusterProfiler)
library(ggthemes)
library(SummarizedExperiment)
library(ComplexHeatmap)

basepath <- file.path("/projectnb/mccall/sbandyadka/mmNCC/")
inputdir <- file.path(basepath,"secretome", "data/")
outputdir <- file.path(basepath,"secretome", "analysis/")

raw_intensities <- read.table(paste0(inputdir,"proteinGroups.txt"),sep="\t",header=TRUE)
raw_intensities <- filter(raw_intensities, Reverse != "+", Potential.contaminant != "+")
raw_intensities$Gene.names %>% duplicated() %>% any()
raw_intensities %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
intensities_uniquegenes <- make_unique(raw_intensities, "Gene.names", "Protein.IDs", delim = ";")
intensities_uniquegenes$name %>% duplicated() %>% any()
LFQ_columns <- grep("Intensity.", colnames(intensities_uniquegenes)) # get LFQ column numbers


exptdesign <- read.table(file.path(basepath,"scripts","secretome", "exptdesign.txt"),sep="\t",header=TRUE,stringsAsFactors = FALSE)


runDE <- function(raw_intensities_uniqueprotein,intentisitycolumns, exptdesigntable,samples_to_include_start,samples_to_include_end,contrastname,imputationmethod, poscolor, negcolor,outputfolder){
  raw_se <- make_se(raw_intensities_uniqueprotein, intentisitycolumns, exptdesigntable)
  raw_se <- raw_se[, samples_to_include_start:samples_to_include_end]
  
  #test <- as.data.frame(assays(raw_se)[[1]])
  #test %>% filter(row.names(test) %in% "Cbp20") %>% print()
  
  plot_frequency(raw_se)
  filtered <- filter_missval(raw_se, thr = 0)
  plot_numbers(filtered)
  plot_coverage(filtered)
  
  normalized <- normalize_vsn(filtered)
  plot_normalization(filtered, normalized)
  plot_missval(filtered)
  plot_detect(filtered)
  
  imp <- impute(normalized, fun = imputationmethod)
  
  #test <- as.data.frame(assays(imp)[[1]])
  #test %>% filter(row.names(test) %in% "Cbp20") %>% print()
  
  de_imp <- test_diff(imp, type = "manual", test = contrastname)
  plot_imputation(normalized, imp)

  
  dep <- add_rejections(de_imp, alpha = 0.05)
  plot_pca(dep, x = 1, y = 2, n = 100, point_size = 4)
  plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "YlGnBu")
  
  
  #png(filename =  paste0(outputdir,contrastname,"_heatmap.png"),width=10,height=10)
  plot_heatmap(dep, type = "centered", kmeans = TRUE, 
               k = 2, col_limit = 4, show_row_names = TRUE,
               indicate = c("condition", "replicate"))
  #draw(hm)
  #dev.off()

  res <- get_results(dep)
  colnames(res) <- c("name","ID","p.val","p.adj","signif","signifcontrast","ratio","cent1","cent2")
  res$genesymbol <-mapIds(org.Dm.eg.db, 
                          keys=res$ID, 
                          column="SYMBOL", 
                          keytype="UNIPROT", 
                          multiVals="first")
  
  res <- res %>% mutate("pointcolors" = ifelse(signifcontrast == "TRUE" & ratio >0 , poscolor,
                                               ifelse(signifcontrast == "TRUE" & ratio < 0,negcolor,"grey"))) %>%
    mutate(genelabels=ifelse(signifcontrast == "TRUE" & abs(ratio)>2,genesymbol,""))
  
  
  volcanoplot <- ggplot(res, aes(x=ratio,
               y=-log10(p.adj)))+
    geom_point(aes(color=pointcolors))+
    geom_text_repel(aes(label=genelabels),size=6,max.overlaps = 3000 )+
    geom_vline(xintercept = 0,linetype="longdash")+
    scale_color_identity()+
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
  
  print(volcanoplot)
 
  ggsave(volcanoplot,filename = paste0(outputdir,contrastname,"_volcano.png"),width=10,height=10)
  
  write.table(res, file = paste0(outputdir,contrastname,".tsv"),sep="\t",row.names = FALSE)
  return(res)
}

phagoptosis_de <- runDE(intensities_uniquegenes,LFQ_columns,exptdesign,1,4,"FedSFCs_vs_FedMBFCs","knn", "#65BADA","#068E8C",outputdir)

phagoptosis_de %>% filter(genesymbol=="Plekhm1")
phagoptosis_de %>% mutate(enrich=ifelse(ratio>0, "SFC","MBFC")) %>%
  group_by(signif, enrich) %>%
  tally() 

starvation_de <- runDE(intensities_uniquegenes,LFQ_columns,exptdesign,5,8,"StarvedMBFCs_vs_FedMBFCs2","knn","#D86F27","#068E8C",outputdir)
starvation_de %>% mutate(enrich=ifelse(ratio>0, "StarvedMBFC","FedMBFC")) %>%
  group_by(signif, enrich) %>%
  tally()

## Are there any significantly DA proteins in phagoptosis and starvation common to both comparisons? 
phagoptosis_siggenes <- phagoptosis_de %>% filter(signifcontrast==TRUE) %>% pull(genesymbol) 
starvation_siggenes <- starvation_de %>% filter(signifcontrast==TRUE) %>% pull(genesymbol)
intersect(phagoptosis_siggenes,starvation_siggenes) ##None found. 

## Fold change phagoptosis vs. starvation 
phagoptosis_de <- phagoptosis_de %>% mutate(gene_peptide = paste0(genesymbol,"-",ID))
starvation_de <- starvation_de %>% mutate(gene_peptide = paste0(genesymbol,"-",ID))

starvation_vs_phagoptosis <- full_join(phagoptosis_de,starvation_de,
                                       by = 'gene_peptide', suffix = c(".phagoptosis", ".starvation"))

starvation_vs_phagoptosis <- starvation_vs_phagoptosis %>% 
  #mutate(peptide = coalesce(ID.phagoptosis,ID.starvation)) %>%
  #mutate(symbol = paste0(genesymbol,"-",peptide)) %>%
  select(ratio.phagoptosis,p.adj.phagoptosis,ratio.starvation, p.adj.starvation, gene_peptide) %>%
  mutate(ratio.phagoptosis = ifelse(is.na(ratio.phagoptosis), 0, ratio.phagoptosis),
         ratio.starvation = ifelse(is.na(ratio.starvation), 0, ratio.starvation))


starvation_vs_phagoptosis_annotated <-starvation_vs_phagoptosis %>% 
  mutate(signifcolors = ifelse(p.adj.phagoptosis < 0.05 & p.adj.starvation < 0.05, "black", 
                               ifelse(p.adj.phagoptosis < 0.05 & p.adj.starvation >= 0.05, "#65BADA",
                                      ifelse(p.adj.phagoptosis >= 0.05 & p.adj.starvation < 0.05, "#D86F27",
                                        ifelse(p.adj.phagoptosis >= 0.05 & p.adj.starvation >= 0.05, "grey","grey"))))) %>%
  mutate(signifcolors2 = ifelse(p.adj.phagoptosis < 0.05 & is.na(p.adj.starvation), "#65BADA","grey")) %>%
  mutate(signifcolors3 = ifelse(p.adj.starvation < 0.05 & is.na(p.adj.phagoptosis), "#D86F27","grey"))

starvation_vs_phagoptosis_annotated <- starvation_vs_phagoptosis_annotated %>%
  mutate(colors = coalesce(signifcolors,signifcolors2)) %>%
  mutate(allcolors = coalesce(colors,signifcolors3)) %>%
  mutate(labels = ifelse(allcolors =="#65BADA" | allcolors== "#D86F27", gene_peptide,"" )) 

#starvation_vs_phagoptosis_annotated <- starvation_vs_phagoptosis_annotated %>% 
#  distinct(labels, .keep_all = TRUE)

head(starvation_vs_phagoptosis_annotated)

starvation_vs_phagoptosis_plot <- ggplot(starvation_vs_phagoptosis_annotated, aes(x=ratio.phagoptosis, y=ratio.starvation,color=allcolors)) +
  geom_point(size=1) +
  scale_color_identity()+
  xlab("Log2 Fold Change (Fed SFCs / Fed MBFCs)") + 
  ylab("Log2 Fold Change (Starved MBFC / Fed MBFCs)") +
  labs(fill='Category') +
  geom_text_repel(label = starvation_vs_phagoptosis_annotated$labels,
                  size=6, max.overlaps = 30, force = 5)+
  geom_hline(yintercept = 0,linetype=2) +
  geom_vline(xintercept = 0,linetype=2) +
  scale_x_continuous(limits=c(-7,7)) +
  scale_y_continuous(limits=c(-7,7)) +
  theme(text = element_text(family="Arial"),
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        legend.position = "none") +
  theme_classic(base_size = 16) 
                  
ggsave(starvation_vs_phagoptosis_plot, filename = paste0(outputdir,"/starvation_vs_phagoptosis_scatter.png"),width = 11,height=10)


## GO-term enrichment


phagoptosis_enrichment <- phagoptosis_de %>% filter(signifcontrast==TRUE) %>% 
  mutate(enrichment = ifelse(ratio>0, "enriched in fed SFCs","enriched in fed MBFCs")) %>% 
  select(genesymbol,enrichment,ratio)


starvation_enrichment <- starvation_de %>% filter(signifcontrast==TRUE) %>% 
  mutate(enrichment = ifelse(ratio>0, "enriched in starved MBFCs","enriched in fed MBFCs")) %>% 
  select(genesymbol,enrichment,ratio)

comparegoenrichment <- rbind(phagoptosis_enrichment,starvation_enrichment)
## Functional enrichment

goenrich <- function(deresult,resultsdir){
  deresult$entrez <- mapIds(org.Dm.eg.db, 
                                  keys=deresult$genesymbol,
                                  column="ENTREZID", 
                                  keytype="SYMBOL",
                                  multiVals="first")

  formula_res_bp <- compareCluster(entrez~enrichment, data=deresult, fun="enrichGO",OrgDb='org.Dm.eg.db',ont="BP")
  formula_res_mf <- compareCluster(entrez~enrichment, data=deresult, fun="enrichGO",OrgDb='org.Dm.eg.db',ont="MF")
  formula_res_cc <- compareCluster(entrez~enrichment, data=deresult, fun="enrichGO",OrgDb='org.Dm.eg.db',ont="CC")
  theme_text_size = 10
  
  
  simplifiedgo <- simplify(formula_res_mf, cutoff=1, by="p.adjust")
  
  #pdf(file.path(resultsdir,paste0(contrast_name,"enrichGO_BP",".pdf")), width = 15, height = 10)
  go_dotplot<- dotplot(simplifiedgo,  x="enrichment") +
    scale_color_continuous()+
    #scale_color_continuous(high="#f2b083",low="#6d2f04") +
    scale_y_discrete(guide = guide_axis(check.overlap = FALSE)) +
    scale_x_discrete(guide = guide_axis(check.overlap = FALSE)) +
    theme_classic(base_size = 15) +
    theme(axis.text = element_text(size= 10 ),
          panel.background = element_rect(fill = NA),
          panel.grid.major = element_line(colour = "#eeeeee",linetype = "dotted"),
          panel.ontop = FALSE)
  
  print(go_dotplot)
  ggsave(go_dotplot,filename=file.path(resultsdir,paste0("gotermenrich",".png")),width = 10, height = 5)
  #dev.off()
  
  cnet_genes <- setReadable(simplifiedgo, 'org.Dm.eg.db', 'ENTREZID')
  go_cnet <- cnetplot(cnet_genes, node_label="all",
           cex_label_gene = 1.5,circular=FALSE,
           cex_label_category=1.5) 
  ggsave(go_cnet,filename=file.path(resultsdir,paste0("gotermenrich_cnet",".png")),width = 7, height = 7)

}

goenrich(comparegoenrichment,outputdir)



