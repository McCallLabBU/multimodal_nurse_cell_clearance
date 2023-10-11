#############################################################################################################################

### Integrate single cell atlases of Drosophila melanogaster ovaries from Rust et al 2020 and Jevitt et al 2020

############################################################################################################################


library(Seurat)
library(tidyverse)
library("glmGamPoi")
library(grDevices)
library(RColorBrewer)
library(tradeSeq)
library("ggsci")
library("ggbeeswarm")
library("ggthemes")
library(viridis)
library(scater)
library(gam)


projectdir <- file.path("/projectnb/mccall/sbandyadka/mmNCC")
source(file.path(projectdir, "scripts", "scRNAseq_integration","integration_utils.R")) 
input_dir <- file.path(projectdir,  "scrnaseq","data/")
output_dir <- file.path(projectdir,  "scrnaseq","analysis/")

# Read in raw data from both Rust2020 and Jevitt2020
jevitt2020_path <- file.path(input_dir,"jevitt2020","GSE146040")
jevitt2020_raw <- Read10X(data.dir = jevitt2020_path)
jevitt2020_raw <- CreateSeuratObject(counts = jevitt2020_raw, 
                           project = "jevitt2020")
                           #min.cells = 3, 
                           #min.features = 200)

jevitt_clusteridents <- readRDS(file=file.path(input_dir, "jevitt2020", "jevitt2020_clusterIdents.RDS"))
jevitt2020_raw$jevitt2020.ident <- jevitt_clusteridents ## Keep the cluster names from Jevitt2020 paper


# Read in raw data from both Rust2020 and Jevitt2020
rust2020_d1_raw <-  Read10X(file.path(input_dir,"rust2020","GSE136162","GSM4615209"))
rust2020_d2_raw <-  Read10X(file.path(input_dir,"rust2020","GSE136162","GSM4615210"))
rust2020_d3_raw <-  Read10X(file.path(input_dir,"rust2020","GSE136162","GSM4615211"))

rust2020_d1_raw <- CreateSeuratObject(counts = rust2020_d1_raw, 
                                     project = "rust2020_dataset1")
rust2020_d2_raw <- CreateSeuratObject(counts = rust2020_d2_raw, 
                                      project = "rust2020_dataset1")
rust2020_d3_raw <- CreateSeuratObject(counts = rust2020_d3_raw, 
                                      project = "rust2020_dataset1")

rust2020_path <- file.path(input_dir,"rust2020","Drosophila_ovary_atlas_Rust_2020.rds")
rust2020_raw <- readRDS(rust2020_path)
DefaultAssay(rust2020_raw) <- "RNA"
rust2020_raw[['integrated']] <- NULL 
rust2020_raw$orig.ident <- "Rust2020"
rust2020_raw$rust2020.ident <- Idents(rust2020_raw) ## Keep the cluster names from Rust2020 paper
Idents(rust2020_raw) <- "Rust2020"

# Pre-integration QC
## Check features overlapping or missing between the 2 datasets. 
jevitt2020_features <- rownames(jevitt2020_raw)
rust2020_features <- rownames(rust2020_raw)
diffgenes <- setdiff(jevitt2020_features,rust2020_features) ##3811 ## Most are regulatory or pseudogenes(CR)
diffgenes <- setdiff(rust2020_features,jevitt2020_features) ##2 ## RFP and eGFP
length(diffgenes) 

## subset using original paper QC parameters 

jevitt2020_raw <- PercentageFeatureSet(jevitt2020_raw, pattern = "-m", col.name = "percent.mt")
jevitt2020_raw_subset <- subset(jevitt2020_raw, 
                   subset = nFeature_RNA > 775 
                   & nFeature_RNA < 55000 
                   & nCount_RNA < 30000 
                   & percent.mt < 0.05)

rust2020_raw <- PercentageFeatureSet(rust2020_raw, pattern = "-m", col.name = "percent.mt")
rust2020_raw_subset <- subset(rust2020_raw, subset = nFeature_RNA > 1000 
                & nFeature_RNA < 3000 
                & nCount_RNA > 1000 
                &  nCount_RNA < 40000
                )


data_list <- c(jevitt2020_raw_subset,rust2020_raw_subset)
#data_list <- c(jevitt2020_raw,rust2020_d1_raw, rust2020_d2_raw, rust2020_d3_raw)


data_list <- lapply(X = data_list, FUN = function(x) {
  #print(VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
  x <- SCTransform(x,method = "glmGamPoi",min_cells=1)
  
})

## Check the number of cells expressing some key genes of interest
sum(GetAssayData(object = jevitt2020_raw_subset, slot = "counts")["Hml",]>0)
sum(GetAssayData(object = rust2020_raw_subset, slot = "counts")["Hml",]>0)
sum(GetAssayData(object = jevitt2020_raw_subset, slot = "counts")["br",]>5)
sum(GetAssayData(object = rust2020_raw_subset, slot = "counts")["br",]>5)

# Select features for integration 

genelists <- loadGeneLists()
features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000)

pub_markers <- append(genelists$rust2020_markers,genelists$jevitt2020_markers)
sfcgenes <- append(genelists$sfc_tfc_markers,genelists$v_atpases)
feature_addons <- append(pub_markers,sfcgenes)
features <- append(features,feature_addons) 
features <- unique(features)

# Integrate data using Seurat RPCA
set.seed(1)

data_list <- PrepSCTIntegration(object.list = data_list, anchor.features = features)

data_list <- lapply(X = data_list, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = TRUE)
})

anchors <- FindIntegrationAnchors(object.list = data_list, 
                                  normalization.method = "SCT",
                                  anchor.features = features, 
                                  dims = 1:30, reduction = "rpca", 
                                  k.anchor = 100)
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)


combined.sct <- RunPCA(combined.sct, verbose = TRUE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30)
combined.clusters <- FindClusters(combined.sct, resolution = 5)
combined.clusters[["ClusterNames_5"]] <- Idents(object = combined.clusters)

allmarkers <- FindAllMarkers(combined.clusters,only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

# Annotate previously identified Stretch FC and Transitional FC populations from both original datasets in the integrated object
## extract  cell barcodes of cells previously identified as SFCs by Jevitt2O20 and Rust2020
jevittsfc1_cells <- FetchData(object = combined.clusters, vars = c("jevitt2020.ident")) %>% filter(jevitt2020.ident=="Stretched Cells 1") %>% rownames()
jevittsfc2_cells <- FetchData(object = combined.clusters, vars = c("jevitt2020.ident")) %>% filter(jevitt2020.ident=="Stretched Cells 2") %>% rownames()
rust67_cells <- FetchData(object = combined.clusters, vars = c("rust2020.ident")) %>% filter(rust2020.ident=="2.1.0 stretch cell 6-7") %>% rownames()
rust8_cells <- FetchData(object = combined.clusters, vars = c("rust2020.ident")) %>% filter(rust2020.ident=="2.1.1 stretch cell 8+") %>% rownames()

transitionalFC_cells <- FetchData(object = combined.clusters, vars = c("jevitt2020.ident")) %>% filter(jevitt2020.ident=="Transitional FCs (Stg. 6-7)") %>% rownames()



## Save new annotations 
combined.clusters$known_sfc_cells<- ifelse(colnames(combined.clusters) %in% jevittsfc1_cells, "JevittSFC 1",
                                           ifelse(colnames(combined.clusters) %in% jevittsfc2_cells,"JevittSFC 2",
                                                  ifelse(colnames(combined.clusters) %in% rust67_cells, "rustSFC 6-7",
                                                         ifelse(colnames(combined.clusters) %in% rust8_cells,"rustSFC 8+",
                                                                ifelse(colnames(combined.clusters) %in% transitionalFC_cells,"Jevitt TFC","NotSFCs"))))) 





saveRDS(anchors,file=file.path(output_dir,"anchors.RDS"))
saveRDS(combined.clusters,file=file.path(output_dir,"integrated.RDS"))
saveRDS(allmarkers,file=file.path(output_dir,"allmarkers.RDS"))

