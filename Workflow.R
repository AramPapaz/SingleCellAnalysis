# Libraries
library(dplyr) 
library(Seurat) 
library(patchwork)
library(DoubletFinder)
library(SingleR)
#library(velocyto.R)
library(enrichR)
library(CellChat)
library(SingleCellExperiment)
library(SeuratWrappers)
library(tidyverse)
library(celldex)
library(monocle3)
library(presto)



# Bone Marrow Mononuclear Cells
bmc1=readRDS(choose.files())
bmc2=readRDS(choose.files())
bmc1=CreateSeuratObject(bmc1)
bmc2=CreateSeuratObject(bmc2)
# CD34+ Enriched Bone Marrow Cells
cd1=readRDS(choose.files())
cd2=readRDS(choose.files())
cd1=CreateSeuratObject(cd1)
cd2=CreateSeuratObject(cd2)



# Filtering
bmc1[["percent.mt"]] <- PercentageFeatureSet(bmc1, pattern = "^MT-")
bmc2[["percent.mt"]] <- PercentageFeatureSet(bmc2, pattern = "^MT-")
cd1[["percent.mt"]] <- PercentageFeatureSet(cd1, pattern = "^MT-")
cd2[["percent.mt"]] <- PercentageFeatureSet(cd2, pattern = "^MT-")

bmc1 <- subset(bmc1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
bmc2 <- subset(bmc2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cd1 <- subset(cd1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cd2 <- subset(cd2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalization
bmc1 <- NormalizeData(bmc1)
bmc2 <- NormalizeData(bmc2)
cd1 <- NormalizeData(cd1)
cd2 <- NormalizeData(cd2)

# Feature Selection
bmc1 <- FindVariableFeatures(bmc1)
bmc2 <- FindVariableFeatures(bmc2)
cd1 <- FindVariableFeatures(cd1)
cd2 <- FindVariableFeatures(cd2)

# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(bmc1), 10)

# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(bmc1)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot2


# Scaling 
bmc1 <- ScaleData(bmc1)
bmc2 <- ScaleData(bmc2)
cd1 <- ScaleData(cd1)
cd2 <- ScaleData(cd2)


# PCA
bmc1=RunPCA(bmc1)
bmc2=RunPCA(bmc2)
cd1=RunPCA(cd1)
cd2=RunPCA(cd2)


# Cluster
bmc1 <- FindNeighbors(bmc1, dims = 1:15)
bmc1 <- FindClusters(bmc1, resolution = 0.5)

bmc2 <- FindNeighbors(bmc2, dims = 1:15)
bmc2 <- FindClusters(bmc2, resolution = 0.5)

cd1 <- FindNeighbors(cd1, dims = 1:15)
cd1 <- FindClusters(cd1, resolution = 0.5)

cd2 <- FindNeighbors(cd2, dims = 1:15)
cd2 <- FindClusters(cd2, resolution = 0.5)

# UMAP
bmc1 <- RunUMAP(bmc1, dims = 1:15)
bmc2 <- RunUMAP(bmc2, dims = 1:15)
cd1 <- RunUMAP(cd1, dims = 1:15)
cd2 <- RunUMAP(cd2, dims = 1:15)
#DimPlot(bmc1,reduction = "umap")


#   BMC1 Doublets

## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.bmc1 <- paramSweep(bmc1, PCs = 1:15, sct = FALSE)
sweep.stats_bmc1 <- summarizeSweep(sweep.res.bmc1, GT = FALSE)
bcmvn_bmc1 <- find.pK(sweep.stats_bmc1)
pK <- bcmvn_bmc1 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

annotations=bmc1@meta.data$seurat_clusters
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*nrow(bmc1@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

bmc1 <- doubletFinder(bmc1, PCs = 1:15, 
                                         pN = 0.25, 
                                         pK = pK, 
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = FALSE, sct = FALSE)



bmc1=subset(bmc1, subset = DF.classifications_0.25_0.005_227 == "Singlet")


#   BMC2 Doublets
sweep.res.bmc2 <- paramSweep(bmc2, PCs = 1:15, sct = FALSE)
sweep.stats_bmc2 <- summarizeSweep(sweep.res.bmc2, GT = FALSE)
bcmvn_bmc2 <- find.pK(sweep.stats_bmc2)
pK <- bcmvn_bmc2 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

annotations=bmc2@meta.data$seurat_clusters
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*nrow(bmc2@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

bmc2 <- doubletFinder(bmc2, PCs = 1:15, 
                      pN = 0.25, 
                      pK = pK, 
                      nExp = nExp_poi.adj,
                      reuse.pANN = FALSE, sct = FALSE)



bmc2=subset(bmc2, subset = DF.classifications_0.25_0.17_231 == "Singlet")

#    CD1 doublets
sweep.res.cd1 <- paramSweep(cd1, PCs = 1:15, sct = FALSE)
sweep.stats_cd1 <- summarizeSweep(sweep.res.cd1, GT = FALSE)
bcmvn_cd1 <- find.pK(sweep.stats_cd1)
pK <- bcmvn_cd1 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

annotations=cd1@meta.data$seurat_clusters
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.016*nrow(cd1@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

cd1 <- doubletFinder(cd1, PCs = 1:15, 
                      pN = 0.25, 
                      pK = pK, 
                      nExp = nExp_poi.adj,
                      reuse.pANN = FALSE, sct = FALSE)




cd1=subset(cd1, subset = DF.classifications_0.25_0.11_22 == "Singlet")



#    CD2 doublets
sweep.res.cd2 <- paramSweep(cd2, PCs = 1:15, sct = FALSE)
sweep.stats_cd2 <- summarizeSweep(sweep.res.cd2, GT = FALSE)
bcmvn_cd2 <- find.pK(sweep.stats_cd2)
pK <- bcmvn_cd2 %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))

annotations=cd2@meta.data$seurat_clusters
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.039*nrow(cd2@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

cd2 <- doubletFinder(cd2, PCs = 1:15, 
                     pN = 0.25, 
                     pK = pK, 
                     nExp = nExp_poi.adj,
                     reuse.pANN = FALSE, sct = FALSE)



cd2=subset(cd2, subset = DF.classifications_0.25_0.005_184 == "Singlet")



# Merge
merged_data=merge(bmc1,c(bmc2,cd1,cd2),add.cell.ids = c("bmc1", "bmc2", "cd1", "cd2"))
merged_data[["RNA"]]=JoinLayers(merged_data[["RNA"]])

# Batch Effect Correction
merged_data$sample <- rownames(merged_data@meta.data)

# split sample column
merged_data@meta.data <- separate(merged_data@meta.data, col = 'sample', into = c('Patient', 'Type', 'Barcode'), 
                                    sep = '_')
# perform standard workflow steps to figure out if we see any batch effects --------
merged_data <- NormalizeData(object = merged_data)
merged_data <- FindVariableFeatures(object = merged_data)
merged_data <- ScaleData(object = merged_data)
merged_data <- RunPCA(object = merged_data)
ElbowPlot(merged_data,ndims = 100)

merged_data <- FindNeighbors(object = merged_data, dims = 1:30)
merged_data <- FindClusters(object = merged_data)
merged_data <- RunUMAP(object = merged_data, dims = 1:30)
merged_data <- RunTSNE(object = merged_data, dims = 1:30)


# plot
DimPlot(merged_data, reduction = 'umap', group.by = 'Patient')
DimPlot(merged_data, reduction = 'umap', group.by = 'Donor')
DimPlot(merged_data, reduction = 'umap', group.by = 'seurat_clusters')
DimPlot(merged_data, reduction = 'tsne', group.by = 'seurat_clusters')

# Correction
# perform integration to correct for batch effects ------
obj.list <- SplitObject(merged_data, split.by = 'Patient')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}


# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
ElbowPlot(seurat.integrated)

seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:20)


DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Patient')

# Annotations
merged_data@meta.data <- separate(merged_data@meta.data, col = 'Barcode', into = c('Info', 'Barcode'), 
                                  sep = ':')
merged_data@meta.data <- merged_data@meta.data %>% mutate(Info=sub("(.{2})(.*)", "\\1_\\2", Info))
merged_data@meta.data <- separate(merged_data@meta.data, col = 'Info', into = c('Donor', 'Replicate'), 
                                  sep = '_')

fun <- function(x){
  if(x=="D2"){
    "M"
  }
  else{
    "F"
  }
}
merged_data@meta.data$Sex <- sapply(merged_data@meta.data$Donor,FUN = fun)

# DE analysis
merged_data[["RNA"]]=JoinLayers(merged_data[["RNA"]])
markers <- FindAllMarkers(merged_data, only.pos = TRUE)

VlnPlot(merged_data,features = "CD14")
FeaturePlot(merged_data,features = "CD14")

# Automatic Annotation
hpca.se <- HumanPrimaryCellAtlasData()

pbmc_counts <- GetAssayData(merged_data, slot = 'counts')

pred <- SingleR(test = pbmc_counts,
                ref = hpca.se,
                labels = hpca.se$label.main)

merged_data$singleR.labels <- pred$labels[match(rownames(merged_data@meta.data), rownames(pred))]
DimPlot(merged_data, reduction = 'umap', group.by = 'singleR.labels',label = TRUE)

# Change Idents
Idents(merged_data) <- merged_data@meta.data$singleR.labels

# Subset data
traj<- subset(merged_data,idents=c("Monocyte","Myelocyte","Pro-Myelocyte","HSC_-G-CSF","Pre-B_cell_CD34-"))
DimPlot(traj, reduction = 'umap', group.by = 'singleR.labels')

# Convert
cds <- as.cell_data_set(traj)
fData(cds)$gene_short_name <- rownames(fData(cds))
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)


cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 

list_cluster <- traj@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- traj@reductions$umap@cell.embeddings

cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = "redefined_cluster",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan')) +
  theme(legend.position = "right")


#  Learn trajectory graph 
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = 'redefined_cluster',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

# Order the cells in pseudotime 

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == "Pro-Myelocyte"]))
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(singleR.labels, monocle3_pseudotime, median), fill = singleR.labels)) +
  geom_boxplot()

# Finding genes that change as a function of pseudotime
deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()
FeaturePlot(traj, features = c('HES4', 'PARK7', 'ENO1'))

traj$pseudotime <- pseudotime(cds)
Idents(traj) <- traj$singleR.labels
FeaturePlot(traj, features = "pseudotime", label = T)

# DE between certain conditions
Idents(merged_data) <- merged_data@meta.data$orig.ident
group_markers <- FindMarkers(merged_data, ident.1 = "BMMC", ident.2 = "CD34")

# Top 5 with p-value
group_markers %>% arrange(p_val_adj) %>% head(5)
# Top 5 fold change
group_markers %>% arrange(abs(avg_log2FC)) %>% tail(5)
sigmarkers <- group_markers %>% filter(p_val_adj <= 0.05)

# Go Enrichment analysis
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()
dbs %>% filter(grepl("GO",libraryName))
dbs <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023")
if (websiteLive) {
  enriched <- enrichr(rownames(sigmarkers), dbs)
}

if (websiteLive) {
  plotEnrich(enriched[[3]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}

# CellChat
BM <- subset(merged_data,ident="BMMC")
CD <- subset(merged_data,ident="CD34")
Idents(BM) <- BM@meta.data$singleR.labels
Idents(CD) <- CD@meta.data$singleR.labels

shared_cell_types <- intersect(unique(BM@meta.data$singleR.labels),unique(CD@meta.data$singleR.labels))
Idents(merged_data) <- merged_data@meta.data$singleR.labels
merged_data <- subset(merged_data, idents=shared_cell_types)
Idents(merged_data) <- merged_data@meta.data$orig.ident
cellchat <- createCellChat(object = CD)

CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)

# set the used database in the object
cellchat@DB <- CellChatDB
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)


execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
#> [1] 12.42671
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:8) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

netVisual_aggregate(cellchat, signaling = "CD99", layout = "circle")
