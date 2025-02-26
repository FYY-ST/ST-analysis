# Day 25 and Day 70 Dorsal-Ventralpatterned forebrain organoids 

# Load libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(glmGamPoi)
library(heatmaply)
library(dplyr)
library(Matrix)
library(topicmodels)
library(devtools)
library(Rfast2)

# Load data from Rds objects

list.files("C:/Users/feiyu/Desktop/SamsungLaptop11-7-2024/GEO submission ST organoid/GEO submission-ST organoid") 
DVday25 <- readRDS("Day25_DV_FB Processed.Rds")
DVday70 <- readRDS("Day70_DV_FB.Rds")
#===================================================================================================================
# https://satijalab.org/seurat/articles/spatial_vignette
# Day25 D-V patterned organoid
# Spatial plot and UMAP
p1<- DimPlot(DVday25, reduction = "umap", pt.size = 2, label = FALSE)+ theme(legend.position = "right")
p2<- SpatialDimPlot(DVday25, label = FALSE, label.size = 10, pt.size.factor = 8)+ theme(legend.position = "top")
p1|p2

colfunc <- colorRampPalette(c('black','#A8186E','#F2541B','#FFCE61','#FFE58A'))

# Dorsal
SpatialFeaturePlot(DVday25, features = c('PAX6'), pt.size.factor = 8, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(DVday25, features = c('TBR1'), pt.size.factor = 8, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(DVday25, features = c('LHX5'), pt.size.factor = 8, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))

# MGE
SpatialFeaturePlot(DVday25, features = c('NKX2-1'), pt.size.factor = 8, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(DVday25, features = c('NKX6-2'), pt.size.factor = 8, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(DVday25, features = c('LHX6'), pt.size.factor = 8, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))

# LGE/CGE
SpatialFeaturePlot(DVday25, features = c('GSX2'), pt.size.factor = 8, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(DVday25, features = c('DLX2'), pt.size.factor = 8, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(DVday25, features = c('MEIS2'), pt.size.factor = 8, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))

# Dot plot
MarkerGenesD25 <- c('FOXG1','VIM','SOX2','NES',
                 'PAX6','TBR1','LHX1','LHX5','WNT7B',
                 'GSX2','MEF2C','SIX3','DLX1','DLX2','DLX5',
                 'NKX6-2','SOX6','NKX2-1','NKX2-2','LHX6','SIX6')
DotPlot(DVday25, feature=MarkerGenesD25) + RotatedAxis()

#===================================================================================================================
# Day70 D-V patterned organoid
# Spatial plot and UMAP
mycol <- c('RGC'="#F8766D", 'VNPC'= "#C49A00", 'GABA'="#53B400",
            'IN'="#00B6EB", 'Unknown'="lightgrey", 'DNPC'="#A58AFF", 
            'EX'="#FB61D7")

p1<- DimPlot(DVday70, reduction = "umap", pt.size = 2, label = FALSE, cols = mycol)+ theme(legend.position = "right")
p2<- SpatialDimPlot(DVday70, label = FALSE, label.size = 10, pt.size.factor = 7, cols = mycol)+ theme(legend.position = "top")
p1|p2


# Dorsal
SpatialFeaturePlot(DVday70, features = c('NEUROD2'), pt.size.factor = 6, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(DVday70, features = c('TBR1'), pt.size.factor = 6, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(DVday70, features = c('HOPX'), pt.size.factor = 6, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))

# Ventral
SpatialFeaturePlot(DVday70, features = c('SIX3'), pt.size.factor = 6, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(DVday70, features = c('GABRB2'), pt.size.factor = 6, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(DVday70, features = c('NES'), pt.size.factor = 6, min.cutoff = 0) + scale_fill_gradientn(colours = colfunc(5))

# Dot plot, remove "Unkown" cell type
DVday701 <- subset(DVday70, idents = c('VNPC', 'DNPC','RGC','EX', 'IN','GABA'))
                                       
MarkerGenesD70  <- c('SOX2','HES1','NES','GLI3','SFRP1','PAX6','LHX2','SLC1A3','EOMES','EMX1','INSM1','NHLH1','HES6',
                   'HOPX','PEA15','LGALS3BP','MOXD1','FABP7','PTN',
                   'NEUROD6','NEUROD2','TBR1','PLXNA4','SATB2','BHLHE22','CACNA2D1',
                   'SP8','SCGN','CALB2','LHX5','GAD1','GAD2',
                   'SST','SLC32A1','NKX2-1','SIX3','PBX3','GABRB3','GABRB2','GABRG1')


DotPlot(DVday701, feature=MarkerGenesD70) + RotatedAxis()

#===================================================================================================================
#===================================================================================================================
#===================================================================================================================
# Umap integrations, please install the Harmony packages
library(tidyverse)
library(Rcpp)
library(harmony)

new.cluster.ids1 <- c("Ventral_VNPC", 'Dorsal_DNPC','Dorsal_RGC','Dorsal_EX','Ventral_IN',
                     'Ventral_GABA', 'Unknown_Unknown')
names(new.cluster.ids1) <- levels(DVday70)
DVday70T <- RenameIdents(DVday70, new.cluster.ids1)
DVday70T@meta.data$CellType <- Idents(DVday70T) 
p1<- DimPlot(DVday70T, reduction = "umap", pt.size = 2, label = FALSE)+ theme(legend.position = "right")
p2<- SpatialDimPlot(DVday70T, label = FALSE, label.size = 10, pt.size.factor = 7)+ theme(legend.position = "top")
p1|p2

new.cluster.ids2<- c("Dorsal_Dorsalday25", 'Ventral_LGE','Ventral_MGE')
names(new.cluster.ids2) <- levels(DVday25)
DVday25T <- RenameIdents(DVday25, new.cluster.ids2)
DVday25T@meta.data$CellType <- Idents(DVday25T ) 
p1<- DimPlot(DVday25T, reduction = "umap", pt.size = 2, label = FALSE)+ theme(legend.position = "right")
p2<- SpatialDimPlot(DVday25T, label = FALSE, label.size = 10, pt.size.factor = 8)+ theme(legend.position = "top")
p1|p2

DVmerged <- merge(DVday25T, y = c(DVday70T),
                   add.cell.ids = c('DVday25', 'DVday70'),
                   project = "DVmerge")

DVmerged$sample <- rownames(DVmerged@meta.data)
DVmerged@meta.data <- separate(DVmerged@meta.data, col = "sample", into = c("OrID", 'Barcode'),
                                sep = '_')
DVmerged@meta.data <- separate(DVmerged@meta.data, col = "CellType", into = c("DV",'CType'),
                                sep = '_')
# sanity check
unique(DVmerged@meta.data$OrID)
unique(DVmerged@meta.data$CType)
unique(DVmerged@meta.data$DV)

DVmerged <- SCTransform(DVmerged, assay = "Spatial", verbose = FALSE)
DVmerged <- RunPCA(DVmerged, assay = "SCT", verbose = FALSE)
DVmerged <- FindNeighbors(DVmerged, reduction = "pca", dims = 1:20)
DVmerged <- RunUMAP(DVmerged, dims = 1:20)
Idents(DVmerged)  <- DVmerged@meta.data$DV
DVmerged1 <- subset(DVmerged, idents = c(
  'Ventral','Dorsal'))   

DVmerge.harmony <- DVmerged1  %>%
  RunHarmony(group.by.vars = 'OrID', plot_convergence = FALSE)

DVmerge.harmony@reductions

DVmerge.harmony.embed <- Embeddings(DVmerge.harmony, "harmony")

DVmerge.harmony <- DVmerge.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1)

mycolorall <- c("#A58AFF",'pink',"#FB61D7","#53B400", "#00B6EB",'cyan','blue',"#F8766D","#C49A00")


p1 <- DimPlot(DVmerge.harmony, reduction = 'umap', group.by = 'OrID')
p2 <- DimPlot(DVmerge.harmony, reduction = 'umap', group.by = 'DV', label = FALSE, pt.size =1)
p3 <- DimPlot(DVmerge.harmony, reduction = 'umap', group.by = 'CType', label = FALSE, cols =mycolorall)
p1|p2|p3

#===================================================================================================================
#===================================================================================================================
# Trajectory analysis, please install the Monocle3 package
library(monocle3)
library(SeuratWrappers)
# ...1 Convert to cell_data_set object ------------------------
cds <- as.cell_data_set(DVmerge.harmony)
cds

# to get cell metadata
colData(cds)
# to gene metdata
fData(cds)
rownames(fData(cds))[1:10]

# since it misses the gene_short_name column, let's add it
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)

cds@clusters$UMAP$partitions

# assign partitions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@colData
reacreate.partition
cds@clusters$UMAP$partitions <- reacreate.partition
unique(cds@colData$CType)

# Manually assign partitions based on CellTypes
re.partitionDorsal <- c(rep(1,length(cds@colData@rownames)))
for (i in 1:802){
  if (cds@colData$CType[i] == c('MGE')){
    re.partitionDorsal[i] <- 2}
  
  if (cds@colData$CType[i] == c('LGE')){
    re.partitionDorsal[i] <- 2}
  
  if (cds@colData$CType[i] == c('IN')){
    re.partitionDorsal[i] <- 2}
  
  if (cds@colData$CType[i] == c('VNPC')){
    re.partitionDorsal[i] <- 2}
  
  if (cds@colData$CType[i] == c('GABA')){
    re.partitionDorsal[i] <- 2}
}

re.partitionDorsal

names(re.partitionDorsal) <- cds@colData@rownames
re.partitionDorsal <- as.factor(re.partitionDorsal)

list_cluster <-  DVmerge.harmony@active.ident
cds@clusters$UMAP$clusters <- list_cluster
cds@int_colData@listData$reducedDims$UMAP <- DVmerge.harmony@reductions$umap@cell.embeddings

# plot
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = 'CType',
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) 

theme(legend.position = "right")

cluster.before.trajectory | cluster.names

# ...3. Learn trajectory graph ------------------------
cds <- learn_graph(cds, use_partition = TRUE)

p1 <- plot_cells(cds,
                 color_cells_by = 'CType',
                 label_groups_by_cluster = FALSE,
                 label_branch_points = FALSE,
                 label_roots = FALSE,
                 label_leaves = FALSE,
                 group_label_size = 0,
                 cell_size = 0.5)+
theme(legend.position = "right") #+ scale_color_manual(values = mycolorVentral)
# ...4. Order the cells in pseudotime -------------------
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == c('MGEDV02','LGEDV02', 'DV02')]))

p2 <- plot_cells(cds,
                 color_cells_by = 'pseudotime',
                 label_groups_by_cluster = FALSE,
                 label_branch_points = FALSE,
                 label_roots = FALSE,
                 label_leaves = FALSE, cell_size = 0.5)+
  theme(legend.position = "right") 

p1|p2

# cells ordered by monocle3 pseudotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(CType, monocle3_pseudotime, median), fill = CType)) +
  geom_boxplot()

# visualizing pseudotime in seurat
DVmerge.harmony$pseudotime <- pseudotime(cds)
Idents(DVmergeDorsal.harmony) <- DVmergeDorsal.harmony$CType
Idents(DVmerge.harmony) <- DVmergeDorsal.harmony$CType
FeaturePlot(DVmergeDorsal.harmony, features = "pseudotime", label = T)
FeaturePlot(DVmerge.harmony, features = "pseudotime", label = F, pt.size = 2)

#===================================================================================================================
#===================================================================================================================
# STdeconvolve Day70 D-V patterned forebrain organoid for modeling single cell resolutions.
# https://jef.works/STdeconvolve/getting_started.html
# Please install STdeconvolve package
#=======================================================================================
# Get positions from Seurat object
posDVraw05 <- GetTissueCoordinates(DVday70)
DV05count <- DVday70@assays$Spatial$counts  # count matrix

#Load library
library(STdeconvolve)

corpusDV05Apha0.05 <- restrictCorpus(DV05count, 
                                     removeAbove=5, 
                                     removeBelow = 0.05, 
                                     nTopOD = NA, 
                                     alpha = 0.05)


## double check genes of interest are kept
genes.test <- c('PAX6', 'DLX2', 'TBR1','NKX2-1','SIX3')
genes.test %in% rownames(corpusDV05Apha0.05)

ldas7Apha <- fitLDA(t(as.matrix(corpusDV05Apha0.05)), Ks = c(7))
optLDA7Apha <- optimalModel(models = ldas7Apha, opt = "7")
results7Apha <- getBetaTheta(optLDA7Apha, perc.filt = 0.05, betaScale = 1000)
deconProp7Apha <- results7Apha$theta
deconGexp7Apha <- results7Apha$beta

## visualize deconvolved cell-type proportions
colnames(posDVraw05) <- c('x', 'y')
vizAllTopics(deconProp7Apha, posDVraw05[rownames(deconProp7Apha),], r=2, lwd=0)	

# plot 
# Loop through all cell-types to visulize together
ps <- lapply(colnames(deconProp7Apha), function(celltype) {
  
  vizTopic(theta = deconProp7Apha, pos =posDVraw05, topic = celltype, plotTitle = paste0("X", celltype),
           size = 2, stroke = 1, alpha = 1,
           low = "white",
           high = "red")
  
})
gridExtra::grid.arrange(grobs = ps, layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6), c(7)))

# Umap for STdeconvolve
#===============================================================

pcs.info <- stats::prcomp(t(log10(as.matrix(DV05count)+1)), center=TRUE)
nPcs <- 20 ## let's take the top 5 PCs
pcs <- pcs.info$x[,1:nPcs]

emb <- Rtsne::Rtsne(pcs,
                    is_distance=FALSE,
                    perplexity=30,
                    num_threads=1,
                    verbose=FALSE)$Y
rownames(emb) <- rownames(pcs)
colnames(emb) <- c("x", "y")

#===============================================================
k <- 20
com <- MERINGUE::getClusters(pcs, k, weight=TRUE, method = igraph::cluster_louvain)

levels(com)

tempCom <- com

# colorManual <- c("grey", 'gold','green',"pink",'purple','red')

colorManual <- c('grey',"grey","pink",'purple','gold', 'grey','blue','cyan')

#colorManual <- c("grey", 'gold','green',"pink",'purple','red')

dat <- data.frame("emb1" = posDVraw05[,"x"],
                  "emb2" = posDVraw05[,"y"],
                  "Cluster" = tempCom)


plt <- ggplot2::ggplot(data = dat) +
  
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 5) +
  
  
  #  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +
  ggplot2::scale_color_manual(values = colorManual) +
  
  # ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +
  
  ggplot2::labs(title = "",
                x = "x",
                y = "y") +
  
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +
  
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +
  
  ggplot2::coord_equal()

plt

# and on the 2D embedding
tempCom <- com

dat <- data.frame("emb1" = emb[,1],
                  "emb2" = emb[,2],
                  "Cluster" = tempCom)

## cluster labels
cent.pos <- do.call(rbind, tapply(1:nrow(emb), tempCom, function(ii) apply(emb[ii,,drop=F],2,median)))
cent.pos <- as.data.frame(cent.pos)
colnames(cent.pos) <- c("x", "y")
cent.pos$cluster <- rownames(cent.pos)
cent.pos <- na.omit(cent.pos)

plt <- ggplot2::ggplot(data = dat) +
  ggplot2::geom_point(ggplot2::aes(x = emb1, y = emb2,
                                   color = Cluster), size = 2) +
  
  ggplot2::scale_color_manual(values = rainbow(n = length(levels(tempCom)))) +
  #  ggplot2::scale_color_manual(values = colorManual) +
  
  ggplot2::scale_y_continuous(expand = c(0, 0), limits = c( min(dat$emb2)-1, max(dat$emb2)+1)) +
  ggplot2::scale_x_continuous(expand = c(0, 0), limits = c( min(dat$emb1)-1, max(dat$emb1)+1) ) +
  
  ggplot2::labs(title = "",
                x = "t-SNE 1",
                y = "t-SNE 2") +
  
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=15, color = "black"),
                 axis.text.y = ggplot2::element_text(size=15, color = "black"),
                 axis.title.y = ggplot2::element_text(size=15),
                 axis.title.x = ggplot2::element_text(size=15),
                 axis.ticks.x = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(size=15),
                 legend.text = ggplot2::element_text(size = 12, colour = "black"),
                 legend.title = ggplot2::element_text(size = 15, colour = "black", angle = 0, hjust = 0.5),
                 panel.background = ggplot2::element_blank(),
                 plot.background = ggplot2::element_blank(),
                 panel.grid.major.y =  ggplot2::element_blank(),
                 axis.line = ggplot2::element_line(size = 1, colour = "black")
                 # legend.position="none"
  ) +
  
  ggplot2::geom_text(data = cent.pos,
                     ggplot2::aes(x = x,
                                  y = y,
                                  label = cluster),
                     fontface = "bold") +
  
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), ncol = 2)
  ) +
  
  ggplot2::coord_equal()

plt

# See the proportions of each deconvolved cell-type across the embedding
ps <- lapply(colnames(deconProp6Apha), function(celltype) {
  
  vizTopic(theta = deconProp6Apha, pos = emb, topic = celltype, plotTitle = paste0("X", celltype),
           size = 2, stroke = 0.5, alpha = 1,
           low = "white",
           high = "red") +
    
    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")
  
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3),
                        c(4, 5, 6)
  )
)

#======================================================================================================
# Correlations
# com_proxyTheta showing each belonging of each barcode
# If one barcode is in one cluster, it show "1" in the matrix, otherwise it shows "0"

# proxy theta for the txn clusters
com_proxyTheta <- model.matrix(~ 0 + com)
rownames(com_proxyTheta) <- names(com)
# fix names
colnames(com_proxyTheta) <- unlist(lapply(colnames(com_proxyTheta), function(x) {
  unlist(strsplit(x, "com"))[2]
}))
com_proxyTheta <- as.data.frame.matrix(com_proxyTheta)
com_proxyTheta[1:5,1:5]

corMat_prop <- STdeconvolve::getCorrMtx(m1 = as.matrix(com_proxyTheta),
                                        m2 = deconProp7Apha,
                                        type = "t")

rownames(corMat_prop) <- paste0("com_", seq(nrow(corMat_prop)))
colnames(corMat_prop) <- paste0("decon_", seq(ncol(corMat_prop)))

## order the cell-types rows based on best match (highest correlation) with each community 
pairs <- STdeconvolve::lsatPairs(corMat_prop)
m <- corMat_prop[pairs$rowix, pairs$colsix]

STdeconvolve::correlationPlot(mat = m,
                              colLabs = "Transcriptional clusters",
                              rowLabs = "STdeconvolve") +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 90)
  )
#=================================================================================================




















