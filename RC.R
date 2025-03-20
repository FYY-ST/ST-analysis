# Day 20 rostal-caudal forebrain-midbrain patterned organoid
# Day 45 rostal-caudal forebrain-midbrain patterned organoid
# Day 45 rostal-caudall fore-mid/hindbrain/SC patterned organoid

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

list.files("~") 
FMday20 <- readRDS("Day20_RC_FM_Processed.Rds")
FMday45 <- readRDS("Day45_RC_FM_Processed.Rds")
FMHday45 <- readRDS("Day45_RC_FMH_Processed.Rds")
#===================================================================================================================
# https://satijalab.org/seurat/articles/spatial_vignette
# Day20 R-C fore-midbrain patterned organoid
# Spatial plot and UMAP
p1<- DimPlot(FMday20, reduction = "umap", pt.size = 2, label = FALSE)+ theme(legend.position = "right")
p2<- SpatialDimPlot(FMday20, label = FALSE, label.size = 10, pt.size.factor = 10)+ theme(legend.position = "top")
p1|p2

colfunc <- colorRampPalette(c('black','#A8186E','#F2541B','#FFCE61','#FFE58A'))

# Forebrain
SpatialFeaturePlot(FMday20, features = "FOXG1", pt.size.factor = 10)+ scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(FMday20, features = "SOX9", pt.size.factor = 10)+ scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(FMday20, features = "LHX2", pt.size.factor = 10)+ scale_fill_gradientn(colours = colfunc(5))

# Midbrain
SpatialFeaturePlot(FMday20, features = "FOXA2", pt.size.factor = 10)+ scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(FMday20, features = "WLS", pt.size.factor = 10)+ scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(FMday20, features = "FOXD1", pt.size.factor = 10)+ scale_fill_gradientn(colours = colfunc(5))

# Progenitors, Scale was adjusted to match the high expression profiles.
SpatialFeaturePlot(FMday20, features = "OTX2", pt.size.factor = 10)+ scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(FMday20, features = "NES", pt.size.factor = 10)+ scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(FMday20, features = "SOX2", pt.size.factor = 10)+ scale_fill_gradientn(colours = colfunc(5))

# Hindbrain/SC, HOX genes were not expressed in this organoid. 
SpatialFeaturePlot(FMday20, features = "IRX3", pt.size.factor = 10)+ scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(FMday20, features = "OLIG2", pt.size.factor = 10)+ scale_fill_gradientn(colours = colfunc(5))
SpatialFeaturePlot(FMday20, features = "HOXA2", pt.size.factor = 10)+ scale_fill_gradientn(colours = colfunc(5))

# Dotplot
MarkerGenesFM <- c('FOXG1', 'LHX2', 'GLI3', 'EMX1','SIX3','SOX9','EOMES','FEZF1', 'FEZF2',
                   'OTX2','OTX1','FOXA1','FOXA2','PAX5', 'LMX1A', 'LMX1B', 'SIM1', 'EN1' ,'TH',
                   'SDC1','WLS','FOXB1','PTX3','NKX6-1','RSPO2','FOXD1','PTCH1','MGST1','POU3F4',
                   'OLIG2','OLIG3','IRX3')
DotPlot(FMday20, features =MarkerGenesFM )+ RotatedAxis()

#===================================================================================================================
# Day45 R-C fore-midbrain patterned organoid
# Spatial plot and UMAP
p1<- DimPlot(FMday45, reduction = "umap", pt.size = 2, label = FALSE)+ theme(legend.position = "right")
p2<- SpatialDimPlot(FMday45, label = FALSE, label.size = 10, pt.size.factor = 10)+ theme(legend.position = "top")
p1|p2

colfunc <- colorRampPalette(c('black','#301D7D', '#A8186E','#F2541B','#FFCE61', '#FFE58A'))
colfunc(6)

SpatialFeaturePlot(FMday45, features = c('SIX3'), pt.size.factor = 10, min.cutoff = 0.5)  + scale_fill_gradientn(colours = colfunc(6))
SpatialFeaturePlot(FMday45, features = c('FOXG1'), pt.size.factor = 10)  + scale_fill_gradientn(colours = colfunc(6))
SpatialFeaturePlot(FMday45, features = c('OTX1'), pt.size.factor = 10)  + scale_fill_gradientn(colours = colfunc(6))
SpatialFeaturePlot(FMday45, features = c('FOXA1'), pt.size.factor = 10) + scale_fill_gradientn(colours = colfunc(6))
SpatialFeaturePlot(FMday45, features = c('NES'), pt.size.factor = 10) + scale_fill_gradientn(colours = colfunc(6))
SpatialFeaturePlot(FMday45, features = c('HES1'), pt.size.factor = 10) + scale_fill_gradientn(colours = colfunc(6))

# Dot plot, remove "Unkown" cell type
FMday45T <- subset(FMday45, idents = c('FB', 'MB','Progenitor'))
MarkerGenesFM45<-  c('SOX2','SOX1','NES','HES5','LEF1','HES1','GLI3',
                     'FOXG1', 'DLX2','DLX1','SST', 'NKX2-1','SIX3','OTX2',
                     'OTX1','FOXA1','FOXA2','FGF8','PAX5', 'LMX1A', 'LMX1B', 'WLS','FOXB1',
                     'NKX6-1','RSPO2','FOXD1','POU3F4')
DotPlot(FMday45T, features =MarkerGenesFM45)+ RotatedAxis()

#===================================================================================================================
# Day45 R-C fore-mid/hind/sc patterned organoid
p1<- DimPlot(FMHday45, reduction = "umap", pt.size = 2, label = FALSE)+ theme(legend.position = "right")
p2<- SpatialDimPlot(FMHday45, label = FALSE, label.size = 10, pt.size.factor = 10)+ theme(legend.position = "top")
p1|p2

#Forebrain
SpatialFeaturePlot(FMHday45, features = c('SST'), pt.size.factor = 8, min.cutoff = 0.5)  + scale_fill_gradientn(colours = colfunc(6))
SpatialFeaturePlot(FMHday45, features = c('GSX2'), pt.size.factor = 8, min.cutoff = 0.5)  + scale_fill_gradientn(colours = colfunc(6))
SpatialFeaturePlot(FMHday45, features = c('SOX9'), pt.size.factor = 8, min.cutoff = 0.5)  + scale_fill_gradientn(colours = colfunc(6))

#Midbrain
SpatialFeaturePlot(FMHday45, features = c('LMX1B'), pt.size.factor = 8, min.cutoff = 0.5)  + scale_fill_gradientn(colours = colfunc(6))
SpatialFeaturePlot(FMHday45, features = c('EN1'), pt.size.factor = 8, min.cutoff = 0.5)  + scale_fill_gradientn(colours = colfunc(6))
SpatialFeaturePlot(FMHday45, features = c('PAX5'), pt.size.factor = 8, min.cutoff = 0.5)  + scale_fill_gradientn(colours = colfunc(6))

#Hindbrain, Pick the top spatial variable genes
SpatialFeaturePlot(FMHday45, features = c('HOXA6'), pt.size.factor = 8, min.cutoff = 0.5)  + scale_fill_gradientn(colours = colfunc(6))
SpatialFeaturePlot(FMHday45, features = c('HOXB9'), pt.size.factor = 8, min.cutoff = 0.5)  + scale_fill_gradientn(colours = colfunc(6))
SpatialFeaturePlot(FMHday45, features = c('IRX4'), pt.size.factor = 8, min.cutoff = 0.5)  + scale_fill_gradientn(colours = colfunc(6))

#Progenitor
SpatialFeaturePlot(FMHday45, features = c('OTX2'), pt.size.factor = 8, min.cutoff = 0.5)  + scale_fill_gradientn(colours = colfunc(6))
SpatialFeaturePlot(FMHday45, features = c('HES1'), pt.size.factor = 8, min.cutoff = 0.5)  + scale_fill_gradientn(colours = colfunc(6))
SpatialFeaturePlot(FMHday45, features = c('GLI3'), pt.size.factor = 8, min.cutoff = 0.5)  + scale_fill_gradientn(colours = colfunc(6))

# Dot plot, remove "Unkown" cell type
FMHday45T <- subset(FMHday45, idents = c('Forebrain', 'Hindbrain'))
MarkerGenesFMH45<-c('SST', 'DLX2','SOX9','PTN','GSX2', 'NEUROD2','SLC2A1','SLC1A3','OTX2',
                    'PAX5','LMX1B','STMN2','PAX8','EN1','EN2',
                    'IRX4','HOXA6', 'HOXB6', 'HOXB8','HOXB9','HOXC6')
DotPlot(FMHday45T, features =MarkerGenesFMH45)+ RotatedAxis()                    
#===================================================================================================================                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
