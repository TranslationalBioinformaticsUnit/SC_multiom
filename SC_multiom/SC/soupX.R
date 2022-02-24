library(Seurat)#, lib.loc="/home/mlagasag/R/x86_64-pc-linux-gnu-library/4.0/")
#library(SoupX)
#library(Signac)


setwd("/home/mlagasag/F_Riney/SC/MM_HR_907052/")
mm<-readRDS("soupx_0.2.rds")

 mm<-NormalizeData(mm)
 mm <- ScaleData(mm, features = rownames(mm))
 
 mm<-FindVariableFeatures(mm)
 
 ## Cell Cycle Analysis
 mm<-CellCycleScoring(mm,s.features = cc.genes$s.genes,g2m.features = cc.genes$g2m.genes)
 
 # Using PCA, we see if there is an effect caussed due to cell cycle phase
 mm<-RunPCA(mm,features = c(cc.genes$s.genes,cc.genes$g2m.genes))
 DimPlot(mm,group.by = "Phase")
 mm$CC.Difference <- mm$S.Score - mm$G2M.Score

 
 # SC Transform (https://satijalab.org/seurat/v3.1/sctransform_vignette.html)
 mm<-SCTransform(mm,vars.to.regress = c("S.Score","G2M.Score","percent.mt","nFeature_RNA"))
 
 
 # PCA
 mm <- RunPCA(mm, npcs = 50, verbose = FALSE)
 P1<-DimPlot(mm,group.by = "Phase")
 P2<-DimPlot(mm)
 P1+P2
 
 # Select components to use in further steps
 ElbowPlot(mm,ndims = 50)
 JackStrawPlot(mm, dims = 1:15)
 ## UMAP
 mm <- RunUMAP(mm, reduction = "pca", dims = 1:20)
 
 # Clustering
 mm <- FindNeighbors(mm, reduction = "pca", dims = 1:20)
 mm <- FindClusters(mm, resolution = 0.8, n.start = 1000) 
 
pdf("UMAP_clusters_soupX0.2.pdf")
DimPlot(mm, group.by="seurat_clusters")
dev.off

saveRDS(mm, "soupx0.2_processed.rds")