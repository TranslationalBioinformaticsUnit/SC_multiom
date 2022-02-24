#.libPaths()
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(Seurat)#, lib.loc= "/home/mlagasag/R/x86_64-pc-linux-gnu-library/4.0")

#/mnt/resultados/oncohematologia/resultados/mlagasag/MM_no_HR_937961_counts/outs/filtered_feature_bc_matrix/


set.seed(1234567)
##set work directory
setwd("/home/mlagasag/F_Riney/SC/MM_no_HR_937961/")
#####
# CD34 Bone Marrow

# Read CellRanger Output (.h5 files are faster to load)
mm <- Read10X(data.dir = "/datos/intercambio/mlagasag/MM_no_HR_937961_counts/outs/filtered_feature_bc_matrix/")

metadata<-read.csv("/datos/intercambio/mlagasag/MM_no_HR_937961_counts/outs/per_barcode_metrics.csv", sep=",")

metadata<-metadata[metadata$is_cell==1,]

rownames(metadata)<-metadata$barcode

metadata$FRIP<-metadata$atac_peak_region_fragments/metadata$atac_fragments

metadata$FRIP_group[metadata$FRIP>0.7]<-"high_FRIP"

metadata$FRIP_group[metadata$FRIP<=0.7]<-"low_FRIP"
#
#
rna_counts <- mm$`Gene Expression`
atac_counts <- mm$Peaks
#
mm <- CreateSeuratObject(counts = rna_counts, meta.data=metadata)
mm[["percent.mt"]] <- PercentageFeatureSet(mm, pattern = "^MT-")
#
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb (ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
#
#
frag.file <- "/datos/intercambio/mlagasag/MM_no_HR_937961_counts/outs/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
mm[["ATAC"]] <- chrom_assay

fresh<-read.table("/home/mlagasag/F_Riney/SC/fresh_multiom/metadata_fresh.txt")

min<-seq(from = 0, to = 29000, by = 100)
max<-seq(from = 1000, to = 30000, by = 100)
empty_list<-vector(mode = "list", length = length(min))
conf <- vector(mode = "list", length = length(min))


 for (j in 1:length(min)){
 for (i in 1:dim(fresh)[1]){
   if(fresh$nCount_RNA[i]>min[j] && (fresh$nCount_RNA[i]<max[j])){
   empty_list[[j]] <- c(empty_list[[j]],fresh$nFeature_RNA[i])}
   }
 conf[[j]][1]<- quantile(empty_list[[j]], 0.025)
 conf[[j]][2]<- quantile(empty_list[[j]], 0.975)
 }
 
 mm@meta.data$sel<-0
 for (j in 1:length(min)){
 for (i in 1:dim(mm@meta.data)[1]){
   if(mm@meta.data$nCount_RNA[i]>min[j] && mm@meta.data$nCount_RNA[i]<max[j] && mm@meta.data$nFeature_RNA[i]>conf[[j]][1] && mm@meta.data$nFeature_RNA[i]<conf[[j]][2]){
   mm@meta.data$sel[i]<- 1 }
   }
 }
 
### FILTER DATA

mm <- subset(mm, subset= sel==1 & FRIP>0.7 & percent.mt<15 & nCount_RNA>quantile(mm$nCount_RNA, 0.1) & nCount_RNA<quantile(mm$nCount_RNA,0.9) & nCount_ATAC>quantile(mm$nCount_ATAC, 0.1) & nCount_ATAC<quantile(mm$nCount_ATAC,0.9) & nFeature_RNA>quantile(mm$nFeature_RNA, 0.1) & nFeature_RNA<quantile(mm$nFeature_RNA, 0.9) & nFeature_ATAC>quantile(mm$nFeature_ATAC, 0.1) & nFeature_ATAC<quantile(mm$nFeature_ATAC, 0.9))
 

#pdf("QC.pdf")
#VlnPlot(mm, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
#       log = TRUE, pt.size = 0) + NoLegend()
#dev.off()
#
## RNA analysis
DefaultAssay(mm) <- "RNA"
mm <- SCTransform(mm, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
#
#
DefaultAssay(mm) <- "ATAC"
mm <- RunTFIDF(mm)
mm <- FindTopFeatures(mm, min.cutoff = 'q0')
mm <- RunSVD(mm)
mm <- RunUMAP(mm, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

mm <- FindMultiModalNeighbors(mm, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
mm <- RunUMAP(mm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
mm <- FindClusters(mm, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

saveRDS(mm, "data_filtered.rds")
#
p1 <- DimPlot(mm, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(mm, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(mm, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
#
pdf("UMAP_fltered.pdf", width=15)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

#
## MT %
#pdf("prediction_MT_ATAC.pdf")
#p0 = DimPlot(object = mm, reduction = "umap.atac", group.by = "percent.mt")
#ggplot(p0$data, aes(p0$data$atacUMAP_1, p0$data$atacUMAP_2)) +
#  geom_point(aes(colour = mm@meta.data$percent.mt)) +
#  scale_colour_gradient2()
#dev.off()

#mm<-readRDS("data.rds")
#
#DefaultAssay(mm)<-"RNA"
#
#Fanconi.Markers<-FindAllMarkers(mm,only.pos = TRUE,min.pct = 0.25, min.diff.pct = 0.2, logfc.threshold = 0.25)
#Fanconi.Markers$Specific<-Fanconi.Markers$pct.1-Fanconi.Markers$pct.2
#
#write.table(Fanconi.Markers, "markers_cluster.txt")
