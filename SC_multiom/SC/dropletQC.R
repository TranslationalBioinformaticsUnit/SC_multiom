####### removing ambient RNA########

#library(Seurat)
#library(Signac)
#library(ggplot2)
#library(DropletQC)
library(diem, lib.loc="/home/mlagasag/R/x86_64-pc-linux-gnu-library/4.0")

#------SET WORK DIRECTORY AND LOAD DATA

setwd("/home/mlagasag/F_Riney/SC/smouldering_594489/")
#mm<-readRDS("data.rds")
#
##------Calculating the nuclear fraction
#
#nf1 <- nuclear_fraction_tags(
#  outs = "/mnt/resultados/oncohematologia/resultados/mlagasag/smouldering_594489_counts/outs/",
#  tiles = 1, cores = 1, verbose = FALSE)
#head(nf1)
#
## Get data frame with the nuclear fraction in the first column and umi counts in
## the second
#
#
#data<-cbind(mm@meta.data,nf1)
#nf_umi <- data.frame(nf=data$nuclear_fraction,
#                     umi=data$nCount_RNA)
#
#
#gbm.ed <- identify_empty_drops(nf_umi=nf_umi)
#head(gbm.ed)
#table(gbm.ed$cell_status)
#
#gbm.ed$celltype<-"MM"
#
#gbm.ed.dc <- identify_damaged_cells(gbm.ed, verbose = FALSE, output_plots = TRUE)
#
#table(gbm.ed.dc [[1]]$cell_status)
#
#pdf("dropletQC.pdf")
#wrap_plots(gbm.ed.dc[[2]], nrow = 1)
#dev.off()

################################
####DIEM##################
###########################

counts <- read_10x("/datos/intercambio/mlagasag/smouldering_594489_counts/outs/raw_feature_bc_matrix") 
counts<-counts[c(1:36600),]
mb_small <- create_SCE(counts, name="prueba")
dim(mb_small)
class(mb_small)

drop_data <- droplet_data(mb_small)
head(drop_data)
summary(drop_data)


mt_genes <- grep(pattern = "^mt-", x = rownames(mb_small@gene_data),
                 ignore.case = TRUE, value = TRUE)
mb_small <- get_gene_pct(x = mb_small, genes = mt_genes, name = "pct.mt")
genes <- grep(pattern = "^malat1$", x = rownames(mb_small@gene_data),
              ignore.case = TRUE, value = TRUE)
mb_small <- get_gene_pct(x = mb_small, genes = genes, name = "MALAT1")
drop_data <- droplet_data(mb_small, type="all")
summary(drop_data)


require(gridExtra)
p1 <- plot_data(mb_small, feat_x = "total_counts", feat_y = "n_genes", 
                log_x = TRUE, log_y = TRUE, ret = TRUE, data_type = "all")
p2 <- plot_data(mb_small, feat_x = "n_genes", feat_y = "pct.mt", 
                log_x = TRUE, log_y = FALSE, ret = TRUE, data_type = "all")
p3 <- plot_data(mb_small, feat_x = "n_genes", feat_y = "MALAT1", 
                log_x = TRUE, log_y = FALSE, ret = TRUE, data_type = "all")
p4 <- plot_data(mb_small, feat_x = "pct.mt", feat_y = "MALAT1", 
                log_x = FALSE, log_y = FALSE, ret = TRUE, data_type = "all")
pdf("info_QC.pdf")
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()
barcode_rank_plot(mb_small, title = "MouseBrain")

mb_small <- set_debris_test_set(mb_small, min_counts = 500)
length(mb_small@test_set)
length(mb_small@bg_set)

mb_small <- filter_genes(mb_small, cpm_thresh = 0)
genes <- gene_data(mb_small)
summary(genes)

mb_small <- get_pcs(mb_small,
                    n_var_genes = 50,
                    n_pcs = 10)
                    
mb_small <- init(mb_small,
                 k_init = 20,
                 nstart_init = 30,
                 min_size_init = 2) 
                 
                 
mb_small <- assign_clusters(mb_small)
drop_data <- droplet_data(mb_small, type="test")
table(drop_data[,"Cluster"])
tapply(drop_data[,"n_genes"],
       drop_data[,"Cluster"],
       mean)
tapply(drop_data[,"pct.mt"],
       drop_data[,"Cluster"],
       mean)

mb_small <- run_em(mb_small)
mb_small <- assign_clusters(mb_small)
drop_data <- droplet_data(mb_small, type="test")
table(drop_data[,"Cluster"])
tapply(drop_data[,"n_genes"],
       drop_data[,"Cluster"],
       mean)
tapply(drop_data[,"pct.mt"],
       drop_data[,"Cluster"],
       mean)
tapply(drop_data[,"MALAT1"],
       drop_data[,"Cluster"],
       mean)
       
       mb_small <- estimate_dbr_score(mb_small,
                               thresh_genes = 50,
                               thresh_p = 0.5)
de_genes <- debris_genes(mb_small, p_adj = 0.5)
head(de_genes)


p1 <- plot_clust(mb_small, feat_x = "n_genes", feat_y = "score.debris", 
                 log_x = TRUE, ret = TRUE)
p2 <- plot_clust(mb_small, feat_x = "total_counts", feat_y = "score.debris",
                 log_x = TRUE, ret = TRUE)
p3 <- plot_clust(mb_small, feat_x = "pct.mt", feat_y = "score.debris", 
                 log_x = FALSE, ret = TRUE)
p4 <- plot_clust(mb_small, feat_x = "MALAT1", feat_y = "score.debris", 
                 log_x = FALSE, ret = TRUE)
grid.arrange(p1, p2, p3, p4, ncol = 2)

sm <- summarize_clusters(mb_small, top_n = 20)
par(mfrow=c(1,2))
plot(sm[,"avg_n_genes"], sm[,"avg_dbr_score"], pch= NA,
     xlab="Average number Genes", ylab="Avergae debris score")
text(sm[,"avg_n_genes"], sm[,"avg_dbr_score"], sm[,"Cluster"])
plot(sm[,"avg_n_counts"], sm[,"avg_dbr_score"], pch= NA,
     xlab="Average total counts", ylab="Avergae debris score")
text(sm[,"avg_n_counts"], sm[,"avg_dbr_score"], sm[,"Cluster"])


clust <- clusters(mb_small)
counts <- raw_counts(mb_small)
deg <- de_ttest_all(counts = counts, labels = clust)
head(deg)


mb_small <- call_targets(mb_small,
                         thresh_score = 0.5,
                         min_genes = 0)

drop_data <- droplet_data(mb_small, type="test")
tapply(drop_data[,"n_genes"],
       drop_data[,"Call"],
       mean)
tapply(drop_data[,"score.debris"],
       drop_data[,"Call"],
       mean)
tapply(drop_data[,"pct.mt"],
       drop_data[,"Call"],
       mean)
tapply(drop_data[,"MALAT1"],
       drop_data[,"Call"],
       mean)
tapply(drop_data[,"Cluster"],
       drop_data[,"Call"],
       table)

mb_small <- call_targets(mb_small,
                         clusters = "debris",
                         thresh_score = NULL,
                         min_genes = 0)
drop_data <- droplet_data(mb_small)
tapply(drop_data[,"n_genes"],
       drop_data[,"Call"],
       mean)
tapply(drop_data[,"score.debris"],
       drop_data[,"Call"],
       mean)
tapply(drop_data[,"pct.mt"],
       drop_data[,"Call"],
       mean)
tapply(drop_data[,"MALAT1"],
       drop_data[,"Call"],
       mean)
tapply(drop_data[,"Cluster"],
       drop_data[,"Call"],
       table)   


seur <- convert_to_seurat(mb_small, min.features = 50)       

save.image("diem_filtered.Rda")


# Data matrix overrepresented genes
#for each cell identify the overrepresented genes (>2% of total counts)
data_no_zeros <-mm@assays$RNA@data[rowSums(mm@assays$RNA@data)!=0,]
prop_counts <- data_no_zeros
total_counts_per_cell <- colSums(data_no_zeros)
for(i in 1:ncol(data_no_zeros)){
  prop_counts[,i] <- data_no_zeros[,i]/total_counts_per_cell[i]
  cat(i,"-")
}
prop_counts_x_100 <- prop_counts*100
top_genes <- vector()
for(i in 1:ncol(prop_counts_x_100)){
  top_genes <- c(top_genes,rownames(prop_counts_x_100)[prop_counts_x_100[,i]>2])
}
unique_top_genes <- unique(top_genes)
pdf(paste(format(Sys.time(),"%Y%m%d"),"_prop_counts_heatmap_",sample_name,".pdf", sep=""))
Heatmap(prop_counts_x_100[unique_top_genes,],
        show_column_names=FALSE, row_names_gp = gpar(fontsize = 3),
        col=jdb_palette("flame_watts"))
dev.off()

my_data <- mm@meta.data
my_data_order <- my_data[order(my_data$nFeature_RNA, decreasing=FALSE),]
my_data_order$rank <- c(1:nrow(my_data_order))
my_data_order$logumi <-log10(my_data_order$nCount_RNA)
for(i in 1:dim(my_data_order)[1]){
 my_data_order$percentage[i] <- my_data_order$nFeature_RNA[i]/my_data_order$nCount_RNA[i]
}
ggplot(mm@meta.data, aes(x=nFeature_RNA, y=nFeature_ATAC)) + 
  geom_point(size=2) +
  scale_fill_value 
  geom_segment(aes(x=nFeature_RNA, 
                   xend=nFeature_RNA, 
                   y=nCount_RNA, 
                   yend=nCount_RNA)) + 
  labs(title="Lollipop Chart", 
       subtitle="Make Vs Avg. Mileage", 
       caption="source: mpg") + 
  theme(axis.text.x = element_text(angle=65, vjust=0.6))


# MT %
pdf("prediction_MT_ATAC.pdf")
p0 = DimPlot(object = mm, reduction = "umap.atac", group.by = "percent.mt")
ggplot(my_data_order, aes(x=rank, y=logumi)) +
  geom_point(aes(colour = percentage)) + 
  scale_fill_gradient()
dev.off()

my_counts  <- as.matrix(mm@assays$RNA@counts)
cumulative_gene_cell <- apply(my_counts, 1, function(c)sum(c!=0))
cells <- c(1:ncol(my_counts))
num_genes <- vector()
for(i in 1:length(cells)){
 num_genes[i] <- sum(cumulative_gene_cell==i)
}
pdf("histogram_cumulative_gene_cell.pdf")
hist(cumulative_gene_cell[cumulative_gene_cell!=0])
dev.off()

  , ylim=c(0,5000)
# Read CellRanger Output (.h5 files are faster to load)
mm<-Read10X("DATA/")

# Create Seurat Object with initial filtering
# A gene must be expressed in at  least 3 cells to be considered
# A cell must have at least 400 detected genes to be considered
mm <- CreateSeuratObject(counts = mm, project = "mm", min.cells = 3, min.features = 400)

