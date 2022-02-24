###############################
######
###### QC of multiom SC data###
######
###############################


library(ggplot2)



setwd("/home/mlagasag/F_Riney/SC/merge_normal_multiome/")

metadata<-read.table("ATAC_MM.txt", row.names=NULL)

metadata$log10GenesPerUMI <- log10(metadata$nFeature_RNA) / log10(metadata$nCount_RNA)

my_data_order <- metadata[order(metadata$nFeature_RNA, decreasing=FALSE),]
my_data_order$rank <- c(1:nrow(my_data_order))
my_data_order$logumi <-log10(my_data_order$nCount_RNA)
for(i in 1:dim(filtered)[1]){
 filtered$percentage[i] <- filtered$nFeature_RNA[i]/filtered$nCount_RNA[i]
}
for(i in 1:dim(my_data_order)[1]){
 my_data_order$gene_peak_ratio[i] <- my_data_order$nFeature_RNA[i]/my_data_order$nFeature_ATAC[i]
}

metadata$Sample<-factor(metadata$Sample, levels=c("fresh_multiom", "MM_HR_295303", "MM_HR_907052", "MM_HR_907052_RE", "MM_no_HR_937961", "smouldering_594489"))


my_data_order$Sample<-factor(my_data_order$Sample,levels=c("fresh_multiom", "MM_HR_295303", "MM_HR_907052", "MM_HR_907052_RE", "MM_no_HR_937961", "smouldering_594489"))
# histogram ratio gene/counts per sample


pdf("histogram_percentage_FIL.pdf", width=12)
ggplot(filtered, aes(x=percentage, fill=Sample)) +
  geom_histogram( color="black", breaks=seq(0, 1, by = 0.025)) +facet_wrap(factor(filtered$Sample)) + theme_bw()
dev.off()
# boxplot of percentage

pdf("boxplot_percentage_FIL.pdf", width=10)  
ggplot(filtered, aes(y=percentage, x=Sample, fill=Sample)) +
  geom_boxplot( color="black")
dev.off()


# density of percentage gene/count per sample

pdf("density_FRIP_FIL.pdf") 
ggplot(filtered, aes(x = FRIP, y = Sample, fill = Sample)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
dev.off()


# scater plot ncounts and nfeatures


pdf("features_counts_all_samples_FIL.pdf")
p <- ggplot(filtered, aes(nFeature_RNA, nCount_RNA),colour = factor(Sample))+
 geom_point(aes(colour = factor(Sample)), alpha=0.2)
 
p+theme_bw()

dev.off()

# scater plot ncounts and nfeatures split by sample


pdf("features_counts_split_sample_FIL.pdf", width=10)
ggplot(mm@meta.data, aes(x=nFeature_RNA, y=nCount_RNA)) +
  geom_point(aes(colour = Sample)) +facet_wrap(as.factor(mm@meta.data$Sample))+
  scale_fill_gradient()
 

dev.off()

pdf("FRIP_peak_counts_FIL.pdf", width=13)
ggplot(filtered, aes(x=nFeature_ATAC, y=nCount_ATAC)) +
  geom_point(aes(colour = FRIP)) +facet_grid(filtered$FRIP_group~filtered$Sample)+theme_bw()+
  scale_fill_gradient()
 

dev.off()
# scater plot nfeatures atac and RNA
pdf("features_peaks_split.pdf", width=10)
p <- ggplot(metadata, aes(nFeature_ATAC, nFeature_RNA),colour = factor(Sample))

p + geom_point(aes(colour = factor(Sample)), alpha=0.1)

dev.off()

pdf("features_peaks_split_sample_RE.pdf", width=10)
p <- ggplot(metadata, aes(nFeature_ATAC, nFeature_RNA),colour = factor(Sample))

p + geom_point(aes(colour = factor(Sample)))+facet_wrap(as.factor(metadata$Sample)) +theme_bw()

dev.off()

pdf("features_peaks_split_sample_percentage.pdf", width=10)
p <- ggplot(my_data_order, aes(nFeature_ATAC, nFeature_RNA),colour = factor(percentage))

p + geom_point(aes(colour = percentage)) +facet_wrap(as.factor(my_data_order$Sample)) +theme_bw()+
  scale_color_gradient(low = "grey", high = "red")

dev.off()


for (var in unique(my_data_order$Sample)) {
    dev.new()
    print( ggplot(my_data_order[my_data_order$Sample==var,], aes(nFeature_ATAC, nFeature_RNA)) + geom_point(aes(colour = percentage)) +scale_color_gradient(low = "grey", high = "red") )
}


# density counts
pdf("density_percent_mt_FIL.pdf")
ggplot(filtered, aes(x = percent.mt, y = Sample, fill = Sample)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
dev.off()

 # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
pdf("density_log10genesumi_FIL.pdf")
ggplot(filtered, aes(x = log10GenesPerUMI, y = Sample, fill = Sample)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
dev.off()

 

  # boxplot

p<-ggplot(metadata, aes(x=Sample, y=nFeature_RNA, fill=Sample)) +
  geom_boxplot()+
boxplot_framework(upper_limit = 6000)
 # density genes

ggplot(metadata, aes(nCount_RNA, colour = Sample)) +
  geom_density()         
 
 
 # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI


  	ggplot(my_data_order, aes(x=log10GenesPerUMI, color = Sample, fill=Sample)) +
  	geom_density() +
  	theme_classic() +
  	geom_vline(xintercept = 0.8) 
  

  	ggplot(my_data_order,aes(x=nCount_RNA, y=nFeature_RNA, color=percentage)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	#geom_vline(xintercept = 500) +
  	#geom_hline(yintercept = 250) +
  	facet_wrap(~Sample)    
    
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
 
 metadata$sel<-0
 for (j in 1:length(min)){
 for (i in 1:dim(metadata)[1]){
   if(metadata$nCount_RNA[i]>min[j] && metadata$nCount_RNA[i]<max[j] && metadata$nFeature_RNA[i]>conf[[j]][1] && metadata$nFeature_RNA[i]<conf[[j]][2]){
   metadata$sel[i]<- 1 }
   }
 }


## number of cells by selection
 pdf("barplot_number_cells_100_RE.pdf", width=10)
ggplot(metadata, aes(x=Sample,fill=factor(sel))) +
  geom_bar(position="dodge") 

dev.off()

# scater plot ncounts and nfeatures split by sample


pdf("features_count_SEL.pdf", width=10)
ggplot(metadata, aes(x=nFeature_RNA, y=nCount_RNA)) +
  geom_point(aes(colour = factor(sel))) +facet_wrap(as.factor(metadata$Sample)) + theme_bw()


dev.off()

 


