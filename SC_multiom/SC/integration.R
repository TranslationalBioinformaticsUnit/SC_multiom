##################################`
# 06/10/2020
#INTEGRATION HEALTHY AND FANCONI##
# MIREN LASAGA
##################################
#, lib.loc="/home/mlagasag/R/x86_64-pc-linux-gnu-library/4.0")
#library(ggpubr,lib.loc="/home/mlagasag/R/x86_64-pc-linux-gnu-library/4.0")
#, lib.loc="/home/mlagasag/R/x86_64-pc-linux-gnu-library/4.0/")

#loadNamespace(Seurat)

library(Rcpp, lib.loc="/home/mlagasag/R/x86_64-pc-linux-gnu-library/4.0")

#library(RcppArmadillo, lib.loc="/home/mlagasag/R/x86_64-pc-linux-gnu-library/4.0")

library(Seurat, lib.loc="/home/mlagasag/R/x86_64-pc-linux-gnu-library/4.0")



#library(ggplot2)



#>> Work directory and seed

setwd("/home/mlagasag/F_Riney/SC/integration_all_multiom/")
set.seed(1234567)



fresh<-readRDS("/home/mlagasag/F_Riney/SC/fresh_multiom/data_filtered.rds")
fresh$Sample<-"fresh"

MM<-readRDS("/home/mlagasag/F_Riney/SC/MM_HR_295303/data_filtered.rds")
MM$Sample<-"MM_HR_295303"

MM_90<-readRDS("/home/mlagasag/F_Riney/SC/MM_HR_907052/data_filtered.rds")
MM_90$Sample<-"MM_HR_907052"


MM_no_HR<-readRDS("/home/mlagasag/F_Riney/SC/MM_no_HR_937961/data_filtered.rds")
MM_no_HR$Sample<-"MM_no_HR_937961"


smouldering<-readRDS("/home/mlagasag/F_Riney/SC/smouldering_594489/data_filtered.rds")
smouldering$Sample<-"smouldering_594489"


#>> INTEGRATION
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 10000 * 1024^2)

combined.list <- c(fresh,MM,MM_90,MM_no_HR,smouldering)

for (i in 1:length(combined.list)) {
  combined.list[[i]] <- SCTransform(combined.list[[i]], verbose = FALSE)
}


#save.image("Integration/SCT.Rdata")
fanconi.features <- SelectIntegrationFeatures(object.list = combined.list, nfeatures = 2000)
#
#
combined.list <- PrepSCTIntegration(object.list = combined.list, anchor.features = fanconi.features, 
                                    verbose = FALSE)

combined.anchors <- FindIntegrationAnchors(object.list = combined.list, normalization.method = "SCT", 
                                           anchor.features = fanconi.features, verbose = FALSE)

#save.image("Integration/SCT.Rdata")
fanconi.integrated <- IntegrateData(anchorset = combined.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)
#                                     
saveRDS(fanconi.integrated, file="integration.rds")

