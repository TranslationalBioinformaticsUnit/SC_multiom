
library(Rcpp, lib.loc="/home/mlagasag/R/x86_64-pc-linux-gnu-library/4.0")

#library(RcppArmadillo, lib.loc="/home/mlagasag/R/x86_64-pc-linux-gnu-library/4.0")

library(Seurat, lib.loc="/home/mlagasag/R/x86_64-pc-linux-gnu-library/4.0")


library(harmony, lib.loc="/home/mlagasag/R/x86_64-pc-linux-gnu-library/4.0")

setwd("/home/mlagasag/F_Riney/SC/integration_all_multiom/")

mm<-readRDS("integrated_def.rds")

mm<-RunHarmony(mm, "Sample", assay.use="integrated")

saveRDS(mm, "harmnoy.rds")