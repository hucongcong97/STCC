library(BASS)
library(Seurat)
library(tidyverse)

#-------------------------------------------------------------------------------
# 1.Mouse mPFC data by STARmap
# https://github.com/zhengli09/BASS-Analysis/blob/master/analysis/STARmap.Rmd
#-------------------------------------------------------------------------------
load('../data/starmap_mpfc.RData')

# a list of gene expression count matrices
cnts <- starmap_cnts 
# a list of spatial coordinates matrices
xys <- lapply(starmap_info, function(info.i){
  as.matrix(info.i[, c("x", "y")])
}) 

set.seed(2024)

C_list = seq(10,20,1)
for(C in C_list){
  # hyper-parameters
  # C <- 15 # number of cell types
  R <- 4 # number of spatial domains
  
  # Set up BASS object
  BASS <- createBASSObject(cnts, xys, C = C, R = R, beta_method = "SW")
  
  # Data pre-processing:
  # 1.Library size normalization followed with a log2 transformation
  # 2.Dimension reduction with PCA after standardizing all the genes
  # 3.Batch effect adjustment using the Harmony package
  
  m1 = n1 = matrix(nrow = dim(starmap_cnts[["20180417_BZ5_control"]])[2],ncol = 10)
  m2= n2 = matrix(nrow = dim(starmap_cnts[["20180419_BZ9_control"]])[2],ncol = 10)
  m3 = n3 = matrix(nrow = dim(starmap_cnts[["20180424_BZ14_control"]])[2],ncol = 10)
  
  pcs = seq(10,55,5)
  for(i in 1:10){
    print(paste0('K = ',C,',PC = ',pcs[i]))
    BASS <- BASS.preprocess(BASS, doLogNormalize = TRUE, 
                            doPCA = TRUE, scaleFeature = TRUE, nPC = pcs[i])
    
    # Run BASS algorithm
    BASS <- BASS.run(BASS)
    
    # Post-process posterior samples:
    # 1.Adjust for label switching with the ECR-1 algorithm
    # 2.Summarize the posterior samples to obtain the cell type labels, spatial 
    #   domain labels, and the cell type proportion matrix
    BASS <- BASS.postprocess(BASS)
    clabels <- BASS@results$c # cell type clusters
    zlabels <- BASS@results$z # spatial domain labels
    
    m1[,i] = clabels[[1]]
    m2[,i] = clabels[[2]]
    m3[,i] = clabels[[3]]
    n1[,i] = zlabels[[1]]
    n2[,i] = zlabels[[2]]
    n3[,i] = zlabels[[3]]
  }
  
  cdata = list("20180417_BZ5_control" = m1,
               "20180419_BZ9_control" = m2,
               "20180424_BZ14_control" = m3)
  
  # zdata = list("20180417_BZ5_control" = n1,
  #              "20180419_BZ9_control" = n2,
  #              "20180424_BZ14_control" = n3)
  
  save(cdata, file=paste0("../result/BASS/clabels_k=",C,".RData"))
  # save(zdata, file=paste0("../result/BASS/zlabels_k=",C,".RData"))
}

