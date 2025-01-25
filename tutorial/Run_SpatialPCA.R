#-------------------------------------------------------------------------------
# This tutorial is for analyzing mouse mPFC data using SpatialPCA
#-------------------------------------------------------------------------------

library(SpatialPCA)

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
pcs = seq(10,55,5)
for(cluster_number in seq(10,20,1)) {
  for(i in names(cnts)){
    cnt = cnts[[i]]
    xy = xys[[i]]
    colnames(xy) = c('row','col')
    
    ST_raw = CreateSpatialPCAObject(counts=as.matrix(cnt), location=as.matrix(xy), 
                                    project = "SpatialPCA",gene.type="spatial",
                                    sparkversion="spark", gene.number=3000,
                                    customGenelist=NULL,min.loctions = 5, min.features=5)
    
    df = data.frame(rep(1,dim(xy)[1]))
    for (pc in pcs) {
      # Estimate spatial PCs
      print('******now start buildkernel!******')
      ST = SpatialPCA_buildKernel(ST_raw, bandwidthtype="SJ")
      ST = SpatialPCA_EstimateLoading(ST,fast=FALSE,SpatialPCnum=pc)
      ST = SpatialPCA_SpatialPCs(ST, fast=FALSE)
      
      # Detect spatial domains
      clusterlabel= walktrap_clustering(cluster_number, ST@SpatialPCs,round(sqrt(dim(ST@location)[1])))
      clusterlabel_refine=refine_cluster_10x(clusterlabel,ST@location,shape="square")
      
      df[paste0('SpatialPCA_',pc)] = clusterlabel_refine
    }
    df = df[,-1]
    write.csv(df,file = paste0('../result/SpatialPCA/',i,'_k=',cluster_number,'.csv'))
  }
}
