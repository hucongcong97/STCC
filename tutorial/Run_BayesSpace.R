library(BayesSpace)
library(SingleCellExperiment)

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

seeds = 1:10
for (cluster_number in seq(10,20,1)) { 
  for(i in names(cnts)){
    cnt = cnts[[i]]
    xy = xys[[i]]
    colnames(xy) = c('row','col')
    
    sce <- SingleCellExperiment(assays=list(counts=as(cnt, "dgCMatrix")),
                                colData=xy)
    
    df = data.frame(rep(1,dim(xy)[1]))
    for (seed in seeds) {
      # Pre-processing data
      set.seed(seed)
      sce <- spatialPreprocess(sce, platform="Visium", 
                               n.PCs=20, n.HVGs=2000, log.normalize=T)
      
      # Clustering with BayesSpace
      set.seed(seed)
      sce = spatialCluster(sce, q=cluster_number, platform="Visium", d=20,
                           init.method="mclust", model="t", gamma=2,
                           nrep=10000, burn.in=100,
                           save.chain=TRUE)
      # q:The number of clusters
      # d:Number of top principal components to use when clustering
      result = colData(sce)
      
      df[paste0('BayesSpace_',seed)] = result$spatial.cluster
    }
    df = df[,-1]
    write.csv(df,file = paste0('../result/BayesSpace\\',i,'_k=',cluster_number,'.csv'))
  }
}
