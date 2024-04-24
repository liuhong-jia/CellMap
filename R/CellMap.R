
##########################################################################################################
#' @title CellMap
#' @param st.data st data Seurat object.
#' @param coord.df The spatial coordinates of spatial transcriptomics images.
#' @param ref.expr Single-cell gene expression profile.
#' @param ref.anno Cell type information of gene expression profile, corresponding to the above `ref.expr`.
#' @param factor.size Factor size for scaling the weight of gene expression. Default: 0.1.
#' @param seed.num Number of seed genes of each cell type for recognizing candidate markers. Default: 10.
#' @param pvalue.cut Threshold for filtering cell type marker genes. Default: 0.1
#' @param knn The number of nearest neighboring single cells for each spot. Default: 0.1
#' @param verbose Show running messages or not. Default: TRUE.
#' @export 
#' @example	
##########################################################################################################

#' @example scMap <- CellMap(st.data, coord.df, ref.expr, ref.anno,factor.size = 0.1,seed.num = 10,pvalue.num = 0.1,verbose = TRUE)

CellMap <- function(st.data = st.data,
                       coord.df = coord.df,
                       ref.expr = ref.expr,
                       ref.anno = ref.anno,
                       factor.size = 0.1,
                       seed.num = 10,
                       pvalue.cut = 0.1,
                       knn =5,
                       verbose = TRUE)
{
  st.data.counts <- st.data@assays$Spatial@counts %>% as.matrix
  st.data.matrix <- getNormalizeData(st.data.counts)
  #gene.anno <- readRDS(system.file('data', 'mouse.gene.anno.rds', package = 'CellMap'))
  gene.anno <- readRDS(system.file('data', 'human.gene.anno.rds', package = 'CellMap'))
  sc.data <- createSeuratObjAndQC(ref.expr, ref.anno)
  sc.data <- sc.data[which(rownames(sc.data) %in% gene.anno[, 'Symbol']), ]
  sc.data <- sc.data %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30)
  print('[INFO] Searching candidate marker genes...',verbose = verbose)
  avg.expr.ref <- AverageExpression(sc.data, slot = "data")$RNA
  seed.genes <- identSeedGenes(avg.expr.ref, factor.size, count = seed.num)
  ref.markers <- searchMarkersByCorr(
    GetAssayData(sc.data) %>% as.matrix, 
    seed.genes, 
    scale.data = TRUE,
    p.cut = pvalue.cut,
    seed.num = seed.num,
    tcga.data.u = NULL
  )
  
  markers <- buildSigMatrix(
    ref.markers, 
    avg.expr.ref, 
    min.size = 10, 
    max.size = 200,
  )
  
  print('[INFO] Estimate the number of single cells in the spot',verbose = verbose)
  num.cells <- getCellNumber(st.data.matrix,mean.cell.num = 5)
 
  print('[INFO] Integrate single-cell and spatial spot data',verbose = verbose)
  sc.st.obj <- getInterdata(sc.data, st.data,coord = c("imagerow","imagecol"),markers)
  
  dist <- getDistMatrix(sc.st.obj)
  #nearCells <- getUniqueCells(dist, k = 3)
  
  nearCells <- getCells(dist,knn = 5)
  
  print('[INFO] Train a random forest model and predict',verbose = verbose)
  
  pred <- trainModel(st.data,sc.data,nearCells,markers)
  sim <- getSimMatrix(pred,st.data,sc.data,markers)
  
  sim.cells <- linearAllocation(sim,num.cells)
  
  print('[INFO] Map single cells onto spatial spots',verbose = verbose)
  cell.coords <- getRandomCoords(coord.df,sim.cells)
  mapping <- mapSctoSpatial(st.data,sc.data,sim.cells,coord.df)
  
  print('[INFO] Construct Seurat object',verbose = verbose)
  sc.out <- createSeuratObj(st.data,sc.data, mapping)
  return(sc.out = sc.out)
}
