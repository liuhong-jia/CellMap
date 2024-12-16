
##########################################################################################################
#' @title CellMap
#' @param st.obj Seurat object of spatial transcriptome data.
#' @param sc.obj Seurat object of scRNA-seq data.
#' @param coord Coordinates column names in ST images slot.coord = c("x","y") or coord = c("imagerow","imagecol").
#' @param celltype.column The column name for cell type in the single-cell Seurat object, with the default value as "idents".
#' @param sc.sub.size Downsampling proportion or number for scRNA-seq data. Default: NULL.
#' @param min.sc.cells The minimum number of cell types in scRNA-seq data.Default: 50
#' @param factor.size Factor size for scaling the weight of gene expression. Default: 0.1.
#' @param seed.num Number of seed genes of each cell type for recognizing candidate markers. Default: 10.
#' @param pvalue.cut Threshold for filtering cell type marker genes. Default: 0.1.
#' @param knn The number of nearest neighboring single cells for each spot. Set to 5 for low-resolution data and 1 for high-resolution data.
#' @param mean.cell.num The average number of single cells in the spot.Set to 5 for low-resolution data and 1 for high-resolution data.
#' @param max.cell.num The maximum number of cells within each spot, if equal to 1, indicates that each spot contains only a single cell.
#' @param n.workers Number of cores to be used for parallel processing. Default: 4.
#' @param verbose Show running messages or not. Default: TRUE.
#' @export Seurat object of assignment result. 
##########################################################################################################

							  


CellMap <- function(st.obj = st.obj,
                    sc.obj = sc.obj,
					coord = c("x","y"),
					celltype.column = "idents",
					sc.sub.size = NULL,
					min.sc.cell = 50,
                    factor.size = 0.1,
                    seed.num = 10,
                    pvalue.cut = 0.1,
                    knn = 5,
					mean.cell.num = 5,
					max.cell.num = 10,
					n.workers = 4,
                    verbose = TRUE)
{ 
	checkInputParams(st.obj, sc.obj, coord, celltype.column, sc.sub.size, min.sc.cell,factor.size, seed.num, pvalue.cut, knn, mean.cell.num, max.cell.num,n.workers, verbose)
	st.obj <- processSpatialData(st.obj,resolution = resolution)
	sc.obj <- processScData(sc.obj,celltype.column = "idents")
	genes <- intersect(rownames(st.obj),rownames(sc.obj))
  st.data.counts <- GetAssayData(st.obj,slot = "counts")
	
	print('[INFO] Searching candidate marker genes...',verbose = verbose)
	avg.expr.ref <- AverageExpression(sc.obj[genes,], assay ="SCT")$SCT 
	seed.genes <- identSeedGenes(avg.expr.ref, factor.size, count = seed.num)
  
	ref.markers <- searchMarkersByCorr(
		GetAssayData(sc.obj[genes,]) %>% as.matrix, 
		seed.genes, 
		scale.data = TRUE,
		p.cut = pvalue.cut,
		seed.num = seed.num,
		n.workers = 4
	)
  
	markers <- buildSigMatrix(
		ref.markers, 
		avg.expr.ref, 
		min.size = 10, 
		max.size = 200,
	)
	print('[INFO] Estimate the number of single cells in the spot',verbose = verbose)
	
	if (max.cell.num == 1) {
		num.cells <- rep(1, ncol(st.data.counts))
		names(num.cells) <- colnames(st.data.counts)
    } else {
		num.cells <- getCellNumber(st.data.counts, mean.cell.num = mean.cell.num)
    } 
	
	print('[INFO] Integrate single-cell and spatial spot data',verbose = verbose)
  
	sc.st.obj <- getInterdata(sc.obj, st.obj,coord = coord ,markers)
  
	dist <- getDistMatrix(sc.st.obj)
	nearCells <- getCells(dist,k = knn)
  
	print('[INFO] Train a random forest model and predict,waiting...',verbose = verbose)


	if (ncol(sc.obj) > ncol(st.obj)) {
    pred <- trainModel(st.obj, sc.obj, nearCells, markers)
    sim <- getSimMatrix(pred$prediction, st.obj, sc.obj, markers)
    } else {
    cell.to.cluster <- getCellClusterInfo(nearCells, st.obj)
    sim <- getSimMatrix(cell.to.cluster, st.obj, sc.obj, markers)
       }
    print('[INFO] Assign single cells to spatial spots,waiting...',verbose = verbose)
	sim.cells <- linearAllocation(sim,num.cells)
	
	print('[INFO] Map single cells onto spatial spots',verbose = verbose)
	
	cell.coords <- getRandomCoords(st.obj,sim.cells,n.workers = 4)
	
	mapping <- mapSctoSpatial(st.obj, sc.obj, sim.cells)  
	
	print('[INFO] Construct Seurat object',verbose = verbose)
	sc.out <- createSeuratObj(st.obj,sc.obj, mapping)
  
	##deconve
	metadata <- sc.out@meta.data[,c("CellType","SpotName")]
	type.counts <- metadata %>% count(SpotName, CellType) 

	celltype.prop <- type.counts %>% 
	group_by(SpotName) %>% 
		mutate(Proportion = n / sum(n)) %>% 
		ungroup() %>% 
		select(-n) %>% 
		tidyr::spread(key = CellType, value = Proportion, fill = 0) %>%
		as.data.frame
  
	rownames(celltype.prop) <- celltype.prop$SpotName 
	celltype.prop <- celltype.prop[,-1]
  
	print('[INFO] Finish!',verbose = verbose)
  
	return(list(sc.out = sc.out,decon = celltype.prop))
}
