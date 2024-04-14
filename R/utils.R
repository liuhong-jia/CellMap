
##########################################################################################################
#' @title getNormalizeData
#' @description Data normalization.
#' @param data Data that needs to be normalized.
#' @return Normalized expression data.
#' @export
##########################################################################################################

getNormalizeData <- function(data) {
  data <- data * (10^4 / colSums(data, na.rm = TRUE))
  data <- log2(data + 1)
  data <- ifelse(is.nan(data), 0, data)
  return(data)
}


##########################################################################################################
#' @title searchMarkersBySeurat
#' @description Search marker genes of each cell type using Seurat method.
#' @param obj Seurat object.
#' @return A list of markers for cell types.
#' @export

##########################################################################################################

searchMarkersBySeurat <- function(obj, ...) {
  options(future.globals.maxSize = 40000 * 1024^2)
  future::plan("multiprocess", workers = parallel::detectCores())
  deg.markers <- FindAllMarkers(obj, only.pos = TRUE, ...)	
  return(split(deg.markers$gene, deg.markers$cluster))
}

##########################################################################################################
#' @title downSamplSeurat
#' @description Sample single cells by cell type.
#' @param obj Seurat object.
#' @param cnt Sample size for each ident, default: 200.
#' @param seed Randon seed, default: 123.
#' @return Subset of seurat object.
#' @export

##########################################################################################################

downSamplSeurat <- function(obj, cluster.col = NULL, cnt = 200, seed = 123, percent = 0.3) {
  set.seed(seed)
  if (!is.null(cluster.col)) Idents(obj) <- obj@meta.data[, cluster.col]
  cells <- Idents(obj) %>% table
  sub.cells <- sapply(names(cells), function(xx) {
    sub.cells <- Idents(obj)[Idents(obj) == xx] %>% names
    cnt <- ifelse(is.null(percent), cnt, length(sub.cells) * percent)
    if (length(sub.cells) > cnt) sub.cells <- sample(sub.cells, cnt, replace = FALSE)
    return(sub.cells)
  }) %>% unlist(use.names = F)
  subset(obj, cells = sub.cells)
}


##########################################################################################################
#' @title getCellNumber
#' @description Estimate the number of cells in each spot using linear regression.
#' @param ST.data The ST expression data has genes as rows and spots as columns.
#' @param mean.cell.num Mean cell number of spot,for 10X visium data,default:5.
#' @return The number of cells of each spot.
#' @export
#' @example
##########################################################################################################

#' @example num.cells <- getCellNumber(ST.data,mean.cell.num = 5)

getCellNumber <- function(ST.data, mean.cell.num = 5){
  data.filtered <- ST.data[rowSums(ST.data) > 0, ]
  gene.var <- apply(data.filtered, 1, var)
  stable.genes <- rownames(data.filtered)[gene.var <= 0.5]
  spot.total <- colSums(ST.data[stable.genes,])
  spot.mean <- sum(ST.data[stable.genes,])/dim(ST.data)[2]
  spot.cell <- round((spot.total/spot.mean)* mean.cell.num) 
  spot.cell[spot.cell == 0] <- 1
  return(spot.cell)
}


##########################################################################################################
#' @title getELDist
#' @description Calculate the Euclidean distance between single cells and spots.
#' @param scRNA.data scRNA-seq data of sampled cells.
#' @param ST.data ST data.
#' @return Euclidean distance of single cells and spots.
#' @export
#' @example
##########################################################################################################

#' @example dist.matrix <- getELDist(scRNA.data,ST.data)

getELDist <- function(scRNA.data, ST.data) {
  over.genes <- intersect(rownames(scRNA.data), rownames(ST.data))
  ST.data <- ST.data[over.genes, , drop = FALSE]
  scRNA.data <- scRNA.data[over.genes, , drop = FALSE]
  
  #euclidean distance
  euclidean.dist <- function(vec1, vec2) {
    distance <- sqrt(sum((vec1 - vec2)^2))
    return(distance)
  }	
  
  num.cells <- ncol(scRNA.data)
  num.spots <- ncol(ST.data)
  dist.matrix <- matrix(0, nrow = num.cells, ncol = num.spots)
  
  #Parallel computing
  num.cores <- detectCores()
  results <- mclapply(1:num.cells, function(i) {
    cell.vec <- scRNA.data[, i]
    dis <- sapply(1:num.spots, function(j) {
      spot.vec <- ST.data[, j]
      euclidean.dist(cell.vec, spot.vec)
    })
    dis
  }, mc.cores = num.cores - 1)
  dist.matrix <- do.call(cbind, results)
  rownames(dist.matrix) <- colnames(ST.data)
  colnames(dist.matrix) <- colnames(scRNA.data)
  return(dist.matrix)
}


##########################################################################################################
#' @title getSimCells
#' @description Get the most similar single cell for each spot.
#' @param sim.matrix Similarity matrix of single cell and spatial spot,row:spot,column:single cell.
#' @param num.cells The number of cells of each spot.
#' @return A list of cells corresponding to each spot.
#' @export getSimCells
#' @example
##########################################################################################################

#' @example sim.cells <- getSimCells(similar.matrix, num.cells %>% as.vector)  

#%>% unlist(., use.names = FALSE)

getSimCells <- function(sim.matrix, num.cells) {
  num.spots <- nrow(sim.matrix)
  sim.cells <- lapply(1:num.spots, function(i) {
    spot.sim <- sim.matrix[i,]
    top.cells <- names(spot.sim[order(spot.sim, decreasing = TRUE)[1:num.cells[i]]])
    return(top.cells)
  }) %>% `names<-`(rownames(sim.matrix))
  return(sim.cells)
} 


##########################################################################################################
#' @title getRandomCoords
#' @description Get random coordinates for mapping single-cells to spatial.
#' @param spot.coords The spatial coordinates of the spot, including two columns, imagerow and imagecol.
#' @param sim.cells A list of single cells corresponding to spots.
#' @return Spatial coordinates of simulated single cells.
#' @export
##########################################################################################################

#' @examples cell.coords <- getRandomCoords(coord.df,sim.cells)


getRandomCoords <- function(spot.coords, sim.cells) {
  #Estimated spot diameter
  kNN.dist <- dbscan::kNN(spot.coords,k= 4)$dist
  spot.diameter <- median(kNN.dist) %>% round
  spot.coords <- as.matrix(spot.coords)
  sc.coords <- lapply(1:nrow(spot.coords), function(i){
    sp.coord <- spot.coords[i,]
    cells <- lengths(sim.cells)[i]
    theta <- runif(cells,0,2*pi) 
    dis <- sqrt(runif(cells)) * (spot.diameter/2)
    cell.x <- sp.coord[1] + dis * cos(theta)
    cell.y <- sp.coord[2] + dis * sin(theta) 
    coords <- data.frame(cell.x,cell.y) 
    coords$SPOT <- paste0(spot.coords[i, 1], "x", spot.coords[i, 2])
    coords$centerX <- spot.coords[i, 1]
    coords$centerY <- spot.coords[i, 2]
    coords$SpotName <- rownames(spot.coords)[i] 
    coords
  }) %>% do.call(rbind, .) %>% as.data.frame
  return (sc.coords)
}

##########################################################################################################
#' @title mapSctoSpatial
#' @description Assign single cells to spatial coordinates.
#' @param st.obj Seurat object of Spatial Transcriptome data.
#' @param sc.obj Seurat object of single cells data.
#' @param sim.cells A list of cells corresponding to each spot.
#' @return Mapping results for single cells and spatial coordinates.
#' @export mapSctoSpatial
#' @example
##########################################################################################################

#' @example mapping <- mapSctoSpatial(st.obj, sc.obj, sim.cells,coord.df)  


mapSctoSpatial <- function(st.obj, sc.obj, sim.cells,coord.df) {
  coord.df <- coord.df %>% .[intersect(rownames(.), names(sim.cells)), ]
  st.obj <- st.obj[,rownames(coord.df)]
  sim.cells <- sim.cells[rownames(coord.df)]
  cord.new <- getRandomCoords(coord.df,sim.cells)
  
  map <- lapply(1:length(sim.cells), function(i) {
    cords <- coord.df[i, ]
    sub.cords <- cord.new[cord.new$SPOT == paste0(cords[1], "x", cords[2]), ]
    #sc.obj$CellType
    sub.cords$CellType <- Idents(sc.obj)[sim.cells[[i]]] %>% as.vector()
    sub.cords$Cell <- sim.cells[[i]]
    sub.cords
  }) %>% do.call(rbind, .) %>% as.data.frame()
  return(map)
}	


#########################################################################################################
#' @title createSeuratObj
#' @description Create a Seurat object.
#' @param st.obj Seurat object of Spatial Transcriptome data.
#' @param sc.obj Seurat object of single cell data.
#' @param mapping single cell mapping spatial results.
#' @return Seurat object.
#' @export
#########################################################################################################

#' @examples sc.out <- createSeuratObj(st.obj ,ref.obj, mapping)


createSeuratObj <- function(st.obj, sc.obj, mapping){
  
  mapping$Cell.new <- make.names(mapping$Cell,unique = T)
  sc.obj$cell <- names(sc.obj$orig.ident)
  obj <- CreateSeuratObject(counts = sc.obj@assays$RNA[, mapping$Cell] %>% set_colnames(mapping$Cell.new),
                            project = 'Cell2Space', assay = "RNA",
                            meta.data = sc.obj@meta.data[mapping$Cell, ] %>%
                              dplyr::rename(Cell.raw = cell) %>%
                              mutate(Cell.new = mapping$Cell.new) %>%
                              set_rownames(mapping$Cell.new))
  
  obj@meta.data <- dplyr::left_join(obj@meta.data, mapping) %>% data.frame %>% set_rownames(obj$Cell.new)
  Idents(obj) <- obj$CellType
  obj[["RNA"]]@data <- obj[["RNA"]]@counts										
  sc.coord.obj <- CreateDimReducObject(embeddings = mapping %>%
                                         dplyr::mutate(coord1 = cell.y, coord2 = max(cell.x)+ min(cell.x)- cell.x) %>%
                                         dplyr::select(c(coord1, coord2)) %>% set_rownames(mapping$Cell.new) %>% as.matrix,
                                       assay = "RNA", key = 'Cell2Space_')
  
  sc.pca.obj <- CreateDimReducObject(embeddings = sc.obj@reductions$pca@cell.embeddings[mapping$Cell, ] %>%
                                       set_rownames(mapping$Cell.new) %>% as.matrix, assay = "RNA", key = 'pca_')
  
  sc.umap.obj <- CreateDimReducObject(embeddings = sc.obj@reductions$umap@cell.embeddings[mapping$Cell, ] %>%
                                        set_rownames(mapping$Cell.new) %>% as.matrix, assay = "RNA", key = 'umap_')
  
  obj@reductions$Cell2Space <- sc.coord.obj
  obj@reductions$pca <- sc.pca.obj
  obj@reductions$umap <- sc.umap.obj
  
  if (!is.null(st.obj)) {
    obj@images <- st.obj@images
    obj@images[[1]]@assay <- DefaultAssay(obj)
    obj@images[[1]]@coordinates <- data.frame(imagerow = mapping$cell.x, imagecol = mapping$cell.y) %>% set_rownames(mapping$Cell.new)
    obj@images[[1]]@scale.factors <- st.obj@images[[1]]@scale.factors
  }
  return(obj)
}



