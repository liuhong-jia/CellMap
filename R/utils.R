

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
#' @title createSpatialObject
#' @description Creates a Spatial Seurat object for spatial transcriptomics (ST) data from various platforms.
#' This function ensures that only matching barcodes between the count matrix and coordinates are included.
#' @param counts A matrix of UMI counts, where rows represent genes and columns represent barcodes.
#' @param coord.df A data frame containing spatial coordinates for the barcodes. The coordinates should be in the first two columns.
#' @param coord.label A character vector specifying the column names in coord.df for the spatial coordinates. Default is c("x", "y").
#' @param meta.data Optional metadata data frame, where rows are barcodes.
#' @param class A character string specifying the class name for the Spatial image slot. Default is "SlideSeq".
#' @return A Seurat object for spatial transcriptomics data with the spatial image slot populated.
#' @export

##########################################################################################################

createSpObj <- function(counts, coord.df, coord.label = c("x", "y"), meta.data = NULL,class = "SlideSeq") {
    matched_spots <- intersect(colnames(counts), rownames(coord.df))
    if (!is.null(meta.data)) {
        matched_spots <- intersect(matched_spots, rownames(meta.data))
    }
    if (length(matched_spots) == 0) {
        stop("No matching barcodes found between the count matrix and the coordinates.")
    }

    obj <- CreateSeuratObject(counts = counts[, matched_spots, drop = FALSE], assay = "Spatial")
    obj@images$image <- new(
	    Class = class,
        assay = "Spatial",
        key = "image_",
        coordinates = coord.df[matched_spots, coord.label, drop = FALSE]
    )
    if (!is.null(meta.data)) {
        obj@meta.data <- cbind(obj@meta.data, meta.data[matched_spots, , drop = FALSE])
    }
    return(obj)
}


##########################################################################################################
#' @title Process Spatial Transcriptomics Data
#' @description A function to perform SCTransform normalization for spatial transcriptomics data.
#' @param st.obj Seurat object of spatial transcriptome data.
#' @param assay The assay to be used for SCTransform normalization. Default is "Spatial".
#' @return A Seurat object with SCTransform normalization.
#' @export

##########################################################################################################
#' @examples st.obj <- processSpatialData(st.obj)


processSpatialData <- function(st.obj, assay = "Spatial") {
  st.obj <- st.obj[, intersect(colnames(st.obj),rownames(GetTissueCoordinates(st.obj)))]
  st.obj <- SCTransform(st.obj,assay = assay)
            return(st.obj)
}


##########################################################################################################
#' @title Process scRNA-seq data.
#' @description A function to perform SCTransform normalization, PCA, UMAP, neighbor finding, and clustering on a Seurat object for spatial transcriptomics data.
#' @param sc.obj Seurat object of scRNA-seq data.
#' @param celltype.column The column name for cell type in the single-cell Seurat object, with the default value as "idents".
#' @param sc.sub.size Downsampling proportion or number for scRNA-seq data. Default: NULL
#' @param min.sc.cells The minimum number of cell types in scRNA-seq data.Default: 50
#' @return A Seurat object with SCTransform normalization.
#' @export
##########################################################################################################

#' @examples sc.obj <- processScData(sc.obj)



processScData <- function(sc.obj, celltype.column = "idents", sc.sub.size = NULL, min.sc.cells = 50) {
  
  if (!is.null(sc.sub.size)) {
    sc.obj <- if (sc.sub.size > 1) {
      downSamplSeurat(sc.obj, cnt = sc.sub.size)
    } else {
      downSamplSeurat(sc.obj, percent = sc.sub.size)
    }
  }
  
  if (celltype.column != "idents") {
    Idents(sc.obj) <- sc.obj@meta.data[, celltype.column]
  }

  cell.types <- intersect(levels(sc.obj), unique(Idents(sc.obj)))
  levels(sc.obj) <- cell.types
  
  cell.counts <- table(Idents(sc.obj))
  filtered.idents <- names(cell.counts[cell.counts >= min.sc.cells])
  sc.obj <- subset(sc.obj, idents = filtered.idents)
 
  sc.obj <- SCTransform(sc.obj)%>% FindVariableFeatures() %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims = 1:30)
  sc.obj <- subset(sc.obj, subset = nFeature_RNA > 0)
  return(sc.obj)
}



##########################################################################################################
#' @title Set Up Parallel Processing
#' @description Configures and opens multiple workers for Seurat or other parallel processing in R.
#' @param n.workers Number of cores to be used for parallel processing. Default: 8.
#' @param max.size Maximum size of global objects allowed in the future framework, in megabytes (MB). Default: 40,000 MB.
#' @export
#' @return NULL

##########################################################################################################


multipleProcess <- function(n.workers = 4, max.size = 40000) {
    options(future.globals.maxSize = max.size * 1024^2, future.seed = TRUE)
    future::plan("multicore", workers = n.workers)
    message("Parallel processing configured with ", n.workers, " workers.")
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
  future::plan("multicore", workers = parallel::detectCores())
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
#' @param st.data The spatial transcriptome expression data has genes as rows and spots as columns.
#' @param mean.cell.num Mean cell number of spot;set to 5 for 10X Visium data; set to 1 for high-resolution data.
#' @return The number of cells of each spot.
#' @export
##########################################################################################################

#' @example num.cells <- getCellNumber(st.data,mean.cell.num = 5)

getCellNumber <- function(st.data, mean.cell.num = mean.cell.num ){
	if (!inherits(st.data, "dgCMatrix")) {
		stop("st.data must be a dgCMatrix sparse matrix")
	}
	
	st.data <- st.data[, Matrix::colSums(st.data, na.rm = TRUE) > 0]  
    data.filtered <- st.data[Matrix::rowSums(st.data, na.rm = TRUE) > 0, ] 

	gene.var <- apply(data.filtered, 1, var)
	stable.genes <- rownames(data.filtered)[gene.var <= 0.5]
	if (length(stable.genes) == 0) {
		warning("No stable genes found with variance <= 0.5")
	}
	if (length(stable.genes) > 0) {
		spot.total <- Matrix::colSums(st.data[stable.genes, ])
		spot.mean <- sum(st.data[stable.genes, ]) / ncol(st.data)
		spot.cell <- round((spot.total / spot.mean) * mean.cell.num)
		spot.cell[spot.cell == 0] <- 1
	} else {
		spot.cell <- rep(mean.cell.num, ncol(st.data))
		}
  return(spot.cell)
}




##########################################################################################################
#' @title getELDist
#' @description Calculate the Euclidean distance between single cells and spots.
#' @param sc.data scRNA-seq data of sampled cells.
#' @param st.data Spatial transcriptomics expression data.
#' @return Euclidean distance of single cells and spots.
#' @export
##########################################################################################################

#' @example dist.matrix <- getELDist(sc.data,st.data)

getELDist <- function(sc.data, st.data) {
  over.genes <- intersect(rownames(sc.data), rownames(st.data))
  st.data <- st.data[over.genes, , drop = FALSE]
  sc.data <- sc.data[over.genes, , drop = FALSE]
  
  #euclidean distance
  euclidean.dist <- function(vec1, vec2) {
    distance <- sqrt(sum((vec1 - vec2)^2))
    return(distance)
  }	
  
  num.cells <- ncol(sc.data)
  num.spots <- ncol(st.data)
  dist.matrix <- matrix(0, nrow = num.cells, ncol = num.spots)
  
  #Parallel computing
  num.cores <- detectCores()
  results <- mclapply(1:num.cells, function(i) {
    cell.vec <- sc.data[, i]
    dis <- sapply(1:num.spots, function(j) {
      spot.vec <- st.data[, j]
      euclidean.dist(cell.vec, spot.vec)
    })
    dis
  }, mc.cores = num.cores - 1)
  dist.matrix <- do.call(cbind, results)
  rownames(dist.matrix) <- colnames(st.data)
  colnames(dist.matrix) <- colnames(sc.data)
  return(dist.matrix)
}


##########################################################################################################
#' @title getSimCells
#' @description Get the most similar single cell for each spot.
#' @param sim.matrix Similarity matrix of single cell and spatial spot,row:spot,column:single cell.
#' @param num.cells The number of cells of each spot.
#' @return A list of cells corresponding to each spot.
#' @export getSimCells
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
#' @param st.obj Seurat object of spatial transcriptome data.
#' @param sim.cells A list of single cells corresponding to spots.
#' @param n.workers Number of cores to be used for parallel processing. Default: 4.
#' @return Spatial coordinates of simulated single cells.
#' @export
##########################################################################################################

#' @examples cell.coords <- getRandomCoords(st.obj,sim.cells,n.workers = 4)


getRandomCoords <- function(st.obj,sim.cells,n.workers = 4) {
  future::plan("multicore", workers = n.workers)
  spot.coords <- GetTissueCoordinates(st.obj,scale = NULL) %>% .[intersect(rownames(.),names(sim.cells)),1:2]
  #Estimated spot diameter
  kNN.dist <- dbscan::kNN(spot.coords,k=4)$dist
  spot.diameter <- median(kNN.dist) %>% round
  spot.coords <- as.matrix(spot.coords)
  sc.coords <- future.apply::future_lapply(1:nrow(spot.coords), function(i){
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
  },future.seed = TRUE) %>% do.call(rbind, .) %>% as.data.frame
  return (sc.coords)
}

##########################################################################################################
#' @title mapSctoSpatial
#' @description Assign single cells to spatial coordinates.
#' @param st.obj Seurat object of Spatial Transcriptome data.
#' @param sc.obj Seurat object of scRNA-seq data.
#' @param sim.cells A list of cells corresponding to each spot.
#' @return Mapping results for single cells and spatial coordinates.
#' @export mapSctoSpatial
##########################################################################################################

#' @example mapping <- mapSctoSpatial(st.obj, sc.obj, sim.cells)  


mapSctoSpatial <- function(st.obj, sc.obj,sim.cells) {
  coord.df <- GetTissueCoordinates(st.obj,scale = NULL) %>% .[intersect(rownames(.),names(sim.cells)),1:2]
  st.obj <- st.obj[,rownames(coord.df)]
  sim.cells <- sim.cells[rownames(coord.df)]
  cord.new <- getRandomCoords(st.obj,sim.cells)
  map <- lapply(1:length(sim.cells), function(i) {
    cords <- coord.df[i, ]
    sub.cords <- cord.new[cord.new$SPOT == paste0(cords[1], "x", cords[2]), ]
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

#' @examples sc.out <- createSeuratObj(st.obj,sc.obj, mapping)


createSeuratObj <- function(st.obj, sc.obj, mapping){
  
  mapping$Cell.new <- make.names(mapping$Cell,unique = T)
  sc.obj$cell <- names(sc.obj$orig.ident)
  obj <- CreateSeuratObject(counts = sc.obj@assays$RNA[, mapping$Cell] %>% set_colnames(mapping$Cell.new),
                            project = 'CellMap', assay = "RNA",
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
                                       assay = "RNA", key = 'CellMap_')
  
  sc.pca.obj <- CreateDimReducObject(embeddings = sc.obj@reductions$pca@cell.embeddings[mapping$Cell, ] %>%
                                       set_rownames(mapping$Cell.new) %>% as.matrix, assay = "RNA", key = 'pca_')
  
  sc.umap.obj <- CreateDimReducObject(embeddings = sc.obj@reductions$umap@cell.embeddings[mapping$Cell, ] %>%
                                        set_rownames(mapping$Cell.new) %>% as.matrix, assay = "RNA", key = 'umap_')
  
  obj@reductions$CellMap <- sc.coord.obj
  obj@reductions$pca <- sc.pca.obj
  obj@reductions$umap <- sc.umap.obj
  
  if (!is.null(st.obj)) {
    obj@images <- st.obj@images
    obj@images[[1]]@assay <- DefaultAssay(obj)
    obj@images[[1]]@coordinates <- data.frame(imagerow = mapping$cell.x, imagecol = mapping$cell.y) %>%
                                    set_rownames(mapping$Cell.new)

   if ("scale.factors" %in% slotNames(st.obj@images[[1]])) {
    if (!is.null(st.obj@images[[1]]@scale.factors)) {
        obj@images[[1]]@scale.factors <- st.obj@images[[1]]@scale.factors
    } else {
        obj@images[[1]]@scale.factors <- NULL
    }
   }
}
  return(obj)
}



