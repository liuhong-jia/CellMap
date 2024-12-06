##########################################################################################################
#' @title getInterdata
#' @description Integrate single-cell(sc) and spatial transcriptomic(st) data.
#' @param sc.obj sc.obj Seurat object of scRNA-seq data.
#' @param st.obj Seurat object of spatial transcriptome data.
#' @param markers Feature genes for intergration.
#' @param coord Coordinates column names in st images slot.
#' @return Seurat object of integrated sc and st data.
#' @export 

##########################################################################################################	

#' @example sc.st.obj <- getInterdata(sc.obj, st.obj,coord = c("imagerow","imagecol"),markers)


getInterdata <- function(sc.obj, st.obj, coord = c("imagerow","imagecol"),markers) {
  
  sc.obj$id <- colnames(sc.obj)
  st.obj$id <- colnames(st.obj)
  
  sc.obj$type <- "sc"
  st.obj$type <- "st"
  
  st.obj <- st.obj[,rownames(st.obj@images[[1]]@coordinates)]
  
  st.obj$image_x <- st.obj@images[[1]]@coordinates[, coord[1]]
  st.obj$image_y <- st.obj@images[[1]]@coordinates[, coord[2]]
  
  DefaultAssay(st.obj) <- "Spatial"
  DefaultAssay(sc.obj) <- "RNA"
  
  sc.st.features <- intersect(markers,rownames(st.obj))
  
  # Find transfer anchors
  sc.st.anchors <- FindTransferAnchors(reference = sc.obj, query = st.obj,
                                       features = sc.st.features, reduction = "cca")
  
  # Transfer data
  st.obj.trans <- TransferData(anchorset = sc.st.anchors, refdata = GetAssayData(sc.obj,slot = "data")[sc.st.features,],weight.reduction = "cca")
  
  st.obj@assays$trans <- st.obj.trans
  
  # Create integrated Seurat object
  sc.st.meta <- bind_rows(st.obj@meta.data, sc.obj@meta.data)
  counts <- cbind(data.frame(st.obj@assays$trans@data), data.frame(sc.obj@assays$RNA@data[sc.st.features,]))
  
  rownames(sc.st.meta) <- sc.st.meta$id
  colnames(counts) <- sc.st.meta$id
  
  sc.st.obj <- CreateSeuratObject(counts = counts, assay = "integrate", meta.data = sc.st.meta)

  sc.st.obj@assays$integrate@data <- sc.st.obj@assays$integrate@counts
  sc.st.obj@assays$integrate@counts <- matrix(nrow = 0,ncol = 0)
  
  # Run PCA and UMAP
  sc.st.obj <- ScaleData(sc.st.obj,features = rownames(sc.st.obj)) %>% RunPCA(features = rownames(sc.st.obj)) %>% RunUMAP(dim = 1:30)
  
  sc.st.obj@images <- st.obj@images
  sc.st.obj@images[[1]]@coordinates <- data.frame(imagerow = sc.st.obj$image_x,
                                                  imagecol = sc.st.obj$image_y)
  return(sc.st.obj)
}


##########################################################################################################
#' @title getDistMatrix
#' @description Calculate the distance matrix between single cells and spatial spots.
#' @param sc.st.obj Seurat object of integrated sc and st data.
#' @return A distance matrix of single cell and spot.
#' @export
##########################################################################################################	

#' @example dist <- getDistMatrix(sc.st.obj)

getDistMatrix <- function(sc.st.obj) {
  
  UMAP <- sc.st.obj@reductions$umap@cell.embeddings
  
  st.meta.data <- sc.st.obj@meta.data[sc.st.obj@meta.data$type == "st",]
  sc.meta.data <- sc.st.obj@meta.data[sc.st.obj@meta.data$type == "sc",]
  
  st.UMAP <- UMAP[rownames(st.meta.data),]
  sc.UMAP <- UMAP[rownames(sc.meta.data),]
  
  #euclidean distance
  euclidean.dist <- function(vec1, vec2) {
    distance <- sqrt(sum((vec1 - vec2)^2))
    return(distance)
  }	
  num.cells <- nrow(sc.UMAP)
  num.spots <- nrow(st.UMAP)
  dist.matrix <- matrix(0, nrow = num.cells, ncol = num.spots)
  
  #Parallel computing
  num.cores <- parallel::detectCores()
  
  dist.matrix <- do.call(cbind, parallel::mclapply(1:num.cells, function(i) {
    cell.vec <- sc.UMAP[i,]
    dis <- sapply(1:num.spots, function(j) {
      spot.vec <- st.UMAP[j,]
      euclidean.dist(cell.vec, spot.vec)
    })
    dis
  }, mc.cores = num.cores - 1))
  
  rownames(dist.matrix) <- rownames(st.UMAP)
  colnames(dist.matrix) <- rownames(sc.UMAP)
  
  return(dist.matrix)
}

##########################################################################################################
#' @title getCells
#' @description Obtaining the k nearest single cells of each spot.
#' @param sc.st.dist A distance matrix of single cell and spot.
#' @param k Number of nearest neighboring single cells for each spot.
#' @return A list of nearest single cells corresponding to the spot.
#' @export
##########################################################################################################	

#' @example nearCells <- getCells(dist, knn = 5)

getCells <- function(sc.st.dist, k = knn) {
  neighbors <- lapply(1:nrow(sc.st.dist), function(i){
    near.index <- order(sc.st.dist[i,])[1:k]
    colnames(sc.st.dist)[near.index]
  })
  names(neighbors) <- rownames(sc.st.dist)
  return(neighbors)
}

##########################################################################################################
#' @title getUniqueCells
#' @description Obtaining the k nearest single cells of each spot,and each single cell is selected only once.
#' @param dist.matrix A distance matrix of single cell and spot.
#' @param k Number of single cells for each spot.Default:3
#' @return A list of nearest single cells corresponding to the spot.
#' @export

##########################################################################################################	
#example nearCells <- getUniqueCells(dist, k = 3)

getUniqueCells <- function(dis.matrix, k = 3) {
  n.spots <- nrow(dis.matrix)
  n.cells <- ncol(dis.matrix)
  neighbors <- vector("list", length = n.spots)
  selected.cells <- c()  
  
  for (i in 1:n.spots) {
    sorted.indices <- order(dis.matrix[i, ])
    j <- 1
    count <- 0
    while (count < k && j <= n.cells) {
      cell.index <- sorted.indices[j]
      cell.name <- colnames(dis.matrix)[cell.index]
      
      if (!cell.name %in% selected.cells) {
        if (is.null(neighbors[[i]])) {
          neighbors[[i]] <- vector("character", length = k)
        }
        neighbors[[i]][count + 1] <- cell.name
        selected.cells <- c(selected.cells, cell.name)  
        count <- count + 1
      }
      j <- j + 1
    }
  }
  
  names(neighbors) <- rownames(dis.matrix)
  return(neighbors)
}



##########################################################################################################
#' @title trainModel
#' @description Train a random forest model and make predictions.
#' @param st.obj Seurat object of spatial transcriptome data.
#' @param sc.obj Seurat object of scRNA-seq data.
#' @param nearCells A list of nearest single cells corresponding to the spot.
#' @param markers Cell type specific marker genes.
#' @param n.workers Number of cores to be used for parallel processing. Default: 4.
#' @return A list of predicted results.
#' @export

##########################################################################################################	
#' @example pred <- trainModel(st.obj,sc.obj,nearCells,markers,n.workers = 4)


trainModel <- function(st.obj,sc.obj,nearCells,markers,n.workers = 4){
  
	st.meta <- data.frame(Spot = rownames(st.obj@meta.data),cluster = st.obj$seurat_clusters)
  
	future::plan("multicore", workers = n.workers)
	spot.df <- future.apply::future_lapply(names(nearCells), function(spot.name) {
        data.frame(Spot = spot.name, Cell = nearCells[[spot.name]], stringsAsFactors = FALSE)
        },future.seed = TRUE) %>%
        do.call(rbind, .)
  
	merged.data <- merge(spot.df,st.meta, by = "Spot")
	genes <- intersect(rownames(st.obj), markers)
  
	sc.matrix <- GetAssayData(sc.obj,slot = "data")
	x.train <- sc.matrix[genes,merged.data$Cell] %>% as.matrix %>% t
	y.train <- merged.data$cluster %>% as.character %>% as.factor
 
	#model training
	rf.model <- randomForest(x.train, y.train, ntree = 500)

	#model prediction
	x.test <- GetAssayData(sc.obj[genes,],slot = "data") %>% as.matrix %>% t
	rf.pred <- predict(rf.model, newdata = x.test, type = "response") 
	sc.to.cluster <- data.frame(cell = rownames(x.test),cluster = rf.pred) 
	return(list(model = rf.model, prediction = sc.to.cluster))
}

##########################################################################################################
#' @title getSimMatrix
#' @description Calculate the similarity of single cells and spots.
#' @param sc.to.cluster The correspondence between cells and clusters.
#' @param st.obj st data Seurat object.
#' @param sc.obj sc data Seurat object.
#' @param sig.mat Marker gene expression matrix.
#' @return A similarity matrix of single cells and spots.
#' @export
##########################################################################################################	

#' @example sim <- getSimMatrix(pred,st.obj,sc.obj,markers)


getSimMatrix <- function(sc.to.cluster,st.obj,sc.obj,markers) {
  
  #cosine similarity
  cosine.sim <- function(vec1, vec2) {
    dot.product <- sum(vec1 * vec2)
    norm1 <- sqrt(sum(vec1^2))
    norm2 <- sqrt(sum(vec2^2))
    similar <- dot.product/(norm1 * norm2)
    return(similar)
  }
  
  #sc.to.cluster <- pred$prediction
  spot.cluster <- data.frame(cluster = st.obj$seurat_clusters)
  unique.clusters <- unique(spot.cluster$cluster) %>% as.vector
  
  genes <- intersect(rownames(st.obj), markers)
  sc.data.matrix <- GetAssayData(sc.obj,slot = "data")
  st.data.matrix <- GetAssayData(st.obj,slot = "data")
  
  sc.data.matrix <- sc.data.matrix[genes,]
  st.data.matrix <- st.data.matrix[genes,]
  
  sim.results <- list()
  for (cluster.id in unique.clusters) {
    cluster.cells <- sc.to.cluster$cell[sc.to.cluster$cluster == cluster.id]
    cluster.spots <- row.names(spot.cluster)[spot.cluster$cluster == cluster.id]
    
    num.cells <- length(cluster.cells)
    num.spots <- length(cluster.spots)
    sim.matrix <- matrix(0, nrow = num.cells, ncol = num.spots)
    
    #Parallel computing
    num.cores <- parallel::detectCores()
    results <- parallel::mclapply(1:num.cells, function(i) {
      cell.vec <- sc.data.matrix[, cluster.cells[i]]
      similarity <- sapply(1:num.spots, function(j) {
        spot.vec <- st.data.matrix[, cluster.spots[j]]
        cosine.sim(cell.vec, spot.vec)
      })
      similarity
    }, mc.cores = num.cores - 1)
    sim.matrix <- do.call(cbind, results)
    rownames(sim.matrix) <- cluster.spots
    colnames(sim.matrix) <- cluster.cells
    sim.results[[cluster.id]] <- sim.matrix
  }
  return(sim.results)
}


##########################################################################################################
#' @title getNearCells
#' @description Compute the nearest single cells for each spot.
#' @param sim Similarity matrix of single cells and spots.
#' @param num.cells The number of single cells contained in each spot.
#' @return A list of single cells corresponding to spots.
#' @export
##########################################################################################################	

#' @example sim.cells <- getNearCells(sim,num.cells) 

getNearCells <- function(sim,num.cells) {
  
  cells <- unlist(lapply(names(sim), function(i) {
    sim.matrix <- sim[[i]]
    spot.name <- rownames(sim.matrix)
    num.cells.clu <- num.cells[spot.name]
    
    near.cells <- sapply(1:nrow(sim.matrix), function(j) {
      spot.sim <- sim.matrix[j,]
      top.cells <- names(spot.sim[order(spot.sim, decreasing = TRUE)[1:num.cells.clu[j]]])
    }) %>% `names<-`(rownames(sim.matrix))
    return(near.cells)
  }),recursive = FALSE)
} 

##########################################################################################################	
#' @title getCellClusterInfo
#' 
#' @description This function retrieves the cluster information for each cell based on the spot-cell correspondence.
#' @param nearCells A list where the names are spot identifiers and the values are cell identifiers.
#' @param st.obj Seurat object of spatial transcriptomics data with metadata including cluster information.
#' @export
##########################################################################################################	

getCellClusterInfo <- function(nearCells, st.obj) {
  meta_data <- st.obj@meta.data
  result <- data.frame(cell = character(), cluster = factor(), stringsAsFactors = FALSE)
  
  for (spot in names(nearCells)) {
    cells <- nearCells[[spot]]
    cluster_info <- meta_data[spot, "seurat_clusters"]
    result <- rbind(result, data.frame(cell = cells, cluster = cluster_info, stringsAsFactors = FALSE))
  }
  
  return(result)
}



##########################################################################################################
#' @title linearAllocation
#' @description Apply linear allocation algorithm to achieve one-to-one mapping of single cells to spots.
#' @param sim Similarity matrix of single cells and spots.
#' @param num.cells The number of single cells contained in each spot.
#' @return A list of single cells corresponding to spots.
#' @export
##########################################################################################################	

#' @example sim.cells <- linearAllocation(sim,num.cells) 


linearAllocation <- function(sim,num.cells){
  
  cells <- unlist(lapply(names(sim),function(i){
    sim.matrix <- sim[[i]]
    cell.num <- num.cells[rownames(sim.matrix)]
    mat <- matrix(0, nrow = sum(cell.num), ncol = ncol(sim.matrix))
    
    if(nrow(mat) < ncol(mat)){
      index <- 1  
      for (spot.idx in 1:length(cell.num)) {
        num <- cell.num[spot.idx]
        mat[index:(index + num - 1), ] <- matrix(rep(sim.matrix[spot.idx, ], each = num), nrow = num, ncol = ncol(sim.matrix))
        index <- index + num
        colnames(mat) <- colnames(sim.matrix)
      }
      ##linear solver
      spot.name <- rep(names(cell.num),cell.num) %>% as.vector
      cell.name <- colnames(sim.matrix) %>% as.vector
      cost.mat <- 1 - mat
      
      solve <- solve_LSAP(cost.mat)
      cell <- cell.name[solve]
      df <- data.frame(SpotName = spot.name,Cell = cell)
      df.list <- split(df$Cell,df$SpotName)
    }
    else{
      df.list <- sapply(1:nrow(sim.matrix), function(j) {
        spot.sim <- sim.matrix[j,]
        top.cells <- names(spot.sim[order(spot.sim, decreasing = TRUE)[1:cell.num[j]]])
      }) %>% `names<-`(rownames(sim.matrix))
    }
    return(df.list)
  }),recursive = FALSE)
}





