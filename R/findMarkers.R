##########################################################################################################
#' @title createSeuratObjAndQC
#' @param ref.expr Reference gene expression profile.
#' @param ref.anno Cell types corrsponding to the single cells in the reference.
#' @return Seurat object.
##########################################################################################################


createSeuratObjAndQC <- function(ref.expr, ref.anno) {
	ref.obj <- CreateSeuratObject(
		counts = ref.expr, 
		project = 'Reference', 
		meta.data = ref.anno %>% as.data.frame(.) %>% `colnames<-`('CellType') %>% `rownames<-`(colnames(ref.expr))
	)
	Idents(ref.obj) <- ref.anno
	ref.obj[["percent.mt"]] <- PercentageFeatureSet(ref.obj, pattern = "^MT-")
	ref.obj[["percent.rb"]] <- PercentageFeatureSet(ref.obj, pattern = "^RP[SL]")
	ref.obj <- subset(
		ref.obj, 
		subset = ( 
			nFeature_RNA > min(quantile(ref.obj@meta.data$nFeature_RNA, 0.05), 300) & 
			nFeature_RNA < max(quantile(ref.obj@meta.data$nFeature_RNA, 0.95), 5000) & 
			nCount_RNA   > min(quantile(ref.obj@meta.data$nCount_RNA, 0.05), 500) & 
			nCount_RNA   < max(quantile(ref.obj@meta.data$nCount_RNA, 0.95), 10000) & 
			percent.mt   < max(quantile(ref.obj@meta.data$percent.mt, 0.95), 15) & 
			percent.rb   < max(quantile(ref.obj@meta.data$percent.rb, 0.95), 35)
		)
	)
	return(ref.obj)
}


##########################################################################################################
#' @title identSeedGenes
#' @param avg.expr Reference gene expression profile.
#' @param factor.size Factor size to scale the weight. Default: 0.1.
#' @param count Number of candidate query genes. Default: 10.
#' @return A list of seed genes for cell types.
#' @export 
##########################################################################################################

identSeedGenes <- function(avg.expr, factor.size = 0.1, count = 10) {
	genes.weight <- apply(avg.expr, 1, median) / median(apply(avg.expr, 1, median))
	genes.weight <- tanh(factor.size * genes.weight)
	seed.genes <- lapply(1 : ncol(avg.expr), FUN = function(idx) {
		xx <- (avg.expr[, idx] / rowMeans(avg.expr[, -idx]) * genes.weight) %>% .[!is.na(.)]
		xx[order(xx) %>% rev][1 : count]
	}) %>% `names<-`(colnames(avg.expr %>% data.frame))
	return(seed.genes)
}

##########################################################################################################
#' @description Search marker genes on the basis of seed genes for each cell type.
#' @param ref.expr Refernce gene expression profile.
#' @param query.genes A list of query genes.
#' @param scale.data Scale gene expression or not, default: TRUE.
#' @param p.cut Threshold for filtering cell type-specific markers. Default: 0.1.
#' @param seed.num Number of seed genes. Default: 10.
#' @return Seurat object 
#' @export
##########################################################################################################


searchMarkersByCorr <- function(ref.expr, query.genes, scale.data = TRUE, p.cut = 0.1, seed.num = 10, tcga.data.u = NULL) {
    options(warn = -1)
	if(is.null(tcga.data.u)){
	if (!inherits(query.genes, 'list')) stop(sprintf('%s must be a list of query genes, exiting...', deparse(substitute(query.genes))))
	if (scale.data) {
		ref.expr <- sweep(ref.expr, 1, apply(ref.expr, 1, mean), "-")
		ref.expr <- sweep(ref.expr, 1, sqrt(apply(ref.expr, 1, var)), "/")
	}
	ref.expr <- ref.expr[rowSums(is.na(ref.expr)) == 0, ]
	p1 <- irlba::prcomp_irlba(ref.expr %>% t, n = min(length(query.genes) * 2, min(dim(ref.expr))))
	data.u <- p1$rotation %>% as.data.frame %>% t
	colnames(data.u) <- rownames(ref.expr)
	}
	else {data.u <- tcga.data.u}
	#ifelse (!is.null(tcga.data.u),data.u <- tcga.data.u, data.u <- data.u)
	
	marker.lst <- lapply(query.genes, function(QGs) {
		QGs <- names(QGs) 
		cov.mat <- t(data.u[, QGs]) %*% data.u
		diags <- sqrt(colSums((data.u)^2))
		cor.mat <- cov.mat / diags[QGs]
		cor.mat <- sweep(cor.mat, 2, diags, '/')
		colnames(cor.mat) <- colnames(data.u)
		rownames(cor.mat) <- QGs
		if (seed.num > 1) {
			cor.vec <- colMeans(atanh(cor.mat[seq_len(length(QGs)), -c(which(colnames(cor.mat) %in% QGs))]))
		} else {
			cor.vec <- colMeans(atanh(cor.mat[seq_len(length(QGs)), -c(which(colnames(cor.mat) %in% QGs))]) %>% as.data.frame %>% t)
		}
		cor.vec[order(cor.vec) %>% rev]		
		pvals <- pnorm(cor.vec, mean = mean(cor.vec), sd = sd(cor.vec), lower.tail = FALSE)
		c(QGs, pvals[pvals < p.cut] %>% names) %>% unique
	})
	return(marker.lst)
}


##########################################################################################################
#' @title buildSigMatrix
#' @param markers.lst A list of markers for cell types.
#' @param avg.expr Average gene expression of reference.
#' @param min.size Minimum size of each subset. Default: 10.
#' @param max.size Maximum size of each subset. Default: 200.
#' @return Signature matrix.
#' @export 
##########################################################################################################

buildSigMatrix <- function(markers.lst, avg.expr, min.size = 10, max.size = 200, subset.size = NULL) {
	if (!is.null(subset.size)) {
		sigmat.res <- lapply(markers.lst, function(genes) genes[1 : min(subset.size, length(genes))]) %>% unlist(.) %>% unique(.) %>% avg.expr[., ]
		return(sigmat.res)
	}
	max.size <- lapply(markers.lst, length) %>% unlist(., use.names = FALSE) %>% min(.) %>% min(., max.size)
	if (max.size < 10) {
		message('Warning: too few marker genes!')
		return(avg.expr[unlist(markers.lst) %>% unique, ])
	}	
	kappa.vals <- lapply(min.size : max.size, function(it.num) {		
		kp.val <- lapply(markers.lst,function(genes) genes[1 : it.num]) %>% 
			unlist(., use.names = FALSE) %>% 
			unique(.) %>% 
			avg.expr[., ] %>% 
			kappa(.) 
		kp.val 
	})
    group.size <- min.size + which.min(kappa.vals) - 1
	sigmat.res <- lapply(markers.lst, function(genes) genes[1 : group.size]) %>% unlist(.) %>% unique(.) %>% avg.expr[., ]
	markers <- rownames(sigmat.res)
	return(markers)
}




