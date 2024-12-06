
##########################################################################################################
#' @title checkInputParams
#' @return TRUE.

##########################################################################################################


checkInputParams <- function(st.obj, sc.obj, coord,celltype.column, sc.sub.size, min.sc.cell,
                              factor.size, seed.num, pvalue.cut, knn, mean.cell.num, max.cell.num,n.workers, verbose) {
  
  if (!inherits(st.obj, 'Seurat')) {
    stop(sprintf('%s must be a Seurat object, exiting...', deparse(substitute(st.obj))))
  }
  
  if (!inherits(sc.obj, 'Seurat')) {
    stop(sprintf('%s must be a Seurat object, exiting...', deparse(substitute(sc.obj))))
  }
  
  if (!is.character(coord) || length(coord) != 2) {
    stop('coord must be a character vector of length 2, exiting...')
  }
  
  if (!is.character(celltype.column)) {
    stop('celltype.column must be a character string, exiting...')
  }

  if (!is.null(sc.sub.size) && !is.numeric(sc.sub.size)) {
    stop('sc.sub.size must be numeric or NULL, exiting...')
  }

  if (!is.numeric(min.sc.cell) || min.sc.cell <= 0 || min.sc.cell %% 1 != 0) {
    stop('min.sc.cell must be a positive integer, exiting...')
  }

  if (!is.numeric(factor.size) || factor.size <= 0) {
    stop('factor.size must be a positive number, exiting...')
  }

  if (!is.numeric(seed.num) || seed.num %% 1 != 0) {
    stop('seed.num must be an integer, exiting...')
  }

  if (!is.numeric(pvalue.cut) || pvalue.cut < 0 || pvalue.cut > 1) {
    stop('pvalue.cut must be a number between 0 and 1, exiting...')
  }

  if (!is.numeric(knn) || knn <= 0 || knn %% 1 != 0) {
    stop('knn must be a positive integer, exiting...')
  }

 if (!(is.numeric(mean.cell.num) && mean.cell.num > 0 && mean.cell.num %% 1 == 0)) {
  stop('mean.cell.num must be a positive integer, exiting...')
}
  
 if (!(is.numeric(max.cell.num) && max.cell.num > 0 && max.cell.num %% 1 == 0)) {
  stop('max.cell.num must be a positive integer, exiting...')
}

  if (!is.numeric(n.workers) || n.workers <= 0 || n.workers %% 1 != 0) {
    stop('n.workers must be a positive integer, exiting...')
  }

  if (!is.logical(verbose)) {
    stop('verbose must be a logical value (TRUE or FALSE), exiting...')
  }

  return(TRUE)
}