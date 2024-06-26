% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{downSamplSeurat}
\alias{downSamplSeurat}
\title{downSamplSeurat}
\usage{
downSamplSeurat(obj, cluster.col = NULL, cnt = 200, seed = 123, percent = 0.3)
}
\arguments{
\item{obj}{Seurat object.}

\item{cnt}{Sample size for each ident, default: 200.}

\item{seed}{Randon seed, default: 123.}
}
\value{
Subset of seurat object.
}
\description{
Sample single cells by cell type.
}
