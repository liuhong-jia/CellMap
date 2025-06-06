% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inter.R
\name{getSimMatrix}
\alias{getSimMatrix}
\title{getSimMatrix}
\usage{
getSimMatrix(sc.to.cluster, st.obj, sc.obj, markers)
}
\arguments{
\item{sc.to.cluster}{The correspondence between cells and clusters.}

\item{st.obj}{st data Seurat object.}

\item{sc.obj}{sc data Seurat object.}

\item{sig.mat}{Marker gene expression matrix.}
}
\value{
A similarity matrix of single cells and spots.
}
\description{
Calculate the similarity of single cells and spots.
}
