# CellMapper 
CellMapper is an innovative tool crafted to precisely map individual cells onto spatial coordinates within tissue slices. Its broad utility lies in unveiling the spatial distribution of cell types, dissecting cellular compositions within tissue section spatial spots, and identifying critical functional structures within biological systems.

![image](https://github.com/liuhong-jia/CellMapper/blob/main/vignettes/workflow.png)

In this tutorial, we will illustrate the installation and usage of CellMapper using the mouse brain cortex dataset as an example.

## Installing the package
To install CellMapper,we recommed using devtools:

```
library(devtools)
devtools::install_github("liuhong-jia/CellMapper")  
```

## Dependencies
● R version >= 4.3.0.

● R packages: Seurat, dplyr, ggplot2, Matrix, jsonlite, magrittr, randomForest, parallel

## Importing packages and preparing input data(scRNA-seq data and spatial transcriptomes data)

```
library(CellMapper)
library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(jsonlite)
library(magrittr)
library(randomForest)
library(parallel)
```

```
sc.data <- readRDS(system.file("data", "cortex.sc.rds", package = "CellMapper"))
st.data <- readRDS(system.file("data", "st.rds",package = "CellMapper"))
ref.repr <- GetAssayData(sc.data, slot = 'counts') %>% as.data.frame
ref.anno <- sc.data$subclass %>% as.vector
coord.df <- st.data@images$anterior1@coordinates[,c(4,5)]
```

## Setting the parameters
