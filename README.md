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
- R version >= 4.3.0.

- R packages: Seurat, dplyr, ggplot2, Matrix, jsonlite, magrittr, randomForest, parallel

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
sc.data <- readRDS(system.file("data", "sc.rds", package = "CellMapper"))
st.data <- readRDS(system.file("data", "st.rds",package = "CellMapper"))
ref.repr <- GetAssayData(sc.data, slot = 'counts') %>% as.data.frame
ref.anno <- sc.data$subclass %>% as.vector
coord.df <- st.data@images$anterior1@coordinates[,c(4,5)]
```

## Setting the parameters
|**Parameters**|**Description**                      |
|----------|-----------------------------------------|
|st.data   |A Seurat object of spatial transcriptomics data.|
|coord.df  |The spatial coordinates of tissue sections in spatial transcriptomic data.|
|ref.expr  | The gene expression profile of scRNA-seq data.|
|ref.anno  |Cell type information of scRNA-seq data, corresponding to the above `ref.expr`.|
|factor.size|Factor size for scaling the weight of gene expression. Default: 0.1.|
|seed.num|Number of seed genes of each cell type for recognizing marker genes. Default: 10.|
|pvalue.cut|Threshold for filtering cell type-specific markers. Default: 0.01|
|knn       |The number of nearest neighboring single cells for each spot. Default: 5|
|verbose   |Show running messages or not. Default: TRUE.|

## Run CellMapper  to assign single cells to spatial locations.
Details of the results is described in the table below.
|**output**|**details**|
|------|-------|
|sc.out|Seurat object of spatial transcriptomic data with single-cell resolution.|
|decon |The cellular composition of each spot in tissue sections.|

    results <- CellMapper(st.data = st.data,
                          coord.df = coord.df,
                          ref.expr = ref.expr,
                          ref.anno = ref.anno,
                          factor.size = 0.1,
                          seed.num = 10,
                          pvalue.cut = 0.1,
		              knn = 5,
                          verbose = TRUE)
   
     [INFO] Identification of cell type-specific genes...
     [INFO] Estimate the number of single cells in the spot
     [INFO] Integrate single-cell and spatial spot data
     [INFO] Train a random forest model and predict,waiting...
     [INFO] Map single cells onto spatial spots
     [INFO] Construct Seurat object
     [INFO] Finish!
  
## Visualization 
```
sc.data <- sc.data %>%
           NormalizeData() %>%
           FindVariableFeatures() %>%
           ScaleData() %>%
           RunPCA() %>%
           RunUMAP(dims = 1:30)
```

```
colors <- c("Astro" = "brown3","L2/3 IT" = "chartreuse3","L4" = "khaki1","L5 IT" = "plum1","L5 PT" = "#06986b","L6 CT" ="lightskyblue",
             "L6 IT" = "darkorange","L6b" ="hotpink","Lamp5" ="seashell3","NP" ="#a75629","Pvalb" ="#986ada","Sst" = "#3973d5","Vip" = "#edad0e")
p1 <- DimPlot(sc.data,group.by= "subclass",label = T,label.size = 6,
        cols = colors, pt.size = 1.5 , repel = T ) + 
        NoLegend() + labs(x = "UMAP1",y = "UMAP2", title = "CellType") +
        theme(panel.border = element_rect(fill=NA,color= "black",size= 1,linetype="solid"))+
        theme(axis.title.x =element_text(size=24), axis.title.y=element_text(size=24))+
        theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),
        axis.text=element_text(size=12,face = "bold"),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14))
p2 <- SpatialDimPlot(results$sc.out, group.by = "CellType",pt.size.factor = 1,label.size = 8,cols = colors)+
      theme(legend.text=element_text(size= 10))
```
