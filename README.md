# CellMap 
CellMap is an innovative tool crafted to precisely map individual cells onto spatial coordinates within tissue slices. Its broad utility lies in unveiling the spatial distribution of cell types, dissecting cellular compositions within tissue section spatial spots, and identifying critical functional structures within biological systems.

![image](https://github.com/liuhong-jia/CellMap/blob/main/vignettes/workflow.png)

In this tutorial, we will illustrate the installation and usage of CellMap using the HER2+ breast cancer dataset as an example.

## Installing the package
To install CellMap,we recommed using devtools:

```
library(devtools)
devtools::install_github("liuhong-jia/CellMap")  
```

## Dependencies
- R version >= 4.3.0.

- R packages: Seurat, dplyr, ggplot2, Matrix, clue, jsonlite, magrittr, randomForest, parallel

## Importing packages and preparing input data(scRNA-seq data and spatial transcriptomes data)

```
library(CellMap)
library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(clue)
library(jsonlite)
library(magrittr)
library(randomForest)
library(parallel)
```
A zip file containing single-cell and spatial transcriptomics data can be downloaded from the following link:

[example data](https://drive.google.com/file/d/1lu0Y8hknGm6aVKogXZmQAUM7PsxAZghX/view?usp=drive_link)

```
sc.data <- readRDS("sc.obj.rds")
st.data <- readRDS("st.obj.rds")
```

## Setting the parameters
|**Parameters**|**Description**                      |
|----------|-----------------------------------------|
|st.obj    |Seurat object of spatial transcriptome data.|
|sc.obj    |Seurat object of scRNA-seq data.|
|coord     |Coordinates column names in ST images slot.coord = c("x","y") or coord = c("imagerow","imagecol").|
|celltype.column|The column name for cell type in the single-cell Seurat object, with the default value as "idents".|
|sc.sub.size|Downsampling proportion or number for scRNA-seq data. Default: NULL.|
|min.sc.cells|The minimum number of cell types in scRNA-seq data.Default: 50.|
|factor.size|Factor size for scaling the weight of gene expression. Default: 0.1.|
|seed.num|Number of seed genes of each cell type for recognizing candidate markers. Default: 10.|
|pvalue.cut|Threshold for filtering cell type marker genes. Default: 0.1.|
|knn|The number of nearest neighboring single cells for each spot. Set to 5 for low-resolution data and 1 for high-resolution data.|
|mean.cell.num|The average number of single cells in the spot.Set to 5 for low-resolution data and 1 for high-resolution data.|
|max.cell.num|The maximum number of cells within each spot, if equal to 1, indicates that each spot contains only a single cell.|
|n.workers|Number of cores to be used for parallel processing. Default: 4.|
|verbose|Show running messages or not. Default: TRUE.|

## Run CellMapper  to assign single cells on 10X Visium human HER2+ breast data.
Details of the results is described in the table below.
|**output**|**details**|
|------|-------|
|sc.out|Seurat object of spatial transcriptomic data with single-cell resolution.|
|decon |The cellular composition of each spot in tissue sections.|

	results <-  CellMap(st.obj = st.obj,
                        sc.obj = sc.obj,
		        coord = c("imagerow","imagecol"),
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
   
     [INFO] Identification of cell type-specific genes...
     [INFO] Estimate the number of single cells in the spot
     [INFO] Integrate single-cell and spatial spot data
     [INFO] Train a random forest model and predict,waiting...
     [INFO] Map single cells onto spatial spots
     [INFO] Construct Seurat object
     [INFO] Finish!
  
## Visualization 

```
colors <-c("B-cells" = "#e68fac","CAFs" = "#a1caf1","Cancer Epithelial" = "#f7b565","Endothelial" = "#875692",
           "Myeloid" = "#d14c6f","Normal Epithelial" = "#894846","Plasmablasts" = "#848482","PVL" = "#56af8f","T-cells" = "#0067a5")

p1 <- DimPlot(sc.data,group.by= "celltype_major",label = T,label.size = 6,
              cols = colors, pt.size = 1.5 , repel = T ) + 
              NoLegend() + labs(x = "UMAP1",y = "UMAP2", title = "CellType") +
              theme(panel.border = element_rect(fill=NA,color= "black",size= 1,linetype="solid"))+
              theme(axis.title.x =element_text(size=24), axis.title.y=element_text(size=24))+
              theme(plot.title = element_text(hjust = 0.5,size = 20, face = "bold"),
              axis.text=element_text(size=12,face = "bold"),
              axis.title.x=element_text(size=14),
              axis.title.y=element_text(size=14))

p2 <- SpatialDimPlot(results$sc.out, group.by = "CellType", pt.size.factor = 1, label.size = 8, cols = colors) + 
                   theme(legend.title = element_text(size = 14),  
                   legend.text = element_text(size = 12))

p1 + p2
```
![image](https://github.com/liuhong-jia/CellMapper/blob/main/vignettes/mapping.png)

## Run CellMap to assign single cells on high-resolution ST data ,such as Slide-seq V2,Stereo-seq,Visium HD and Imaging-based ST platform .
