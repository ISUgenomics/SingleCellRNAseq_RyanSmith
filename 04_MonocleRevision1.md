#  Need to make some modified plots of fed and unfed cells

  I would like for you to redo some iterations of the pseudotime analysis (Figures 4A and B in the paper) to break down the samples by feeding status with/without the inclusion of Cluster 1. I have added the data files to “revision” folder in the  “Ryan_Smith_Monocle” folder on Box. There are five similar data sets included in the folder with the desired breakdowns, and would like the paired UMAP/pseudotime plots for each file.


##### Versions used
```
ml monocle3/0.2.0
singularity shell  /opt/rit/singularity/images/monocle3/0.2.0/monocle3.simg

```
### Cluster 1 naive

##### CREATE TABLES
```
/work/gif/remkv6/Smith_Ryan/01_Monocle/01_Cluster1Naive
#Some cells did not have expression in some genes, so modify AllGeneMetadata for each run.

#formatted in excel, hrmm. didnt format well in excel, fixing
awk '{print $1,$2,$3,$4}' Cl1NaiveCellMetadata.txt |tr " " "\t" >FixedCl1NaiveCellMetadata.txt
awk '{print $1}' Cl1NaiveGeneExpressionMatrix.txt |grep -f - AllGeneMetadata |cat <(awk 'NR==1' AllGeneMetadata) - >Cl1NaiveGeneMetadata
```

### Create plot
```
library(monocle3)
library(dplyr)
expression_matrix2<- round(as.matrix(read.table("Cl1NaiveGeneExpressionMatrix.txt",header = T)))
gene_metadata2 <- as.matrix(read.table("Cl1NaiveGeneMetadata",header = T,row.names = 1))
cell_metadata2<- as.matrix(read.table("FixedCl1NaiveCellMetadata.txt",header = T,row.names = 1))
library(Matrix)
M1 <- as(expression_matrix2, "dgCMatrix")
cds <- new_cell_data_set(M1,cell_metadata = cell_metadata2,gene_metadata = gene_metadata2)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds,root_cells=row.names(colData(cds[,3],on=5)))
plot_cells(cds, reduction_method="UMAP", color_cells_by="group",cell_size=1.3,label_cell_groups=FALSE,label_leaves=TRUE,graph_label_size=2,label_branch_points=TRUE)


#set the root for the gradient
get_earliest_principal_node <- function(cds,on=5){
  cell_ids <- which(colData(cds)[, "group"] == 5)

  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds, reduction_method="UMAP", color_cells_by="pseudotime",cell_size=1.3,label_cell_groups=FALSE,label_leaves=FALSE,graph_label_size=2,label_branch_points=TRUE)

```


### Cluster1 BLOOOOD!


##### prepare tables
```
#/work/gif/remkv6/Smith_Ryan/01_Monocle/02_Cluster1Blood


awk '{print $1}' Cl1BloodGeneExpressionMatrix.txt |grep -f - AllGeneMetadata |cat <(awk 'NR==1' AllGeneMetadata) - >Cl1BloodGeneMetadata

#replace cluster header with group
Cl1BloodCellMetadata.txt

#made in excel with transpose
Cl1BloodGeneExpressionMatrix.txt
```

##### Cluster cells and create plots
```
library(monocle3)
library(dplyr)
expression_matrix2<- round(as.matrix(read.table("Cl1BloodGeneExpressionMatrix.txt",header = T)))
gene_metadata2 <- as.matrix(read.table("Cl1BloodGeneMetadata",header = T,row.names = 1))
cell_metadata2<- as.matrix(read.table("Cl1BloodCellMetadata.txt",header = T,row.names = 1))
library(Matrix)
M1 <- as(expression_matrix2, "dgCMatrix")
cds <- new_cell_data_set(M1,cell_metadata = cell_metadata2,gene_metadata = gene_metadata2)

cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds,root_cells=row.names(colData(cds[,3],on=5)))
plot_cells(cds, reduction_method="UMAP", color_cells_by="group",cell_size=1.3,label_cell_groups=FALSE,label_leaves=TRUE,graph_label_size=2,label_branch_points=TRUE)


#set the root for the gradient
get_earliest_principal_node <- function(cds,on=5){
  cell_ids <- which(colData(cds)[, "group"] == 5)

  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds, reduction_method="UMAP", color_cells_by="pseudotime",cell_size=1.3,label_cell_groups=FALSE,label_leaves=FALSE,graph_label_size=2,label_branch_points=TRUE)


```


### All naive cells


##### prepare tables
```
#/work/gif/remkv6/Smith_Ryan/01_Monocle/03_Allnaive


awk '{print $1}' AllnaiveGeneExpressionMatrix.txt |grep -f - AllGeneMetadata |cat <(awk 'NR==1' AllGeneMetadata) - >AllnaiveGeneMetadata

#replace cluster header with group
AllnaiveCellMetadata.txt

#made in excel with transpose
AllnaiveGeneExpressionMatrix.txt
```

##### Cluster cells and create plots
```
library(monocle3)
library(dplyr)
expression_matrix2<- round(as.matrix(read.table("AllnaiveGeneExpressionMatrix.txt",header = T)))
gene_metadata2 <- as.matrix(read.table("AllnaiveGeneMetadata",header = T,row.names = 1))
cell_metadata2<- as.matrix(read.table("AllnaiveCellMetadata.txt",header = T,row.names = 1))
library(Matrix)
M1 <- as(expression_matrix2, "dgCMatrix")
cds <- new_cell_data_set(M1,cell_metadata = cell_metadata2,gene_metadata = gene_metadata2)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds,root_cells=row.names(colData(cds[,3],on=5)))
plot_cells(cds, reduction_method="UMAP", color_cells_by="group",cell_size=1.3,label_cell_groups=FALSE,label_leaves=TRUE,graph_label_size=2,label_branch_points=TRUE)


#set the root for the gradient
get_earliest_principal_node <- function(cds,on=5){
  cell_ids <- which(colData(cds)[, "group"] == 5)

  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds, reduction_method="UMAP", color_cells_by="pseudotime",cell_size=1.3,label_cell_groups=FALSE,label_leaves=FALSE,graph_label_size=2,label_branch_points=TRUE)


```

### All BLOOOOOD! cells


##### prepare tables
```
#/work/gif/remkv6/Smith_Ryan/01_Monocle/04_Allblood


awk '{print $1}' AllbloodGeneExpressionMatrix.txt |grep -f - AllGeneMetadata |cat <(awk 'NR==1' AllGeneMetadata) - >AllbloodGeneMetadata.txt

#replace cluster header with group
AllbloodCellMetadata.txt

#made in excel with transpose
AllbloodGeneExpressionMatrix.txt
```

##### Cluster cells and create plots
```
library(monocle3)
library(dplyr)
expression_matrix2<- round(as.matrix(read.table("AllbloodGeneExpressionMatrix.txt",header = T)))
gene_metadata2 <- as.matrix(read.table("AllbloodGeneMetadata.txt",header = T,row.names = 1))
cell_metadata2<- as.matrix(read.table("AllbloodCellMetadata.txt",header = T,row.names = 1))
library(Matrix)
M1 <- as(expression_matrix2, "dgCMatrix")
cds <- new_cell_data_set(M1,cell_metadata = cell_metadata2,gene_metadata = gene_metadata2)

cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds,root_cells=row.names(colData(cds[,3],on=5)))
plot_cells(cds, reduction_method="UMAP", color_cells_by="group",cell_size=1.3,label_cell_groups=FALSE,label_leaves=TRUE,graph_label_size=2,label_branch_points=TRUE)


#set the root for the gradient
get_earliest_principal_node <- function(cds,on=5){
  cell_ids <- which(colData(cds)[, "group"] == 5)

  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds, reduction_method="UMAP", color_cells_by="pseudotime",cell_size=1.3,label_cell_groups=FALSE,label_leaves=FALSE,graph_label_size=2,label_branch_points=TRUE)


```



### All samples blood and naive

##### prepare tables
```
#/work/gif/remkv6/Smith_Ryan/01_Monocle/04_Allblood


awk '{print $1}' EverythingGeneExpressionMatrix.txt |grep -f - AllGeneMetadata |cat <(awk 'NR==1' AllGeneMetadata) - >EverythingGeneMetadata.txt

#replace cluster header with group
EverythingCellMetadata.txt

#made in excel with transpose
EverythingGeneExpressionMatrix.txt
```

##### Cluster cells and create plots
```
library(monocle3)
library(dplyr)
expression_matrix2<- round(as.matrix(read.table("EverythingGeneExpressionMatrix.txt",header = T)))
gene_metadata2 <- as.matrix(read.table("EverythingGeneMetadata.txt",header = T,row.names = 1))
cell_metadata2<- as.matrix(read.table("EverythingCellMetadata.txt",header = T,row.names = 1))
library(Matrix)
M1 <- as(expression_matrix2, "dgCMatrix")
cds <- new_cell_data_set(M1,cell_metadata = cell_metadata2,gene_metadata = gene_metadata2)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds,root_cells=row.names(colData(cds[,3],on=5)))
plot_cells(cds, reduction_method="UMAP", color_cells_by="group",cell_size=1.3,label_cell_groups=FALSE,label_leaves=TRUE,graph_label_size=2,label_branch_points=TRUE)


#set the root for the gradient
get_earliest_principal_node <- function(cds,on=5){
  cell_ids <- which(colData(cds)[, "group"] == 5)

  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds, reduction_method="UMAP", color_cells_by="pseudotime",cell_size=1.3,label_cell_groups=FALSE,label_leaves=FALSE,graph_label_size=2,label_branch_points=TRUE)


```
