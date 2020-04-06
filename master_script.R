#start Seurat tutorial for Guided Clustering 
#load packages make sure that Seurat is loaded or else you get errors with 10X command
install.packages("Seurat")
library(dplyr)
library(Seurat)
library(patchwork)
#Read the data from the desktop folder and apply the Read10X function which reads in the pipeline from the 10X database--need to tell R that Read10X comes from Seurat
pbmc.data <- Seurat::Read10X(data.dir ="/Users/paigehalas/Desktop/filtered_gene_bc_matrices/hg19/Seraut_tutorial/")
#creates a Seurat object where counts is your data set, the project is next, then you add in the minimum amount of cells and features you want to detect--this could be dictated by the project itself) 
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbm3ck", min.cells = 3, min.features = 200)
#print 
pbmc
#put at symbol after object name to see all of the objects associated with it 
#with Ryan--can add a metadata column by putting your dataset name [["column.name"]] <-
pbmc[["column_name"]] <- test
pbmc(#start QC process-first is looking at percentage of features that are mitochondrial genes 
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#violin plots of  where we look at RNA features and the count of RNA, and then total percent mitochondrial
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Make feature scatter plots to look at the relationships between features: first is RNA versus the amount that are mitochondrial
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
#make feature scatter plot to look between features: count and features 
plot2 <- FeatureScatter(pbmc, feature1= "nCount_RNA", feature2 = "nFeature_RNA")
plot2
#display both plots together
plot1 + plot2
#filter the cells so that we eliminate ones with less than 200 and greater than 2500 features and remove cells with greater than 5% mitochondrial genes--this was based on violin plots from before
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc
#plot the filtered data to see if it is filter correctly
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
plot2 <- FeatureScatter(pbmc, feature1= "nCount_RNA", feature2 = "nFeature_RNA")
plot2
plot1 + plot2
#normalize filtered data on a log scale (pbmc has been modified) 
pbmcn <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#Identifying variable features that have greatest variability (vst) (highly or lowly expressed) 
pbmcvf <- FindVariableFeatures(pbmcn, selection.method = "vst", nfeatures = 2000)
#Identify the top 10 highly variable genes 
top10 <- head(VariableFeatures(pbmcvf), 10)
top10
#plot without variable names 
plot1 <- VariableFeaturePlot(pbmcvf)
plot1
#plot with variable names
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
plot1 + plot2
#scale the data to shift expression so that expression is 0 (rather than the giant negative number given before) and variance is 1 
all.genes <- rownames(pbmcvf)
pbmc <- ScaleData(pbmcvf, features = all.genes)
pbmc
#perform PCA (principle component analysis) can reduce dimensions of the variables to see which are contributing more than others
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#visualize what the PCA sees which is being contributed and which is not contributing as much (each PC is like a pie chart that needs to add up to 100)
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
#visualize plot of PCAs 
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
#visualize using DimPlot 
DimPlot(pbmc, reduction = "pca")
#visualize using DimHeatmap where it is plotted for 500 cells for one PC (dimension) and is balanced
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
#visualize using DimHeadmap where it is plotted for 500 cells for multiple PCs (1-15)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
#next we determine the dimensionality of the dataset because there are PCs that are picked out that do not have components that arent as critical-want to filter to just the ones that are critical (p values that are low) 
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
#Visualize p values within each PC, significant features have low pvalue (curves closest to black dashed line) for PCs 1-15 
JackStrawPlot(pbmc, dims = 1:15)
#Can also visualize these using an Elbow plot where we can see where the true signal is based on the Elbow of the plot
ElbowPlot(pbmc)
#begin to cluser the cells--first use K Nearest Neighbor use the 10 dimentionsthat was found to be the elbow of the PCA 
#can continue to find more clusters off of this object and you can continue to store this in metadata
pbmc <- FindNeighbors(pbmc, dims = 1:10)
#next identify the clusters from this KNN 
pbmc <- FindClusters(pbmc, resolution = 0.5)
#play around to see what the clusters with look at with a resolution of 3
pbmc <- FindClusters(pbmc, resolution = 3)
#look at the cluster IDs of the first 5 cells 
head(Idents(pbmc), 5)
#Install UMAP 
install.packages('umap')
library(umap)
#Run UMAP analysis of the data with the relevent PCs (1-10)
pbmc <- RunUMAP(pbmc, dims = 1:10)
#Plot the PCs so that you get clusters 
DimPlot(pbmc, reduction = "umap")
#Now we are doing to define the clusters by looking at what genes are differentially expressed (what are the biomarkers) looking at it with the minimum.pct argument which is the minimul percentage of two groups of cells 
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
#look at the top 5 markers 
head(cluster1.markers, n = 5)
#can also look at the markers for clusters between eachother (here we are looking at differences between cluster 5 with cluser 0 and 3)
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
#can find markers for every cluster compared to all the remaining cells (will calculate for all clusters one by one) and report only positive ones 
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
#Using ROC to classify differential expression of individual markers
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#violin plot to visualize gene expression
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
#can also plot raw counts 
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
#can look multiple genes at once 
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))
#can use a heat map to look at the top markers (in this case 10)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
#use cononical markers to match unbiased clustering to known cell types 
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#start tutorial with Ryan 
#Look at the barcodes that are associated with each cluster
pbmc@active.ident
head(pbmc)
head(pbmc@meta.data)
#look at random genes that are expressed 
rownames(pbmc@assays$RNA)[1:5]
#made a list of genes that are expressed in the blood samples 
genes <- list("random" = c("AL627309.1","AP006222.2", "RP11-206L10.2"))
#score genes based on how they are expressed
pbmc <- AddModuleScore(pbmc, genes, name = "Module")
#check to see if the genes are scores 
head(pbmc@meta.data)
#plot the genes in multiple ways and group by the RNA resolution 0.5 where applicable 
VlnPlot(pbmc, feature = "Module1", group.by = "RNA_snn_res.0.5")
FeaturePlot(pbmc, feature = "Module1")
DotPlot(pbmc, feature = "Module1", group.by = "RNA_snn_res.0.5")
#not as commonly used because not as informative 
RidgePlot(pbmc, feature = "Module1", group.by = "RNA_snn_res.0.5")

#make two lists of genes because you might want to categorize genes based on module to plot them 
genes <- list("random" = c("AL627309.1","AP006222.2", "RP11-206L10.2"),"Random2" = c("AL627309.1"))
pbmc <- AddModuleScore(pbmc, genes, name = "Module")
head(pbmc@meta.data)
pbmc@meta.data
#old columns are still here so we had to remove them, we used the - function and removed the columns 
pbmc@meta.data <- pbmc@meta.data[,-c(9:11)]
head(pbmc@meta.data)]
#try it again to see if it will print both lists 
pbmc <- AddModuleScore(pbmc, genes, name = "Module")
head(pbmc@meta.data)
#try to print it again with both lists and see if the columns are named after the genes if you dont name them (they arent)
pbmc <- AddModuleScore(pbmc, genes)
head(pbmc@meta.data)

