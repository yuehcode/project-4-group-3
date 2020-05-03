


library(Seurat)
library(tximport)
library(dplyr)
library(Matrix)
library(ggplot2)
library(biomaRt)
library(stringr)

raw_counts<- read.table(file = "final_mtx.csv", header = T, row.names=1, sep=",",  stringsAsFactors = FALSE)
rownames(raw_counts) 
head(raw_counts)

dim(raw_counts)

# check if there is any duplicated feactures
sum(duplicated(rownames(raw_counts)))

# get the rownames
id<-rownames(raw_counts)

# remove digits &.
geneid<- str_replace(id, pattern = ".[0-9]+$",replacement = "")
# check value
sum(duplicated(geneid))

###transfer to  gene symbols
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
transtogene <- getBM(attributes=c("ensembl_gene_id","external_gene_name"),
                     filters= "ensembl_gene_id", values = geneid, mart = ensembl)


# get a match list
matchedgene <- match(geneid,transtogene$ensembl_gene_id)
sum(duplicated(matchedgene))
matchedgene <- transtogene$external_gene_name[matchedgene]
sum(duplicated((matchedgene)))#33

#change datatype and add rownames
raw_counts_ma<-as.matrix(raw_counts)
rownames(raw_counts_ma) <- matchedgene

# check dupicated value & remove
sum(duplicated(rownames(raw_counts_ma)))
noduplicated=raw_counts_ma[which(!duplicated(rownames(raw_counts_ma))),]

sum(duplicated(rownames(noduplicated)))


mydata <- CreateSeuratObject(counts = noduplicated, min.cells = 3, min.features = 200, project = "project4")
dim(mydata@meta.data)
head(mydata@meta.data)
mydata
mydata[["percent.mt"]] <- PercentageFeatureSet(mydata, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(mydata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

mydata <- subset(mydata, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
mydata 

# Normalize the data
mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 10000)
head(mydata@meta.data)

# feature selesction & filter out low variance genes
mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mydata), 10)
top10
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mydata)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

#Scaling the data
all.genes <- rownames(mydata)
mydata <- ScaleData(mydata, features = all.genes)

# Perform linear dimensional reduction
mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata))
print(mydata[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mydata, dims = 1:2, reduction = "pca")
DimPlot(mydata, reduction = "pca")
DimHeatmap(mydata, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(mydata, dims = 1:9, cells = 500, balanced = TRUE)

# computation time
mydata <- JackStraw(mydata, num.replicate = 100)
mydata <- ScoreJackStraw(mydata, dims = 1:20)
JackStrawPlot(mydata, dims = 1:15)
ElbowPlot(mydata)

mydata <- FindNeighbors(mydata, dims = 1:10)
mydata <- FindClusters(mydata, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(mydata), 5)

# Non-Linear Dimensional Reduction umap
mydata <- RunUMAP(mydata, dims = 1:10)

# make umap & save result
# individual clusters
DimPlot(mydata, reduction = "umap", label = T)
saveRDS(mydata, file = "10final.rds")

# make grapgh
makegraph <- (table(Idents(mydata)))
