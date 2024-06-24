#!/usr/bin/env Rscript

# v 0.6

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check the number of arguments
if (length(args) > 1) {
  stop("Zero or one argument must be supplied (optional CSV file).", call. = FALSE)
}

optional_csv_file <- ifelse(length(args) == 1, args[1], NA)

# If optional CSV file is provided, load it
if (!is.na(optional_csv_file)) {
  optional_data <- read.csv(optional_csv_file, row.names=1)
  # cat(optional_data)
  # Perform additional analysis with optional_data
  cat("Optional CSV file provided: ", optional_csv_file, "\n")
}

library(Seurat)
library(dplyr)
library(SoupX)
source("/users/ds286q/project0001/Dom/pipeline/commonFunctions.R")

# Search the pipeline folder for mapped data
list.data = c() 
to.skip = c() 
path.to.read = "/users/ds286q/project0001/Dom/pipeline/mapped" 
runs = list.dirs(path = path.to.read, full.names = FALSE, recursive = FALSE) 

date=Sys.Date()
pdf.filename = paste0(date,"output.pdf") 
pdf(paste0("/users/ds286q/project0001/Dom/pipeline/mapped/",pdf.filename))

cat("Reading 10X folders with SoupX\n")

#Read the 10X folder, with SoupX
for (name in runs) {
  if (name %in% to.skip) {
  } else {
    cat(name)
    cat("\n")
    path_name=paste(path.to.read, name, sep="/")
    path_name=paste0(path_name,"/outs")
    #print(path_name)
    #list.data[[name]]=Read10X(path_name)
    list.data[[name]]=SoupX.clean.from.CellRanger(path_name)
    cat("\n")
  }
}
rm(runs)
gc()
#length(list.data)

#Create seurat object
Seurat.list = c() 
for (i in 1:length(list.data)) {
  #print(i)
  name = names(list.data[i])
  #print(name)
  Seurat.object = CreateSeuratObject(counts = list.data[[i]], project = name)
  Seurat.object = Calc.Perc.Features(Seurat.object)

  Seurat.list[[name]] = Seurat.object
  rm(Seurat.object)
  gc()
}
length(Seurat.list)==length(list.data)

merged.Seurat = merge(Seurat.list[[1]], y=Seurat.list[2:length(Seurat.list)])
cat("Before QC\n")
plot(VlnPlot(merged.Seurat, features = "nCount_RNA", split.by = "orig.ident"))
plot(VlnPlot(merged.Seurat, features = "nFeature_RNA", split.by = "orig.ident"))
plot(VlnPlot(merged.Seurat, features = "percent.mt", split.by = "orig.ident"))
plot(VlnPlot(merged.Seurat, features = "percent.hb", split.by = "orig.ident"))
plot(VlnPlot(merged.Seurat, features = "percent.ribo", split.by = "orig.ident"))
plot(VlnPlot(merged.Seurat, features = "percent.MALAT1", split.by = "orig.ident"))
rm(merged.Seurat)
gc()

# QC, Norm, Scale, RunPCA
for (i in 1:length(Seurat.list)) {
  #print(i)
  Seurat.list[[i]] = QC.n.mad(Seurat.list[[i]], n.mad = 5)
  
  #TO DO: save plots?
  #plot1 <- FeatureScatter(GSK1.list[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  #plot2 <- FeatureScatter(GSK1.list[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #plot3 <- FeatureScatter(GSK1.list[[i]], feature1 = "percent.MALAT1", feature2 = "percent.mt")
  #plot4 <- VlnPlot(GSK1.list[[i]], features = c("percent.ribo", "percent.hb"), ncol = 2)
  #plot <- ((plot1 + plot2) / (plot3 + plot4))
  #plot = plot + plot_annotation(title = names(GSK1.list[[i]]), theme = theme(plot.title = element_text(hjust = 0.5)))
  #print(plot)

  Seurat.list[[i]] <- NormalizeData(Seurat.list[[i]])
  Seurat.list[[i]] <- FindVariableFeatures(Seurat.list[[i]], selection.method = "vst", nfeatures = 2000)
  Seurat.list[[i]] <- ScaleData(Seurat.list[[i]])
  Seurat.list[[i]] <- RunPCA(object = Seurat.list[[i]], npcs = 50)

  #Find number of PCA to use explainig 99% variability (assuming npcs used is 100%)
  min.pca=find.significant.PCs(Seurat.list[[i]], 0.95)
  print(paste("PCA to be used", min.pca , sep=" "))

  Seurat.list[[i]] <- RunUMAP(object = Seurat.list[[i]], dims = 1:min.pca)
  Seurat.list[[i]] <- FindNeighbors(Seurat.list[[i]], dims = 1:min.pca) %>% FindClusters(resolution = 0.1)

  #plot = DimPlot(Seurat.list[[i]])
  #plot = plot + plot_annotation(title = names(Seurat.list[i]), theme = theme(plot.title = element_text(hjust = 0.5)))
  #print(plot)

  #Find doublets
  Seurat.list[[i]] = DoubletMark(Seurat.list[[i]], n.cell.recovered = Seurat.list[[i]]@misc$cell.recovered, dim=min.pca)
  #plot1 <- DimPlot(GSK1.list[[i]], group.by = "Doublets Low stringency")
  #plot2 <- DimPlot(GSK1.list[[i]], group.by = "Doublets High stringency")
  #plot <- (plot1 + plot2)}
  #plot = plot + plot_annotation(title = names(GSK1.list[[i]]), theme = theme(plot.title = 
  #element_text(hjust = 0.5)))
  #print(plot)
}

length(Seurat.list)==length(list.data)
rm(list.data)
gc()

merged.Seurat = merge(Seurat.list[[1]], y=Seurat.list[2:length(Seurat.list)])
cat("After QC\n")
plot(VlnPlot(merged.Seurat, features = "nCount_RNA", split.by = "orig.ident"))
plot(VlnPlot(merged.Seurat, features = "nFeature_RNA", split.by = "orig.ident"))
plot(VlnPlot(merged.Seurat, features = "percent.mt", split.by = "orig.ident"))
plot(VlnPlot(merged.Seurat, features = "percent.hb", split.by = "orig.ident"))
plot(VlnPlot(merged.Seurat, features = "percent.ribo", split.by = "orig.ident"))
plot(VlnPlot(merged.Seurat, features = "percent.MALAT1", split.by = "orig.ident"))
#rm(merged.Seurat)
gc()

dev.off()

if (!is.na(optional_csv_file)) {
	Idents(merged.Seurat) <- "orig.ident"
	merged.Seurat <- make.add.meta(merged.Seurat, optional_data)
}

#save
filename=paste0(date,"_Seurat.merged.rds")
#BUG: this path need to be relative
saveRDS(merged.Seurat, paste0("/users/ds286q/project0001/Dom/pipeline/mapped/",filename))
