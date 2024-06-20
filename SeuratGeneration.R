#v 0.2
library(Seurat)
library(dplyr)
library(SoupX)
source("/users/ds286q/project0001/Dom/pipeline/commonFunctions.R")

# Search the pipeline folder for mapped data
list.data = c() 
to.skip = c() 
path.to.read = "/users/ds286q/project0001/Dom/pipeline/mapped" 
runs = list.dirs(path = path.to.read, full.names = FALSE, recursive = FALSE) 

#Read the 10X folder
for (name in runs) {
  if (name %in% to.skip) {
  } else {
    print(name)
    path_name=paste(path.to.read, name, sep="/")
    path_name=paste0(path_name,"/outs")
    #print(path_name)
    #list.data[[name]]=Read10X(path_name)
    list.data[[name]]=SoupX.clean.from.CellRanger(path_name)
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
}

length(Seurat.list)==length(list.data) 
rm(list.data) 
gc()

#save
saveRDS(Seurat.list, "/users/ds286q/project0001/Dom/pipeline/mapped/Seurat.list.rds")
