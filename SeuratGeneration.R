#!/usr/bin/env Rscript

# v 0.7

### Define functions ###

SoupX.clean.from.CellRanger <- function(cellranger.folder) {
  #Clean your scRNAseq data from ambient contamination
  #Use:
  #clean.data = SoupX.clean(cellranger.folder="/my/cellranger/folder/")

  # if (dir.exists(paste(cellranger.folder,"filtered_feature_bc_matrix",sep = "/"))==FALSE){
  #   stop("filtered_feature_bc_matrix folder didn't find")
  # }
  #
  # if (dir.exists(paste(cellranger.folder,"raw_feature_bc_matrix",sep = "/"))==FALSE){
  #   stop("raw_feature_bc_matrix folder didn't find")
  # }
  #
  # if (dir.exists(paste(cellranger.folder,"analysis",sep = "/"))==FALSE){
  #   stop("analysis folder didn't find")
  # }

  #Load data
  sc = SoupX::load10X(cellranger.folder)
  #Estimante the contamination
  sc = SoupX::autoEstCont(sc)

  print("Genes with highest expression in background:")
  print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))

  #Filterout the contamination
  out = SoupX::adjustCounts(sc)
  return(out)
}

SoupX.on.Seurat <- function(Seurat.object, cellranger.folder, min.cells = 5, min.features = 50) {
  #Clean your scRNAseq data from ambient contamination
  #Seurat object should be UNFILTERED (MT genes and doublets), and seurat_clusters should be present
  #Use:
  #clean.data = SoupX.clean(Seurat.object, cellranger.folder="/my/cellranger/folder/")

  if (dir.exists(paste(cellranger.folder,"filtered_feature_bc_matrix",sep = "/"))==FALSE){
    stop("filtered_feature_bc_matrix folder didn't find")
  }

  if (dir.exists(paste(cellranger.folder,"raw_feature_bc_matrix",sep = "/"))==FALSE){
    stop("raw_feature_bc_matrix folder didn't find")
  }

  if (any(colnames(Seurat.object@meta.data) == "seurat_clusters")) {
    print("seurat_clusters present...")
  } else {
    stop("seurat_clusters not present")
  }

  raw.matrix = Seurat::Read10X(paste(cellranger.folder,"raw_feature_bc_matrix",sep = "/"))
  filt.matrix = Seurat::Read10X(paste(cellranger.folder,"filtered_feature_bc_matrix",sep = "/"))

  sc  <- SoupX::SoupChannel(raw.matrix, filt.matrix)

  meta    <- Seurat.object@meta.data
  umap    <- Seurat.object@reductions$umap@cell.embeddings
  sc  <- SoupX::setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc  <- SoupX::setDR(sc, umap)
  head(meta)

  #With defined clusters, run the main SoupX function, calculating ambient RNA profile.
  sc  <- SoupX::autoEstCont(sc)

  print("Genes with highest expression in background:")
  print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))

  adj.matrix = SoupX::adjustCounts(sc, roundToInt = T)

  New.Seurat.object = Seurat::CreateSeuratObject(counts = adj.matrix, min.cells = 5, min.features = 50, meta.data = meta)
  return(New.Seurat.object)
}

make.add.meta <- function(Seurat.Object, metadata, return.only.table=FALSE, verbose=FALSE) {
  #From a metadata table and a Seurat.Object,
  #it creates a proper metadata table with cell barcodeID and add it to Seurat.object
  # Use:
  # Seurat.object = make.add.meta(Seurat.Object, metadata)

  if (verbose) {
    print(head(metadata))
  }

  #if there is only 1 row, add it to the whole Seurat.Object
  if (nrow(metadata)==1) {
    if (verbose) {
      print("nrow(metadata)==1")
    }
    df.cells <- data.frame(row.names = colnames(Seurat.Object))
    for (name in colnames(metadata)) {
      #print(name)
      df.cells[name]=metadata[name]
    }
  #if there is only 1 column it's a list
  } else if (ncol(metadata)==1) {
    if (verbose){
      print(colnames(metadata))
      print(names(metadata))
    }

    #check if length of idents is different of lenght metadata
    if (length(setdiff(Idents(Seurat.Object), rownames(metadata)))!=0 || length(setdiff(rownames(metadata),Idents(Seurat.Object)))!=0) {
      stop("Seurat object Idents and metadata rows are not matching.")
    }

    df.cells <- data.frame()

    #Select the Idents matching “Name” in metadata
    #Idents(Seurat.Object) <- "Condition_name"

    #For each idents within the seurat object (cluster or sample you want)
    for (name in unique(Idents(Seurat.Object)))
    {
      print(name)

      #Select the cells and put them in another df
      new_df <- data.frame(row.names = WhichCells(Seurat.Object, idents = name))

      #Select the row of interest from metadata, corresponding to the metadata to add to those cells
      meta_row=metadata[rownames(metadata) == name,]
      new_df[colnames(metadata)]=meta_row
      df.cells=dplyr::bind_rows(df.cells, new_df)
      }
      #merge in a big df
  #if else is a dataframe
  } else {
    if (verbose) {
      print("nrow(metadata)!=1")
    }
    #check if length of idents is different of lenght metadata
    if (length(setdiff(Idents(Seurat.Object), rownames(metadata)))!=0 || length(setdiff(rownames(metadata),Idents(Seurat.Object)))!=0) {
      stop("Seurat object Idents and metadata rows are not matching.")
    }

    df.cells <- data.frame()

    #Select the Idents matching “Name” in metadata
    #Idents(Seurat.Object) <- "Condition_name"

    #For each idents within the seurat object (cluster or sample you want)
    for (name in unique(Idents(Seurat.Object)))
    {
      print(name)

      #Select the cells and put them in another df
      new_df <- data.frame(row.names = WhichCells(Seurat.Object, idents = name))

      #Select the row of interest from metadata, corresponding to the metadata to add to those cells
      meta_row=metadata[rownames(metadata) == name,]
      if (verbose) {
        #print(metadata)
        #print(rownames(metadata))
        print("meta_row:")
        print(meta_row)
        #print(nrow(meta_row))
        print(colnames(meta_row))
        print(ncol(meta_row))
      }
      if (is.null(colnames(meta_row))) {
        stop("Contact Dom")
      }
      if (ncol(meta_row)==1) {
        stop("There is only 1 column...contact dom, need to be fixed")
      }
      for (col_name in colnames(meta_row)){
        #add the specific metadata you need
        if (verbose) {
          print("col_name is:")
          print(col_name)
        }
        new_df[col_name]=meta_row[col_name]
      }
      #merge in a big df
      df.cells=dplyr::bind_rows(df.cells, new_df)
    }
  }

  #Finally add to the object

  if (return.only.table==TRUE) {
    return(df.cells)
  } else {
    if (verbose) {
      print("Adding to Seurat...")
    }
    Seurat.Object <- AddMetaData(Seurat.Object, metadata = df.cells)
    return(Seurat.Object)
  }
}

plot.depth.seq <- function(Seurat.object, metadata.col=NULL) {
  #metadata.col need to be a metadata column in Seurat.object@meta.data

  if (is.null(metadata.col)) {
    metadata.col = "orig.ident"
  }
  Idents(Seurat.object) <- metadata.col
  libs = unique(Idents(Seurat.object))

  for (i in 1:length(libs)) {
    lib.seurat = subset(Seurat.object, idents = libs[i])

    colSums.counts = as.data.frame(colSums(lib.seurat@assays$RNA@counts))
    colSums.counts$colSums <- colSums.counts$`colSums(lib.seurat@assays$RNA@counts)`
    colSums.counts$`colSums(lib.seurat@assays$RNA@counts)` <- NULL
    #print(colnames(colSums.counts))
    a=ggplot(colSums.counts, aes(x=colSums))+ geom_histogram(color="black", fill="white", bins=50) + ggtitle(libs[[i]])
    b=ggplot(colSums.counts, aes(x=colSums))+ geom_histogram(color="black", fill="white", bins=50, aes(y=..density..))+geom_density(alpha=.2, fill="#FF6666") + ggtitle(libs[[i]])
    print(a+b)
    #print(ggplot(colSums.counts, aes(x=colSums)) + geom_histogram(color="black", fill="white", bins=50, aes(y=..density..)) + geom_density(alpha=.2, fill="#FF6666") + ggtitle(unique(Idents(libs[[i]]))))
  }
}

Mark.cells <- function(Seurat.object, nFeature.low = 250, nFeature.high = 5000, mt.high = 25, nCount.high = 18000) {
  #Mark cells with low/high nFeature, nCount.high, mt genes,
  #Example use: Seurat.Object <- Mark.cells(Seurat.Object)

  # Known issues, to improve:
  # this function works on "percent.mt", "nFeature_RNA" and "nCount_RNA" metadata columns...need to be more flexible

  poscells.high.mt <- WhichCells(Seurat.object, expression = percent.mt > mt.high)
  Seurat.object$high.mt<- ifelse(colnames(Seurat.object) %in% poscells.high.mt, "Pos", "Neg")
  print("high mt cells:")
  print(table(Seurat.object$high.mt))

  poscells.high.ft <- WhichCells(Seurat.object, expression = nFeature_RNA > nFeature.high)
  Seurat.object$high.nFeature <- ifelse(colnames(Seurat.object) %in% poscells.high.ft, "Pos", "Neg")
  print("high nFeature cells:")
  print(table(Seurat.object$high.nFeature))

  poscells.low.ft <- WhichCells(Seurat.object, expression = nFeature_RNA < nFeature.low)
  Seurat.object$low.nFeature <- ifelse(colnames(Seurat.object) %in% poscells.low.ft, "Pos", "Neg")
  print("low nFeature cells:")
  print(table(Seurat.object$low.nFeature))

  poscells.high.nc <- WhichCells(Seurat.object, expression = nCount_RNA > nCount.high)
  Seurat.object$high.nCount <- ifelse(colnames(Seurat.object) %in% poscells.high.nc, "Pos", "Neg")
  print("high UMIs:")
  print(table(Seurat.object$high.nCount))

  return(Seurat.object)
}

find.significant.PCs <- function(Seurat.Object, variance=0.9, st.dev=0.05) {
  ##Find a min PCA significant (90% variance & st.dev<5%) to use
  ##The Seurat.Object need "pca" slot filleds
  ##Example use:
  ## min.pca = find.significant.PCs(Seurat.Object)

  if (variance>=1) {
    simpleError("variance over 100%")
  }
  variance=variance*100
  st.dev=st.dev*100

  stdv <- Seurat.Object[["pca"]]@stdev
  sum.stdv <- sum(Seurat.Object[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > variance & percent.stdv < st.dev)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] -
                       percent.stdv[2:length(percent.stdv)]) > 1-(variance/100)),
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  print(min.pc)
  return(min.pc)
}

Calc.Perc.Features <- function(Seurat.object, mt.pattern = "^MT-", hb.pattern = "^HB[^(P)]", ribo.pattern = ("RPS|RPL"), MALAT1.name="MALAT1", plot.name="") {
  #alternative version: Mod for not plotting
  Seurat.object@misc$cell.recovered = ncol(Seurat.object)
  Seurat.object[["percent.mt"]] <- PercentageFeatureSet(Seurat.object, pattern = mt.pattern)
  Seurat.object[["percent.hb"]] <- PercentageFeatureSet(Seurat.object, pattern = hb.pattern)
  Seurat.object[["percent.ribo"]] <- PercentageFeatureSet(Seurat.object, pattern = ribo.pattern)
  Seurat.object[["percent.MALAT1"]] <- PercentageFeatureSet(Seurat.object, pattern = MALAT1.name)
  #plot1 <- FeatureScatter(Seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
  #plot2 <- FeatureScatter(Seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #plot3 <- FeatureScatter(Seurat.object, feature1 = "percent.MALAT1", feature2 = "percent.mt")
  #plot4 <- VlnPlot(Seurat.object, features = c("percent.ribo", "percent.hb"), ncol = 2)
  #plot <- ((plot1 + plot2) / (plot3 + plot4))
  #plot = plot + plot_annotation(title = plot.name, theme = theme(plot.title = element_text(hjust = 0.5)))
  #print(plot)
  return(Seurat.object)
}

QC.n.mad <- function(Seurat.object, n.mad=4) {
  #mod version for not plotting
  #Based on https://matthieuxmoreau.github.io/EarlyPallialNeurogenesis/html-Reports/Quality_Control.html code
  Cell.QC.Stat <- Seurat.object@meta.data
  print(nrow(Cell.QC.Stat))
  # high and low median absolute deviation (mad) thresholds to exclude outlier cells
  max.mito.thr <- median(Cell.QC.Stat$percent.mt) + n.mad*mad(Cell.QC.Stat$percent.mt)
  min.mito.thr <- median(Cell.QC.Stat$percent.mt) - n.mad*mad(Cell.QC.Stat$percent.mt)

  #Plot
  #p1 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mt)) +
  #  geom_point() +
  #  geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
  #  geom_hline(aes(yintercept = min.mito.thr), colour = "red", linetype = 2) +
  #  annotate(geom = "text", label = paste0(as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr | Cell.QC.Stat$percent.mt < min.mito.thr)[2])," cells removed\n", as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr | Cell.QC.Stat$percent.mt < min.mito.thr)[1])," cells remain"), x = 6000, y = 0.1)

  Cell.QC.Stat <- Cell.QC.Stat %>% dplyr::filter(percent.mt < max.mito.thr) %>% dplyr::filter(percent.mt > min.mito.thr)

  # Set low and hight thresholds on the number of detected genes
  min.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) - n.mad*mad(log10(Cell.QC.Stat$nFeature_RNA))
  max.Genes.thr <- median(log10(Cell.QC.Stat$nFeature_RNA)) + n.mad*mad(log10(Cell.QC.Stat$nFeature_RNA))
  # Set hight threshold on the number of transcripts
  max.nUMI.thr <- median(log10(Cell.QC.Stat$nCount_RNA)) + n.mad*mad(log10(Cell.QC.Stat$nCount_RNA))

  #p2 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
  #  geom_point() +
  #  geom_smooth(method="lm") +
  #  geom_hline(aes(yintercept = min.Genes.thr), colour = "green", linetype = 2) +
  #  geom_hline(aes(yintercept = max.Genes.thr), colour = "green", linetype = 2) +
  #  geom_vline(aes(xintercept = max.nUMI.thr), colour = "red", linetype = 2)

  #p3 <- p1 / p2
  #p1=ggExtra::ggMarginal(p1, type = "histogram", fill="lightgrey", bins=100)
  #p2=ggExtra::ggMarginal(p2, type = "histogram", fill="lightgrey")
  #print(typeof(p1))
  #p3 <- p1 / p2
  #p3 <- p3 + plot_annotation(title = names(Seurat.object), theme = theme(plot.title = element_text(hjust = 0.5)))
  #print(p3)

  # Filter cells based on these thresholds

  Cell.QC.Stat <- Cell.QC.Stat %>% dplyr::filter(log10(nFeature_RNA) > min.Genes.thr) %>% dplyr::filter(log10(nFeature_RNA) < max.Genes.thr) %>% dplyr::filter(log10(nCount_RNA) < max.nUMI.thr)
  print(nrow(Cell.QC.Stat))
  print("########")
  Seurat.object <- subset(Seurat.object, cells = rownames(Cell.QC.Stat))
  return(Seurat.object)
}

#Function to plot PCA of samples, return an avarageExpression Seurat object
PCA.sample <- function(Seurat.Object, sample.name = "orig.ident", nfeatures = 2000, show.label=FALSE, legend=TRUE){
        ##Thif function plot PCA of samples, return an avarageExpression Seurat object
        ##Example use:
  ##My.new.object = PCA.sample(Seurat.Object, sample.name = "Patient", nfeatures = 1000)

        ##Describe the code
  Idents(object = Seurat.Object) <- sample.name
  sample.object <- AverageExpression(object = Seurat.Object, return.seurat = TRUE)
  sample.object <- NormalizeData(sample.object)
  sample.object <- FindVariableFeatures(sample.object, selection.method = "vst", nfeatures = nfeatures)
  sample.object <- ScaleData(sample.object)
  sample.object <- RunPCA(object = sample.object, npcs = 2)
  if (legend == TRUE) {
    print(DimPlot(sample.object, reduction = "pca", pt.size = 5, label = show.label))
  } else {
    print(DimPlot(sample.object, reduction = "pca", pt.size = 5, label = show.label)+ NoLegend())
  }

  return(sample.object)
  #You can load metadata to this seurat object, and plot factors with group.by
}


##############
### Script ###
##############


# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check the number of arguments
if (length(args) > 2) {
  stop("one or two argument must be supplied (folder path containing mapped data and optional CSV metadata file).", call. = FALSE)
}

path.to.read <- args[1]
optional_csv_file <- ifelse(length(args) == 2, args[2], NA)

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

# Search the pipeline folder for mapped data
list.data = c() 
to.skip = c() 
#path.to.read = "/users/ds286q/project0001/Dom/pipeline/mapped" 
runs = list.dirs(path = path.to.read, full.names = FALSE, recursive = FALSE) 

date=Sys.Date()
pdf.filename = paste0(date,"_output.pdf") 
pdf(normalizePath(file.path(path.to.read,pdf.filename)))

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

  ###Find doublets
  # DOESN'T WORK AT THE MOMENT
  ##Seurat.list[[i]] = DoubletMark(Seurat.list[[i]], n.cell.recovered = Seurat.list[[i]]@misc$cell.recovered, dim=min.pca)
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
#print(paste0(path.to.read,filename))
print(normalizePath(file.path(path.to.read, filename)))
saveRDS(merged.Seurat, normalizePath(file.path(path.to.read, filename)))
