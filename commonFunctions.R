#This is a collection of functions useful for Seurat objects.
#To load it run:
#source("commonFunctions.R")

#This is part of VlnPlot.median function
median.stat <- function(x){
	out <- quantile(x, probs = c(0.5))
	names(out) <- c("ymed")
	return(out) 
	}

VlnPlot.median <- function(object, features, colour = "black", pt.size = 0, ncol = NULL, legend = TRUE, vln.colors=NULL) {
	#Tris function creates VlnPlot plus a median 
	#VlnPlot.median(object, c("IL1B","CLEC10A"))
	
	myplots <- vector("list")
	
	#Create a plot for each gene and add median
	for (gene in features) {
		#print(gene)
		if (legend == TRUE) {
		myplots[[gene]] <- Seurat::VlnPlot(object = object, features = gene, pt.size = pt.size, cols = vln.colors) +
		stat_summary(fun = median.stat, geom='point', size = 1, colour = colour)
		} else {
		myplots[[gene]] <- Seurat::VlnPlot(object = object, features = gene, pt.size = pt.size, cols = vln.colors) +
		stat_summary(fun = median.stat, geom='point', size = 1, colour = colour) + NoLegend()
		}
	}
	#patchwork function to combine multiple plots
	patchwork::wrap_plots(myplots, ncol=ncol)
	}

Add.SNPs.HT <- function(Seurat.Object, souporcell.file, verbose=FALSE) {
  #Tris function creates will add your souporcell cluster.tsv file as metadata to a Seurat object 
  #PBMC = Add.SNPs.HT(PBMC,"Mapped/clusters.tsv")
  
  SNPs = read.csv(souporcell.file, sep = "\t")
  common_barcode = intersect(colnames(Seurat.Object), SNPs$barcode)
  print(paste("Number of cells barcoded: ",length(common_barcode)))
  
  if (length(common_barcode)/length(colnames(Seurat.Object))<0.5) {
    warning("Less than 50% barcode found")
  }
  rownames(SNPs) <- SNPs$barcode
  SNPs <- SNPs[, c('status', 'assignment')]
  SNPs$SNP_cluster <- SNPs$assignment
  SNPs$SNP_status <- SNPs$status
  SNPs$assignment <- NULL
  SNPs$status <- NULL
  if (verbose){
    head(SNPs)  
  }
  Seurat.Object <- SeuratObject::AddMetaData(Seurat.Object, SNPs)
  return(Seurat.Object)
}

#Generic form
Add.ADT <- function(Seurat.Object, ADT.folder.path, replace.any=FALSE, verbose=FALSE){
  #Add ADT table to your Seurat object
  #Example use: Seurat.Object =  Add.ADT(Seurat.Object, ADT.folder.path="your/path/umi_count/",replace.any=c("HLA"="MHCII")) 
  
  #Read in the ADT library
  ADT.data <- Seurat::Read10X(ADT.folder.path, gene.column=1)
  #Add "-1"
  colnames(ADT.data) <- paste0(colnames(ADT.data),"-1")
  
  #Just select the common barcodes
  common.barcodes <- intersect(colnames(Seurat.Object), colnames(ADT.data))
  print(paste0("Number of cells in Seurat with ADT: ", length(common.barcodes)))
  
  if (length(common.barcodes)==0) {
    stop("Stop! No cells have ADT")
  }

  no.ADT <- setdiff(colnames(Seurat.Object), colnames(ADT.data))
  extra.ADT <- setdiff(colnames(ADT.data), colnames(Seurat.Object))
  
  if (length(extra.ADT)!=0) {
    ADT.data <- ADT.data[, !colnames(ADT.data) %in% extra.ADT]
  }
  
  if (verbose==TRUE) {
    print(paste0("Number of cells in Seurat: ",ncol(Seurat.Object)))
    print(paste0("Number of cells in ADT: ",ncol(ADT.data)))
    print(head(ADT.data))
    print(paste0("Number of cells in without any ADT: ",length(no.ADT)))
  }
  
  if (length(no.ADT)!=0) {
  
    #create matrix with 4 columns
    no.ADT.table <- matrix(rep(0, times=nrow(ADT.data)*length(no.ADT)), ncol=length(no.ADT), byrow=TRUE)
    
    if (verbose==TRUE) {
      print(head(no.ADT.table))
    }
    
    #define column names and row names of matrix
    colnames(no.ADT.table) <- no.ADT
    rownames(no.ADT.table) <- rownames(ADT.data)
  
    if (verbose==TRUE) {
      print(ncol(no.ADT.table))
      print(head(no.ADT.table))
      print(typeof(no.ADT.table))
    }
    
    ADT.final <- cbind(ADT.data, no.ADT.table)
    
    if (verbose==TRUE) {
      print(ncol(ADT.final))
    }
  } else {
    ADT.final=ADT.data
  }
  
  ## ADT preparation
  #For ADT remove unmapped
  ADT.data.filtered = ADT.final[c(1:(length(rownames(ADT.final))-1)) ,]
  #Remove barcode from the name (they are separated with an "-")
  rownames(ADT.data.filtered) = sapply(strsplit(rownames(ADT.data.filtered),"-"), '[', 1)
  
  #Replace a ADT name with another
  if (replace.any!=FALSE) {
    print("Replacing ADT names...")
    for (i in 1:length(replace.any))
    if (any(rownames(ADT.data.filtered)==names(replace.any[i]))) {
      pos = grep(names(replace.any[i]), rownames(ADT.data.filtered))
      rownames(ADT.data.filtered)[pos] <- replace.any[[i]]
    }
  }
  
  print(barplot(rowSums(ADT.data.filtered), main = "ADT total reads", las=2))
  print(head(ADT.data.filtered))
  
  #Add the ADT as an independent assay, before to subset
  Seurat.Object[["ADT"]] <- SeuratObject::CreateAssayObject(counts = ADT.data.filtered)
  # Normalize ADT data, here we use centered log-ratio (CLR) transformation
  Seurat.Object <- Seurat::NormalizeData(Seurat.Object, assay = "ADT", normalization.method = "CLR", margin = 2)
  
  return(Seurat.Object)
}

extract.HTO <- function(path, barcodeWhitelist, minCountPerCell = 5, methods = c("bff_cluster", "multiseq","dropletutils"), datatypeName = NULL) {
	#Use:
	#> my.HTO.table <- extract.HTO("/mypath/", c("HTO1","HTO6")))
  
  #once Scott update the package version on todata3 I can uncomment this
  #if (utils::packageVersion("cellhashR")<"1.0.4") {
  #  stop("cellhashR version = or > 1.0.4 needed")
  #}
  print("Process Count Matrix...")
	barcodeData <- cellhashR::ProcessCountMatrix(rawCountData = path, minCountPerCell = minCountPerCell, barcodeWhitelist = barcodeWhitelist, datatypeName = datatypeName)
	print("Generate Cell Hashing Calls...")
	calls.HTO <- cellhashR::GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = methods)

	HTOtable <- data.frame(row.names = calls.HTO$cellbarcode)
	HTOtable$HTO <- calls.HTO$consensuscall
	HTOtable$HTO_status <- calls.HTO$consensuscall.global
	rownames(HTOtable) <- paste0(rownames(HTOtable),"-1")

	# Inspect negative cells:
	print("Summarize Cells By Classification...")
	cellhashR::SummarizeCellsByClassification(calls = calls.HTO, barcodeMatrix = barcodeData)
	return(HTOtable)
}

Add.HTO <- function(Seurat.Object, path, barcodeWhitelist, minCountPerCell = 5, methods = c("bff_cluster", "multiseq","dropletutils"), datatypeName = NULL) {
	#Add HTO table to your Seurat object
	#Use:
	#> Seurat.Object <- (Seurat.Object, "/mypath/HTO_folder/", c("HTO1","HTO6")))
	HTOtable <- extract.HTO(path, barcodeWhitelist, minCountPerCell = minCountPerCell, methods = methods, datatypeName = datatypeName)
	Seurat.Object <- SeuratObject::AddMetaData(Seurat.Object, HTOtable)
	return(Seurat.Object)
}

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

Read.BD.Rhap <- function(path, project.name="BD Rhapsody", remove.multiplets = TRUE, remove.undetermined=TRUE,  min.cells = 3, min.features = 20, show.plots=TRUE, verbose=FALSE) {
  # path = folder containing DBEC_MolsPerCell and SampleTag files from 7B
  # TO DO: implement RSEC
  # See https://scomix.bd.com/hc/en-us/articles/360044971032-Bioinformatics
  
  fileList <- list.files(path=path)
  
  #Keep only files of interest (i.e. those which contain UMInormalized counts and are sample tagged)
  print("Reading files...")
  DBEC.file <- fileList[grepl("DBEC_MolsPerCell", fileList)]
  #RSEC.file <- fileList[grepl("RSEC_MolsPerCell", fileList)]  
  sampleTag.file<- fileList[grepl("Sample_Tag", fileList)]
  
  count <- read.csv(file = paste(path, DBEC.file, sep="/"), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, comment.char = "#")
  
  #traspose gene/cells
  count <- t(count)
  if (verbose) {
    print("Counts rownames:")
    print(head(rownames(count)))
    print("Counts colnames:")
    print(head(colnames(count)))
  }
  
  #Read sampleTag
  ST <- read.csv(file = paste(path, sampleTag.file, sep="/"), sep = ',', header = TRUE, row.names = 1, check.names = FALSE, comment.char = "#")
  if (verbose) {
    print(sampleTag.file)
    print("SampleTag:")
    print(head(ST))
  }
  
  #Create Seurat
  Seurat.object <- CreateSeuratObject(counts = count, min.cells = min.cells, min.features = min.features, project = project.name)
  Seurat.object = AddMetaData(Seurat.object, ST)
  Seurat.object@meta.data$file <- sampleTag.file
  print(Seurat.object)
  
  Idents(Seurat.object) <- "Sample_Name"
  if (show.plots==TRUE){
    print(qplot(Seurat.object$Sample_Name, main=project.name))
  }
  
  if (remove.multiplets==TRUE){
    print("Removing Multiplets...")
    Seurat.object = subset(Seurat.object, idents = "Multiplet", invert = TRUE)
  }
  if (remove.multiplets==TRUE){
    print("Removing Undetermined")
    Seurat.object = subset(Seurat.object, idents ="Undetermined", invert = TRUE)
  }
  if (show.plots==TRUE){
    hist(colSums(Seurat.object@assays$RNA@counts), main=project.name, breaks = 50)
  }
  return(Seurat.object)
}

Read.BD.Rhap.simple <- function(MolsPerCell.file, project.name="BD Rhapsody", min.cells = 3, min.features = 20, show.plots=TRUE, verbose=FALSE) {
  # MolsPerCell.file = file containing MolsPerCell and SampleTag files from 7B
  
  #Keep only files of interest (i.e. those which contain UMInormalized counts and are sample tagged)
  print("Reading file...")
  
  count <- read.csv(file = MolsPerCell.file, sep = ',', header = TRUE, row.names = 1, check.names = FALSE, comment.char = "#")
  
  #traspose gene/cells
  count <- t(count)
  if (verbose) {
    print("Counts rownames:")
    print(head(rownames(count)))
    print("Counts colnames:")
    print(head(colnames(count)))
  }
  
  #Create Seurat
  Seurat.object <- CreateSeuratObject(counts = count, min.cells = min.cells, min.features = min.features, project = project.name)
  print(Seurat.object)

  if (show.plots==TRUE){
    hist(colSums(Seurat.object@assays$RNA@counts), main=project.name, breaks = 50)
  }
  return(Seurat.object)
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

DoubletMark <- function(Seurat.Object, dim, sct = FALSE, num.cores = 5, prop.doublets = NULL, n.cell.recovered=NULL, verbose = FALSE) {
  # v. 0.5
  # can provide either prop.doublets expected or n.cell.recovered are mandatory
  # improved % doublet version calculation
  # bug dim fixed
  
  #if (is.null(dim)) {
  #  dim = find.significant.PCs(Seurat.Object)
  #}
  if (is.null(prop.doublets)) {
    if (is.null(n.cell.recovered)) {
      stop("prop.doublets or n.cell.recovered need to be provided")
    }
    #Theoretical prop.doublets calculation
    #from https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/76
    #Improved version (thanks Zuzanna for finding the issue), working better at high cell number 
    prop.doublets = ((0.00076*n.cell.recovered)+0.0527)/100
    
    print(paste("Proportion of cell doublets is (in %):",prop.doublets*100, sep = " "))
  }
  
  if (prop.doublets>0.5) {
    stop("prop.doublets is too high")
  }
  
  print("paramSweep calculation...")
  sweep.res <- DoubletFinder::paramSweep_v3(Seurat.Object, PCs = 1:dim, sct = sct, num.cores = num.cores)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT = FALSE)
  if (verbose){
    print(sweep.stats)
  }
  
  bcmvn <- DoubletFinder::find.pK(sweep.stats)
  if (verbose){
    print(bcmvn)
  }
  
  ##From here:
  ##https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/62
  pK=as.numeric(as.character(bcmvn$pK))
  BCmetric=bcmvn$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  
  #plot pk:
  #par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
  #plot(x = pK, y = BCmetric, pch = 16,type="b", col = "blue",lty=1)
  #abline(v=pK_choose,lwd=2,col='red',lty=2)
  #title("The BCmvn distributions")
  #text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- Seurat.Object@meta.data$seurat_clusters
  homotypic.prop <- DoubletFinder::modelHomotypic(annotations)          
  
  nExp_poi <- round(prop.doublets*nrow(Seurat.Object@meta.data))  ## Assuming prop.doublets doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  if (verbose) {
    print(names(Seurat.Object@meta.data))
  }
  Seurat.Object <- doubletFinder_v3(Seurat.Object, PCs = 1:dim, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = sct)
  to_delete = names(Seurat.Object@meta.data)[length(colnames(Seurat.Object@meta.data))-1]
  
  #DF.cl.name = unique(colnames(Seurat.Object@meta.data))[length(colnames(Seurat.Object@meta.data))]
  names(Seurat.Object@meta.data)[length(colnames(Seurat.Object@meta.data))] = "Doublets Low stringency"
  
  if (verbose){
    print(names(Seurat.Object@meta.data))
    print("####")
    print(to_delete)
  }
  
  #Seurat.Object@meta.data$to_delete <- NULL
  Seurat.Object <- doubletFinder_v3(Seurat.Object, PCs = 1:dim, pK = pK_choose, nExp = nExp_poi.adj, reuse.pANN = "Doublets Low stringency", sct = sct)
  #DF.cl.name.2 = unique(colnames(Seurat.Object@meta.data))[length(colnames(Seurat.Object@meta.data))]
  
  #to_delete = names(Seurat.Object@meta.data)[length(colnames(Seurat.Object@meta.data))-1]
  names(Seurat.Object@meta.data)[length(colnames(Seurat.Object@meta.data))] = "Doublets High stringency"
  #print(to_delete)
  Seurat.Object@meta.data[to_delete] <- NULL
  
  if (verbose){
    print(names(Seurat.Object@meta.data))
    print(DimPlot(Seurat.Object, group.by = "Doublets Low stringency"))
    print(DimPlot(Seurat.Object, group.by = "Doublets High stringency"))
  }
  
  return(Seurat.Object)
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

#Generic form
#Unique.Name.Of.Function <- function(Seurat.Object, var1, varN){
##Describe your function
##Example use: Unique.Name.Of.Function(Seurat.Object, "var1", "var2") 

##Describe the code
#Code
#}







