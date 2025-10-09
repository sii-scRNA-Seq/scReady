#!/usr/bin/env Rscript

################################################################################
#
#
# 			FUNCTIONS DEFITION SECTION
#
#
# This section contains the main functions for the scReady pipeline.
################################################################################

################################################################################
#' SoupX.clean.from.CellRanger
#'
#' Remove ambient RNA contamination from scRNA-seq data using Cell Ranger output folders
#'
#' @param cellranger.folder Path to Cell Ranger output directory containing raw and filtered matrices
#' @return SoupX object with adjusted counts after ambient RNA removal
#' @details
#' 1. Loads data from Cell Ranger output using SoupX::load10X
#' 2. Estimates contamination profile with SoupX::autoEstCont
#' 3. Adjusts counts to remove ambient RNA using SoupX::adjustCounts
#' @examples
#' clean.data <- SoupX.clean.from.CellRanger(cellranger.folder = "/my/cellranger/folder/")
################################################################################
SoupX.clean.from.CellRanger <- function(cellranger.folder) {
  # Load data from Cell Ranger output directory
  sc = SoupX::load10X(cellranger.folder)
  
  # Estimate ambient RNA contamination profile
  sc = SoupX::autoEstCont(sc)

  # Print top 20 genes with highest background expression for QC
  print("Genes with highest expression in background:")
  print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))

  # Adjust counts to remove ambient RNA contamination
  out = SoupX::adjustCounts(sc)
  return(out)
}

################################################################################
#' SoupX.from.h5
#'
#' Remove ambient RNA contamination from scRNA-seq data using .h5 files
#'
#' @param cellranger.folder Path to directory containing raw and filtered .h5 matrices
#' @return SoupX object with adjusted counts after ambient RNA removal
#' @details
#' 1. Validates presence of required .h5 files
#' 2. Loads raw and filtered matrices using Seurat::Read10X_h5
#' 3. Creates SoupChannel object and performs clustering for contamination estimation
#' 4. Estimates contamination and adjusts counts
#' @note Requires both raw_feature_bc_matrix.h5 and filtered_feature_bc_matrix.h5
#' @examples
#' clean.data <- SoupX.from.h5("/folder/containing/h5/files")
################################################################################
SoupX.from.h5 <- function(cellranger.folder) {
  # Validate input files exist
  raw_file <- file.path(cellranger.folder, "raw_feature_bc_matrix.h5")
  filt_file <- file.path(cellranger.folder, "filtered_feature_bc_matrix.h5")
  
  if (!file.exists(raw_file)) {
    stop("raw_feature_bc_matrix.h5 not found in specified directory")
  }
  if (!file.exists(filt_file)) {
    stop("filtered_feature_bc_matrix.h5 not found in specified directory")
  }

  # Load matrices
  raw.matrix <- Seurat::Read10X_h5(raw_file)
  filt.matrix <- Seurat::Read10X_h5(filt_file)

  # Create SoupChannel object
  sc  <- SoupX::SoupChannel(raw.matrix, filt.matrix)

  # Create Seurat object for clustering
  Seurat.object <- CreateSeuratObject(filt.matrix)
  Seurat.object <- NormalizeData(Seurat.object)
  Seurat.object <- FindVariableFeatures(Seurat.object, selection.method = "vst", nfeatures = 750)
  Seurat.object <- ScaleData(Seurat.object)
  Seurat.object <- RunPCA(object = Seurat.object, npcs = 15)  
  Seurat.object <- RunUMAP(object = Seurat.object, dims = 1:15)
  Seurat.object <- FindNeighbors(Seurat.object, dims = 1:15) 
  Seurat.object <- FindClusters(Seurat.object, resolution = 0.1)

  # Extract metadata and UMAP coordinates for SoupX
  meta    <- Seurat.object@meta.data
  umap    <- Seurat.object@reductions$umap@cell.embeddings
  
  # Set clusters and dimensionality reduction in SoupX object
  sc  <- SoupX::setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc  <- SoupX::setDR(sc, umap)

  #With defined clusters, run the main SoupX function, calculating ambient RNA profile.
  sc  <- SoupX::autoEstCont(sc)

  # Print top 20 background genes for QC
  print("Genes with highest expression in background:")
  print(head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20))

  # Adjust counts to remove ambient RNA
  out = SoupX::adjustCounts(sc)
  return(out)
}

################################################################################
#' make.add.meta
#'
#' Add metadata to a Seurat object from an external metadata table
#'
#' @param Seurat.Object A Seurat object to which metadata will be added
#' @param metadata A data frame containing metadata to add to the Seurat object
#' @param return.only.table Logical, if TRUE returns only the metadata table (default: TRUE)
#' @param verbose Logical, if TRUE prints debugging information (default: FALSE)
#'
#' @return A Seurat object with added metadata or just the metadata table if return.only.table=TRUE
#'
#' @details
#' This function handles three metadata scenarios:
#' 1. Single row metadata (applied to all cells)
#' 2. Single column metadata (treated as a list matching Seurat object identities)
#' 3. Multi-row/column metadata (matched to Seurat object identities)
#'
#' @examples
#' # Add metadata to Seurat object
#' Seurat.object <- make.add.meta(Seurat.Object, metadata)
#'
#' # Return only the metadata table
#' metadata.table <- make.add.meta(Seurat.Object, metadata, return.only.table=TRUE)
################################################################################
make.add.meta <- function(Seurat.Object, metadata.tab, return.only.table=FALSE, verbose=FALSE) {
  # Display metadata.tab preview if verbose mode is enabled
  if (verbose) {
    print(head(metadata.tab))
  }
  
  if (is.null(metadata.tab) || nrow(metadata.tab) == 0) {
	stop("Metadata is empty or NULL")
  }

  # Case 1: Single row metadata.tab - apply to all cells in the Seurat object
  if (nrow(metadata.tab)==1) {
    if (verbose) {
      print("nrow(metadata.tab)==1")
    }
	
	# Create empty dataframe with cell barcodes as row names
    df.cells <- data.frame(row.names = colnames(Seurat.Object))
	
	# Add each metadata.tab column to the dataframe
    for (name in colnames(metadata.tab)) {
      #print(name)
      df.cells[name]=metadata.tab[name]
    }
	
  # Case 2: Single column metadata.tab - treat as a list matching Seurat object identities
  } else if (ncol(metadata.tab)==1) {
    if (verbose){
      print(colnames(metadata.tab))
      print(names(metadata.tab))
    }

    # Validate that Seurat object identities match metadata.tab rows
    if (length(setdiff(Idents(Seurat.Object), rownames(metadata.tab)))!=0 || length(setdiff(rownames(metadata.tab),Idents(Seurat.Object)))!=0) {
      stop("Seurat object Idents and metadata.tab rows are not matching.")
    }

	# Initialize empty dataframe
    df.cells <- data.frame()

	# For each identity in the Seurat object
    for (name in unique(Idents(Seurat.Object))) {
      print(name)

	  # Create dataframe for cells with current identity
      new_df <- data.frame(row.names = WhichCells(Seurat.Object, idents = name))

	  # Get corresponding metadata.tab row
      meta_row=metadata.tab[rownames(metadata.tab) == name,]
	  
	  # Add metadata.tab to the dataframe
      new_df[colnames(metadata.tab)]=meta_row
	  
	  # Combine with main dataframe
      df.cells=dplyr::bind_rows(df.cells, new_df)
      }
  } 
  # Case 3: Multi-row/column metadata.tab - match to Seurat object identities
  else {
    if (verbose) {
      message("Multi-row/column metadata.tab detected")
    }
    
	# Validate that Seurat object identities match metadata.tab rows
    if (length(setdiff(Idents(Seurat.Object), rownames(metadata.tab)))!=0 || length(setdiff(rownames(metadata.tab),Idents(Seurat.Object)))!=0) {
      stop("Seurat object Idents and metadata.tab rows are not matching.")
    }

	# Initialize empty dataframe
    df.cells <- data.frame()

	# For each identity in the Seurat object
    for (name in unique(Idents(Seurat.Object))) {
      print(name)

	  # Create dataframe for cells with current identity
      new_df <- data.frame(row.names = WhichCells(Seurat.Object, idents = name))

	  # Get corresponding metadata.tab row
      meta_row=metadata.tab[rownames(metadata.tab) == name,]
      if (verbose) {
        print("meta_row:")
        print(meta_row)
        print(colnames(meta_row))
        print(ncol(meta_row))
      }
	  
	  # Check for potential issues with metadata.tab
      if (is.null(colnames(meta_row))) {
        stop("Contact Dom")
      }
      if (ncol(meta_row)==1) {
        stop("There is only 1 column...contact dom, need to be fixed")
      }
	  
	  # Add each metadata.tab column to the dataframe
      for (col_name in colnames(meta_row)){
        #add the specific metadata.tab you need
        if (verbose) {
          print("col_name is:")
          print(col_name)
        }
        new_df[col_name]=meta_row[col_name]
      }
	  
	  # Combine with main dataframe
      df.cells=dplyr::bind_rows(df.cells, new_df)
    }
  }

  # Return either just the metadata.tab table or the Seurat object with added metadata.tab
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

################################################################################
#' find.significant.PCs
#'
#' Determine the minimum number of significant principal components (PCs) to use
#' based on variance explained and standard deviation thresholds
#'
#' @param Seurat.Object A Seurat object with PCA already computed (must have "pca" slot)
#' @param variance Numeric value (0-1) representing the minimum cumulative variance to explain (default: 0.9)
#' @param st.dev Numeric value (0-1) representing the maximum standard deviation threshold (default: 0.05)
#'
#' @return Integer representing the minimum number of significant PCs
#'
#' @details
#' This function uses two criteria to determine significant PCs:
#' 1. Cumulative variance explained exceeds the specified threshold
#' 2. Individual PC standard deviation falls below the specified threshold
#' The function returns the minimum PC number satisfying either criterion
#'
#' @examples
#' # Find significant PCs explaining 90% of variance with <5% standard deviation
#' min.pc <- find.significant.PCs(Seurat.Object, variance=0.9, st.dev=0.05)
################################################################################
find.significant.PCs <- function(Seurat.Object, variance=0.9, st.dev=0.05) {
  # Validate variance input
  if (variance>=1) {
    simpleError("variance over 100%")
  }
  
  # Convert percentages to decimal fractions
  variance=variance*100
  st.dev=st.dev*100

  # Extract standard deviations from PCA
  stdv <- Seurat.Object[["pca"]]@stdev
  sum.stdv <- sum(Seurat.Object[["pca"]]@stdev)
  
  # Calculate percentage of variance explained by each PC
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  
  # Find first PC where cumulative variance exceeds threshold AND individual SD is below threshold
  co1 <- which(cumulative > variance & percent.stdv < st.dev)[1]
  
  # Find first PC where drop in variance explained is significant (elbow method)
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] -
                       percent.stdv[2:length(percent.stdv)]) > 1-(variance/100)),
              decreasing = T)[1] + 1
			  
  # Return the minimum significant PC from both methods
  min.pc <- min(co1, co2, na.rm = TRUE)
  message("Selected ", min.pc, " significant principal components")
  return(min.pc)
}

################################################################################
#' Calc.Perc.Features
#'
#' Calculate percentage metrics for mitochondrial, hemoglobin, ribosomal genes,
#' and MALAT1 in a Seurat object
#'
#' @param Seurat.object A Seurat object containing gene expression data
#' @param mt.pattern Regular expression pattern for mitochondrial genes (default: "^MT-")
#' @param hb.pattern Regular expression pattern for hemoglobin genes (default: "^HB[^(P)]")
#' @param ribo.pattern Regular expression pattern for ribosomal genes (default: "RPS|RPL")
#' @param MALAT1.name Name of MALAT1 gene (default: "MALAT1")
#' @param plot.name Optional name for plots (not currently used)
#'
#' @return Seurat object with added percentage metrics in metadata
#'
#' @details
#' This function calculates the percentage of counts mapping to:
#' - Mitochondrial genes (percent.mt)
#' - Hemoglobin genes (percent.hb)
#' - Ribosomal genes (percent.ribo)
#' - MALAT1 (percent.MALAT1)
#' Note: Currently only works with gene names, not gene IDs
#'
#' @examples
#' # Calculate quality control metrics
#' Seurat.object <- Calc.Perc.Features(Seurat.object)
################################################################################
Calc.Perc.Features <- function(Seurat.object, mt.pattern = "^MT-", hb.pattern = "^HB[^(P)]", ribo.pattern = ("RPS|RPL"), MALAT1.name="MALAT1", plot.name="") {
  #Known issue: works only with gene names, not gene IDs
  
  # Store number of cells before QC
  Seurat.object@misc$cell.recovered = ncol(Seurat.object)
  
  # Calculate percentage metrics for different gene sets
  Seurat.object[["percent.mt"]] <- PercentageFeatureSet(Seurat.object, pattern = mt.pattern)
  Seurat.object[["percent.hb"]] <- PercentageFeatureSet(Seurat.object, pattern = hb.pattern)
  Seurat.object[["percent.ribo"]] <- PercentageFeatureSet(Seurat.object, pattern = ribo.pattern)
  Seurat.object[["percent.MALAT1"]] <- PercentageFeatureSet(Seurat.object, pattern = MALAT1.name)
  return(Seurat.object)
}

################################################################################
#' QC.n.mad
#'
#' Perform quality control filtering using median absolute deviation (MAD)
#' on mitochondrial content, gene counts, and UMI counts
#'
#' @param Seurat.object A Seurat object with calculated QC metrics
#' @param n.mad Number of median absolute deviations to use for outlier detection (default: 4)
#'
#' @return Filtered Seurat object with low-quality cells removed
#'
#' @details
#' This function filters cells based on:
#' 1. Mitochondrial content (percent.mt)
#' 2. Number of detected genes (nFeature_RNA)
#' 3. Number of transcripts (nCount_RNA)
#' Uses MAD-based thresholds to identify and remove outlier cells.
#' Based on code from: https://matthieuxmoreau.github.io/EarlyPallialNeurogenesis/
#'
#' @examples
#' # Perform QC filtering with 4 MAD threshold
#' Seurat.object <- QC.n.mad(Seurat.object, n.mad = 4)
################################################################################
QC.n.mad <- function(Seurat.object, n.mad=4) {
  Cell.QC.Stat <- Seurat.object@meta.data
  message("Starting with ", nrow(Cell.QC.Stat), " cells")

  # Set MAD-based thresholds for mitochondrial content
  max.mito.thr <- median(Cell.QC.Stat$percent.mt) + n.mad*mad(Cell.QC.Stat$percent.mt)
  min.mito.thr <- median(Cell.QC.Stat$percent.mt) - n.mad*mad(Cell.QC.Stat$percent.mt)

  # Plot mitochondrial content vs. feature count
  p1 <- ggplot(Cell.QC.Stat, aes(x=nFeature_RNA, y=percent.mt)) +
    geom_point() +
    geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
    geom_hline(aes(yintercept = min.mito.thr), colour = "red", linetype = 2) +
    annotate(geom = "text",
             label = paste0(as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr | Cell.QC.Stat$percent.mt < min.mito.thr)[2]),
                           " cells removed\n",
                           as.numeric(table(Cell.QC.Stat$percent.mt > max.mito.thr | Cell.QC.Stat$percent.mt < min.mito.thr)[1]),
                           " cells remain"),
             x = 6000, y = 0.1)

  # Filter cells based on mitochondrial content
  Cell.QC.Stat <- Cell.QC.Stat %>%
    dplyr::filter(percent.mt < max.mito.thr) %>%
    dplyr::filter(percent.mt > min.mito.thr)

  message("After mitochondrial filtering: ", nrow(Cell.QC.Stat), " cells remain")

  # Set thresholds for gene counts and UMI counts using the filtered data
  log10.nFeature <- log10(Cell.QC.Stat$nFeature_RNA)
  min.Genes.thr <- median(log10.nFeature) - n.mad*mad(log10.nFeature)
  max.Genes.thr <- median(log10.nFeature) + n.mad*mad(log10.nFeature)

  log10.nCount <- log10(Cell.QC.Stat$nCount_RNA)
  max.nUMI.thr <- median(log10.nCount) + n.mad*mad(log10.nCount)

  # Plot gene counts vs. UMI counts
  p2 <- ggplot(Cell.QC.Stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_hline(aes(yintercept = min.Genes.thr), colour = "green", linetype = 2) +
    geom_hline(aes(yintercept = max.Genes.thr), colour = "green", linetype = 2) +
    geom_vline(aes(xintercept = max.nUMI.thr), colour = "red", linetype = 2)

  # Combine plots using patchwork::wrap_plots
  combined_plot <- patchwork::wrap_plots(
    plotlist = list(p1, p2),
    ncol = 1,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = names(Seurat.object),
      theme = theme(plot.title = element_text(hjust = 0.5))
    )

  # Print the combined plot
  print(combined_plot)

  # Filter cells based on gene count and UMI thresholds
  Cell.QC.Stat <- Cell.QC.Stat %>%
    dplyr::filter(log10(nFeature_RNA) > min.Genes.thr) %>%
    dplyr::filter(log10(nFeature_RNA) < max.Genes.thr) %>%
    dplyr::filter(log10(nCount_RNA) < max.nUMI.thr)

  message("After all filtering: ", nrow(Cell.QC.Stat), " cells remain")
  message("----------------------------------------")

  # Subset Seurat object to keep only passing cells
  Seurat.object <- subset(Seurat.object, cells = rownames(Cell.QC.Stat))
  return(Seurat.object)
}

################################################################################
#' PCA.sample
#'
#' Perform PCA analysis on sample-level averaged expression data
#'
#' @param Seurat.Object A Seurat object containing single-cell data
#' @param sample.name Column name in metadata to group cells by (default: "orig.ident")
#' @param nfeatures Number of variable features to select (default: 2000)
#' @param show.label Logical, whether to show labels on the plot (default: FALSE)
#' @param legend Logical, whether to show legend on the plot (default: TRUE)
#'
#' @return A Seurat object containing sample-level averaged expression and PCA results
#'
#' @details
#' This function:
#' 1. Groups cells by the specified sample.name
#' 2. Calculates average expression per sample group
#' 3. Performs normalization, variable feature selection, scaling, and PCA
#' 4. Generates a PCA plot of the sample-level data
#'
#' The returned object can be used for further analysis or visualization.
#' You can add metadata to this object and plot factors using group.by.
#'
#' @examples
#' # Perform PCA on sample averages
#' sample.pca <- PCA.sample(Seurat.Object, sample.name = "Patient", nfeatures = 1000)
################################################################################
PCA.sample <- function(Seurat.Object, sample.name = "orig.ident", nfeatures = 2000, show.label=FALSE, legend=TRUE){

  # Set cell identities to the specified sample name
  Idents(object = Seurat.Object) <- sample.name
  
  # Calculate average expression per sample group
  sample.object <- AggregateExpression(object = Seurat.Object, return.seurat = TRUE)
  
  # Perform standard preprocessing steps
  sample.object <- NormalizeData(sample.object)
  sample.object <- FindVariableFeatures(sample.object, selection.method = "vst", nfeatures = nfeatures)
  sample.object <- ScaleData(sample.object)
  sample.object <- RunPCA(object = sample.object, npcs = 2)
  
  # Generate PCA plot with or without legend
  if (legend == TRUE) {
    print(DimPlot(sample.object, reduction = "pca", pt.size = 5, label = show.label))
  } else {
    print(DimPlot(sample.object, reduction = "pca", pt.size = 5, label = show.label)+ NoLegend())
  }

  return(sample.object)
}

################################################################################
#' Add.SNPs.HT
#'
#' Add SNP-based clustering results from Souporcell to a Seurat object
#'
#' @param Seurat.Object A Seurat object to which SNP data will be added
#' @param souporcell.file Path to Souporcell cluster.tsv output file
#' @param verbose Logical, whether to print verbose output (default: FALSE)
#' @param barcode.suffix Optional suffix to append to barcodes for matching (default: NULL)
#'
#' @return Seurat object with added SNP cluster and status metadata
#'
#' @details
#' This function:
#' 1. Reads the Souporcell cluster assignments
#' 2. Optionally appends a suffix to barcodes for matching
#' 3. Verifies barcode matching between Seurat object and Souporcell results
#' 4. Adds SNP cluster assignments and status as metadata to the Seurat object
#'
#' The function warns if less than 50% of barcodes are matched.
#'
#' @examples
#' # Add SNP cluster assignments to Seurat object
#' Seurat.Object <- Add.SNPs.HT(Seurat.Object, "Mapped/clusters.tsv")
################################################################################
Add.SNPs.HT <- function(Seurat.Object, souporcell.file, verbose=FALSE, barcode.suffix=NULL) {
  # Read Souporcell cluster assignments
  SNPs = read.csv(souporcell.file, sep = "\t")
  
  # Optionally append suffix to barcodes
  if (is.null(barcode.suffix)) {
  } else {
    SNPs$barcode = paste0(SNPs$barcode, barcode.suffix)
  }
  
  # Display barcode information for verification
  message("Seurat object barcode preview:")
  print(head(colnames(Seurat.Object)))
  message("SNPs file barcode preview:")
  print(head(SNPs$barcode))
  
  # Find common barcodes between Seurat object and SNPs file
  common_barcode = intersect(colnames(Seurat.Object), SNPs$barcode)
  message(paste("Number of matching barcodes found:", length(common_barcode)))
  
  # Warn if less than 50% of barcodes match
  if (length(common_barcode)/length(colnames(Seurat.Object))<0.5) {
    warning("Less than 50% of barcodes matched between Seurat object and SNPs file")
  }
  
  # Prepare SNP metadata for addition to Seurat object
  rownames(SNPs) <- SNPs$barcode
  SNPs <- SNPs[, c('status', 'assignment')]
  SNPs$SNP_cluster <- SNPs$assignment
  SNPs$SNP_status <- SNPs$status
  SNPs$assignment <- NULL
  SNPs$status <- NULL
  
  # Display preview if verbose
  if (verbose){
	message("SNP metadata preview:")
    print(head(SNPs))
  }
  
  # Add SNP metadata to Seurat object
  Seurat.Object <- SeuratObject::AddMetaData(Seurat.Object, SNPs)
  return(Seurat.Object)
}


################################################################################
#' extract.HTO
#'
#' Extract and process HTO (Hashtag Oligonucleotide) data using cellhashR
#'
#' @param path Path to directory containing HTO count data
#' @param barcodeWhitelist Optional vector of expected HTO barcodes (default: NULL)
#' @param minCountPerCell Minimum count threshold per cell (default: 5)
#' @param methods Vector of methods to use for HTO calling (default: c("bff_cluster", "multiseq","dropletutils"))
#' @param datatypeName Optional name for the datatype (default: NULL)
#'
#' @return Data frame containing HTO assignments and status for each cell
#'
#' @details
#' This function:
#' 1. Processes the HTO count matrix using cellhashR
#' 2. Generates HTO calls using specified methods
#' 3. Creates a table with consensus HTO assignments and status
#' 4. Returns a data frame with cell barcodes and HTO information
#'
#' Note: Requires cellhashR package (version 1.0.4 or higher recommended).
#'
#' @examples
#' # Extract HTO information
#' HTO.table <- extract.HTO("/path/to/HTO/data/", c("HTO1","HTO6"))
################################################################################
extract.HTO <- function(path, barcodeWhitelist = NULL, minCountPerCell = 5, methods = c("bff_cluster", "multiseq","dropletutils"), datatypeName = NULL) {

	if (utils::packageVersion("cellhashR")<"1.0.4") {
	  stop("cellhashR version = or > 1.0.4 needed")
	}
	
	message("Processing HTO count matrix...")
    # Process the HTO count matrix
	barcodeData <- cellhashR::ProcessCountMatrix(rawCountData = path, minCountPerCell = minCountPerCell, barcodeWhitelist = barcodeWhitelist, datatypeName = datatypeName)
	
	message("Generating HTO calls...")
    # Generate HTO calls using specified methods
	calls.HTO <- cellhashR::GenerateCellHashingCalls(barcodeMatrix = barcodeData, methods = methods)

	# Create HTO table with consensus calls
	HTOtable <- data.frame(row.names = calls.HTO$cellbarcode)
	HTOtable$HTO <- calls.HTO$consensuscall
	HTOtable$HTO_status <- calls.HTO$consensuscall.global
	
	# Append suffix to row names to match Seurat object format
	rownames(HTOtable) <- paste0(rownames(HTOtable),"-1")

	# Display cell classification summary
    message("HTO cell classification summary:")
	cellhashR::SummarizeCellsByClassification(calls = calls.HTO, barcodeMatrix = barcodeData)
	return(HTOtable)
}

################################################################################
#' Save a Seurat Object in Various Formats
#'
#' This function saves a Seurat object in one of three supported formats:
#' - Seurat (.rds)
#' - AnnData (.h5ad via SeuratDisk)
#' - SingleCellExperiment (.rds)
#'
#' @param seurat.object A Seurat object to be saved.
#' @param path Character string specifying the directory where the file will be saved.
#' @param date Character string representing the date (or any prefix) to be included in the filename.
#'             If missing, the current date in `YYYY-MM-DD` format will be used.
#' @param file.type Character string specifying the output format. Must be one of:
#'             `"Seurat"`, `"AnnData"`, or `"SingleCellExperiment"`.
#' @param extra Character string to be inserted between the file type and extension (e.g., "merged").
#'             Default is NA (no extra text).
#'
#' @return Character string with the full path to the saved file.
#'
#' @examples
#' # Save a Seurat object as an .rds file
#' WritObject(seurat.object = my_seurat_obj, path = "results/", date = "2025-10-09", file.type = "Seurat", extra = "merged")
#'
#' # Save as AnnData (.h5ad)
#' WritObject(seurat.object = my_seurat_obj, path = "results/", date = "2025-10-09", file.type = "AnnData", extra = "integrated")
#'
#' # Save as SingleCellExperiment (.rds)
#' WritObject(seurat.object = my_seurat_obj, path = "results/", date = "2025-10-09", file.type = "SingleCellExperiment", extra = "processed")
#'
#' @importFrom SeuratDisk SaveH5Seurat Convert
#' @export
################################################################################
WritObject <- function(seurat.object, path, date, file.type, extra = NA) {

    # Validate file.type
    if (!(file.type %in% c("Seurat", "AnnData", "SingleCellExperiment"))) {
        stop("Unsupported file.type. Must be one of: 'Seurat', 'AnnData', or 'SingleCellExperiment'.")
    }
    
    # Build filename components
    filename <- date
    
    # Add extra text if provided
    if (!is.na(extra)) {
        filename <- paste(filename, extra, sep = "_")
    }
    
    # Add file type and extension
    if (file.type == "Seurat") {
        filename <- paste0(filename, ".Seurat.rds")
        saveRDS(seurat.object, file.path(path, filename))
    } else if (file.type == "AnnData") {
        filename <- paste0(filename, "_AnnData.h5Seurat")
        SeuratDisk::SaveH5Seurat(seurat.object, filename = file.path(path, filename), overwrite = TRUE)
		SeuratDisk::Convert(file.path(path, filename), dest = "h5ad", overwrite = TRUE)
    } else if (file.type == "SingleCellExperiment") {
        filename <- paste0(filename, "_SCE.rds")
        sce <- as.SingleCellExperiment(seurat.object)
        saveRDS(sce, file.path(path, filename))
    }
    
    message("Saved ", file.type, " object to: ", file.path(path, filename))
    return(file.path(path, filename))
}

################################################################################
#
#
# 			SCRIPT SECTION
#
#
# This section contains the main workflow execution code for the scReady pipeline.
# It handles configuration loading, data processing, and quality control steps.
################################################################################

################################################################################
# CONFIGURATION LOADING AND VALIDATION
#
# Load and validate the configuration file containing all pipeline parameters.
# The configuration file should define all required variables for the pipeline.
################################################################################

# Define list of required configuration variables
required_vars <- c(
  "min.cells",           # Minimum number of cells required per sample
  "min.features",        # Minimum number of features required per cell
  "ambient.removal",     # Whether to perform ambient RNA removal
  "HTO.features",        # HTO features for demultiplexing
  "demultiplex.algorims",# Algorithms to use for demultiplexing
  "ADT.normalize.scale", # Whether to normalize and scale ADT data
  "souporcell_folder",   # Path to Souporcell output folder
  "save.pre_QC",         # Whether to save data before QC filtering
  "n.mad",               # Number of median absolute deviations for outlier detection
  "var.features",        # Number of variable features to select
  "max.pca",             # Maximum number of principal components to calculate
  "var.explained",       # Proportion of variance to explain with PCA
  "integration",         # Whether to perform data integration
  "integration.method",  # Method to use for data integration
  "file.type"			 # File datatype to save
)

# Attempt to load config file from current directory. If missing, load default config
config_file <- "scReady.config"
if (!file.exists(config_file)) {
  message("scReady.config not found in current directory. Loading default config instead.")
  config_file <- file.path("/opt", "app", "config", config_file)
} else {
  message(paste("Loading local", config_file, sep=" "))
}
source(config_file)

# Check for missing required variables
missing_vars <- setdiff(required_vars, ls())
if (length(missing_vars) > 0) {
  stop("Missing variables in scReady.config: ", paste(missing_vars, collapse = ", "))
}

################################################################################
# COMMAND-LINE ARGUMENT PROCESSING
#
# This section handles command-line arguments passed to the script.
# The script accepts:
# 1. A required path to the folder containing mapped data
# 2. An optional path to a CSV metadata file
################################################################################

# Get command-line arguments passed to the script
args <- commandArgs(trailingOnly = TRUE)

# Validate the number of arguments provided
if (length(args) > 2) {
  stop("Incorrect number of arguments. ",
       "Please provide one or two arguments: ",
       "(1) folder path containing mapped data, and ",
       "(2) optional CSV metadata file path.",
       call. = FALSE)
}

# Extract the required path to the data folder
path.to.read <- args[1]

# Check if an optional CSV metadata file was provided
optional_csv_file <- ifelse(length(args) == 2, args[2], NA)

# If an optional CSV file is provided, load and process it
if (!is.na(optional_csv_file)) {
  # Verify the file exists before attempting to read it
  if (!file.exists(optional_csv_file)) {
    stop("Optional metadata file not found: ", optional_csv_file)
  }

  # Read the CSV file, using the first column as row names
  optional_data <- read.csv(optional_csv_file, row.names = 1)

  # Inform the user that the optional file was successfully loaded
  message("Optional CSV metadata file loaded: ", optional_csv_file)
  message("Metadata contains ", nrow(optional_data), " rows and ",
          ncol(optional_data), " columns.")
}

################################################################################
# LOAD REQUIRED LIBRARIES
#
# Load all R packages required for the pipeline execution.
# These packages are installed in the Docker container.
################################################################################

# Load Seurat for single-cell data analysis
library(Seurat)

# Load SingleCellExperiment for handling single-cell data structures
library(SingleCellExperiment)

# Load dplyr for data manipulation
library(dplyr)
library(hdf5r)

# Load SoupX for ambient RNA removal
library(SoupX)

# Load libraries for plots
library(ggplot2) 
library(ggExtra)
library(patchwork)

# Note: Additional packages are loaded as needed in specific functions
# to minimize memory usage and startup time.

################################################################################
# DATA DISCOVERY AND IMPORT
#
# This section searches the pipeline folder for mapped data and imports it
# into the appropriate data structures for downstream processing.
# It handles different data formats:
# - Cell Ranger multi (per_sample_outs structure)
# - Cell Ranger count (filtered_feature_bc_matrix)
# - H5 files (filtered_feature_bc_matrix.h5)
# - CITE-seq data (Antibody Capture)
# - HTO data (Hashtag Oligonucleotides)
# - SNP data from Souporcell
################################################################################

# Initialize lists to store different types of data
list.data = list()            # Main gene expression data
list.data.citeseq = list()    # CITE-seq antibody data
list.snp = list()             # SNP data from Souporcell
list.hto = list()             # HTO data
list.ambient.result = list()  # Results of ambient RNA removal
to.skip = c()                 # Samples to skip processing

# Remove trailing slash from path if present
path.to.read = sub("/$", "", path.to.read) 

# Get list of sample directories in the input path
runs = list.dirs(path = path.to.read, full.names = FALSE, recursive = FALSE) 

# Initialize PDF for output
date=Sys.Date()
pdf.filename = paste0(date,"_output.pdf") 
pdf(normalizePath(file.path(path.to.read, pdf.filename)))

message("Reading 10X folders with SoupX...")

# Process each sample directory
for (name in runs) {
  # Skip samples in the to.skip list
  if (name %in% to.skip) {
  } else {
    cat(name)
    cat("\n")
	
	# Construct full path to the sample directory
	path_name = file.path(path.to.read, name)
    
	# Check if 'outs' directory exists, if not use the sample directory itself
	path_name_outs = if (dir.exists(file.path(path_name, "outs"))) {
		file.path(path_name, "outs")
	} else {
		message("'outs' directory not found in ", path_name, ". Using sample directory directly.")
		path_name
	}

	# Check for Souporcell results in the sample directory
	souporcell_path = file.path(path_name, souporcell_folder, "clusters.tsv")
	if (file.exists(souporcell_path)) {
		message("Found Souporcell SNP data")	
		list.snp[[name]] = souporcell_path
	}
	
    # Define paths for cellranger multi or count
    sample_id = basename(name)
    multi_path_filtered = file.path(path_name_outs, "per_sample_outs", sample_id, "count", "sample_filtered_feature_bc_matrix")
    
	# Process data based on the directory structure
    if (dir.exists(multi_path_filtered)) {
	  # Cell Ranger multi structure
      message("Found Cell Ranger multi structure")
      path_name_filtered = multi_path_filtered
	  data=Read10X(path_name_filtered)
	  list.data[[name]]=try(SoupX.clean.from.CellRanger(path_name_outs))
	  
	} else if (dir.exists(file.path(path_name_outs, "filtered_feature_bc_matrix")) &&
			   dir.exists(file.path(path_name_outs, "raw_feature_bc_matrix"))) {
	# Cell Ranger count structure with both filtered and raw matrices      print("Found cellranger count structure")
      path_name_filtered = file.path(path_name_outs, "filtered_feature_bc_matrix")
      path_name_raw = file.path(path_name_outs, "raw_feature_bc_matrix")
	  print("Removing ambient RNA")
      list.data[[name]]=try(SoupX.clean.from.CellRanger(path_name_outs))
	  data=Read10X(path_name_filtered)
	  
	} else if (file.exists(file.path(path_name_outs, "filtered_feature_bc_matrix.h5")) &&
			   file.exists(file.path(path_name_outs, "raw_feature_bc_matrix.h5"))) {
	  # H5 file structure with both filtered and raw matrices
	  path_name_filtered = file.path(path_name_outs,"filtered_feature_bc_matrix.h5")
	  path_name_raw = file.path(path_name_outs, "raw_feature_bc_matrix")
	  data=Read10X_h5(path_name_filtered)
	  list.data[[name]]=try(SoupX.from.h5(path_name_outs))
	} else {
		stop("Required files not found in: ", path_name_outs, "\n",
       "Expected either:\n",
       "1. filtered_feature_bc_matrix and raw_feature_bc_matrix directories, or\n",
       "2. filtered_feature_bc_matrix.h5 and raw_feature_bc_matrix.h5 files, or\n",
       "3. Cell Ranger multi structure (per_sample_outs/.../sample_filtered_feature_bc_matrix)")
	}

	# Process CITE-seq data if present
    if (typeof(data)=="list") {
      print("Found CITEseq data")
      list.data.citeseq[[name]]=data$`Antibody Capture`

      print("trying hashtag...")
	  # Extract HTO data if present
      HTO = try(extract.HTO(path_name_filtered, barcodeWhitelist = HTO.features, datatypeName = "Antibody Capture"))  
    if (inherits(HTO, "try-error")) { 
	  print("HTO failed")
      } else {
	  list.hto[[name]]=HTO 
	  print("HTO added")
	  print(names(list.hto))
	  }
    }

	# Handle cases where SoupX fails (e.g., low cell diversity)
    if (inherits(list.data[[name]], "try-error")) {
	  list.ambient.result[[name]] = "Ambient removal failed"
      if (typeof(data)=="list"){
         list.data[[name]]=data$`Gene Expression`
      } else {
         list.data[[name]]=data
      }
    } else {
	  list.ambient.result[[name]] = "Ambient RNA removed"
	}
    cat("\n")
  }
}

# Clean up and free memory
rm(runs)
gc()

# Report summary of loaded data
message("Successfully loaded data from ", length(list.data), " samples")

################################################################################
# SEURAT OBJECT CREATION AND PREPROCESSING
#
# This section creates Seurat objects from the imported data and performs:
# 1. Basic object creation
# 2. Quality control metric calculation
# 3. Doublet detection
# 4. Integration of CITE-seq, SNP, and HTO data
# 5. Merging of individual sample objects
# 6. Pre-QC visualization
################################################################################

# Initialize list to store Seurat objects
Seurat.list = list()

# Process each sample to create a Seurat object
for (i in 1:length(list.data)) {
  name = names(list.data[i])
  
  # Create Seurat object from the processed data
  Seurat.object = CreateSeuratObject(counts = list.data[[i]], project = name, min.cells = min.cells, min.features = min.features)
  
  # Calculate quality control metrics
  Seurat.object = Calc.Perc.Features(Seurat.object)

  # Doublet detection using scDblFinder
  sce <- as.SingleCellExperiment(Seurat.object)
  sce <- scDblFinder::scDblFinder(sce)

  # Extract doublet detection results and add to metadata
  results <- data.frame("Barcode" = rownames(colData(sce)), "scDblFinder_DropletType" = sce$scDblFinder.class, "scDblFinder_Score" = sce$scDblFinder.score) 
  rownames(results) <- results$Barcode
  results$Barcode <- NULL
  Seurat.object <- AddMetaData(Seurat.object, results)
  
  # Add ambient RNA removal results to metadata
  Seurat.object$ambient.results = list.ambient.result[[name]]
  
  # Clean up SingleCellExperiment object
  rm(sce)
  gc()

  # Add CITE-seq data if available
  if (name %in% names(list.data.citeseq)) {
    Seurat.object[['Protein']]=CreateAssayObject(counts = list.data.citeseq[[name]])
  }

  # Add SNP data if available
  if (name %in% names(list.snp)) {
    Seurat.object = try(Add.SNPs.HT(Seurat.object, list.snp[[name]]))
  }
  
  # Add HTO data if available
  if (name %in% names(list.hto)) {
	print("if ok, trying to add HTO")
    Seurat.object = try(SeuratObject::AddMetaData(Seurat.object, list.hto[[name]]))
	print(colnames(Seurat.object@meta.data))
  }

  # Store the Seurat object in the list
  Seurat.list[[name]] = Seurat.object
  
  # Clean up
  rm(Seurat.object)
  gc()
  
  message("Completed processing sample: ", name, "\n")
}

# Verify all samples were processed
if (length(Seurat.list) != length(list.data)) {
  stop("Number of processed samples (", length(Seurat.list),
          ") does not match input samples (", length(list.data), ")")
}

# Save Seurat objects before QC if configured
if (save.pre_QC) {
	filename=paste0(date,"_Seurat.list.preQC.rds")
	save_path = normalizePath(file.path(path.to.read, filename))
	message("Saving pre-QC Seurat objects to: ", save_path)
	saveRDS(Seurat.list, save_path)
}

# Merge individual Seurat objects into a single object
if(length(Seurat.list)>1) {
  message("Merging ", length(Seurat.list), " Seurat objects")
  merged.Seurat = merge(Seurat.list[[1]], y=Seurat.list[2:length(Seurat.list)])
} else {
  message("Using single Seurat object (no merging needed)")
  merged.Seurat = Seurat.list[[1]]
}

# Generate pre-QC visualization plots
message("Generating pre-QC visualization plots")
plot(VlnPlot(merged.Seurat, features = "nCount_RNA", split.by = "orig.ident"),
     main = "UMI Counts by Sample")
plot(VlnPlot(merged.Seurat, features = "nFeature_RNA", split.by = "orig.ident"),
     main = "Feature Counts by Sample")
plot(VlnPlot(merged.Seurat, features = "percent.mt", split.by = "orig.ident"),
     main = "Mitochondrial Content by Sample")
plot(VlnPlot(merged.Seurat, features = "percent.hb", split.by = "orig.ident"),
     main = "Hemoglobin Content by Sample")
plot(VlnPlot(merged.Seurat, features = "percent.ribo", split.by = "orig.ident"),
     main = "Ribosomal Content by Sample")
plot(VlnPlot(merged.Seurat, features = "percent.MALAT1", split.by = "orig.ident"),
     main = "MALAT1 Content by Sample")

#clean up
rm(merged.Seurat)
gc()

################################################################################
# QUALITY CONTROL, NORMALIZATION, AND DIMENSIONALITY REDUCTION
#
# This section performs the following steps for each sample:
# 1. Quality control filtering using MAD thresholds
# 2. Data normalization
# 3. Variable feature selection
# 4. Data scaling
# 5. Principal Component Analysis (PCA)
# 6. Determination of significant PCs
# 7. UMAP dimensionality reduction
# 8. Clustering
# 9. Doublet visualization
################################################################################

# Process each Seurat object for QC, normalization, and dimensionality reduction
for (i in 1:length(Seurat.list)) {
  # 1. Apply quality control filtering using MAD thresholds
  Seurat.list[[i]] = QC.n.mad(Seurat.list[[i]], n.mad = n.mad)

  # 2. Normalize data
  Seurat.list[[i]] <- NormalizeData(Seurat.list[[i]])
  
  # 3. Find variable features
  Seurat.list[[i]] <- FindVariableFeatures(Seurat.list[[i]], selection.method = "vst", nfeatures = var.features)
  
  # 4. Scale data
  Seurat.list[[i]] <- ScaleData(Seurat.list[[i]])
  
  # 5. Run PCA
  Seurat.list[[i]] <- RunPCA(object = Seurat.list[[i]], npcs = max.pca)

  # 6. Determine number of significant PCs explaining X% of variance
  min.pca=find.significant.PCs(Seurat.list[[i]], var.explained) 
  message("Using ", min.pca, " significant PCs explaining ~",
          round(var.explained * 100), "% of variance")

  # 7. Run UMAP using significant PCs  
  Seurat.list[[i]] <- RunUMAP(object = Seurat.list[[i]], dims = 1:min.pca)
  
  # 8. Find neighbors and clusters
  Seurat.list[[i]] <- FindNeighbors(Seurat.list[[i]], dims = 1:min.pca) %>% FindClusters(resolution = 0.1)

  # 9. Visualize doublets
  plot <- DimPlot(Seurat.list[[i]], group.by = "scDblFinder_DropletType")
  plot = plot + patchwork::plot_annotation(title = names(Seurat.list[[i]]), theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  print(plot)
 }

# Verify all samples were processed
if (length(Seurat.list) != length(list.data)) {
  warning("Number of processed samples (", length(Seurat.list),
          ") does not match input samples (", length(list.data), ")")
}

#clean up
rm(list.data)
gc()

# Merge individual Seurat objects into a single object
if(length(Seurat.list)>1) { 
  merged.Seurat = merge(Seurat.list[[1]], y=Seurat.list[2:length(Seurat.list)])
  # Join RNA layers to combine data from multiple samples
  merged.Seurat[["RNA"]] <- JoinLayers(merged.Seurat[["RNA"]])
} else {
  merged.Seurat = Seurat.list[[1]]
}

# Generate post-QC visualization plots
message("Generating post-QC visualization plots...")
Idents(merged.Seurat) <- "orig.ident"
plot(VlnPlot(merged.Seurat, features = "nCount_RNA", split.by = "orig.ident"),
     main = "Post-QC: UMI Counts by Sample")
plot(VlnPlot(merged.Seurat, features = "nFeature_RNA", split.by = "orig.ident"),
     main = "Post-QC: Feature Counts by Sample")
plot(VlnPlot(merged.Seurat, features = "percent.mt", split.by = "orig.ident"),
     main = "Post-QC: Mitochondrial Content by Sample")
plot(VlnPlot(merged.Seurat, features = "percent.hb", split.by = "orig.ident"),
     main = "Post-QC: Hemoglobin Content by Sample")
plot(VlnPlot(merged.Seurat, features = "percent.ribo", split.by = "orig.ident"),
     main = "Post-QC: Ribosomal Content by Sample")
plot(VlnPlot(merged.Seurat, features = "percent.MALAT1", split.by = "orig.ident"),
     main = "Post-QC: MALAT1 Content by Sample")

# Clean up memory
gc()

################################################################################
# FINAL PROCESSING, INTEGRATION, AND OUTPUT
#
# This section performs:
# 1. Optional metadata addition from CSV file
# 2. Variable feature selection on merged data
# 3. Data scaling
# 4. Principal Component Analysis (PCA)
# 5. Dimensionality reduction and clustering
# 6. Optional data integration
# 7. Saving final results
################################################################################

# Add metadata from optional CSV file if provided
if (!is.na(optional_csv_file)) {
    Idents(merged.Seurat) <- "orig.ident"

    # Get the unique sample IDs from the Seurat object
    sample_ids <- unique(merged.Seurat$orig.ident)

    # Check if the sample IDs exist in the metadata
    matching_samples <- intersect(sample_ids, rownames(optional_data))

    if (length(matching_samples) > 0) {
        # Subset to keep only matching samples and ensure we get a data frame
        metadata.tab <- optional_data[matching_samples, , drop = FALSE]

        # Check if metadata.tab is valid
        if (!is.null(metadata.tab) && nrow(metadata.tab) > 0) {
            # Add metadata to the Seurat object
            merged.Seurat <- make.add.meta(merged.Seurat, metadata.tab)
        } else {
            warning("Metadata subset is empty or invalid")
        }
    } else {
        warning("No matching samples found between Seurat object and metadata")
    }
}

# Find variable features in the merged dataset
merged.Seurat <- FindVariableFeatures(merged.Seurat, selection.method = "vst", nfeatures = var.features)

# Scale the merged data
merged.Seurat <- ScaleData(merged.Seurat)

# Run PCA on the merged dataset
merged.Seurat <- RunPCA(object = merged.Seurat, npcs = max.pca)

# Determine number of significant PCs explaining the specified variance
min.pca=find.significant.PCs(merged.Seurat, var.explained)
message("Using ", min.pca, " significant PCs explaining ~",
        round(var.explained * 100), "% of variance")

# Run UMAP using significant PCs - preintegration
merged.Seurat <- RunUMAP(object = merged.Seurat, dims = 1:min.pca)

# Find neighbors and clusters
merged.Seurat <- FindNeighbors(merged.Seurat, dims = 1:min.pca) %>% FindClusters(resolution = 0.1)

# Save the merged Seurat object
WritObject(merged.Seurat, path.to.read, date, file.type = file.type, extra = "merged")

# Perform data integration if enabled
if (integration) {
	if(integration.method == "harmony") {
		merged.Seurat <- harmony::RunHarmony(merged.Seurat, group.by.vars = "orig.ident", plot_convergence = TRUE)
		
		# Run UMAP and clustering on Harmony results
		merged.Seurat <- RunUMAP(merged.Seurat, reduction = "harmony", dims=1:min.pca)
		merged.Seurat <- FindNeighbors(merged.Seurat, reduction = "harmony", dims=1:min.pca) %>% FindClusters(resolution = 0.1)
	} else if(integration.method == "seurat") {
		#TO DO
	} else {
    warning("Unknown integration method: ", integration.method,
            ". Skipping integration.")
	}
	
	# Save the integrated object
	WritObject(merged.Seurat, path.to.read, date, file.type = file.type, extra = "integrated")

	message("Object saved")
}

# Close the PDF device
message("Closing PDF output device...")
dev.off()

message("Pipeline completed successfully!")
