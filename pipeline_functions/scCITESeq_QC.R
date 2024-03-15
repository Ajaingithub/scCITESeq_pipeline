### Arguments 
# Dir: Directory for the CellRanger output
# Sample: SampleName which has been put during the cellranger command
# saveDir: The directory where all the files will be saved
# min_cells, min_genes, max_genes, mitocpercent, ADT_UMI are the quality control for the single cell CITEseq experiment
scCITE_QC <- function(Dir,Sample,saveDir, min_cells=3, min_genes=200, max_genes=4000, mitopercent=10, ADT_UMI=20000){
  message(paste("\n######################################## Processsing",Sample,"#####################################################\n"))
  counts <- Read10X(paste(Dir,Sample,"/count/sample_feature_bc_matrix/",sep = ""))
  GEX <- counts$`Gene Expression`
  ADT <- counts$`Antibody Capture`
  HTO <- counts$`Multiplexing Capture`
  control_index <- grep("control",rownames(ADT),ignore.case = T) ## Not considering the control antibody as a variable feature
  ADT<- ADT[-control_index,]
  dir.create(saveDir,showWarnings = FALSE)
  
  # Creating the seurat Object
  message(paste("\ncreating GEX Seurat Object\n"))
  GEX_obj <- CreateSeuratObject(counts = GEX, project = Sample, min.cells = min_cells)
  GEX_obj[["percent.mt"]] <- PercentageFeatureSet(GEX_obj, pattern = "^MT-")
  dir.create(paste(saveDir,"QC_Vln",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"QC_Vln/",Sample,"_RNA_before.pdf",sep = ""),width = 8, height = 6)
  print(VlnPlot(GEX_obj, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3))
  dev.off()
  
  message(paste("\ncreating ADT Seurat Object\n"))
  ADT_obj <- CreateSeuratObject(counts = ADT, project = Sample, min.cells = min_cells, assay = "ADT")
  pdf(paste(saveDir,"QC_Vln/",Sample,"_ADT_before.pdf",sep = ""),width = 6, height = 6)
  print(VlnPlot(ADT_obj, features = c( "nFeature_ADT", "nCount_ADT"), ncol = 2))
  dev.off()
  
  #### Subset #######
  ### Subsetting the object based on the quality matrix
  ## starting with the RNA matrices
  message("\nsubsetting ADT and GEX object\n")
  GEX_obj_2 <- subset(GEX_obj, subset = nFeature_RNA > min_genes & nFeature_RNA < max_genes &  percent.mt < mitopercent)
  pdf(paste(saveDir,"QC_Vln/",Sample,"_RNA_after.pdf",sep = ""),width = 6, height = 6)
  print(VlnPlot(GEX_obj_2, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3))
  dev.off()
  a <- nrow(GEX_obj@meta.data)
  b <- nrow(GEX_obj_2@meta.data)
  c <- a-b
  print(paste("Total GEX cells in",Sample,a,"subsetted",min_genes,"to",max_genes,"genes",b,"removed cells",c))
  
  ADT_obj_2 <- subset(ADT_obj, subset = nCount_ADT < ADT_UMI)
  pdf(paste(savedir,"QC_Vln/",Sample,"_ADT_after.pdf",sep = ""),width = 6, height = 6)
  print(VlnPlot(ADT_obj_2, features = c( "nFeature_ADT", "nCount_ADT"), ncol = 2))
  dev.off()
  a <- nrow(ADT_obj@meta.data)
  b <- nrow(ADT_obj_2@meta.data)
  c <- a-b
  print(paste("Total ADT cells in",Sample,a,"subsetted <",ADT_UMI,"nUMI",b,"removed cells",c))
  
  #### Same cells in both the modality ######
  common_cells <- colnames(GEX_obj_2)[match(colnames(ADT_obj_2),colnames(GEX_obj_2), nomatch = 0)]
  GEX_obj_3 <- subset(GEX_obj_2, cells = common_cells)
  ADT_obj_3 <- subset(ADT_obj_2, cells = common_cells)
  print(paste("Final cells in ADT",nrow(GEX_obj_3@meta.data),"RNA",nrow(ADT_obj_3@meta.data)))
  condition <- all(colnames(GEX_obj_3)==colnames(ADT_obj_3))
  print(condition)
  stopifnot(condition == "TRUE")
  
  #### cell cycle and mitochondria Effect for each GEX object #######
  message("\nAdding cell phase and mito Ratio to cells \n")
  s_genes <- cc.genes$s.genes
  g2m_genes <- cc.genes$g2m.genes
  
  DefaultAssay(GEX_obj_3) <- 'RNA'
  GEX_obj_3 <- NormalizeData(GEX_obj_3, verbose = TRUE)
  GEX_obj_3 <- CellCycleScoring(GEX_obj_3, g2m.features = g2m_genes, s.features = s_genes)
  
  ## Determing whether the cell cycle is the  major source of artifact or have any effect on the UMAP
  GEX_obj_3 <- FindVariableFeatures(GEX_obj_3,selection.method = "vst",nfeatures = max_genes, verbose = TRUE)
  
  # Scale the counts
  all.genes = rownames(GEX_obj_3)
  GEX_obj_3 <- ScaleData(GEX_obj_3, features=all.genes)
  
  # Perform PCA
  GEX_obj_3 <- RunPCA(GEX_obj_3, npcs = 50, verbose = TRUE)
  
  ## Elbow plot
  dir.create(paste(saveDir,"elbow_plots",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"elbow_plots/",Sample,"_after_filter.pdf", sep = ""), width = 8, height = 7)
  print(ElbowPlot(GEX_obj_3, ndims = 50, reduction = "pca"))
  dev.off()
  
  ### Plotting the PCA to see if cell cycle has any effect on that
  # Plot the PCA colored by cell cycle phase
  dir.create(paste(saveDir,"cell_cycle_phase",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"cell_cycle_phase/",Sample,"_after_filter_splitted.pdf", sep = ""), width = 8, height = 7)
  print(DimPlot(GEX_obj_3, reduction = "pca", group.by= "Phase", split.by = "Phase"))
  dev.off()
  
  ### not spliting
  pdf(paste(saveDir,"cell_cycle_phase/",Sample,"_after_filter.pdf", sep = ""), width = 5, height = 5)
  print(DimPlot(GEX_obj_3, reduction = "pca", group.by= "Phase"))
  dev.off()
  
  ### Finding out the mitochondrial effect on the dataset 
  GEX_obj_3$mitoRatio <- GEX_obj_3@meta.data$percent.mt / 100
  summary_mito <- summary(GEX_obj_3@meta.data$mitoRatio) ## based on the summary need to choose which quartile it is present
  summary_df <- as.data.frame(as.matrix(summary_mito[2:5]))
  mid_values <- summary_df[c(1,2,4),]
  
  # Turn mitoRatio into categorical factor vector based on quartile values
  GEX_obj_3@meta.data$mitoFr <- cut(GEX_obj_3@meta.data$mitoRatio, 
                                    breaks=c(-Inf, mid_values, Inf), 
                                    labels=c("Low","Medium","Medium high", "High"))
  
  dir.create(paste(saveDir,"mitoeffect",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"mitoeffect/",Sample,"_after_filter_splitted.pdf", sep = ""), width = 8, height = 7)
  print(DimPlot(GEX_obj_3,reduction = "pca",group.by= "mitoFr",split.by = "mitoFr"))
  dev.off()
  
  pdf(paste(saveDir,"mitoeffect/",Sample,"_after_filter.pdf", sep = ""), width = 5, height = 5)
  print(DimPlot(GEX_obj_3,reduction = "pca", group.by= "mitoFr"))
  dev.off()
  
  #### saving objects #####
  dir.create(paste(saveDir,"saveRDS_obj",sep = ""), showWarnings = FALSE)
  message(paste("Saving GEX and ADT objects for",Sample))
  saveRDS(GEX_obj,paste(saveDir,"saveRDS_obj/",Sample,"_RNA_before_filtering.RDS",sep = ""))
  saveRDS(ADT_obj,paste(saveDir,"saveRDS_obj/",Sample,"_ADT_before_filtering.RDS",sep = ""))
  saveRDS(GEX_obj_3,paste(saveDir,"saveRDS_obj/",Sample,"_RNA.RDS",sep = ""))
  saveRDS(ADT_obj_3,paste(saveDir,"saveRDS_obj/",Sample,"_ADT.RDS",sep = ""))
  
  message(paste(Sample,"\nInitial Filtering is done\n"))
  message(paste("\n Please check for the filtering Violin QC at this location: ",saveDir,"QC_Vln/\n",sep = ""))
  message(paste("\n Please check for the Unwanted cellcycle Effect on your data at this location: ",saveDir,"cell_cycle_phase/\nand for mitochondria at this location: ",saveDir,"mitoeffect/\n",sep = ""))
  message(paste("\n Please check for elbow plot, to decide the PCs or dims for the nearest neigbour at this location:",saveDir,"elbow_plots/\n",sep = ""))
  message(paste("\n Also you can find how many cells has been removed during the filtering process at this location: ",saveDir,"counts/\n",sep = ""))
  message(paste("\n if there is any Effect that could lead to bias in your result, please regress it out at the sctransform\n"))
  message("\nFurther Remove the doublet Finder\n")
  
  a <- nrow(GEX_obj@meta.data)
  b <- nrow(GEX_obj_3@meta.data)
  c <- a-b
  d=c(a,b,c)
  return(d)
}














