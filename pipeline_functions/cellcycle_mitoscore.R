cellcycle_mito <- function(obj, saveDir,ngenes = 4000, process, Assay, samplename){
  dir.create(paste(saveDir,"QC_Vln/",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"QC_Vln/",samplename,".pdf",sep = ""),width = 12, height = 9)
  print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3, group.by = "orig.ident"))
  dev.off()
  
  message(paste("\n",samplename,"Normalization and adding the cellcycle and Mito score\n"))
  DefaultAssay(obj) <- Assay
  obj <- NormalizeData(obj)
  obj <- CellCycleScoring(obj, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes, set.ident = TRUE)
  
  ## Determing whether the cell cycle is the  major source of artifact or have any effect on the UMAP
  obj <- FindVariableFeatures(obj,selection.method = "vst",nfeatures = ngenes,verbose = FALSE)
  
  message("\nRunning PCA\n")
  # Scale the counts
  all.genes=rownames(obj)
  obj <- ScaleData(obj, features=all.genes)
  obj <- RunPCA(obj, npcs=50)
  
  ### Plotting the PCA to see if cell cycle has any effect on that
  # Plot the PCA colored by cell cycle phase
  dir.create(paste(saveDir,"cell_cycle_phase/",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"cell_cycle_phase/",samplename,"_",process,"_",Assay,"_splitted.pdf",sep = ""),width = 8, height = 8)
  print(DimPlot(obj, reduction = "pca", group.by= "Phase", split.by = "Phase"))
  dev.off()
  
  ### not spliting
  pdf(paste(saveDir,"cell_cycle_phase/",samplename,"_",process,"_",Assay,".pdf",sep = ""),width = 8, height = 8)
  print(DimPlot(obj,reduction = "pca",group.by= "Phase"))
  dev.off()
  
  message("Adding the mitochondrial effect on the dataset\n")
  obj$mitoRatio <- obj@meta.data$percent.mt / 100
  summary_mito <- summary(obj@meta.data$mitoRatio) ## based on the summary need to choose which quartile it is present
  summary_df <- as.data.frame(as.matrix(summary_mito[2:5]))
  mid_values <- summary_df[c(1,2,4),]
  
  # Turn mitoRatio into categorical factor vector based on quartile values
  obj@meta.data$mitoFr <- cut(obj@meta.data$mitoRatio,
                              breaks=c(-Inf, mid_values, Inf),
                              labels=c("Low","Medium","Medium high", "High"))
  
  dir.create(paste(saveDir,"mitoeffect/",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"mitoeffect/",samplename,"_",process,"_",Assay,"_splitted.pdf",sep = ""),width = 8, height = 6)
  print(DimPlot(obj, reduction = "pca", group.by= "mitoFr", split.by = "mitoFr"))
  dev.off()
  
  pdf(paste(saveDir,"mitoeffect/",samplename,"_",process,"_",Assay,".pdf",sep = ""),width = 8, height = 6)
  print(DimPlot(obj, reduction = "pca", group.by= "mitoFr"))
  dev.off()
  
  dir.create(paste(saveDir,"saveRDS_obj/",sep = ""), showWarnings = FALSE)
  message(paste("Saving the Object at this path",saveDir,"saveRDS_obj/"))
  #saveRDS(obj,paste(saveDir,"saveRDS_obj/",samplename,"_",process,"_",Assay,".RDS",sep = ""))
  message(paste("\n Please check for the combined sample QC at this location: ",saveDir,"QC_Vln/\n",sep = ""))
  message(paste("\n Please check for the Unwanted cellcycle Effect on your data at this location: ",saveDir,"cell_cycle_phase/\nand for mitochondria at this location: ",saveDir,"mitoeffect/\n",sep = ""))
  message("We have added the cellcycle and Mito Effect to the merged object\n Further needs to normalize and integrate it.\nWe will use the sctransform V2 for RNA, while CLR for ADT")
  return(obj)
}

# savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/usftp21.novogene.com/raw_data/combined_raw_data/Analysis/preprocessing_with_CMO_library/multiconfigHTO_in_ref/outs/per_sample_outs/regressing_scTransform_doublet/ADT_clustering/Pipeline/"
# obj <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "_RNA_singlet.RDS",full.names = TRUE)
# samplename <- gsub("_RNA_singlet.RDS","",basename(obj))
# 
# combined <- merge(readRDS(obj[1]), y = c(readRDS(obj[2]), readRDS(obj[3]),readRDS(obj[4])), 
#                   add.cell.ids = samplename, project = "combined")
# 
# cellcycle_mito(combined, savedir, ngenes = 4000, process = "merged", Assay = "RNA", samplename = "combined")
