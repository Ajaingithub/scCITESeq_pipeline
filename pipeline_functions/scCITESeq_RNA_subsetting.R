# ADT_main is an object from where we will subset the RNA and ADT object
RNA_subsetting <- function(RNA_obj_path, ADT_obj_path, ADT_main, subset_cluster=NULL, cellnames=NULL, ngenes=4000, saveDir, objname, Assay, process){
  message("\n loading both the objects \n")
  RNA_integrated_impute <- RNA_obj_path
  ADT_integrated_impute <- ADT_obj_path
  
  message("\n subsetting for the ADT ",objname," cells \n")
  DefaultAssay(ADT_integrated_impute) <- "ADT"
  CD4_ADT_main <- subset(ADT_main, idents = subset_cluster, cells = cellnames)
  CD4_cell_names <- colnames(CD4_ADT_main)
  
  # CD4_cell_names <- gsub("^DMSO_|^zostavax_|^shingrix_","",CD4_cell_names)
  
  message(paste("\n subsetting for the RNA",objname,"cells \n"))
  CD4_RNA_integrated <- subset(RNA_integrated_impute, cells = CD4_cell_names)
  print(CD4_RNA_integrated)
  
  message(paste("\n subsetting for the ADT",objname,"cells \n"))
  CD4_ADT_integrated_impute <- subset(ADT_integrated_impute, cells = CD4_cell_names)
  print(CD4_ADT_integrated_impute)
  
  DefaultAssay(CD4_RNA_integrated) <- Assay
  message("\n Normalizing, Scaling, and PCA \n")
  CD4_RNA_integrated <- NormalizeData(CD4_RNA_integrated)
  CD4_RNA_integrated <- FindVariableFeatures(CD4_RNA_integrated, selection.method = "vst", nfeatures = ngenes, verbose = FALSE)
  
  all.genes=rownames(CD4_RNA_integrated)
  CD4_RNA_integrated <- ScaleData(CD4_RNA_integrated, features=all.genes)
  CD4_RNA_integrated <- RunPCA(CD4_RNA_integrated, npcs=50, 
                               reduction.name = "pca", assay = "RNA")
  
  message("\n CD4 cell cycle phase and mitochondria Effect \n")
  DefaultAssay(CD4_RNA_integrated) <- "RNA"
  
  dir.create(paste(saveDir,"cell_cycle_phase/",sep = ""),showWarnings = F)
  pdf(paste(saveDir,"cell_cycle_phase/",objname,"_",Assay,"_",process,".pdf",sep = ""), width = 10, height = 7)
  print(DimPlot(CD4_RNA_integrated, reduction = "pca", group.by= "Phase", split.by = "Phase"))
  dev.off()
  
  pdf(paste(saveDir,"cell_cycle_phase/",objname,"_",Assay,"_",process,"_spltted.pdf",sep = ""), width = 6, height = 6)
  print(DimPlot(CD4_RNA_integrated, reduction = "pca", group.by= "Phase"))
  dev.off()
  
  CD4_RNA_integrated$mitoRatio <- CD4_RNA_integrated@meta.data$percent.mt / 100
  summary_mito <- summary(CD4_RNA_integrated@meta.data$mitoRatio) ## based on the summary need to choose which quartile it is present
  
  summary_df <- as.data.frame(as.matrix(summary_mito[2:5]))
  mid_values <- summary_df[c(1,2,4),]
  # Turn mitoRatio into categorical factor vector based on quartile values
  CD4_RNA_integrated@meta.data$mitoFr <- cut(CD4_RNA_integrated@meta.data$mitoRatio,
                                             breaks=c(-Inf, mid_values, Inf),
                                             labels=c("Low","Medium","Medium high", "High"))
  
  dir.create(paste(saveDir,"mitoeffect/",sep = ""),showWarnings = F)
  pdf(paste(saveDir,"mitoeffect/",objname,"_",Assay,"_",process,"_spltted.pdf",sep = ""), width = 8, height = 6)
  print(DimPlot(CD4_RNA_integrated, reduction = "pca", group.by= "mitoFr", split.by = "mitoFr"))
  dev.off()
  
  pdf(paste(saveDir,"mitoeffect/",objname,"_",Assay,"_",process,".pdf",sep = ""), width = 6, height = 6)
  print(DimPlot(CD4_RNA_integrated, reduction = "pca", group.by= "mitoFr"))
  dev.off()
  
  message(paste("\nMaking Elbow Plot\n"))
  dir.create(paste(saveDir,"elbow_plots/",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"elbow_plots/",objname,"_",process,"_",Assay,".pdf",sep = ""),width = 8, height = 8)
  print(ElbowPlot(CD4_RNA_integrated, ndims = 50))
  dev.off()
  
  dir.create(paste(saveDir,"saveRDS_obj/",sep = ""),showWarnings = F)
  message(paste("\n Saving",objname," ADT subsetted \n"))
  #saveRDS(CD4_RNA_integrated,paste(saveDir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = ""))
  saveRDS(CD4_ADT_integrated_impute,paste(saveDir,"saveRDS_obj/",objname,"_",process,"_ADT.RDS",sep = ""))
  
  message("\nWow Please go to this location ",saveDir,"cell_cycle_phase/ for cellcycle effect\n and for mitoeffect go to",saveDir,"mitoeffect/ to identify which effect needs to be regressed out from the object\n")
  return(CD4_RNA_integrated)
}

# # please add the object path for RNA imputed integrated 
# RNA_obj_path <- paste(savedir,"saveRDS_obj/combined_integrated_imputating_RNA.RDS",sep = "")
# ADT_obj_path <- paste(savedir,"saveRDS_obj/combined_ADT_UMAP_QC_integrated.RDS",sep = "")
# CD4_cluster <- c(0,5,8,9) ## Need to see the UMAP to find the which cluster belongs to the CD4
# 
# objname = "CD4"
# Assay = "RNA"
# process = "subset"
# 
# RNA_subsetting(RNA_obj_path = RNA_obj_path, 
#                ADT_obj_path = ADT_obj_path,
#                subset_cluster = CD4_cluster,
#                ngenes = 4000, 
#                saveDir = savedir,
#                objname = "CD4",
#                Assay = "RNA",
#                process = "subset")