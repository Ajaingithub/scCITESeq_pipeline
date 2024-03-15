###### UMAP object 
# obj_path: "RNA"
# saveDir: Saved directory path
# dims: PCA dimension to be used
# RNA_features: Genes featureplot
# Assay: RNA or ADT used
# objname: Name of the object
# ncol: splitted UMAP to be divided into number of cols
# ndims: For Elbow plots number of dims
RNA_integration <- function(obj,saveDir,dims=20,RNA_features = c("CD4","CD8A"), 
                            Assay="RNA", process="integration", objname, ncol=2, ndims = 50){
  message("\n loading the RNA/ADT integrated Object \n")
  RNA_integrated <- obj
  RNA_integrated <- RunUMAP(RNA_integrated, dims = 1:dims, reduction = "pca")
  
  message("\n Saving the UMAP \n")
  # filename <- gsub("_combined.RDS","",basename(obj_path))
  dir.create(paste(saveDir,"UMAP/",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"UMAP/",objname,"_",process,"_",Assay,".pdf",sep = ""),width = 6, height = 5)
  print(DimPlot(RNA_integrated, group.by = "orig.ident"))
  dev.off()
  
  pdf(paste(saveDir,"UMAP/",objname,"_",process,"_",Assay,"_splitted.pdf",sep = ""),width = 10, height = 8)
  print(DimPlot(RNA_integrated, split.by = "orig.ident", group.by="orig.ident", ncol = ncol))
  dev.off()
  
  pdf(paste(saveDir,"UMAP/",objname,"_",process,"_",Assay,"_vaccine.pdf",sep = ""),width = 6, height = 5)
  print(DimPlot(RNA_integrated, group.by = "vaccine"))
  dev.off()
  
  pdf(paste(saveDir,"UMAP/",objname,"_",process,"_",Assay,"_vaccine_splitted.pdf",sep = ""),width = 10, height = 8)
  print(DimPlot(RNA_integrated, split.by = "vaccine", group.by="vaccine", ncol = 3))
  dev.off()
  
  message("\n Making Featureplot \n")
  dir.create(paste(saveDir,"featureplot",sep = ""), showWarnings = FALSE)
  DefaultAssay(RNA_integrated) <- Assay
  RNA_features_name = gsub("-protein","",paste(RNA_features, collapse = "_"))
  pdf(paste(saveDir,"featureplot/",objname,"_",process,"_",Assay,"_",RNA_features_name,".pdf",sep = ""),width = 10, height = 5)
  print(FeaturePlot(RNA_integrated,RNA_features))
  dev.off()
  
  ### Imputing the matrix
  message("\nImputing the matrix\n")
  DefaultAssay(RNA_integrated) <- Assay
  # seurat_integrated <- subset(integrated, features = rownames(integrated[-gene_zero,])) # Since only one row is zero no need of subsetting it
  library(reticulate)
  use_python("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/python/python-3.9.2/bin/python3.9")
  library(Rmagic)
  py_discover_config("magic") # to check
  RNA_integrated_impute <- magic(RNA_integrated, npca=dims) ## imputing the RNA data as for RNA PCs are 20
  DefaultAssay(RNA_integrated_impute) <- paste("MAGIC",Assay,sep = "_")
  RNA_integrated_impute <- ScaleData(RNA_integrated_impute,features=rownames(RNA_integrated_impute)) # Sacling data on all the features as we need it during the feature plot scaled data for gene that are not variable
  
  message("\n Making Imputed featureplot\n")
  DefaultAssay(RNA_integrated_impute) <- paste("MAGIC",Assay,sep = "_")
  pdf(paste(saveDir,"featureplot/",objname,"_",process,"_",Assay,"_",RNA_features_name,"_imputed.pdf",sep = ""),width = 10, height = 5)
  print(FeaturePlot(RNA_integrated_impute,RNA_features))
  dev.off()  
  
  DefaultAssay(RNA_integrated_impute) <- "integrated"
  dir.create(paste(saveDir,"elbow_plots/",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"elbow_plots/",objname,"_",process,"_",Assay,".pdf",sep = ""),width = 8, height = 8)
  print(ElbowPlot(RNA_integrated_impute, ndims = ndims))
  dev.off()
  
  #message("Saving the integrated Imputed Object")
  #dir.create(paste(saveDir,"saveRDS_obj/",sep = ""), showWarnings = FALSE)
  #saveRDS(RNA_integrated_impute,paste(saveDir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = ""))
  message(paste("Please go through the UMAP at ",saveDir,"UMAP/ to check if it is properly ran or not", sep = ""))
  message(paste("Please go through the Elbow Plot at ",saveDir,"elbow_plots/ to check for the Nearest Neighour", sep = ""))
  message(paste("Please go through the feature plot of ",RNA_features," at ",saveDir,"featureplot/ to check if it is properly ran or not",sep = ""))
  return(RNA_integrated_impute)
}

# obj <- paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")
# ## For each sample in RNA assay they had to be regress out for the Mitochondria, however the cell cycle phase looks good
# ## Remember we are not integrating we are merging the samples
# ## By looking at the graph, we add the mitoratio or cellcycle (G2M.Score,S.Score)
# 
# objname = "combined_integrated"
# process = "imputating"
# Assay = "RNA"
# 
# 
# RNA_integration(obj_path = obj,saveDir = savedir, dims = 30, RNA_features = c("CD4","CD8A"), 
#                 Assay = Assay, process = process, objname = objname)
