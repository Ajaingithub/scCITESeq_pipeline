modality_integration <- function(RNA_obj_path, ADT_obj_path, RNA_dims, ADT_dims, saveDir, Assay, process, objname){
  message("\n Reading ADT and RNA object \n")
  RNA_integration <- RNA_obj_path
  ADT_integration <- ADT_obj_path
  
  ### adding ADT counts and Feature to the integrated object
  condition = all(rownames(RNA_integration@meta.data)==rownames(ADT_integration@meta.data))
  print(condition)
  stopifnot(condition == "TRUE")
  
  ## We will add the ADT assay to the RNA object
  RNA_integration[["IADT"]] <- ADT_integration[["integrated"]]
  RNA_integration[["MAGIC_ADT"]] <- ADT_integration[["MAGIC_ADT"]]
  RNA_integration[["ADT"]] <- ADT_integration[["ADT"]]
  RNA_integration[["pcaadt"]] <- ADT_integration[["pca"]]
  RNA_integration[["umapadt"]] <- ADT_integration[["umap"]]
  
  RNA_integration@meta.data$nFeature_ADT <- ADT_integration@meta.data$nFeature_ADT
  RNA_integration@meta.data$nCount_ADT <- ADT_integration@meta.data$nCount_ADT
  
  # PC for RNA is RNA_dims and ADT also ADT_dims
  ### the integration corrects normalized data, so there is no need to do another log-normalization after integration.
  message("\n finding multimodal neighbours \n")
  DefaultAssay(RNA_integration) <- Assay
  RNA_integration <- FindMultiModalNeighbors(RNA_integration,
                                             reduction.list = list("pca", "pcaadt"),
                                             dims.list = list(1:RNA_dims, c(1:ADT_dims)),
                                             modality.weight.name = c("RNA.weight", "ADT.weight"))
  
  message("\n Running Integration \n")
  RNA_integration <- RunUMAP(RNA_integration, nn.name = "weighted.nn",
                             reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  
  pdf(paste(saveDir,"UMAP/",objname,"_",process,"_",Assay,".pdf",sep = ""),width = 5, height = 5)
  print(DimPlot(RNA_integration, reduction = "wnn.umap", group.by = "orig.ident"))
  dev.off()
  
  pdf(paste(saveDir,"UMAP/",objname,"_",process,"_",Assay,"_split.pdf",sep = ""),width = 12, height = 8)
  print(DimPlot(RNA_integration, reduction = "wnn.umap", group.by = "orig.ident", split.by = "orig.ident", ncol = 5))
  dev.off()
  
  pdf(paste(saveDir,"UMAP/",objname,"_",process,"_",Assay,"_vaccine.pdf",sep = ""),width = 5, height = 5)
  print(DimPlot(RNA_integration, reduction = "wnn.umap", group.by = "vaccine"))
  dev.off()

  pdf(paste(saveDir,"UMAP/",objname,"_",process,"_",Assay,"_vaccine_split.pdf",sep = ""),width = 9, height = 6)
  print(DimPlot(RNA_integration, reduction = "wnn.umap", group.by = "vaccine", split.by = "vaccine", ncol = 3))
  dev.off()
  
  ### Indrntfying the resolution in the data using cluster
  ## It does not depends on the Default assay as we are using the weighted from both the modalities
  message("\n Running Cluster Tree \n")
  rm(forres)
  forres <- RNA_integration
  resolution <- c(0.001,0.01,0.03,0.05,0.07,0.09,0.1,0.125,0.15,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2)
  for(i in 1:length(resolution)){
    forres <- FindClusters(forres, resolution = resolution[i], graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  }
  message("\n Saving Tree \n")
  library(clustree)
  pdf(paste(saveDir,"clustertree/",objname,"_",Assay,"_",process,".pdf",sep = ""), width = 7, height = 12)
  print(clustree(forres, prefix = "wsnn_res."))
  dev.off()
  
  #saveRDS(RNA_integration,paste(saveDir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = ""))
  
  message(paste("\n we have ran the modality integration and UMAP, please check the UMAP and cluster tree resolution to be used for the FindCluster at this location \n",saveDir,"UMAP/\n and ",saveDir,"clustertree/",sep = ""))
  return(RNA_integration)
}

# RNA_obj_path <- paste(savedir,"saveRDS_obj/CD4_UMAP_and_QC_integrated.RDS",sep = "")
# ADT_obj_path <- paste(savedir,"saveRDS_obj/CD4_ADT_UMAP_QC_integrated.RDS",sep = "")
# 
# objname <- "CD4"
# process <- "modality_integrate"
# Assay <- "integrated"
# 
# modality_integration(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path, RNA_dims = 20, ADT_dims = 10,
#                      saveDir = savedir, Assay = Assay, process = process, objname = objname)
