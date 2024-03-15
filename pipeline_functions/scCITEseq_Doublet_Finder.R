### Args
# Obj = The RNA object created after QC using scCITE_QC.R function
# dims = PCA dims to include for UMAP
# res = Resolution to be used for Louvian clustering 
# Obj2 = The ADT object created after QC using scCITE_QC.R function
# samplename = The name of the sample for saving the file
# process = To save the output we will use this name to identify at which process this file is created
# Assay = Which assay to be considered for doublet filtering
doublet_scCITEseq <- function(Obj,dims,res,saveDir,Obj2,samplename,process,Assay){
  message(paste("\n\n######################################## Processsing",samplename,"#####################################################\n\n"))
  GEX_obj_3 <- Obj
  DefaultAssay(GEX_obj_3) <- Assay
  dir.create(saveDir,showWarnings = FALSE)
  
  message("\nFinding Cluster and Running UMAP\n")
  GEX_obj_3 <- FindClusters(GEX_obj_3, resolution = res, verbose = FALSE)
  dir.create(paste(saveDir,"UMAP/",sep = ""), showWarnings = FALSE)
  
  message(paste("\nSaving UMAP at this location:",saveDir,"UMAP/ for",samplename,"\n"))
  GEX_obj_3 <- RunUMAP(GEX_obj_3,dims = 1:dims, reduction = "pca")
  pdf(paste(saveDir,"UMAP/",samplename,"_",Assay,".pdf", sep = ""), width = 5, height = 5)
  print(DimPlot(GEX_obj_3, reduction = "umap"))
  dev.off()
  
  ##### Doublet Finder ######
  message(paste("\nRunning Doublet Finder for",samplename,"\n"))
  message("\nIt will takes some time, Take some time Off \n You are doing great\n")
  library(DoubletFinder)
  sweep.res.list <- paramSweep_v3(GEX_obj_3, PCs = 1:dims, sct = FALSE) # sctranform = FALSE
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  pK=as.numeric(as.character(bcmvn$pK))
  BCmetric=bcmvn$BCmetric
  pK_choose = pK[which(BCmetric %in% max(BCmetric))]
  
  dir.create(paste(saveDir,"bcmvn/",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"bcmvn/",samplename,"_",process,".pdf",sep = ""))
  par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
  plot(x = pK, y = BCmetric, pch = 16,type="b",col = "blue",lty=1)
  abline(v=pK_choose,lwd=2,col='red',lty=2)
  title("The BCmvn distributions")
  text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")
  dev.off()
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(GEX_obj_3@meta.data$seurat_clusters) ## ex: annotations <- spleen_subset@meta.data$seurat_clusters
  nExp_poi <- round(0.075*nrow(GEX_obj_3@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  GEX_obj_3 <- doubletFinder_v3(GEX_obj_3, PCs = 1:dims, pN = 0.25, pK = pK_choose, 
                                nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  pANN <- colnames(GEX_obj_3@meta.data[grep("DF",colnames(GEX_obj_3@meta.data))])[1]
  
  GEX_obj_3 <- doubletFinder_v3(GEX_obj_3, PCs = 1:dims, pN = 0.25, pK = pK_choose,
                                nExp = nExp_poi.adj, reuse.pANN = pANN, 
                                sct = FALSE)
  
  pdf(paste(saveDir,"UMAP/",samplename,"_",process,".pdf",sep = ""))
  print(DimPlot(GEX_obj_3, pt.size = 1,reduction = "umap",group.by = pANN)+theme(aspect.ratio = 1))
  dev.off()
  
  #Doublet removal
  message("\nRemoving Doublets\n")
  pANN_index <- grep("DF",colnames(GEX_obj_3@meta.data))[1]
  singlet_idx <- which(GEX_obj_3@meta.data[pANN_index] == "Singlet")
  GEX_obj_4 <- subset(GEX_obj_3,cells=rownames(GEX_obj_3@meta.data[singlet_idx,]))
  
  message("\n subsetting the ADT project also\n")
  ADT_obj_3 <- readRDS(Obj2)
  ADT_obj_4 <- subset(ADT_obj_3, cells = rownames(GEX_obj_4@meta.data))
  
  condition <- all(colnames(ADT_obj_4) == colnames(GEX_obj_4))
  print(condition)
  stopifnot(condition == "TRUE")
  
  ### Saving Singlet samples
  message("\n Saving the GEX  and ADT Singlet RDS\n")
  saveRDS(GEX_obj_4,paste(saveDir,"saveRDS_obj/",samplename,"_RNA_singlet.RDS",sep = ""))
  saveRDS(ADT_obj_4,paste(saveDir,"saveRDS_obj/",samplename,"_ADT_singlet.RDS",sep = ""))
  
  message(paste(samplename,"Doublets are removed and both RNA and ADT has same number of cells\n"))
  message(paste(samplename,"UMAP is at this location: ",saveDir,"UMAP/\n"))
  message(paste("Please check the distribution of cells after getting filtered and doublets removed at this location ",saveDir,"UMAP/\n"))
  
  a <- nrow(GEX_obj_3@meta.data)
  b <- nrow(GEX_obj_4@meta.data)
  c <- a-b
  d=c(a,b,c)
  return(d)
}
