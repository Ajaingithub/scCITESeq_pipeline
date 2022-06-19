## should have vaccine column in the meta data of the object
Modal_Cluster_and_UMAP_QC <- function(modal_integrate_obj_path, saveDir , res, RNA_features, ADT_features, Assay, process, objname){
  library(reshape2)
  message("Running Clustering on RNA and ADT integration\nloading the object")
  
  modal_integrated <- modal_integrate_obj_path
  modal_integrated <- FindClusters(modal_integrated, resolution = res,
                                   graph.name = "wsnn", algorithm = 3, verbose = FALSE)
  
  message("Saving the Cluster UMAP")
  pdf(paste(saveDir,"UMAP/",objname,"_",process,"_",Assay,".pdf",sep = ""), width = 5, height = 5)
  print(DimPlot(modal_integrated, reduction = "wnn.umap", label = TRUE))
  dev.off()
  
  message("Performing the cluster QC")
  ## Since we have already removed the doublet when we ran the CellRanger no Need to run it again
  ### Remove those cluster which has low nFeature RNA or nFeature count and mt.percent -----------------------------------------
  pdf(paste(saveDir,"QC_Vln/",objname,"_",Assay,"_",process,".pdf",sep = ""), width = 20, height = 10)
  print(VlnPlot(modal_integrated, features = c("nFeature_RNA", "nCount_RNA","nFeature_ADT", "nCount_ADT", "percent.mt"), ncol = 5))
  dev.off()
  
  ### Comparing the clusters based on samples and vaccine ######
  # Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
  dir.create(paste(saveDir,"Table",sep = ""), showWarnings = FALSE)
  n_cells <- FetchData(modal_integrated, 
                       vars = c("ident", "orig.ident")) %>%
    dplyr::count(ident, orig.ident) %>%
    tidyr::spread(ident, n)
  
  write.table(n_cells, paste(saveDir,"Table/",objname,"_",process,"_",Assay,".txt", sep = ""), quote = F, row.names = F, col.names = T)
  
  rownames(n_cells) <- n_cells$orig.ident
  n_cells <- n_cells[,-1] # to remove the first colummn
  n_cells[is.na(n_cells)] <- 0
  n_cells_sum <- as.vector(rowSums(n_cells))
  
  ### Making an empty dataframe 
  df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
  rownames(df) <- rownames(n_cells)
  colnames(df) <- colnames(n_cells)
  
  for (i in 1:nrow(n_cells)) {
    df[i,] <- (n_cells[i,]/n_cells_sum[i])*100
  }
  
  df$SampleID <- rownames(df)
  library(reshape2)
  df_melted <- melt(df)
  colnames(df_melted) <- c("SampleID","Cluster","Cell_percentage")
  p <- ggplot(df_melted, aes(x=Cluster, y=Cell_percentage, fill=SampleID)) +
    geom_bar(stat = "identity",color="black", position = "dodge") +
    theme_bw() + ggtitle(paste(objname," normalizing across the samples"))
  
  pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_sample_histogram.pdf",sep = ""), width = 12, height = 8)
  print(p)
  dev.off()
  
  p <- ggplot(df_melted, aes(x=Cluster, y=SampleID, size=Cell_percentage, color = SampleID)) +
    geom_point(alpha=.75) +
    scale_size(range = c(1,15), breaks=seq(0,30,by=10)) +
    theme_bw() + ggtitle(paste(objname,"Cluster Sample"))
  
  pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_sample_bubbleplot.pdf",sep = ""), width = 12, height = 8)
  print(p)
  dev.off()
  
  dir.create(paste(saveDir,"Table/",sep = ""),showWarnings = F)
  write.table(n_cells,
              paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_cluster_and_samples.txt",sep = ""),
              quote = F,
              row.names = T,
              col.names = T)
  
  ### Calculating the percentage in each cluster without normalizing it
  n_cells_sum <- as.vector(colSums(n_cells))
  
  ### Making an empty dataframe 
  df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
  rownames(df) <- rownames(n_cells)
  colnames(df) <- colnames(n_cells)
  
  for (i in 1:ncol(n_cells)) {
    df[,i] <- (n_cells[,i]/n_cells_sum[i])*100
  }
  
  df$SampleID <- rownames(df)
  df_melted <- melt(df)
  colnames(df_melted) <- c("SampleID","Cluster","Cell_percentage")
  p <- ggplot(df_melted, aes(x=Cluster, y=Cell_percentage, fill=SampleID)) +
    geom_bar(stat = "identity",color="black", position = "dodge") +
    theme_bw() + ggtitle(paste(objname," within cluster percentage calculation"))
  
  pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_sample_histogram_within_cluster.pdf",sep = ""), width = 12, height = 8)
  print(p)
  dev.off()
  
  p <- ggplot(df_melted, aes(x=Cluster, y=SampleID, size=Cell_percentage, color = SampleID)) +
    geom_point(alpha=.75) +
    scale_size(range = c(1,15), breaks=seq(0,90,by=15)) +
    theme_bw() + ggtitle(paste(objname,"within cluster Sample"))
  
  pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_sample_bubbleplot_within_cluster.pdf",sep = ""), width = 12, height = 8)
  print(p)
  dev.off()
  
  dir.create(paste(saveDir,"Table/",sep = ""),showWarnings = F)
  write.table(n_cells,
              paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_cluster_and_samples.txt",sep = ""),
              quote = F,
              row.names = T,
              col.names = T)
  
  #### Based on the vaccine
  n_cells <- FetchData(modal_integrated,
                       vars = c("ident", "vaccine")) %>%
    dplyr::count(ident, vaccine) %>%
    tidyr::spread(ident, n)

  rownames(n_cells) <- n_cells$vaccine
  n_cells <- n_cells[,-1] # to remove the first colummn
  n_cells[is.na(n_cells)] <- 0
  n_cells_sum <- as.vector(rowSums(n_cells))

  ### Making an empty dataframe
  df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
  rownames(df) <- rownames(n_cells)
  colnames(df) <- colnames(n_cells)

  for (i in 1:nrow(n_cells)) {
    df[i,] <- (n_cells[i,]/n_cells_sum[i])*100
  }

  df$Vaccine <- rownames(df)
  df_melted <- melt(df)
  colnames(df_melted) <- c("Vaccine","Cluster","Cell_percentage")
  p <- ggplot(df_melted, aes(x=Cluster, y=Cell_percentage, fill=Vaccine)) +
    geom_bar(stat = "identity",color="black", position = "dodge") +
    theme_bw() + ggtitle(paste(objname," normalizing across the samples"))

  pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_vaccine_histogram.pdf",sep = ""), width = 10, height = 6)
  print(p)
  dev.off()


  write.table(n_cells,
              paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_cluster_and_vaccine.txt",sep = ""),
              quote = F,
              row.names = T,
              col.names = T)


  ### Percentage calculation based on the cluster not normalizing
  n_cells_sum <- as.vector(colSums(n_cells))

  ### Making an empty dataframe
  df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
  rownames(df) <- rownames(n_cells)
  colnames(df) <- colnames(n_cells)

  for (i in 1:ncol(n_cells)) {
    df[,i] <- (n_cells[,i]/n_cells_sum[i])*100
  }

  df$Vaccine <- rownames(df)
  df_melted <- melt(df)
  colnames(df_melted) <- c("Vaccine","Cluster","Cell_percentage")
  p <- ggplot(df_melted, aes(x=Cluster, y=Cell_percentage, fill=Vaccine)) +
    geom_bar(stat = "identity",color="black", position = "dodge") +
    theme_bw() + ggtitle(paste(objname," within cluster calculation"))

  pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_vaccine_histogram_percentage_calculation.pdf",sep = ""), width = 10, height = 6)
  print(p)
  dev.off()
  
  
  ## https://github.com/satijalab/seurat/issues/2137 For feature plot we need to use the RNA assay rather than integrated or other assay
  DefaultAssay(modal_integrated) <- "RNA"
  
  RNA_features_name <- paste(RNA_features, collapse = "_")
  pdf(paste(saveDir,"featureplot/",objname,"_",process,"_",Assay,"_",RNA_features_name,".pdf",sep = ""),width = 10, height = 5)
  print(FeaturePlot(modal_integrated, 
                    reduction = "wnn.umap", 
                    features = RNA_features, 
                    order = TRUE,
                    label = TRUE))
  dev.off()
  
  # Select the RNA counts slot to be the default assay
  DefaultAssay(modal_integrated) <- "IADT"
  ADT_features_name <- paste(ADT_features, collapse = "_")
  pdf(paste(saveDir,"featureplot/",objname,"_",process,"_",Assay,"_",ADT_features_name,".pdf",sep = ""),width = 10, height = 5)
  print(FeaturePlot(modal_integrated, 
                    reduction = "wnn.umap", 
                    features = ADT_features, 
                    order = TRUE,
                    label = TRUE))
  dev.off()
  
  ## With imputation
  DefaultAssay(modal_integrated) <- "MAGIC_RNA"
  pdf(paste(saveDir,"featureplot/",objname,"_",process,"_",Assay,"_",RNA_features_name,"_imputed.pdf",sep = ""),width = 10, height = 5)
  print(FeaturePlot(modal_integrated, 
                    reduction = "wnn.umap", 
                    features = RNA_features, 
                    order = TRUE,
                    label = TRUE))
  dev.off()
  
  # Select the RNA counts slot to be the default assay
  DefaultAssay(modal_integrated) <- "MAGIC_ADT"
  pdf(paste(saveDir,"featureplot/",objname,"_",process,"_",Assay,"_",ADT_features_name,"._imputed.pdf",sep = ""),width = 10, height = 5)
  print(FeaturePlot(modal_integrated, 
                    reduction = "wnn.umap", 
                    features = ADT_features, 
                    order = TRUE,
                    label = TRUE))
  dev.off()
  
  message("Please go to this location for the clustering ",saveDir,"/UMAP/\n")
  message("Please go to this location for the featurePlot ",saveDir,"/featureplot/\n")
  message("Please go to this location for the featurePlot ",saveDir,"/QC_Vln/\n")
  message("Please go to this location for the cell number in each cluster with sample and vaccine ",saveDir,"/Table/\n")
  # saveRDS(modal_integrated, paste(saveDir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = ""))
  return(modal_integrated)
}

# obj_path <- paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")
# RNA_features = c("CD4","CD8A")
# ADT_features = c("CD4-protein","CD8a-protein")
# 
# objname = "CD4"
# Assay = "integrated"
# process = "modality_cluster_and_QC"
# Modal_Cluster_and_UMAP_QC(modal_integrate_obj_path = obj_path, saveDir = savedir, res = 0.3, RNA_features = RNA_features, 
#                           ADT_features = ADT_features,Assay = Assay, process = process, objname=objname)