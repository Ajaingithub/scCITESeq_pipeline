cluster_UMAP_and_QC <- function(obj_path,dims,res,saveDir, Assay = "RNA", QC_features ,objname, process){
  library(dplyr)
  library(reshape2)
  message("\nLoading the object\n")
  ADT_integrated_impute <- obj_path
  
  DefaultAssay(ADT_integrated_impute) <- Assay
  ADT_integrated_impute <- FindClusters(ADT_integrated_impute, resolution = res, algorithm = 1, verbose = FALSE)
  
  dir.create(paste(saveDir,"UMAP/",sep = ""),showWarnings = FALSE)
  pdf(paste(saveDir,"UMAP/",objname,"_",Assay,"_",process,".pdf",sep = ""), width = 6, height = 6)
  print(DimPlot(ADT_integrated_impute, reduction = "umap", label = TRUE))
  dev.off()
  
  dir.create(paste(saveDir,"QC_Vln/",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"QC_Vln/",objname,"_",Assay,"_",process,".pdf",sep = ""), width = 9, height = 6)
  print(VlnPlot(ADT_integrated_impute, QC_features))
  dev.off()
  
  message("saving the confusion matrix for sample, vaccine and cluster")
  n_cells <- FetchData(ADT_integrated_impute, 
                       vars = c("ident", "orig.ident")) %>%
    dplyr::count(ident, orig.ident) %>%
    tidyr::spread(ident, n)
  
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
  df_melted <- melt(df)
  colnames(df_melted) <- c("SampleID","Cluster","Cell_percentage")
  p <- ggplot(df_melted, aes(x=Cluster, y=Cell_percentage, fill=SampleID)) +
    geom_bar(stat = "identity",color="black", position = "dodge") +
    theme_bw()
  
  dir.create(paste(saveDir,"Table/",sep = ""),showWarnings = F)

  pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_sample_histogram.pdf",sep = ""), width = 12, height = 8)
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
  
  
  ##### Vaccine
  n_cells <- FetchData(ADT_integrated_impute, 
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
    theme_bw()

  pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_vaccine_histogram.pdf",sep = ""), width = 10, height = 6)
  print(p)
  dev.off()

  
  write.table(n_cells,
              paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_cluster_and_vaccine.txt",sep = ""),
              quote = F,
              row.names = T,
              col.names = T)
  
  #### Percentage calculation based on the cluster not normalizing
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
  
  
  dir.create(paste(saveDir,"saveRDS_obj/",sep = ""), showWarnings = F)
  message("\nSaving the object\n")
  # saveRDS(ADT_integrated_impute, paste(saveDir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = ""))
  
  message(paste("\n Please check for cluster and VlnPlot at this location ",saveDir,"UMAP/\nVln QC at this location ",saveDir,"QC_Vln/",sep = ""))
  message(paste("\n Check for the Table with Shingrix Zostavax and DMSO for each cluster"))
  return(ADT_integrated_impute)
}
