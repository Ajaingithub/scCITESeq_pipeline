# obj = RNA or ADT object
# res = choose the resolution for clustering
# saveDir = save directory
# QC_features = nCountRNA, nFeature_RNA
# objname = samplename
# process = "clustering"
# col_sel = need to make the bar graph for the columns like Old and Young, and different Vaccine 

cluster_UMAP_and_QC <- function(obj_path,dims,res,saveDir, Assay = "RNA", QC_features ,objname, process, col_sel = "orig.ident"){
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  message("\nLoading the object\n")
  ADT_integrated_impute <- obj_path

  DefaultAssay(ADT_integrated_impute) <- Assay
  ADT_integrated_impute <- FindClusters(ADT_integrated_impute, resolution = res, algorithm = 1, verbose = FALSE)

  dir.create(paste(saveDir,"UMAP/",sep = ""),showWarnings = FALSE)
  pdf(paste(saveDir,"UMAP/",objname,"_",Assay,"_",process,".pdf",sep = ""), width = 6, height = 6)
  print(DimPlot(ADT_integrated_impute, reduction = "umap", label = TRUE, raster = TRUE))
  dev.off()

  dir.create(paste(saveDir,"QC_Vln/",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"QC_Vln/",objname,"_",Assay,"_",process,".pdf",sep = ""), width = 9, height = 6)
  print(VlnPlot(ADT_integrated_impute, QC_features))
  dev.off()

  for (i in 1:length(col_sel)) {
    message("Running for ",col_sel[i])
    n_cells_2 <- FetchData(ADT_integrated_impute,vars = c("ident", col_sel[i]))
    n_cells <- table(n_cells_2[,col_sel[i]],n_cells_2[,"ident"])
    n_cells_sum <- as.vector(rowSums(n_cells))

    ### Making an empty dataframe
    df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
    rownames(df) <- rownames(n_cells)
    colnames(df) <- colnames(n_cells)

    for (j in 1:nrow(n_cells)) {
      df[j,] <- (n_cells[j,]/n_cells_sum[j])*100
    }

    df$col <- rownames(df)
    df_melted <- melt(df)
    colnames(df_melted) <- c(col_sel[i],"Cluster","Cell_percentage")
    p <- ggplot(df_melted, aes_string(x="Cluster", y="Cell_percentage", fill=col_sel[i])) +
      geom_bar(stat = "identity",color="black", position = "dodge") +
      theme_bw()

    dir.create(paste(saveDir,"Table/",sep = ""),showWarnings = F)

    pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_",col_sel[i],"_histogram.pdf",sep = ""), width = 20, height = 8)
    print(p)
    dev.off()

    p <- ggplot(df_melted, aes_string(x="Cluster", y=col_sel[i], size="Cell_percentage", color=col_sel[i])) +
      geom_point(alpha=.75) +
      scale_size(range = c(1,15), breaks=seq(0,45,by=15)) +
      theme_bw()

    pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_",col_sel[i],"_bubbleplot.pdf",sep = ""), width = 20, height = 8)
    print(p)
    dev.off()

    dir.create(paste(saveDir,"Table/",sep = ""),showWarnings = F)
    write.table(n_cells,
                paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_cluster_and_",col_sel[i],".txt",sep = ""),
                quote = F,
                row.names = T,
                col.names = T)

    ### Calculating the percentage in each cluster without normalizing it
    n_cells_sum <- as.vector(colSums(n_cells))

    ### Making an empty dataframe
    df <- data.frame(matrix(nrow = nrow(n_cells), ncol = ncol(n_cells)))
    rownames(df) <- rownames(n_cells)
    colnames(df) <- colnames(n_cells)

    j=0
    for (j in 1:ncol(n_cells)) {
      df[,j] <- (n_cells[,j]/n_cells_sum[j])*100
    }

    df$col <- rownames(df)
    df_melted <- melt(df)
    colnames(df_melted) <- c(col_sel[i],"Cluster","Cell_percentage")
    p <- ggplot(df_melted, aes_string(x="Cluster", y="Cell_percentage", fill=col_sel[i])) +
      geom_bar(stat = "identity",color="black", position = "dodge") +
      theme_bw() + ggtitle(paste(objname," within cluster percentage calculation"))

    pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_",col_sel[i],"_histogram_within_cluster.pdf",sep = ""), width = 20, height = 8)
    print(p)
    dev.off()

    p <- ggplot(df_melted, aes_string(x="Cluster", y=col_sel[i], size="Cell_percentage", color=col_sel[i])) +
      geom_point(alpha=.75) +
      scale_size(range = c(1,15), breaks=seq(0,45,by=15)) +
      theme_bw()

    pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_",col_sel[i],"_bubbleplot_within_cluster.pdf",sep = ""), width = 20, height = 8)
    print(p)
    dev.off()
  }

  ### Making a heatmap
  source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/Heatmap.R")
  AJ_heatmap(ADT_integrated_impute, saveDir, paste(objname,"_",process))

  message(paste("\n Please check for cluster and VlnPlot at this location ",saveDir,"UMAP/\nVln QC at this location ",saveDir,"QC_Vln/",sep = ""))
  message(paste("\n Check for the Table with Shingrix Zostavax and DMSO for each cluster"))
  return(ADT_integrated_impute)
}
