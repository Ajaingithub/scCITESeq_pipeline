## This function is used to find the cluster in the modality (RNA and ADT) and batch integrated object
# Args:
# modal_integrate_obj: Batch and Modality integrated object
# saveDir: save directory
# res: resolution to perform louvian clustering
# RNA_features: Genes to make featureplots like CD4 and CD8A
# ADT_features: proteins to make featuresplots, like CD4_protein or CD8A_protein
# Assay: which assay to use, like RNA, integrated
# protein_assay: "ADT"
# process: "Modality integration" In order to save the files
# objname: name of the samples to save files
# col_sel: metadata column selected for performing QC and making a bar graph for different samples within that column like Old VS Young
# scalingdata: while makeing the heatmap do you need to perform the scaling or it is already done?

Modal_Cluster_and_UMAP_QC <- function(modal_integrate_obj, saveDir , res, RNA_features, ADT_features,
                                      Assay, protein_assay="ADT", process, objname,col_sel="orig.ident",scalingdata=TRUE){
  library(Seurat)
  library(ggplot2)
  library(reshape2)
  message("Running Clustering on RNA and ADT integration\nloading the object")

  modal_integrated <- modal_integrate_obj
  modal_integrated <- FindClusters(modal_integrated, resolution = res,
                                   graph.name = "wsnn", algorithm = 3, verbose = FALSE)

  
  source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scRNASeq/pipeline_functions/Heatmap.R")
  AJ_heatmap(modal_integrated,savedir = saveDir, objname = paste(objname,res,sep = "_"), assay=protein_assay,scalingdata=scalingdata)

  CD8_genes <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Hirohisa/scCITESeq/SARS_CoV_2/Run01/Analysis/downstream/Table/CD8_genes")
  CD8_genes_vec <- CD8_genes[,1]
  CD8_genes_vec <- CD8_genes_vec[grep("^MT",CD8_genes_vec, invert = TRUE)]

  dir.create(paste(savedir,"dotplot",sep = ""), showWarnings = FALSE)
  pdf(paste(savedir,"dotplot/",objname,"_",process,"_",Assay,".pdf",sep = ""), width = 17, height = 14)
  print(DotPlot(modal_integrated, assay="RNA",features = CD8_genes_vec) + coord_flip())
  dev.off()

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

  for (i in 1:length(col_sel)) {
    message("Running for ",col_sel[i])
    n_cells_2 <- FetchData(modal_integrated,vars = c("ident", col_sel[i]))
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

    pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_",col_sel[i],"_histogram.pdf",sep = ""), width = 20, height = 10)
    print(p)
    dev.off()

    p <- ggplot(df_melted, aes_string(x="Cluster", y=col_sel[i], size="Cell_percentage", color=col_sel[i])) +
      geom_point(alpha=.75) +
      scale_size(range = c(1,15), breaks=seq(0,45,by=15)) +
      theme_bw()

    pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_",col_sel[i],"_bubbleplot.pdf",sep = ""), width = 14, height = 8)
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

    pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_",col_sel[i],"_histogram_within_cluster.pdf",sep = ""), width = 18, height = 9)
    print(p)
    dev.off()

    p <- ggplot(df_melted, aes_string(x="Cluster", y=col_sel[i], size="Cell_percentage", color=col_sel[i])) +
      geom_point(alpha=.75) +
      scale_size(range = c(1,15), breaks=seq(0,45,by=15)) +
      theme_bw()

    pdf(paste(saveDir,"Table/",objname,"_",Assay,"_",process,"_",col_sel[i],"_bubbleplot_within_cluster.pdf",sep = ""), width = 14, height = 8)
    print(p)
    dev.off()
  }

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
  DefaultAssay(modal_integrated) <- "ADT"
  ADT_features_name <- paste(ADT_features, collapse = "_")
  pdf(paste(saveDir,"featureplot/",objname,"_",process,"_",Assay,"_",ADT_features_name,".pdf",sep = ""),width = 10, height = 5)
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

## To Run the Code
# obj_path <- paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")
# RNA_features = c("CD4","CD8A")
# ADT_features = c("CD4-protein","CD8a-protein")
#
# objname = "CD4"
# Assay = "integrated"
# process = "modality_cluster_and_QC"
# Modal_Cluster_and_UMAP_QC(modal_integrate_obj_path = obj_path, saveDir = savedir, res = 0.3, RNA_features = RNA_features,
#                           ADT_features = ADT_features,Assay = Assay, process = process, objname=objname)
