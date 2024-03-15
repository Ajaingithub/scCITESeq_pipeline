## Making a Pseudobulk. Identifying the differential genes from the single cell experiments.
## Remove samples the samples which are an outlier and which you want to remove like DMSO
## Cluster 1 is first cluster Vs Cluster2 is the second cluster which you want to perform the differential
## You can have multiple cluster in both the cluster 1 and cluster 2
## column_name_split all the samples when splitted with _ or -. what will be the name of the column names when splitted Sprotein_O_BB22031_F_Run01 
## column_name_split = c("vaccine","age","sampleid","gender","run")
## splitting how one wants to split the column_name_split could be "_" or "-"
## Keep the batch column and sample splitting same like Run

COVID_pseudobulk_within_cluster_AJ <- function(obj,savedir, group1, group2, grouping_by, cluster="all",
                                         cluster_group = "seurat_clusters", cell_freq=10,
                                         remove_samples=NULL, gene_min_counts = 5, sample_col = "orig.ident",
                                         batch_col = "Run", column_name_split, splitting = "_"){
  message("Loading Packages")
  set.seed(123)
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(DESeq2)
  library(edgeR)
  library(RColorBrewer)
  # library(Glimma)
  library(ggplot2)
  library(ggrepel)
  library(session)
  library(DEFormats)
  # library(EnhancedVolcano)
  library(limma)
  
  
  if(cluster=="all"){
    print("all")
    dir.create(savedir,showWarnings = FALSE)
    savedir2 <- paste(savedir,"pseudobulk/",sep = "")
    dir.create(savedir2, showWarnings = FALSE)
    remove_samples_name <- paste(remove_samples,collapse="_and_")
    savedir3 <- paste(savedir2,"clus_",cluster,"_",group1,"_vs_",group2,"/",sep = "")
    dir.create(savedir3,showWarnings = FALSE)
    obj@meta.data[,paste(sample_col,grouping_by,sep = "_")] <- paste(obj@meta.data[,sample_col], obj@meta.data[,grouping_by], sep="_")
    
    DefaultAssay(obj) <- "RNA"
    cell_samples <- table(obj@meta.data[,paste(sample_col,grouping_by,sep = "_")]) %>% as.data.frame()
    sample_low <- cell_samples[cell_samples$Freq < cell_freq,1]
    sample_low <- gsub("_","_",sample_low)
    sample_low_2 <- c(sample_low,remove_samples)
    message(paste("Removing this sample",sample_low_2,"\n",sep = " "))
    Treg_cts <- AggregateExpression(obj, group.by = paste(sample_col,grouping_by,sep = "_"),
                                    assays = 'RNA', slot = "counts", return.seurat = FALSE)
    Treg_cts_2 <- as.data.frame(Treg_cts$RNA)
    if(length(sample_low_2) != 0){
      index <- grep(paste(sample_low_2,collapse="|"),colnames(Treg_cts_2), invert=TRUE)
      Treg_cts_3 <- Treg_cts_2[,index]
    } else {
      Treg_cts_3 <- Treg_cts_2
    }
    
    
    Treg_metadata <- data.frame(samples = colnames(Treg_cts_3))
    # Function to split and convert to dataframe
    split_to_dataframe <- function(x) {
      split_elements <- strsplit(x, splitting)[[1]]
      data.frame(t(split_elements))
    }
    
    # Convert split elements to dataframe
    Treg_metadata_2 <- do.call(rbind, lapply(Treg_metadata$samples, split_to_dataframe))
    colnames(Treg_metadata_2) <- c(column_name_split,grouping_by)
    Treg_metadata_2$orig.ident <- Treg_metadata$samples
    # Treg_metadata_2$Run <- Treg_metadata_2$orig.ident
    
    obj@meta.data[,paste(sample_col,batch_col,sep = "_")] <- paste(obj@meta.data[,sample_col],obj@meta.data[,batch_col],sep = "_")
    batch_group <- unique(obj@meta.data[,batch_col])
    # for (i in 1:length(batch_group)) {
    #   donor_change <- unique(obj@meta.data[,paste(sample_col,batch_col,sep = "_")]) %>% grep(batch_group[i],.,value=TRUE) %>% 
    #     gsub(paste(batch_group[i],sep = ""),"",.) %>% gsub("_","-",.) %>% paste(.,collapse=".*.|")
    #   donor_change <- paste(donor_change,".*.",sep="")
    #   Treg_metadata_2$Run <- gsub(donor_change,batch_group[i],Treg_metadata_2$Run)
    # }
    # Treg_metadata_2$Run <- gsub("-.*.","",Treg_metadata_2$Run)
    
    stopifnot(all(Treg_metadata_2$orig.ident == colnames(Treg_cts_3))) # if TRUE move forward
    # Treg_metadata$samplename <- gsub("gE_|DMSO_|_DMSO|_gE|B22015","",Treg_metadata$samples)
    ## Run will be treated as a covariate in the regression model,
    
    design_AJ <- as.formula(paste("~",paste(colnames(Treg_metadata_2[,c(batch_col,grouping_by)]),collapse = "+")))
    Treg_dds <- DESeqDataSetFromMatrix(countData = Treg_cts_3,
                                       colData = Treg_metadata_2,
                                       design = design_AJ)
    
    keep = filterByExpr(Treg_dds, group=colData(Treg_dds)[,grouping_by], min.count=gene_min_counts)
    table(keep)
    Treg_dds_2 <- Treg_dds[keep,]
    
    rld<-vst(Treg_dds_2)
    
    PCA <- DESeq2::plotPCA(object = rld,
                           intgroup = colnames(Treg_metadata_2),
                           returnData=TRUE,ntop = 5000)
    percentVar <- round(100 * attr(PCA, "percentVar"))
    
    savedir4 <- paste(savedir3,"PCA/",sep = "")
    dir.create(savedir4, showWarnings = FALSE)
    
    
    PC1 <-ggplot(PCA, aes_string("PC1", "PC2", label="name", color=batch_col, shape=grouping_by)) +
      geom_point(size=5, alpha=0.7) +
      geom_text_repel(size=2) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      ggtitle("PCA") +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    pdf(paste(savedir4,"PCA_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".pdf",sep = ""),
        width = 9, height = 8)
    print(PC1)
    dev.off()
    
    #### Further we are interested in all the PCs PC1 to all
    pc <- prcomp(t(assay(rld)))
    # attributes(pc)
    loading <- pc$rotation 
    summary <- summary(pc)
    
    rownames(Treg_metadata_2) <- Treg_metadata_2$orig.ident
    pc_with_metadata <- merge(pc$x,Treg_metadata_2,by = 'row.names', all = TRUE)
    
    write.table(pc_with_metadata, paste(savedir4,"PCA_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,"_all_PCs.txt",sep = ""),
                quote = F, row.names = T, col.names = T, sep = "\t")
    
    rm(plot_list)
    plot_list <- list()
    
    plot_list[[1]] <-ggplot(pc_with_metadata, aes_string("PC1", "PC2", label="orig.ident", color=grouping_by, shape=batch_col)) +
      geom_point(size=5, alpha=0.7) +
      geom_text_repel(size=2) +
      xlab(paste0("PC1: ",(round(as.vector(summary$importance[2,][1]),4)*100),"% variance")) +
      ylab(paste0("PC2: ",(round(as.vector(summary$importance[2,][2]),4)*100),"% variance")) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    plot_list[[2]] <-ggplot(pc_with_metadata, aes_string("PC3", "PC4", label="orig.ident", color=grouping_by, shape=batch_col)) +
      geom_point(size=5, alpha=0.7) +
      geom_text_repel(size=2) +
      xlab(paste0("PC1: ",(round(as.vector(summary$importance[2,][3]),4)*100),"% variance")) +
      ylab(paste0("PC2: ",(round(as.vector(summary$importance[2,][4]),4)*100),"% variance")) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    plot_list[[3]] <-ggplot(pc_with_metadata, aes_string("PC5", "PC6", label="orig.ident", color=grouping_by, shape=batch_col)) +
      geom_point(size=5, alpha=0.7) +
      geom_text_repel(size=2) +
      xlab(paste0("PC1: ",(round(as.vector(summary$importance[2,][5]),4)*100),"% variance")) +
      ylab(paste0("PC2: ",(round(as.vector(summary$importance[2,][6]),4)*100),"% variance")) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    
    pdf(paste(savedir4,"PC1_to_PC6_batch_effect_not_removed_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".pdf",sep = ""),
        width = 9, height = 8)
    print(plot_list)
    dev.off()
    
    saveRDS(pc,paste(savedir4,"pc_prcomp_batch_effect_not_removed.RDS",sep = ""))
    
    assay(rld) <- limma::removeBatchEffect(assay(rld),
                                           batch=rld$Run)
    
    # PCA_batch <- DESeq2::plotPCA(object = rld,
    #                              intgroup = colnames(Treg_metadata_2),
    #                              returnData=TRUE,ntop = 5000)
    # 
    # percentVar <- round(100 * attr(PCA_batch, "percentVar"))
    # 
    # PC1 <-ggplot(PCA_batch, aes_string("PC1", "PC2", label="name", color=grouping_by, shape=batch_col)) +
    #   geom_point(size=5, alpha=0.7) +
    #   geom_text_repel(size=2) +
    #   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    #   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    #   theme(plot.title = element_text(hjust = 0.5),
    #         panel.background = element_rect(fill = 'white', colour = 'black'),
    #         panel.grid.minor = element_line(colour = "grey"),
    #         panel.grid.major = element_line(colour = "grey"))
    # 
    # pdf(paste(savedir4,"PCA_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".pdf",sep = ""),
    #     width = 9, height = 8)
    # print(PC1)
    # dev.off()
    
    # write.table(PCA_batch, paste(savedir4,"PCA_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".txt",sep = ""),
    #             quote = F, row.names = T, col.names = T, sep = "\t")
    
    
    #### Further we are interested in all the PCs PC1 to all
    pc <- prcomp(t(assay(rld)))
    # attributes(pc)
    loading <- pc$rotation 
    summary <- summary(pc)
    
    rownames(Treg_metadata_2) <- Treg_metadata_2$orig.ident
    pc_with_metadata <- merge(pc$x,Treg_metadata_2,by = 'row.names', all = TRUE)
    
    write.table(pc_with_metadata, paste(savedir4,"PCA_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,"_all_PCs.txt",sep = ""),
                quote = F, row.names = T, col.names = T, sep = "\t")
    
    rm(plot_list)
    plot_list <- list()
    
    plot_list[[1]] <-ggplot(pc_with_metadata, aes_string("PC1", "PC2", label="orig.ident", color=grouping_by, shape=batch_col)) +
      geom_point(size=5, alpha=0.7) +
      geom_text_repel(size=2) +
      xlab(paste0("PC1: ",(round(as.vector(summary$importance[2,][1]),4)*100),"% variance")) +
      ylab(paste0("PC2: ",(round(as.vector(summary$importance[2,][2]),4)*100),"% variance")) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    plot_list[[2]] <-ggplot(pc_with_metadata, aes_string("PC3", "PC4", label="orig.ident", color=grouping_by, shape=batch_col)) +
      geom_point(size=5, alpha=0.7) +
      geom_text_repel(size=2) +
      xlab(paste0("PC1: ",(round(as.vector(summary$importance[2,][3]),4)*100),"% variance")) +
      ylab(paste0("PC2: ",(round(as.vector(summary$importance[2,][4]),4)*100),"% variance")) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    plot_list[[3]] <-ggplot(pc_with_metadata, aes_string("PC5", "PC6", label="orig.ident", color=grouping_by, shape=batch_col)) +
      geom_point(size=5, alpha=0.7) +
      geom_text_repel(size=2) +
      xlab(paste0("PC1: ",(round(as.vector(summary$importance[2,][5]),4)*100),"% variance")) +
      ylab(paste0("PC2: ",(round(as.vector(summary$importance[2,][6]),4)*100),"% variance")) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    
    pdf(paste(savedir4,"PC1_to_PC6_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".pdf",sep = ""),
        width = 9, height = 8)
    print(plot_list)
    dev.off()
    
    saveRDS(pc,paste(savedir4,"pc_prcomp_batch_effect_removed.RDS",sep = ""))
    
    
    message("Please find the PCA at this location ",savedir4)
    # return(Treg_dds_2)
    return(rld)
  } else {
    message(paste("cluster ", cluster))
    dir.create(savedir,showWarnings = FALSE)
    cluster_name <- paste(cluster,collapse = "_")
    savedir2 <- paste(savedir,"pseudobulk/",sep = "")
    dir.create(savedir2, showWarnings = FALSE)
    remove_samples_name <- paste(remove_samples,collapse="_and_")
    # savedir3 <- paste(savedir2,"clus_",cluster_name,"_removed_",remove_samples_name,"_",group1,"_vs_",group2,"/",sep = "")
    savedir3 <- paste(savedir2,"clus_",cluster_name,"_",group1,"_vs_",group2,"/",sep = "")
    dir.create(savedir3,showWarnings = FALSE)
    cluster_num <- paste("^",cluster,"$",sep = "")
    if (cluster_group == cluster_group) {
      obj@meta.data$merge_cluster <- as.integer(gsub(paste(cluster_num, collapse="|"),100,obj@meta.data[,cluster_group]))
    } else{
      obj@meta.data$merge_cluster <- gsub(paste(cluster_num, collapse="|"),100,obj@meta.data[,cluster_group])
    }
    
    Idents(obj) <- obj@meta.data$merge_cluster
    obj_subset <- subset(obj,idents = c(100))
    obj_subset@meta.data[,paste(sample_col,grouping_by,sep = "_")] <- paste(obj_subset@meta.data[,sample_col], obj_subset@meta.data[,grouping_by], sep="_")
    
    DefaultAssay(obj_subset) <- "RNA"
    cell_samples <- table(obj_subset@meta.data[,paste(sample_col,grouping_by,sep = "_")]) %>% as.data.frame()
    sample_low <- cell_samples[cell_samples$Freq < cell_freq,1]
    sample_low <- gsub("_","_",sample_low)
    sample_low_2 <- c(sample_low,remove_samples)
    message(paste("Removing this sample",sample_low_2,"\n",sep = " "))
    Treg_cts <- AggregateExpression(obj_subset, group.by = paste(sample_col,grouping_by,sep = "_"),
                                    assays = 'RNA', slot = "counts", return.seurat = FALSE)
    Treg_cts_2 <- as.data.frame(Treg_cts$RNA)
    if(length(sample_low_2) != 0){
      index <- grep(paste(sample_low_2,collapse="|"),colnames(Treg_cts_2), invert=TRUE)
      Treg_cts_3 <- Treg_cts_2[,index]
    } else {
      Treg_cts_3 <- Treg_cts_2
    }
    
    Treg_metadata <- data.frame(samples = colnames(Treg_cts_3))
    # Function to split and convert to dataframe
    split_to_dataframe <- function(x) {
      split_elements <- strsplit(x, splitting)[[1]]
      data.frame(t(split_elements))
    }
    
    # Convert split elements to dataframe
    Treg_metadata_2 <- do.call(rbind, lapply(Treg_metadata$samples, split_to_dataframe))
    colnames(Treg_metadata_2) <- c(column_name_split,grouping_by)
    Treg_metadata_2$orig.ident <- Treg_metadata$samples
    # Treg_metadata_2$Run <- Treg_metadata_2$orig.ident
    
    obj@meta.data[,paste(sample_col,batch_col,sep = "_")] <- paste(obj@meta.data[,sample_col],obj@meta.data[,batch_col],sep = "_")
    batch_group <- unique(obj@meta.data[,batch_col])
    # for (i in 1:length(batch_group)) {
    #   donor_change <- unique(obj@meta.data[,paste(sample_col,batch_col,sep = "_")]) %>% grep(batch_group[i],.,value=TRUE) %>% 
    #     gsub(paste(batch_group[i],sep = ""),"",.) %>% gsub("_","-",.) %>% paste(.,collapse=".*.|")
    #   donor_change <- paste(donor_change,".*.",sep="")
    #   Treg_metadata_2$Run <- gsub(donor_change,batch_group[i],Treg_metadata_2$Run)
    # }
    # Treg_metadata_2$Run <- gsub("-.*.","",Treg_metadata_2$Run)
    
    stopifnot(all(Treg_metadata_2$orig.ident == colnames(Treg_cts_3))) # if TRUE move forward
    # Treg_metadata$samplename <- gsub("gE_|DMSO_|_DMSO|_gE|B22015","",Treg_metadata$samples)
    ## Run will be treated as a covariate in the regression model,
    
    design_AJ <- as.formula(paste("~",paste(colnames(Treg_metadata_2[,c(batch_col,grouping_by)]),collapse = "+")))
    Treg_dds <- DESeqDataSetFromMatrix(countData = Treg_cts_3,
                                       colData = Treg_metadata_2,
                                       design = design_AJ)
    
    keep = filterByExpr(Treg_dds, group=colData(Treg_dds)[,grouping_by], min.count=gene_min_counts)
    table(keep)
    Treg_dds_2 <- Treg_dds[keep,]
    
    rld<-rlog(Treg_dds_2)
    PCA <- DESeq2::plotPCA(object = rld,
                           intgroup = colnames(Treg_metadata_2),
                           returnData=TRUE,ntop = 5000)
    percentVar <- round(100 * attr(PCA, "percentVar"))
    
    savedir4 <- paste(savedir3,"PCA/",sep = "")
    dir.create(savedir4, showWarnings = FALSE)
    
    PC1 <-ggplot(PCA, aes_string("PC1", "PC2", label="name", color=batch_col, shape=grouping_by)) +
      geom_point(size=5, alpha=0.7) +
      geom_text_repel(size=2) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      ggtitle("PCA") +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    pdf(paste(savedir4,"PCA_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".pdf",sep = ""),
        width = 9, height = 8)
    print(PC1)
    dev.off()
    
    assay(rld) <- limma::removeBatchEffect(assay(rld),
                                           batch=rld$Run)
    
    PCA_batch <- DESeq2::plotPCA(object = rld,
                                 intgroup = colnames(Treg_metadata_2),
                                 returnData=TRUE,ntop = 5000)
    
    percentVar <- round(100 * attr(PCA_batch, "percentVar"))
    
    PC1 <-ggplot(PCA_batch, aes_string("PC1", "PC2", label="name", color=grouping_by, shape=batch_col)) +
      geom_point(size=5, alpha=0.7) +
      geom_text_repel(size=2) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    pdf(paste(savedir4,"PCA_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".pdf",sep = ""),
        width = 9, height = 8)
    print(PC1)
    dev.off()
    
    write.table(PCA_batch, paste(savedir4,"PCA_removed_batch_effect_cell_gt_",cell_freq,"_cluster_",group1,"_vs_",group2,".txt",sep = ""),
                quote = F, row.names = T, col.names = T, sep = "\t")
    
    message("Please find the PCA at this location ",savedir4)
    return(Treg_dds_2)
  }
}
