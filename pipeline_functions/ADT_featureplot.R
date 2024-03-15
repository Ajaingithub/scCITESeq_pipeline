# This function will create the featureplots for all the protein in the seurat object both imputed and  non-imputed datasets
# Args:
# obj: seurat_obj
# savedir: save directory
# objname: name of the sample
# process: "featureplots"
# reduction: dimensionality reduction like pca, umap, wnn.umap
# x = PC1, TSNE1, UMAP_1"
# y = PC2, TSNE2, UMAP_2

ADT_featureplot <- function(obj, savedir, objname, process, 
                            reduction = "umap",x="UMAP_1", y="UMAP_2"){
  source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/single_cell/express_cell_front.R")
  DefaultAssay(obj) <- "MAGIC_ADT"
  p1  <- featureplot_front(obj, "CD4-protein")
  p2 <- featureplot_front(obj, "CD8a-protein")
  
  pdf(paste(savedir,"featureplot/",objname,"_MAGIC_ADT_",process,"_CD4_CD8a_front.pdf", sep = ""), width = 10, height = 6)
  print(p1|p2)
  dev.off()
  
  DefaultAssay(obj) <- "ADT"
  source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/single_cell/express_cell_front.R")
  protein_names <- rownames(obj@assays$ADT)
  rm(plot_list)
  plot_list <- list()
  for (i in 1:length(protein_names)) {
    p <- print(featureplot_front(obj, 
                                 protein_names[i],
                                 reduction = reduction,
                                 x=x, y=y, size=0.01))
    plot_list[[i]] <- p
  }
  
  require(grid)
  require(gridExtra)
  pdf(paste(savedir,"featureplot/",objname,"_ADT_",process,"_first20_protein.pdf", sep = ""), width = 25, height = 25)
  grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],
               plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],
               plot_list[[11]],plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],
               plot_list[[16]],plot_list[[17]],plot_list[[18]],plot_list[[19]],plot_list[[20]], ncol = 5)
  dev.off()
  
  pdf(paste(savedir,"featureplot/",objname,"_ADT_",process,"_last20_protein.pdf", sep = ""), width = 25, height = 25)
  grid.arrange(plot_list[[21]],plot_list[[22]],plot_list[[23]],plot_list[[24]],plot_list[[25]],
               plot_list[[26]],plot_list[[27]],plot_list[[28]],plot_list[[29]],plot_list[[30]],
               plot_list[[31]],plot_list[[32]],plot_list[[33]],plot_list[[34]],plot_list[[35]],
               plot_list[[36]],plot_list[[37]],plot_list[[38]],plot_list[[39]],plot_list[[40]], ncol = 5)
  dev.off()
  
  DefaultAssay(obj) <- "MAGIC_ADT"
  source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/single_cell/express_cell_front.R")
  protein_names <- rownames(obj@assays$ADT)
  rm(plot_list)
  plot_list <- list()
  for (i in 1:length(protein_names)) {
    p <- print(featureplot_front(obj, 
                                 protein_names[i],
                                 reduction = reduction,
                                 x=x, y=y, size=0.01))
    plot_list[[i]] <- p
  }
  
  require(grid)
  require(gridExtra)
  pdf(paste(savedir,"featureplot/",objname,"_MAGIC_ADT_",process,"_first20_protein.pdf", sep = ""), width = 25, height = 25)
  grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],
               plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],
               plot_list[[11]],plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],
               plot_list[[16]],plot_list[[17]],plot_list[[18]],plot_list[[19]],plot_list[[20]], ncol = 5)
  dev.off()
  
  pdf(paste(savedir,"featureplot/",objname,"_MAGIC_ADT_",process,"_last20_protein.pdf", sep = ""), width = 25, height = 25)
  grid.arrange(plot_list[[21]],plot_list[[22]],plot_list[[23]],plot_list[[24]],plot_list[[25]],
               plot_list[[26]],plot_list[[27]],plot_list[[28]],plot_list[[29]],plot_list[[30]],
               plot_list[[31]],plot_list[[32]],plot_list[[33]],plot_list[[34]],plot_list[[35]],
               plot_list[[36]],plot_list[[37]],plot_list[[38]],plot_list[[39]],plot_list[[40]], ncol = 5)
  dev.off()
  
  message(paste("Please check the ADT Feature plot both imputed and unimputed antibody in the UMAP, Please go to the ",savedir,"featureplot/",sep = ""))
}
