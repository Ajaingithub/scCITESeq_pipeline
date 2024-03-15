## Args
# ADT_obj = ADT_obj
# saveDir = save directory
# dims = PCA dims for ADT
# numfeatures = number of variable features 
# Assay = "ADT"
# process = "Integration"
# objname = name of the sample
# split_by = column which to be used for integration like batch
# reference = Which samples from the split_by need to be considered as reference onto which all other batches which be integrated
# dims = Number of PCA dims for making the UMAP
# sample_tree = Order for the integration
# k_weight = number of cells to be used for integration

# reference set you have to provide the index number that is present in your list
# reference_dataset <- which(names(pbmc.list) == "10x Chromium (v3)")
# sample_tree is a matrix that decide how to integrate the samples
# 1st we integrate the sample for each run then we integrate the both the batches for which we make the matrix
# how to make check the Mayo commands doc
# # k_weight on reintegration if the number of cell in the following is low than we have to adjust the k_weight so that based on the those cell the nearest neighbour will be considered
ADT_merging <- function(ADT_obj, saveDir, dims=10, numfeatures=NULL, Assay="ADT", process, objname, split_by="orig.ident", 
                        reference=NULL, sample_tree = NULL, k_weight=100){
   DefaultAssay(ADT_obj) <- "ADT"
  split_combined_ADT <- SplitObject(ADT_obj, split.by = split_by)
  
  message("Normalizing the object ",objname)
  # normalize and identify variable features for each dataset independently
  split_combined_ADT <- lapply(X = split_combined_ADT, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = 'CLR', margin = 2)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = numfeatures)
  })
  
  message("Scaling and Running PCA ",objname)
  features <- SelectIntegrationFeatures(object.list = split_combined_ADT,
                                        nfeatures = numfeatures)
  split_combined_ADT <- lapply(X = split_combined_ADT, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, npcs = numfeatures, approx=FALSE)
  })
  
  message("Integrating the samples ",objname)
  ## We then identify anchors using the FindIntegrationAnchors() function, which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData().
  
  immune.anchors <- FindIntegrationAnchors(object.list = split_combined_ADT,
                                           anchor.features = features,
                                           reference = reference,
                                           reduction = "rpca")
  
  # rpca is reciprocal pca represents a more conservative approach where cells in different biological states are less likely 
  # to ‘align’ after integration. We therefore,recommend RPCA during integrative analysis where: * A substantial fraction of 
  # cells in one dataset have no matching type in the other 
  ## How to make a integration sample.tree matrix
  # 5 into 3
  # 1 into 3 5
  # 2 into 6
  # 4 into 6 2
  # 3 5 1 into 6 2 4
  # https://stackoverflow.com/questions/58582881/creating-dendrograms-manually-how-to-fix-merge-matrix-has-invalid-contents
  # It shows the warning Hi, this is a warning not an error, just saying that you're running partial SVD but computing most (or all) of the singular values.
  # Integrate across conditions
  # you can think of k.weight as a smoothing parameter. A value of 100 means each cell will be transformed by a weighted combination 
  # of the nearest 100 anchors. The assumption is that there is some amount of randomness/sparsity in the data, making it desirable
  # to combine anchors in the same neighborhood. If you set k.weight very low, you will have less smoothing and are basically 
  # assuming the information in each of your cells is more reliable/complete. This might be the case for your data, but for most
  # single-cell RNA assays there is a good deal of sparsity.
  # https://github.com/satijalab/seurat/issues/3930
  ADT_integrated <- IntegrateData(anchorset = immune.anchors, 
                                  dims = 1:dims, 
                                  sample.tree = sample_tree,
                                  k.weight = k_weight)
  
  # the integration corrects normalized data, so there is no need to do another log-normalization after integration.
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(ADT_integrated) <- "integrated"
  
  message("Scaling and Running PCA on the ",objname," integrated data")
  # Run the standard workflow for visualization and clustering
  ADT_integrated <- ScaleData(ADT_integrated, verbose = FALSE)
  ADT_integrated <- RunPCA(ADT_integrated, npcs = numfeatures, verbose = TRUE, approx=FALSE)
  
  dir.create(paste(saveDir,"elbow_plots/",sep = ""))
  pdf(paste(saveDir,"elbow_plots/",objname,"_",process,"_",Assay,".pdf",sep = ""),width = 8, height = 8)
  print(ElbowPlot(ADT_integrated, ndims = numfeatures, reduction = "pca"))
  dev.off()
  
  # message("\nSaving ADT integrated Object\n")
  # saveRDS(ADT_integrated, paste(saveDir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = ""))
  message(paste("\n Check for elbow plot at this location ",saveDir,"elbow_plots/\n In order to FindCluster and making UMAP",sep = ""))
  return(ADT_integrated)
}

# ADT_obj <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "_ADT_singlet.RDS",full.names = TRUE)
# ADT_samplename <- gsub("_ADT_singlet.RDS","",basename(ADT_obj))
# 
# ## For each sample in ADT assay 
# ## Remember we are not integrating we are merging the samples
# ## This has to be modified everytime when merging different object based on the number of the object
# combined_ADT <- merge(readRDS(ADT_obj[1]), y = c(readRDS(ADT_obj[2]), readRDS(ADT_obj[3]),readRDS(ADT_obj[4])), 
#                       add.cell.ids = ADT_samplename, project = "ADT_combined")
# 
# ## By looking at the graph, we add the mitoratio or cellcycle (G2M.Score,S.Score)
# # numFeatures depends on how many antibody that has been used for ADT features
# ADT_merging(combined_ADT, savedir, dims = 10, numfeatures=40, Assay = "ADT", process = "ADT_merging", objname = "combined")
