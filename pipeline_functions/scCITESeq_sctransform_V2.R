#### Args
# Obj = RNA object created after the removing the doublets
# saveDir = Directory to save the output files
# ngenes = Total number of genes to be considered as variable genes and also the integration of different batches
# regress = To regressed out the genes like cell cycle genes not to be considered for UMAP
# Assay = "RNA"
# process = Process for integration and normalization
# objname = While saving to identify which sample or object is saved
# split_by = In meta data of seurat object which column need to be consider for integration
# reference = Which samples from the split_by need to be considered as reference onto which all other batches which be integrated
# dims = Number of PCA dims for making the UMAP
# sample_tree = Order for the integration
# k_weight = number of cells to be used for integration

# k_weight on reintegration if the number of cell in the following is low than we have to adjust the k_weight so that based on the those cell the nearest neighbour will be considered
sctransform_V2_integration <- function(obj,saveDir,ngenes=4000,regress,Assay,process,objname, 
                                       split_by="orig.ident", reference = NULL, dims=30, sample_tree=NULL, k_weight=100){
  split_obj <- SplitObject(obj, split.by = split_by)
  samplename <- unique(obj@meta.data[,split_by])
  ## providing more memory
  options(future.globals.maxSize = 4000 * 1024^2)
  
  message("\n Performing sctransform and initiated Integration \n")
  # In the SCTransform workflow, we perform a better version of variance stabilization, so we do not scale in this case.
  # We fix the slope parameter of the GLM to ln(10) with log10(total UMI) used as the predictor as proposed by Lause et al.
  # We utilize an improved parameter estimation procedure that alleviates uncertainty and bias that result from fitting GLM models for very lowly expressed genes.
  # We place a lower bound on gene-level standard deviation when calculating Pearson residuals. This prevents genes with extremely low expression (only 1-2 detected UMIs) from having a high pearson residual.
  # https://github.com/satijalab/seurat/issues/3003
  for (i in 1:length(split_obj)) {
    split_obj[[i]] <- SCTransform(split_obj[[i]],
                                  vst.flavor = "v2", # using version 2
                                  vars.to.regress = regress,
                                  assay = "RNA",
                                  variable.features.n = ngenes,
                                  ncells = nrow(split_obj[[i]]@meta.data))
    file <- unique(split_obj[[i]]@meta.data$orig.ident)
    message(paste("\n################## Normalization",file[i],"######################\n"))
    sample <- split_obj[[i]]
    sample <- RunPCA(sample, npcs=50) # This is just to check the whether cellcycle or mito effect were regressed out. So, kept the PCs = 50

    dir.create(paste(saveDir,"cell_cycle_phase/",sep = ""), showWarnings = FALSE)
    pdf(paste(saveDir,"cell_cycle_phase/",objname,"_",samplename[i],"_",Assay,"_",process,"_splitted.pdf",sep = ""),width = 8, height = 8)
    print(DimPlot(sample,reduction = "pca",group.by= "Phase",split.by = "Phase"))
    dev.off()

    pdf(paste(saveDir,"cell_cycle_phase/",objname,"_",samplename[i],"_",Assay,"_",process,".pdf",sep = ""),width = 8, height = 6)
    print(DimPlot(sample,reduction = "pca",group.by= "Phase"))
    dev.off()

    dir.create(paste(saveDir,"mitoeffect/",sep = ""), showWarnings = FALSE)
    pdf(paste(saveDir,"mitoeffect/",objname,"_",samplename[i],"_",Assay,"_",process,"_splitted.pdf",sep = ""),width = 8, height = 8)
    print(DimPlot(sample,reduction = "pca",group.by= "mitoFr",split.by = "mitoFr"))
    dev.off()

    pdf(paste(saveDir,"mitoeffect/",objname,"_",samplename[i],"_",Assay,"_",process,".pdf",sep = ""),width = 8, height = 8)
    print(DimPlot(sample,reduction = "pca",group.by= "mitoFr"))
    dev.off()
  }
  
  print(split_obj)
  
  message(paste("Performing Integration of",objname))
  integ_features <- SelectIntegrationFeatures(object.list = split_obj, nfeatures = ngenes) 
  
  # Prepare the SCT list object for integration
  split_obj <- PrepSCTIntegration(object.list = split_obj, anchor.features = integ_features)
  
  ## perform CCA, find the best buddies or anchors and filter incorrect anchors.
  # Find best buddies - can take a while to run
  integ_anchors <- FindIntegrationAnchors(object.list = split_obj, 
                                          normalization.method = "SCT", 
                                          anchor.features = integ_features,
                                          reference = reference)
  
  # Integrate across conditions
  # you can think of k.weight as a smoothing parameter. A value of 100 means each cell will be transformed by a weighted combination 
  # of the nearest 100 anchors. The assumption is that there is some amount of randomness/sparsity in the data, making it desirable
  # to combine anchors in the same neighborhood. If you set k.weight very low, you will have less smoothing and are basically 
  # assuming the information in each of your cells is more reliable/complete. This might be the case for your data, but for most
  # single-cell RNA assays there is a good deal of sparsity.
  # https://github.com/satijalab/seurat/issues/3930
  integrated <- IntegrateData(anchorset = integ_anchors,
                              normalization.method = "SCT",
                              dims = 1:dims, 
                              sample.tree = sample_tree,
                              k.weight = k_weight)
  
  # PCA or UMAP, PCA will take into consideration only 2 PCs, while UMAP will from the top number of PCs.
  message(paste("\nRunning PCA\n"))
  integrated <- RunPCA(object = integrated)
  
  message(paste("\nMaking Elbow Plot\n"))
  dir.create(paste(saveDir,"elbow_plots/",sep = ""), showWarnings = FALSE)
  pdf(paste(saveDir,"elbow_plots/",objname,"_",process,"_",Assay,".pdf",sep = ""),width = 8, height = 8)
  print(ElbowPlot(integrated, ndims = 50))
  dev.off()
  
  # message(paste("Saving the ",objname," integrated Object"))
  # saveRDS(integrated,paste(saveDir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = ""))
  message(paste("Please check the Elbow Plots for further Integration. Please go to this location ",saveDir,"elbow_plots/",sep = ""))
  message(paste("Please check the regress out Cell cycle or mitoeffect for further Integration. Please go to this location for cellcycle ",saveDir,"cell_cycle_phase/\n and for mitoeffect at ",saveDir,"mitoeffect/",sep = ""))
  message(paste("After normalization,we will perform the sample integration using the function at this location:\n/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R"))
  return(integrated)
}

# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")
# combined <- readRDS(paste(savedir,"saveRDS_obj/combined_merged_RNA.RDS",sep = ""))
# ## For each sample in RNA assay they had to be regress out for the Mitochondria, however the cell cycle phase looks good
# ## Remember we are not integrating we are merging the samples
# ## By looking at the graph, we add the mitoratio or cellcycle (G2M.Score,S.Score)
# sctransform_V2_integration(combined, savedir, regress = c("mitoRatio"), objname = "combined_integrated",
#                            process = "sctransform", Assay = "RNA")
