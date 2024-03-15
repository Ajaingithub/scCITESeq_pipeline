# This function is used for identify the markers, generating featureplots, heatmaps as well as perform Gene Ontology for the differential genes
# obj: Seurat_object
# saveDir: saving directory
# ident_1: Cluster to be used for identify markers and GO
# min_pct: minimum percentage of the features expressing within the clusters while calculating the markers
# objname: Name of the samples
# process: During saving what name to use
# Assay: Which seurat assay to use to calcuate the markers

cluster_markers_GSEA = function(obj, saveDir, ident_1, min_pct = 0.25, objname,process,Assay){
  library(Seurat)
  message(paste("Calculating the RNA and ADT markers for the specific cluster",ident_1, sep = ""))
  CM_EM_marker_ADT <- FindMarkers(obj, assay = "ADT", ident.1 = ident_1, min.pct = min_pct)
  CM_EM_marker_RNA <- FindMarkers(obj, assay = "RNA", ident.1 = ident_1, min.pct = min_pct)
  
  message(paste("Making the heatmap for RNA and FeaturePlots of the ADT"))
  ident_name = paste(ident_1,collapse = "_")
  dir.create(paste(saveDir,"featureplot/",sep = ""),showWarnings = F)
  DefaultAssay(obj) = "ADT"
  pdf(paste(saveDir,"featureplot/",objname,"_",process,"_",Assay,"_",ident_name,".pdf",sep = ""), width = 12, height = 12)
  print(FeaturePlot(obj,
                    reduction = "wnn.umap",
                    features = rownames(CM_EM_marker_ADT),
                    order = TRUE,
                    label = TRUE))
  dev.off()
  
  dir.create(paste(saveDir,"heatmap",sep = ""),showWarnings = F)
  
  RNA_features = rownames(CM_EM_marker_RNA) %>% head(100)
  pdf(paste(savedir,"heatmap/",objname,"_",process,"_",Assay,"_",ident_name,".pdf",sep = ""), width = 12, height = 15)
  print(DoHeatmap(CD4, assay = "RNA", features = RNA_features, group.by = "ident"))
  dev.off()
  
  message("Making the gene list for the Gene Set Enrichment Analysis we cannot process it in this session Please run on the local system")
  message("\n Use the code from this function that has been commented. The location is /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/heatmap_and_pathway_analysis.R")
  
  message("saving the files to be used on the local machines")
  dir.create(paste(saveDir,"GSEA",sep = ""),showWarnings = FALSE)
  write.csv(CM_EM_marker_RNA,paste(saveDir,"GSEA/",objname,"_",process,"_",Assay,"_",ident_name,".csv",sep = ""),quote = FALSE)
  
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(org.Hs.eg.db)

  genename_symbol = rownames(CM_EM_marker_RNA)
  genename_entrez <- bitr(genename_symbol, fromType = "SYMBOL",
                          toType = c("ENTREZID"), OrgDb = "org.Hs.eg.db")
  CM_EM_marker_RNA_entrezid = CM_EM_marker_RNA[match(genename_entrez$SYMBOL,genename_symbol),]
  original_gene_list = CM_EM_marker_RNA_entrezid$avg_log2FC
  
  
  names(original_gene_list) <- genename_entrez$ENTREZID
  gene_list<-na.omit(original_gene_list)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  gse <- gseGO(geneList=gene_list,
               ont ="BP",
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               verbose = TRUE,
               OrgDb = org.Hs.eg.db,
               pAdjustMethod = "none")
  
  require(DOSE)
  dir.create(paste(savedir,"GSEA",sep = ""),showWarnings = F)
  pdf(paste(saveDir,"GSEA/gsea.pdf",sep = ""),showWarnings = FALSE)
  dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  dev.off()  
  
}
