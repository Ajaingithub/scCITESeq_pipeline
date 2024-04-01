# scCITESeq
 
The functions can be used to perform the downstream analysis of the scCITESeq data. I have used all these function in **pipeline_functions/scCITESeq_Execution.R** 

1. **scCITESeq_QC function** - This function will perform the quality check scRNA and scADT dataset before and after fitering the cells based on the cell number, number of genes, mitochondir percent, and ADT UMI.

   #### To Run the Code
       source("./pipeline_functions/CR7_scCITESeq_QC.R")

       CR_outdir ="./post_processing/genotype_demux/Exp01_geno/outs/per_sample_outs/Exp01_geno"

       scCITE_QC(samplepath = CR_outdir,
                 samplename = "BB22003",
                 saveDir = "./output/",
                 min_cells = 3,
                 min_genes = 200,
                 max_genes = 5000,
                 mitopercent = 10,
                 ADT_UMI = 100000)

   This will create the filtered RNA and ADT (protein) object

![Screenshot 2024-04-01 at 5 20 15 PM](https://github.com/Ajaingithub/scCITESeq_pipeline/assets/37553954/354b791f-64ea-4927-8d76-843fd646f628)

   
2. **scCITEseq_Doublet_Finder** - This function will remove the doublets from the object.

   #### To Run the Code
          source("./pipeline_functions/scCITEseq_Doublet_Finder.R")

          doublet_scCITEseq(Obj= RNA_obj,
                            dims = 30,
                            res = 0.8,
                            saveDir = "./output/",
                            Obj2 = ADT_obj,
                            samplename = "BB22003",
                            process = "Doublet_finder",
                            Assay = "RNA")

   This will filter out doublets from both the RNA and ADT
   
   ![Screenshot 2024-04-01 at 5 23 52 PM](https://github.com/Ajaingithub/scCITESeq_pipeline/assets/37553954/597e8948-3f76-43f5-8912-57bf256b8907)

4. **scCITESeq_sctransform_V2** - The object is normalized using scTransform V2 and samples integreated using CCA.

   #### To Run the Code
          source("./pipeline_functions/scCITESeq_sctransform_V2.R")

          sctransform_V2_integration(obj = RNA_obj,
                                     saveDir = "./output/" ,
                                     ngenes=4000,
                                     regress = "mito_score",
                                     Assay = "RNA",
                                     process ="norm",
                                     objname = "BB22003",
                                     split_by="orig.ident",
                                     reference = NULL,
                                     dims=30,
                                     sample_tree=NULL,
                                     k_weight=100)

5. **scCITESeq_RNA_intergation** - The RNA / ADT object integration function.

   ##### To Run the Code
         source("pipeline_functions/scCITESeq_RNA_intergation.R")

         RNA_integration(obj = RNA_obj/ADT_obj,
                         saveDir = "./output/",
                         dims=30,
                         RNA_features = c("CD4","CD8A"), 
                         Assay="RNA",
                         process="integration",
                         objname = "BB22003",
                         ncol=2,
                         ndims = 50)
   
6. **scCITESeq_ADT_merging** - Performing the ADT (protein) Sample Integration using CCA.
   #### To Run the Code
             source("pipeline_functions/scCITESeq_ADT_merging.R")
   
             ADT_merging(ADT_obj = ADT_obj,
                         saveDir= "./output/",
                         dims=10,
                         numfeatures=NULL,
                         Assay="ADT",
                         process="ADT_integration",
                         objname="BB22003",
                         split_by="orig.ident", 
                         reference=NULL,
                         sample_tree = NULL,
                         k_weight=100)

   ![Screenshot 2024-04-01 at 5 25 28 PM](https://github.com/Ajaingithub/scCITESeq_pipeline/assets/37553954/3a2ad526-fd3d-400e-b7a1-caba184ac866)

7. **cluster_UMAP_QC** - Performing clustering and checking the quality for each cluster.
   #### To Run the Code
          source("pipeline_functions/cluster_UMAP_QC.R")

          cluster_UMAP_and_QC(obj_path=RNA_obj/ADT_obj,
                              dims = 30,
                              res = 0.8,
                              saveDir = "./output/",
                              Assay = "RNA",
                              QC_features = c("nCount_RNA","nFeature_RNA"),
                              objname = "BB22003",
                              process = "UMAP_QC",
                              col_sel = c("Age","Run","orig.ident","gender"))
    
    ![Screenshot 2024-04-01 at 5 26 28 PM](https://github.com/Ajaingithub/scCITESeq_pipeline/assets/37553954/5daba60c-9168-4165-bc80-027b0bfe8436)

8. **scCITEseq_modality_integration** -  Modality RNA and ADT integration.
   #### To Run the Code
              source("pipeline_functions/scCITEseq_modality_integration.R")

              modality_integration(RNA_obj = RNA_obj,
                                   ADT_obj = ADT_obj,
                                   RNA_dims = 30,
                                   ADT_dims = 10,
                                   saveDir = savedir,
                                   Assay = "integrated",
                                   process = "RNA_ADT",
                                   objname= "BB22003")

These 7 step will perform the batch as well as modality (RNA and protein) integration for the scCITESeq experiment

![Screenshot 2024-04-01 at 5 27 05 PM](https://github.com/Ajaingithub/scCITESeq_pipeline/assets/37553954/f9437690-1c86-4ce4-becd-150db6005f77)


### Downstream Analysis

**Pseudobulk PCA and Differential**:

           source("pipeline_functions/Psuedobulk_PCA.R")
           
           CD8_subset_dds <- COVID_pseudobulk_within_cluster_AJ(obj = CD8_mem, 
                                                               savedir = savedir, 
                                                               group1 = group1, 
                                                               group2 = group2,
                                                               grouping_by = "Age", 
                                                               cluster = "all", 
                                                               cell_freq = 20, 
                                                               remove_samples = remove_samples,
                                                               cluster_group = "seurat_clusters",
                                                               sample_col = "orig.ident", 
                                                               batch_col = "Run",
                                                               gene_min_counts =  5, 
                                                               column_name_split = c("sampleid","virus","age","age_number","gender","Run"))
                                                               


         source("pipeline_functions/Pseudobulk_differential.R")
         
         design0 <- model.matrix(~ 0 + Age + Run, data = colData(CD8_subset_dds))
         colnames(design0) <- c("O","Y",paste("Run",1:(ncol(design0)-2),sep = ""))
         cm <- makeContrasts(O_VS_Y = O-Y,levels = design0)
         
         desl_clus <- LimmaEdgeR_differential(dds = CD8_subset_dds,
                                              design0 = design0,
                                              cm = cm, 
                                              savedir = savedir2,
                                              logfc = 0.5,
                                              p_value_adj = 0.05)

![Screenshot 2024-04-01 at 5 31 25 PM](https://github.com/Ajaingithub/scCITESeq_pipeline/assets/37553954/591be590-97f7-42a7-969e-d971ef2061f5)

![Screenshot 2024-04-01 at 5 32 06 PM](https://github.com/Ajaingithub/scCITESeq_pipeline/assets/37553954/220691b9-d267-4b84-be68-10bcd8a69cdb)

                                                               




