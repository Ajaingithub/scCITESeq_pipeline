# scCITESeq
 
The functions can be used to perform the downstream analysis of the scCITESeq data. I have used all these function in **pipeline_functions/scCITESeq_Execution.R** 

1. **scCITESeq_QC function** - This function will perform the quality check scRNA and scADT dataset before and after fitering the cells based on the cell number, number of genes, mitochondir percent, and ADT UMI.

   #### To Run the Code
       source("./pipeline_functions/CR7_scCITESeq_QC.R")

       CR_outdir ="./post_processing/genotype_demux/Exp01_geno/outs/per_sample_outs/Exp01_geno

       scCITE_QC(samplepath = CR_outdir,
                 samplename = "BB22003",
                 saveDir = "./output/",
                 min_cells = 3,
                 min_genes = 200,
                 max_genes = 5000,
                 mitopercent = 10,
                 ADT_UMI = 100000)

   This will create the filtered RNA and ADT (protein) object
   
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
   
3. **scCITESeq_sctransform_V2** - The object is normalized using scTransform V2 and samples integreated using CCA.

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

4. **scCITESeq_RNA_intergation** - The RNA / ADT object integration function.

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
   
5. **scCITESeq_ADT_merging** - Performing the ADT (protein) Sample Integration using CCA.
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
   
6. **cluster_UMAP_QC** - Performing clustering and checking the quality for each cluster.
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
    
    
7. **scCITEseq_modality_integration** -  Modality RNA and ADT integration.
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

