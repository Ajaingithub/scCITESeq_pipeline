# scCITESeq
 
The functions can be used to perform the downstream analysis of the scCITESeq data.
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
   
3. **scCITEseq_Doublet_Finder** - This function will remove the doublets from the object.

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
   
4. **scCITESeq_sctransform_V2** - The object is normalized using scTransform V2 and integreate using CCA.

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
                         saveDir = "./object/",
                         dims=30,
                         RNA_features = c("CD4","CD8A"), 
                         Assay="RNA",
                         process="integration",
                         objname = "BB22003",
                         ncol=2,
                         ndims = 50)
   
9. **scCITESeq_ADT_merging** - Performing the ADT merging.
10. **cluster_UMAP_QC** - Performing clustering and checking the quality for each cluster.
11. **scCITEseq_modality_integration** -  Modality RNA and ADT integration.
