# scCITESeq
 
The functions can be used to perform the downstream analysis of the scCITESeq data.
1. **scCITESeq_QC function ** - This function will perform the quality check scRNA and scADT dataset before and after fitering the cells based on the cell number, number of genes, mitochondir percent, and ADT UMI. 
2. **scCITEseq_Doublet_Finder** - This function will remove the doublets from the object.
3. **scCITESeq_sctransform_V2** - The object is normalized using scTransform V2. 
4. **scCITESeq_RNA_intergation** - The RNA / ADT object integration function.
5. **scCITESeq_ADT_merging** - Performing the ADT merging.
6. **cluster_UMAP_QC** - Performing clustering and checking the quality for each cluster.
7. **scCITEseq_modality_integration** - Finally modality integration.
