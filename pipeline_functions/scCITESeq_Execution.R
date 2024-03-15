#### This script is scCITEseq Ines Run01 Demultiplexing using Cell Ranger regressing unwanted Effect and integrating Modalities and samples keeping Feature 4000 #

#### Before Command ########
# https://satijalab.org/seurat/articles/visualization_vignette.html
# Download the fastq data from the Novogene on the screen. This will take 2-2.5 days. Since data is not demultiplexed. Try to download it in parallel but keep track everything is downloaded
# Check the bytes as well as MD5 for your files using md5sum filename with the MD5.txt they have provided 
# ls *.fastq.gz | awk '{print "md5sum "$0" >> checking_MD5.txt"}' > MD5_check.sh
# Merge the data for both RNA and ADT. all R1, R2, and I1. Since we are using Antibody tagged oligos, GEX. We will get
# 2 fastqs. 1. GEX, 2. CSP. Need to make a new folder with the fastq from CSP. Remove the first read from all R1, R2, and I1
# Read this https://kb.10xgenomics.com/hc/en-us/articles/4407386498957-I-used-antibody-tags-for-cell-surface-protein-capture-and-cell-hashing-with-Single-Cell-3-chemistry-How-can-I-use-Cell-Ranger-to-analyze-my-data-
# removing 1st four lines from the CSP from R1, R2, and I1. In order to have separate identity from the CSP for multiplex capture.
# zcat Exp01_CSP_combined_CKDL210029696-1a-SI_GA_G8_H7VK3DSX3_S1_L001_I1_001.fastq.gz | sed '1,4d' > Exp01_CSP_combined_del_4_CKDL210029696-1a-SI_GA_G8_H7VK3DSX3_S1_L001_I1_001.fastq
# ls *.fastq.gz | sed 's/_CSP_/\t/g' | awk '{print "zcat "$1"_CSP_"$2" | sed \x271,4d\x27 > /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run02_and_Run03/combined_raw_data/Exp02_and_03_CSP_multiplex_capture/"$1"_CSP_4_del"$2" &"}' > multiplex_fastq.sh
# please make change \x271 has put q make it '1 in all the files
# gzip Exp01_CSP_combined_del_4_CKDL210029696-1a-SI_GA_G8_H7VK3DSX3_S1_L001_I1_001.fastq
# ls *.fastq | awk '{print "gzip "$0" &"}' > zipping.sh
# Now you have three folder 1. GEX, 2. CSP, 3. multiplexing capture
# Need to make CMO reference for demultiplexing your sample
# id,name,read,pattern,sequence,feature_type
# HTO1,HTO1,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,GTCAACTCTTTAGCG,Multiplexing Capture
# HTO2,HTO2,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TGATGGCCTATTGGG,Multiplexing Capture
# HTO3,HTO3,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,AAGTATCGTTTCGCA,Multiplexing Capture
# HTO4,HTO4,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,GGTTGCCAGATGTCA,Multiplexing Capture

# Also need to make a protein reference keep the negative control also to find if any cluster has high negative control in that
# id,name,read,pattern,sequence,feature_type
# CD4_protein_Total_SeqB,CD4_protein,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TGTTCCCGCTCAACT,Antibody Capture
# CD8a_protein_Total_SeqB,CD8a_protein,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,GCTGCGCTTTCCATT,Antibody Capture
# CD45RA_protein_Total_SeqB,CD45RA_protein,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,TCAATCCTTCCGCTT,Antibody Capture
# CD45RO_protein_Total_SeqB,CD45RO_protein,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,CTCCGAATCATGTTG,Antibody Capture
# CD197_CCR7_protein_Total_SeqB,CD197_CCR7_protein,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,AGTTCAGTCAACCGA,Antibody Capture
# CD95_Fas_protein_Total_SeqB,CD95_Fas_protein,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,CCAGCTCATTAGAGC,Antibody Capture

# Make the multiconfig.csv file, This include 
# [gene-expression]
# reference,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/refdata-gex-GRCh38-2020-A/
# expect-cells,12000
# cmo-set,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/usftp21.novogene.com/raw_data/combined_raw_data/required_files/Ines_CMO_reference.csv
# 
# [feature]
# reference,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/usftp21.novogene.com/raw_data/combined_raw_data/required_files/Ines_CITESeq_feature_ref_Total_seqB_2_with_HTO_reference.csv
# 
# [libraries]
# fastq_id,fastqs,lanes,feature_types,subsample_rate
# Exp01_GEX_combined_CKDL210029696-1a-SI_GA_G8_H7VK3DSX3,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/usftp21.novogene.com/raw_data/combined_raw_data/Exp01_GEX_combined,any,Gene Expression,
# Exp01_CSP_combined_CKDL210029696-1a-SI_GA_G8_H7VK3DSX3,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/usftp21.novogene.com/raw_data/combined_raw_data/Exp01_CSP_combined,any,Antibody Capture,
# Exp01_CSP_combined_del_4_CKDL210029696-1a-SI_GA_G8_H7VK3DSX3,/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/usftp21.novogene.com/raw_data/combined_raw_data/Exp01_CSP_multiplexing_capture,any,Multiplexing Capture
# 
# [samples]
# sample_id,cmo_ids,description
# Sample1,HTO1,Sample1_M_27
# Sample2,HTO2,Sampl2_F_24
# Sample3,HTO3,Sample3_F_23
# Sample4,HTO4,Sample4_M_25

# source ~/.bashrc
# cellranger multi --id=multiconfigHTO --csv=/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/usftp21.novogene.com/raw_data/combined_raw_data/required_files/Ines_multi_config_scCITEseq.csv
# qsub -l h_vmem=10G -pe threaded 32 -N multiHTO -q 1-day -o multiHTO.log -e multiHTO.err -m ae -M jain.abhinav@mayo.edu -cwd multiconfig_HTO.sh
# CITESeq count does not work as it will fetch out the data from Undetermined fastqs https://www.biostars.org/p/429121/

# To check the total number of line that needs to be divided by 4 afterwards please go throught the command before running it
# mkdir read_numbers/
# ls Exp0[0-9]_*/*.fastq.gz | sed 's/.fastq.gz//g' | awk '{print "zcat "$0".fastq.gz | wc -l > read_numbers/"$0" &"}' | sed 's:/Exp0[0-9]_GEX_4::3' | sed 's:/Exp0[0-9]_CSP_[0-9]::3' > read_number.sh
# ls Exp0[0-9]_*/*.fastq.gz | sed 's/.fastq.gz//g' | awk '{print "zcat "$0".fastq.gz | wc -l > read_numbers/"$0" &"}' | sed 's:/Exp0[0-9]_GEX::3' | sed 's:/Exp0[0-9]_CSP::3' > read_numbers.sh
# grep [0-9] * | sed 's/:/\t/g' | sed 's/_CKDL/\t/g' | sed 's/_L0[0-9][0-9]_/\t/g'  | sed 's/_001//g' | awk '{print $1"\t"$3"\t"$4/4}' > final_reads.txt
# then we have to cast to get in proper shape
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/read_number_casting.R")
read_number_casting(tablename = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/Run02/usftp21.novogene.com/raw_data/read_numbers/final_reads.txt")
# Integration of multiple samples together with different modalities
## https://github.com/satijalab/seurat/issues/5346  
## https://github.com/satijalab/seurat/issues/3890
## https://github.com/satijalab/seurat/issues/2207

##### Starting with the QC ####
library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_QC.R")

## Please keep the "/" at the end of the directories
savedir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/usftp21.novogene.com/raw_data/combined_raw_data/Analysis/preprocessing_with_CMO_library/multiconfigHTO_in_ref/outs/per_sample_outs/regressing_scTransform_doublet/ADT_clustering/Pipeline/executing_pipeline_2/"
dir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/scCITESeq/Ines/usftp21.novogene.com/raw_data/combined_raw_data/Analysis/preprocessing_with_CMO_library/multiconfigHTO_in_ref/outs/per_sample_outs/"
files <- list.files(dir,pattern = "^Sample[1-4]$")
table <- data.frame(Filter = integer(), Doublets = integer()) # Making an empty dataframe to put the count before and after filtering
i=1
for (i in 1:length(files)) {
  table[i,1] <- paste(scCITE_QC(dir, files[i], saveDir = savedir),collapse = ",")
}

dir.create(paste(savedir,"counts/",sep = ""), showWarnings = FALSE)
write.table(table, paste(savedir,"counts/cellfiltered.txt",sep = ""), sep = "\t", quote = F, row.names = F, col.names = T)

### Identify the resolution of the object for each sample with the cluster will be used for the doublet finder ####
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/NN_clustertree.R")
object <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "[0-9]_RNA.RDS",full.names = TRUE)
dims <- c(20,20,20,20)
samplename <- c("Sample1","Sample2","Sample3","Sample4")
Assay = "RNA"
process = "nn_clustertree"

for (i in 1:length(object)) {
  nearest_neigbour(Obj =  object[i], Dims = dims[i], saveDir = savedir, samplename = samplename[i], Assay = Assay, process = process)
}

### Running and Filtering the Doublet Finder from the RNA and ADT ####
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITEseq_Doublet_Finder.R")
dims <- c(20,20,20,20)
res <- c(1,1.2,0.4,0.4)
obj <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "nn_clustertree_RNA.RDS",full.names = TRUE)
obj2 <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "[0-9]_ADT.RDS",full.names = TRUE)
samplename <- c("Sample1","Sample2","Sample3","Sample4")
table2 <- read.table(paste(savedir,"counts/cellfiltered.txt",sep = ""),sep = "\t",  header = TRUE)

for (i in 1:length(obj)) {
  table2[i,2]<- paste(doublet_scCITEseq(Obj = obj[i], dims = dims[i],res = res,saveDir = savedir,Obj2 = obj2[i],
                    samplename = samplename[i], process = "doublet_finder", Assay = "RNA"), collapse =",") 
}
dir.create(paste(savedir,"counts/",sep = ""), showWarnings = FALSE)
write.table(table2, paste(savedir,"counts/Doublets_Removed.txt",sep = ""), sep = "\t", quote = F, row.names = F, col.names = T)


#### Adding the unwanted effect i.e. cell cycle phase and mitochondria effect to each samples ####
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cellcycle_mitoscore.R")
obj <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "_RNA_singlet.RDS",full.names = TRUE)
samplename <- gsub("_RNA_singlet.RDS","",basename(obj))

## We have to combine the RNA singlet object to perform the normalization for the unwanted the effects.
combined <- merge(readRDS(obj[1]), y = c(readRDS(obj[2]), readRDS(obj[3]),readRDS(obj[4])), 
                  add.cell.ids = samplename, project = "combined")

samplename = "combined"
process = "merged"
Assay = "RNA"
cellcycle_mito(combined, savedir, ngenes = 4000, process = process, Assay = Assay, samplename = samplename)


### Normalizing and Integrating the RNA dataset ######
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")

combined <- readRDS(paste(savedir,"saveRDS_obj/",samplename,"_",process,"_",Assay,".RDS",sep = ""))
## For each sample in RNA assay they had to be regress out for the Mitochondria, however the cell cycle phase looks good
## Remember we are not integrating we are merging the samples
## By looking at the graph, we add the mitoratio or cellcycle (G2M.Score,S.Score)
objname = "combined"
process = "sctransform"
Assay = "RNA"
sctransform_V2_integration(combined, savedir, regress = c("mitoRatio"), objname = objname,
                           process = process, Assay = Assay)

#### RNA UMAP and FeaturePlot  #######
# After performing the RNA integration we will make the Umap and Feature plot for the certain genes mainly CD4 and CD8 to subset the cells expressing these RNA or protein
# Please check for the emblow plot after the inetgration

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
obj <- paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "") 

objname = "combined_sctransformed"
process = "integration"
Assay = "RNA"
RNA_integration(obj,savedir,dims = 30,RNA_features = c("CD4","CD8A"), Assay=Assay, process=process, objname=objname)


### Performing ADT merging, integration #######
# Based on the ADT UMAP and Feature Plot we will subset the CD4 and CD8 cells
# https://satijalab.org/seurat/articles/integration_rpca.html
# https://github.com/satijalab/seurat/issues/3890
# For protein integration we will use the reciprocal PCA instead of canonical correlation analysis (CCA) by Hao to identify
# the anchors. Normalize using CLR due to less number of ADT.
# In RPCA When determining anchors between any two datasets using RPCA, we project each dataset into the others 
# PCA space and constrain the anchors by the same mutual neighborhood requirement. 
# CCA-based integration therefore enables integrative analysis when experimental conditions or disease states 
# introduce very strong expression shifts, or when integrating datasets across modalities and species.
# However, CCA-based integration may also lead to overcorrection, especially when a large proportion of cells are 
# non-overlapping across datasets.

# RPCA-based integration runs significantly faster, and also represents a more conservative approach where cells 
# in different biological states are less likely to ‘align’ after integration. We therefore,recommend RPCA during 
# integrative analysis where: * A substantial fraction of cells in one dataset have no matching type in the other *

#### Merging and clustering the ADT dataset ######
## Since for ADT we have different object that we have created above, will use them for the integration using RPCA

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
ADT_obj <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "_ADT_singlet.RDS",full.names = TRUE)
ADT_samplename <- gsub("_ADT_singlet.RDS","",basename(ADT_obj))

## For each sample in ADT assay 
## Remember we are not integrating we are merging the samples
## This has to be modified everytime when merging different object based on the number of the object
combined_ADT <- merge(readRDS(ADT_obj[1]), y = c(readRDS(ADT_obj[2]), readRDS(ADT_obj[3]),readRDS(ADT_obj[4])), 
                      add.cell.ids = ADT_samplename, project = "ADT_combined")

## By looking at the graph, we add the mitoratio or cellcycle (G2M.Score,S.Score)
# numFeatures depends on how many antibody that has been used for ADT features
# dims has been decided based on the previous analysis, however after running if you need to change you can make changes in this
objname = "combined"
process = "merging"
Assay = "ADT"
ADT_merging(combined_ADT, savedir, dims = 10, numfeatures=40, Assay=Assay, process=process, objname=objname)

#### ADT Integration ####
# Even though it is named RNA integration then same will work for ADT integration just put all the arguments
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
obj_path= paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")
objname = "combined"
process = "integration"
Assay = "ADT"
RNA_integration(obj_path = obj_path, saveDir = savedir, dims = 10, RNA_features = c("CD4-protein","CD8a-protein"),
                Assay = Assay, process = process, objname = objname)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/NN_clustertree.R")
obj_path = paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")

objname = "combined"
process = "ADT_NN"
Assay = "integrated"
nearest_neigbour(Obj_path = obj_path, Dims = 10, saveDir = savedir, samplename = objname, process = process, Assay = Assay)

#### Running the UMAP and feature pltot for ADT integration ####
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
obj_path = paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")

objname = "combined"
process = "ADT_UMAP_QC"
Assay = "integrated"
cluster_UMAP_and_QC(obj_path = obj_path, dims = 10, res = 0.8, saveDir=savedir, Assay = Assay,
                    QC_features = c("nCount_ADT","nFeature_ADT"), objname = objname, process = process)


## Since ADT has much better clustering, we will adopt with this method for separating both of them.
# we will separate CD4 and CD8, then we will integrate RNA and ADT modalities

######## Separating CD4 and CD8 cluster #########
######### CD4 ##########
# we are using this issue to perform subsetting. Go to the end for atakanekiz,
# Since there is no official comment from Satija lab
# Whose Approach 1 and Approach 2 work perfectly. Showed with the UMAP also. 
# Applying Approach 2 for RNA integration and Approach 2 for ADT integration
#### https://github.com/satijalab/seurat/issues/1883

# Based on the ADT clustering we will subset CD4 clusters look at this location
# For CD4 we have selected cluster 0,5,8,9
########## CD4 ADT #########
set.seed=42
# We need to tak out the cell names from the RNA and ADT integrated object and apply it to the RNA and ADT integrated object
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_subsetting.R")

# In here you to put the final RNA ADT which you have executed. Please go to your <savedir>/saveRDS/
RNA_obj_path <- paste(savedir,"saveRDS_obj/combined_sctransformed_integration_RNA.RDS",sep = "")
ADT_obj_path <- paste(savedir,"saveRDS_obj/combined_ADT_UMAP_QC_integrated.RDS",sep = "")
CD4_cluster <- c(0,5,8,9) ## Need to see the UMAP to find the which cluster belongs to the CD4

objname = "CD4"
Assay = "RNA"
process = "subset"

RNA_subsetting(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path, subset_cluster = CD4_cluster, ngenes = 4000, 
               saveDir = savedir, objname = objname, Assay = Assay, process = process)


# Since both the mito and cell cycle looks good we donot have regress out, however regressing out Mitochondria Effect will be nice
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")
obj <- readRDS(paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = ""))

objname = "CD4"
Assay = "RNA"
process = "sctransform"
sctransform_V2_integration(obj = obj, saveDir = savedir, ngenes = 4000, regress = c("mitoRatio"),
                           Assay = Assay, process = process, objname = objname)

### Checking the Elbow plot, dims = 20

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
obj_path <- paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")

objname = "CD4"
Assay = "RNA"
process = "integration"
RNA_integration(obj_path = obj_path, saveDir = savedir, dims = 20, RNA_features =  c("CD4","CD8A"), 
                Assay=Assay, process=process, objname=objname)

### Checking the Elbow plot, dims = 20

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/NN_clustertree.R")
obj_path <- paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")

objname = "CD4"
process = "nearest_neighbour"
Assay = "integrated"
nearest_neigbour(Obj_path = obj_path, Dims = 20, saveDir = savedir, samplename = objname, process = process, Assay = Assay)

# looking at the cluster tree res = 0.4
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
obj_path <- paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")

objname = "CD4"
process = "UMAP_and_QC"
Assay = "integrated"

cluster_UMAP_and_QC(obj_path = obj_path, dims = 20, res = 0.4, saveDir = savedir, Assay = Assay, 
                    QC_features = c("nCount_RNA","nFeature_RNA","percent.mt"), objname = objname, process = process)


#### CD4 ADT #######
## Since we have already subsetted the RDS file for ADT in the RNA subsetting we will use that one
ADT_obj <- readRDS(paste(savedir,"saveRDS_obj/CD4_subset_ADT.RDS",sep = ""))
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")

objname = "CD4"
Assay = "ADT"
process = "merging"
ADT_merging(ADT_obj = ADT_obj, saveDir = savedir, dims = 10, numfeatures = 40, 
            Assay = Assay, process = process, objname = objname)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
obj_path= paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")
objname = "CD4"
process = "integration"
Assay = "ADT"
RNA_integration(obj_path = obj_path, saveDir = savedir, dims = 10, RNA_features = c("CD4-protein","CD8a-protein"),
                Assay = Assay, process = process, objname = objname)


source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/NN_clustertree.R")
obj_path = paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")

objname = "CD4"
process = "ADT_NN"
Assay = "integrated"
nearest_neigbour(Obj_path = obj_path, Dims = 10, saveDir = savedir, samplename = objname, process = process, Assay = Assay)

# Based on the cluster tree we will take the resolution as 0.2
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
obj_path = paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")

objname = "CD4"
process = "ADT_UMAP_QC"
Assay = "integrated"
cluster_UMAP_and_QC(obj_path = obj_path, dims = 10, res = 0.2, saveDir=savedir, Assay = Assay,
                    QC_features = c("nCount_ADT","nFeature_ADT"), objname = objname, process = process)

#### Integrating the modalities 
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITEseq_modality_integration.R")
# Find the last save object for RNA and ADT
RNA_obj_path <- paste(savedir,"saveRDS_obj/CD4_UMAP_and_QC_integrated.RDS",sep = "")
ADT_obj_path <- paste(savedir,"saveRDS_obj/CD4_ADT_UMAP_QC_integrated.RDS",sep = "")

objname <- "CD4"
process <- "modality_integrate"
Assay <- "integrated"

modality_integration(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path, RNA_dims = 20, ADT_dims = 10,
                     saveDir = savedir, Assay = Assay, process = process, objname = objname)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/Modal_cluster_and_UMAP.R")

obj_path <- paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")
RNA_features = c("CD4","CD8A")
ADT_features = c("CD4-protein","CD8a-protein")

objname = "CD4"
Assay = "integrated"
process = "modality_cluster_and_QC"
Modal_Cluster_and_UMAP_QC(modal_integrate_obj_path = obj_path, saveDir = savedir, res = 0.4, RNA_features = RNA_features, 
                          ADT_features = ADT_features,Assay = Assay, process = process, objname=objname)


###### CD8A ########
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_subsetting.R")

# In here you to put the final RNA ADT which you have executed. Please go to your <savedir>/saveRDS/
RNA_obj_path <- paste(savedir,"saveRDS_obj/combined_sctransformed_integration_RNA.RDS",sep = "")
ADT_obj_path <- paste(savedir,"saveRDS_obj/combined_ADT_UMAP_QC_integrated.RDS",sep = "")

CD8_cluster <- c(1,2,3,4,6,7) ## Need to see the UMAP to find the which cluster belongs to the CD8

objname = "CD8"
Assay = "RNA"
process = "subset"

RNA_subsetting(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path, subset_cluster = CD8_cluster, ngenes = 4000, 
               saveDir = savedir, objname = objname, Assay = Assay, process = process)

# Since both the mito and cell cycle looks good we donot have regress out, however regressing out Mitochondria Effect will be nice
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")
obj <- readRDS(paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = ""))

objname = "CD8"
Assay = "RNA"
process = "sctransform"
sctransform_V2_integration(obj = obj, saveDir = savedir, ngenes = 4000, regress = c("mitoRatio"),
                           Assay = Assay, process = process, objname = objname)

### Checking the Elbow plot, dims = 20

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
obj_path <- paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")

objname = "CD8"
Assay = "RNA"
process = "integration"
RNA_integration(obj_path = obj_path, saveDir = savedir, dims = 20, RNA_features =  c("CD4","CD8A"), 
                Assay=Assay, process=process, objname=objname)

### Checking the Elbow plot, dims = 20
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/NN_clustertree.R")
obj_path <- paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")

objname = "CD8"
process = "nearest_neighbour"
Assay = "integrated"
nearest_neigbour(Obj_path = obj_path, Dims = 20, saveDir = savedir, samplename = objname, process = process, Assay = Assay)

# looking at the cluster tree res = 0.4
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
obj_path <- paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")

objname = "CD8"
process = "UMAP_and_QC"
Assay = "integrated"

cluster_UMAP_and_QC(obj_path = obj_path, dims = 20, res = 0.4, saveDir = savedir, Assay = Assay, 
                    QC_features = c("nCount_RNA","nFeature_RNA","percent.mt"), objname = objname, process = process)

#### CD8 ADT #######
## Since we have already subsetted the RDS file for ADT in the RNA subsetting we will use that one
ADT_obj <- readRDS(paste(savedir,"saveRDS_obj/CD8_subset_ADT.RDS",sep = ""))
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")

objname = "CD8"
Assay = "ADT"
process = "merging"
ADT_merging(ADT_obj = ADT_obj, saveDir = savedir, dims = 10, numfeatures = 40, 
            Assay = Assay, process = process, objname = objname)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
obj_path= paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")
objname = "CD8"
process = "integration"
Assay = "ADT"
RNA_integration(obj_path = obj_path, saveDir = savedir, dims = 10, RNA_features = c("CD4-protein","CD8a-protein"),
                Assay = Assay, process = process, objname = objname)


source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/NN_clustertree.R")
obj_path = paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")

objname = "CD8"
process = "ADT_NN"
Assay = "integrated"
nearest_neigbour(Obj_path = obj_path, Dims = 10, saveDir = savedir, samplename = objname, process = process, Assay = Assay)

# Based on the cluster tree we will take the resolution as 0.8
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
obj_path = paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")

objname = "CD8"
process = "ADT_UMAP_QC"
Assay = "integrated"
cluster_UMAP_and_QC(obj_path = obj_path, dims = 10, res = 0.6, saveDir=savedir, Assay = Assay,
                    QC_features = c("nCount_ADT","nFeature_ADT"), objname = objname, process = process)

#### Integrating the modalities 
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITEseq_modality_integration.R")
# Find the last save object for RNA and ADT
RNA_obj_path <- paste(savedir,"saveRDS_obj/CD8_UMAP_and_QC_integrated.RDS",sep = "")
ADT_obj_path <- paste(savedir,"saveRDS_obj/CD8_ADT_UMAP_QC_integrated.RDS",sep = "")

objname <- "CD8"
process <- "modality_integrate"
Assay <- "integrated"

modality_integration(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path, RNA_dims = 20, ADT_dims = 10,
                     saveDir = savedir, Assay = Assay, process = process, objname = objname)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/Modal_cluster_and_UMAP.R")

obj_path <- paste(savedir,"saveRDS_obj/",objname,"_",process,"_",Assay,".RDS",sep = "")
RNA_features = c("CD4","CD8A")
ADT_features = c("CD4-protein","CD8a-protein")

objname = "CD8"
Assay = "integrated"
process = "modality_cluster_and_QC"
Modal_Cluster_and_UMAP_QC(modal_integrate_obj_path = obj_path, saveDir = savedir, res = 0.8, RNA_features = RNA_features, 
                          ADT_features = ADT_features,Assay = Assay, process = process, objname=objname)
