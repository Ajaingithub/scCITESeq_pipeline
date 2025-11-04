# Processed it # /mnt/data/projects/.immune/Mayo/Ines/preprocessing.R and now running the T cells.

library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.
savedir = "/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Run06/"

#### Generating the singlet file for each run #####
## In directory Run02 by mistake I wrote Run07 so ignore it
dir_run01 = "/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/preprocessing/Run06_CR9/outs/per_sample_outs/"
dir_run01_files = list.files(dir_run01, pattern = "T_",full.names = T)

files <- c(dir_run01_files)
id =  basename(files)
unique(id)

source("/mnt/data/projects/pipeline_functions/CR7_scCITESeq_QC.R")
samplepath = files[which(basename(files) %in% id)]
Sample = basename(samplepath)
table <- data.frame(Filter = integer(), Doublets = integer()) # Making an empty dataframe to put the count before and after filtering
i=1
for (i in 1:length(samplepath)) {
  table[i,1] <- paste(scCITE_QC(samplepath[i],Sample[i], saveDir = savedir, min_cells=3, min_genes=700, max_genes=6000),collapse = ",")
}

rownames(table) <- Sample
dir.create(paste(savedir,"counts/",sep = ""), showWarnings = FALSE)
write.table(table, paste(savedir,"counts/cellfiltered.txt",sep = ""), sep = "\t", quote = F, row.names = T, col.names = T)

## Identify the resolution of the object for each sample with the cluster will be used for the doublet finder ####
## Before running this please change the source file, try to save the RDS files
source("/mnt/data/projects/pipeline_functions/NN_clustertree.R")
object <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "_RNA.RDS",full.names = TRUE)
object <- object[grep("clustertree",object,invert = TRUE)]
dims =20
samplename <- gsub("_RNA.RDS","",basename(object))
Assay = "RNA"
process = "nn_clustertree"

for (i in 1:length(object)) {
  obj = readRDS(object[i])
  obj = nearest_neigbour(Obj =  obj, Dims = dims, saveDir = savedir, samplename = samplename[i], Assay = Assay, process = process)
  assign(paste(samplename[i],"_nn",sep = ""), obj)
}

### Running and Filtering the Doublet Finder from the RNA and ADT ####
# https://github.com/satijalab/seurat/issues/1565 resolution
# Keeping minimum resolution to 0.4 and max 1.2 based on the Seurat essence, there is no correct clustering parameter,
# either you will over or under cluster your data. To compensate for what makes biological sense in the context of your experiment, you can merge certain clusters togethe
source("/mnt/data/projects/pipeline_functions/scCITEseq_Doublet_Finder.R")
dims <- 20
res <- 0.4
# obj <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "_RNA.RDS",full.names = TRUE)
# obj <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "_nn_clustertree_RNA.RDS",full.names = TRUE)
obj = ls(pattern="_nn")
obj2 <- list.files(paste(savedir,"saveRDS_obj",sep = ""),pattern = "ADT.RDS",full.names = TRUE)

samplename <- gsub("_nn","",basename(obj))
table2 <- read.table(paste(savedir,"counts/cellfiltered.txt",sep = ""),sep = "\t",  header = TRUE)

i=1
for (i in 1:length(obj)) {
  table2[i,2]<- paste(doublet_scCITEseq(Obj = get(obj[i]), dims = dims,res = res, saveDir = savedir,Obj2 = obj2[i],
                                        samplename = samplename[i], process = "doublet_finder", Assay = "RNA"), collapse =",")
}

table2_1 <- data.frame(samplename, table2$Doublets)
rownames(table2_1) <- table2_1$samplename
data_frame_merge <- merge(table, table2_1,
                          by = 'row.names', all = TRUE)

rownames(data_frame_merge) <- data_frame_merge$Row.names
data_frame_merge_filtered <- data_frame_merge[,c("Filter","table2.Doublets")]

dir.create(paste(savedir,"counts/",sep = ""), showWarnings = FALSE)
write.table(data_frame_merge_filtered, paste(savedir,"counts/Doublets_Removed.txt",sep = ""), sep = "\t", quote = F, row.names = T, col.names = T)


