### Downstream analysis for the samples between 65 to 80 years old

library("Seurat")
library("dplyr")
library("sctransform") ## we will use the sctransform version 2
library("glmGamPoi") # substantially improves the speed of the learning procedure.
savedir = "/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Tcell_combined/"

### SaveRDS file singlet already geenrated
run01 = list.files("/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Run01/saveRDS_obj",pattern = "_RNA_singlet.RDS", full.names = TRUE)
run02 = list.files("/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Run02/saveRDS_obj",pattern = "_RNA_singlet.RDS", full.names = TRUE)
run03 = list.files("/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Run03/saveRDS_obj",pattern = "_RNA_singlet.RDS", full.names = TRUE)
run04 = list.files("/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Run04/saveRDS_obj",pattern = "_RNA_singlet.RDS", full.names = TRUE)
run05 = list.files("/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Run05/saveRDS_obj",pattern = "_RNA_singlet.RDS", full.names = TRUE)
run06 = list.files("/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Run06/saveRDS_obj",pattern = "_RNA_singlet.RDS", full.names = TRUE)
RNA_RDS <- c(run01,run02,run03,run04,run05,run06)

samplename <- gsub("_RNA_singlet.RDS","",basename(RNA_RDS))

## We have to combine the RNA singlet object to perform the normalization for the unwanted the effects.
combined <- merge(x=readRDS(RNA_RDS[1]), y = c(readRDS(RNA_RDS[2]),readRDS(RNA_RDS[3]),readRDS(RNA_RDS[4]),readRDS(RNA_RDS[5]),readRDS(RNA_RDS[6]),
                                         readRDS(RNA_RDS[7]),readRDS(RNA_RDS[8]),readRDS(RNA_RDS[9]),readRDS(RNA_RDS[10]),readRDS(RNA_RDS[11]),
                                         readRDS(RNA_RDS[12]),readRDS(RNA_RDS[13]),readRDS(RNA_RDS[14]),readRDS(RNA_RDS[15]),readRDS(RNA_RDS[16]),
                                         readRDS(RNA_RDS[17]),readRDS(RNA_RDS[18]),readRDS(RNA_RDS[19]),readRDS(RNA_RDS[20]),readRDS(RNA_RDS[21]),
                                         readRDS(RNA_RDS[22]),readRDS(RNA_RDS[23]),readRDS(RNA_RDS[24]),readRDS(RNA_RDS[25]),readRDS(RNA_RDS[26]),
                                         readRDS(RNA_RDS[27]), readRDS(RNA_RDS[28]),readRDS(RNA_RDS[29]),readRDS(RNA_RDS[30]),readRDS(RNA_RDS[31]),
                                         readRDS(RNA_RDS[32]),readRDS(RNA_RDS[33]),readRDS(RNA_RDS[34]),readRDS(RNA_RDS[35]),readRDS(RNA_RDS[36]),
                                         readRDS(RNA_RDS[37]),readRDS(RNA_RDS[38]),readRDS(RNA_RDS[39]),readRDS(RNA_RDS[40]),readRDS(RNA_RDS[41]),
                                         readRDS(RNA_RDS[42]),readRDS(RNA_RDS[43]),readRDS(RNA_RDS[44]),readRDS(RNA_RDS[45]),readRDS(RNA_RDS[46]),
                                         readRDS(RNA_RDS[47]),readRDS(RNA_RDS[48]),readRDS(RNA_RDS[49]),readRDS(RNA_RDS[50]),readRDS(RNA_RDS[51]),
                                         readRDS(RNA_RDS[52]), readRDS(RNA_RDS[53]),readRDS(RNA_RDS[54]),readRDS(RNA_RDS[55]),readRDS(RNA_RDS[56]),
                                         readRDS(RNA_RDS[57]),readRDS(RNA_RDS[58]),readRDS(RNA_RDS[59]),readRDS(RNA_RDS[60])),
                  add.cell.ids = samplename, project = "combined")

run01 = list.files("/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Run01/saveRDS_obj",pattern = "_ADT_singlet.RDS", full.names = TRUE)
run02 = list.files("/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Run02/saveRDS_obj",pattern = "_ADT_singlet.RDS", full.names = TRUE)
run03 = list.files("/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Run03/saveRDS_obj",pattern = "_ADT_singlet.RDS", full.names = TRUE)
run04 = list.files("/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Run04/saveRDS_obj",pattern = "_ADT_singlet.RDS", full.names = TRUE)
run05 = list.files("/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Run05/saveRDS_obj",pattern = "_ADT_singlet.RDS", full.names = TRUE)
run06 = list.files("/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Run06/saveRDS_obj",pattern = "_ADT_singlet.RDS", full.names = TRUE)
ADT_RDS <- c(run01,run02,run03,run04,run05,run06)
samplename_ADT <- paste(samplename,"ADT",sep = "_")

## We have to combine the RNA singlet object to perform the normalization for the unwanted the effects.
combined_ADT <- merge(x=readRDS(ADT_RDS[1]), y = c(readRDS(ADT_RDS[2]),readRDS(ADT_RDS[3]),readRDS(ADT_RDS[4]),readRDS(ADT_RDS[5]),readRDS(ADT_RDS[6]),
                                         readRDS(ADT_RDS[7]),readRDS(ADT_RDS[8]),readRDS(ADT_RDS[9]),readRDS(ADT_RDS[10]),readRDS(ADT_RDS[11]),
                                         readRDS(ADT_RDS[12]),readRDS(ADT_RDS[13]),readRDS(ADT_RDS[14]),readRDS(ADT_RDS[15]),readRDS(ADT_RDS[16]),
                                         readRDS(ADT_RDS[17]),readRDS(ADT_RDS[18]),readRDS(ADT_RDS[19]),readRDS(ADT_RDS[20]),readRDS(ADT_RDS[21]),
                                         readRDS(ADT_RDS[22]),readRDS(ADT_RDS[23]),readRDS(ADT_RDS[24]),readRDS(ADT_RDS[25]),readRDS(ADT_RDS[26]),
                                         readRDS(ADT_RDS[27]), readRDS(ADT_RDS[28]),readRDS(ADT_RDS[29]),readRDS(ADT_RDS[30]),readRDS(ADT_RDS[31]),
                                         readRDS(ADT_RDS[32]),readRDS(ADT_RDS[33]),readRDS(ADT_RDS[34]),readRDS(ADT_RDS[35]),readRDS(ADT_RDS[36]),
                                         readRDS(ADT_RDS[37]),readRDS(ADT_RDS[38]),readRDS(ADT_RDS[39]),readRDS(ADT_RDS[40]),readRDS(ADT_RDS[41]),
                                         readRDS(ADT_RDS[42]),readRDS(ADT_RDS[43]),readRDS(ADT_RDS[44]),readRDS(ADT_RDS[45]),readRDS(ADT_RDS[46]),
                                         readRDS(ADT_RDS[47]),readRDS(ADT_RDS[48]),readRDS(ADT_RDS[49]),readRDS(ADT_RDS[50]),readRDS(ADT_RDS[51]),
                                         readRDS(ADT_RDS[52]), readRDS(ADT_RDS[53]),readRDS(ADT_RDS[54]),readRDS(ADT_RDS[55]),readRDS(ADT_RDS[56]),
                                         readRDS(ADT_RDS[57]),readRDS(ADT_RDS[58]),readRDS(ADT_RDS[59]),readRDS(ADT_RDS[60])),
                  add.cell.ids = samplename, project = "combined_ADT")

combined_ADT <- NormalizeData(combined_ADT, normalization.method = 'CLR', margin = 2)
combined_ADT <- FindVariableFeatures(combined_ADT, selection.method = "vst", nfeatures = 100)

dir.create(paste0(savedir,"saveRDS_obj"),showWarnings = FALSE)
# saveRDS(combined,paste(savedir,"saveRDS_obj/combined.RDS",sep = ""))
# saveRDS(combined_ADT,paste(savedir,"saveRDS_obj/combined_ADT.RDS",sep = ""))

# combined = readRDS(paste0(savedir,"saveRDS_obj/combined.RDS"))
# combined_ADT = readRDS(paste0(savedir,"saveRDS_obj/combined_ADT.RDS"))

source("/mnt/data/projects/pipeline_functions/cellcycle_mitoscore.R")
samplename = "combined"
process = "merged"
Assay = "RNA"
combined <- cellcycle_mito(combined, savedir, ngenes = 4000, process = process, Assay = Assay, samplename = samplename)

combined@meta.data$sampleID <- gsub("T_","",combined@meta.data$orig.ident) %>% gsub("_.*.","",.)
combined@meta.data$condition <- gsub(".*._","",combined@meta.data$orig.ident)
combined@meta.data$AgeSex <- gsub("T_","",combined@meta.data$orig.ident) %>% gsub("Y[0-9]_|O[0-9]_","",.) %>% gsub("_.*.","",.)

# saveRDS(combined,paste(savedir,"saveRDS_obj/combined.RDS",sep = ""))
# saveRDS(combined_ADT,paste(savedir,"saveRDS_obj/combined_ADT.RDS",sep = ""))

# combined = readRDS(paste(savedir,"saveRDS_obj/combined.RDS",sep = ""))
# combined_ADT = readRDS(paste(savedir,"saveRDS_obj/combined_ADT.RDS",sep = ""))

# GEX_CD4_CD8 <- as.data.frame(t(as.data.frame(combined@assays$RNA@data[c("CD4","CD8A"),])))
combined_ADT <- JoinLayers(combined_ADT)
ADT_TCR_alpha_delta <- as.data.frame(t(as.data.frame(combined_ADT@assays$ADT@layers$data[grep("TCRa-b-protein|TCRVd2-protein",rownames(combined_ADT)),])))
colnames(ADT_TCR_alpha_delta) <- c("TCRa_b_protein","TCRVd2_protein")
grep("TCRa-b-protein|TCRVd2-protein",rownames(combined_ADT),value=TRUE)
rownames(ADT_TCR_alpha_delta) <- colnames(combined_ADT)

p1 <- ggplot(ADT_TCR_alpha_delta, aes_string(x="TCRa-b-protein",y="TCRVd2-protein")) +
  geom_point() + theme_bw() +
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  geom_vline(xintercept=0.75, linetype="dashed", color = "red") +
  ggtitle("ADT TCR Alpha and delta")

dir.create(paste(savedir,"scatter_plot",sep = ""),showWarnings = FALSE)

pdf(paste(savedir,"scatter_plot/ADT_TCR_alpha_and_delta.pdf",sep = ""))
p1
dev.off()

# separating out the CD8 and CD4 T cells
library(Seurat)
savedir = "/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Tcell_combined/"
combined = readRDS(paste0(savedir,"saveRDS_obj/combined.RDS"))
combined_ADT = readRDS(paste0(savedir,"saveRDS_obj/combined_ADT.RDS"))

## First extract out the TCR alpha cells removing the delta cells
# Removing TCR delta cells
TCRalphacells=rownames(ADT_TCR_alpha_delta[ADT_TCR_alpha_delta$TCRVd2_protein < 1,])

### combined_ADT
combined_ADT_TCR_alpha = subset(combined_ADT, cells = TCRalphacells)

ADT_CD4_CD8 <- as.data.frame(t(as.data.frame(combined_ADT_TCR_alpha@assays$ADT@layers$data[grep("CD4-protein|CD8-protein",rownames(combined_ADT_TCR_alpha)),])))
colnames(ADT_CD4_CD8) <- c("CD8_protein","CD4_protein")
grep("CD4-protein|CD8-protein",rownames(combined_ADT),value=TRUE)
rownames(ADT_CD4_CD8) <- colnames(combined_ADT_TCR_alpha)

p1 <- ggplot(ADT_CD4_CD8, aes_string(x="CD4_protein",y="CD8_protein")) +
  geom_point() + theme_bw() +
  geom_hline(yintercept=0.8, linetype="dashed", color = "red") +
  geom_vline(xintercept=0.8, linetype="dashed", color = "red") +
  ggtitle("CD4_CD8 protein")

dir.create(paste(savedir,"scatter_plot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"scatter_plot/CD4_CD8_scatter_plot.pdf",sep = ""))
p1
dev.off()

DN <- subset(ADT_CD4_CD8, CD4_protein<0.8 & CD8_protein<0.8)
DP <- subset(ADT_CD4_CD8, CD4_protein>0.8 & CD8_protein>0.8)

CD8_TP <- subset(ADT_CD4_CD8, CD4_protein<0.8 & CD8_protein>0.8)
CD4_TP <- subset(ADT_CD4_CD8, CD4_protein>0.8 & CD8_protein<0.8)

## extracting out GEX CD4, CD8A, and CD8B normalized gene expression
combined_GEX_DN_ADT = subset(combined, cells = rownames(DN))

combined_GEX_DN_ADT = JoinLayers(combined_GEX_DN_ADT)
GEX_CD4_CD8 <- as.data.frame(t(as.data.frame(combined_GEX_DN_ADT@assays$RNA@layers$data[grep("^CD4$|^CD8A$|^CD8B$",rownames(combined_GEX_DN_ADT)),])))
colnames(GEX_CD4_CD8) <- grep("^CD4$|^CD8A$|^CD8B$",rownames(combined_GEX_DN_ADT),value=TRUE)
rownames(GEX_CD4_CD8) <- colnames(combined_GEX_DN_ADT)

p1 <- ggplot(GEX_CD4_CD8, aes_string(x="CD4",y="CD8A")) +
  geom_point() + theme_bw() +
  geom_hline(yintercept=0.2, linetype="dashed", color = "red") +
  geom_vline(xintercept=0.2, linetype="dashed", color = "red") +
  ggtitle("CD4_CD8 GEX")

dir.create(paste(savedir,"scatter_plot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"scatter_plot/CD4_CD8_GEX_scatter_plot.pdf",sep = ""))
p1
dev.off()

CD8_TP_RNA <- subset(GEX_CD4_CD8, CD4<0.2 & CD8A>0.2)
CD4_TP_RNA <- subset(GEX_CD4_CD8, CD4>0.2 & CD8A<0.2)

DN_RNA <- subset(GEX_CD4_CD8, CD4<0.2 & CD8A<0.2)
DP_RNA <- subset(GEX_CD4_CD8, CD4>0.2 & CD8A>0.2)

### Calculating CD8B from DN CD4 and CD8A
p1 <- ggplot(DN_RNA, aes_string(x="CD4",y="CD8B")) +
  geom_point() + theme_bw() +
  geom_hline(yintercept=0.2, linetype="dashed", color = "red") +
  geom_vline(xintercept=0.2, linetype="dashed", color = "red") +
  ggtitle("CD4_CD8B GEX")

dir.create(paste(savedir,"scatter_plot",sep = ""),showWarnings = FALSE)
pdf(paste(savedir,"scatter_plot/CD4_CD8B_GEX_scatter_plot.pdf",sep = ""))
p1
dev.off()

CD8B_RNA = DN_RNA[DN_RNA$CD8B > 0,]

dir.create(paste(savedir,"Table",sep = ""),showWarnings = FALSE)
write.table(rownames(CD8_TP),paste(savedir,"Table/CD8_TP.txt",sep = ""),quote = F, sep = "\t",row.names = F, col.names = F)
write.table(rownames(CD8_TP_RNA),paste(savedir,"Table/CD8_TP_RNA.txt",sep = ""),quote = F, sep = "\t",row.names = F, col.names = F)
write.table(rownames(DN_RNA),paste(savedir,"Table/DN_RNA.txt",sep = ""),quote = F, sep = "\t",row.names = F, col.names = F)
write.table(rownames(DP_RNA),paste(savedir,"/Table/DP_RNA.txt",sep = ""),quote = F, sep = "\t",row.names = F, col.names = F)

### loading the cells
CD4_TP = read.table(paste0(savedir,"Table/CD4_TP.txt"), header = TRUE, sep = "\t")
CD4_TP_RNA = read.table(paste0(savedir,"Table/CD4_TP_RNA.txt"), header = TRUE, sep = "\t")
CD4_cellnames = c(rownames(CD4_TP), rownames(CD4_TP_RNA))

CD8_TP = read.table(paste0(savedir,"Table/CD8_TP.txt"), header = TRUE, sep = "\t")
CD8A_TP_RNA = read.table(paste0(savedir,"Table/CD8_TP_RNA.txt"), header = TRUE, sep = "\t")
CD8B_TP_RNA = read.table(paste0(savedir,"Table/CD8B_RNA.txt"), header = TRUE, sep = "\t")
CD8_cellnames = c(rownames(CD8_TP), rownames(CD8A_TP_RNA), rownames(CD8B_TP_RNA))

write.table(CD4_cellnames,paste(savedir,"Table/CD4_cellnames.txt",sep = ""), quote = F, col.names = F, row.names = F)
write.table(CD8_cellnames,paste(savedir,"Table/CD8_cellnames.txt",sep = ""), quote = F, col.names = F, row.names = F)

### Since everything is ran in the one batch. Not performing any integration
#### subsetting CD4 and CD8 separately
#### CD4 subsetting ######
savedir = "/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Tcell_combined/"
CD4_cellnames <- read.table(paste(savedir,"Table/CD4_cellnames.txt",sep = ""))[,1]
CD8_cellnames <- read.table(paste(savedir,"Table/CD8_cellnames.txt",sep = ""))[,1]

CD4_combined = subset(combined, cells = CD4_cellnames)
CD4_combined = NormalizeData(CD4_combined)
CD4_combined = FindVariableFeatures(CD4_combined, selection.method = "vst", nfeatures = 4000)
CD4_combined = ScaleData(CD4_combined)
CD4_combined = RunPCA(CD4_combined)
CD4_combined = FindNeighbors(CD4_combined, dims = 1:30)
CD4_combined = RunUMAP(CD4_combined, dims = 1:30)

savedir = "/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Tcell_combined/CD4/"
dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE)
pdf(paste0(savedir,"UMAP/CD4_combined_UMAP.pdf"))
DimPlot(CD4_combined, group.by = "orig.ident") + NoLegend()
dev.off()

pdf(paste0(savedir,"UMAP/CD4_combined_UMAP_condition.pdf"))
DimPlot(CD4_combined, group.by = "condition", label = TRUE)
dev.off()

pdf(paste0(savedir,"UMAP/CD4_combined_UMAP_condition_splitted.pdf"))
DimPlot(CD4_combined, group.by = "condition", label = FALSE, split.by = "condition", ncol = 3)
dev.off()

CD4_combined <- FindClusters(CD4_combined, resolution = 0.8)

pdf(paste0(savedir,"UMAP/CD4_combined_UMAP_seurat_clus.pdf"))
DimPlot(CD4_combined, label = TRUE)
dev.off()

dir.create(paste0(savedir,"saveRDS_obj"), showWarnings = FALSE)
saveRDS(CD4_combined, paste0(savedir,"saveRDS_obj/CD4_combined.RDS"))

CD4_combined = JoinLayers(CD4_combined)
CD4_markers = FindAllMarkers(CD4_combined)

write.table(CD4_markers, paste0(savedir,"Table/CD4_all_marker.txt"), sep = "\t", row.names = T, col.names = T, quote = F)

### Combined ADT
CD4_combined_ADT = subset(combined_ADT, cells = CD4_cellnames)
CD4_combined_ADT = NormalizeData(CD4_combined_ADT, normalization.method = 'CLR', margin = 2)
CD4_combined_ADT = FindVariableFeatures(CD4_combined_ADT, selection.method = "vst", nfeatures = 70)
CD4_combined_ADT = ScaleData(CD4_combined_ADT)
CD4_combined_ADT = RunPCA(CD4_combined_ADT)
savedir = "/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Tcell_combined/CD4/"
dir.create(paste(savedir, "elbow_plots/", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "elbow_plots/CD4_combined_ADT.pdf", sep = ""), width = 8, height = 8)
print(ElbowPlot(CD4_combined_ADT, ndims = 30, reduction = "pca"))
dev.off()
CD4_combined_ADT = FindNeighbors(CD4_combined_ADT, dims = 1:12)
CD4_combined_ADT = RunUMAP(CD4_combined_ADT, dims = 1:12)

savedir = "/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Tcell_combined/CD4/"
dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE)
pdf(paste0(savedir,"UMAP/CD4_combined_ADT_UMAP.pdf"))
DimPlot(CD4_combined_ADT, group.by = "orig.ident") + NoLegend()
dev.off()

CD4_combined_ADT@meta.data$condition= gsub(".*._","",CD4_combined_ADT@meta.data$orig.ident)
pdf(paste0(savedir,"UMAP/CD4_combined_ADT_UMAP_condition.pdf"))
DimPlot(CD4_combined_ADT, group.by = "condition", label = TRUE)
dev.off()

pdf(paste0(savedir,"UMAP/CD4_combined_ADT_UMAP_condition_splitted.pdf"))
DimPlot(CD4_combined_ADT, group.by = "condition", label = FALSE, split.by = "condition", ncol = 3)
dev.off()

CD4_combined_ADT <- FindClusters(CD4_combined_ADT, resolution = 0.8)

pdf(paste0(savedir,"UMAP/CD4_combined_ADT_UMAP_seurat_clus.pdf"))
DimPlot(CD4_combined_ADT, label = TRUE)
dev.off()

saveRDS(CD4_combined_ADT, paste0(savedir,"saveRDS_obj/CD4_combined_ADT.RDS"))

CD4_combined_ADT = JoinLayers(CD4_combined_ADT)
CD4_ADT_markers = FindAllMarkers(CD4_combined_ADT)

write.table(CD4_ADT_markers, paste0(savedir,"Table/CD4_all_marker_ADT.txt"), sep = "\t", row.names = T, col.names = T, quote = F)

#region CD8
#### Performing the same with CD8
### Since everything is ran in the one batch. Not performing any integration
#### CD8 subsetting ######
savedir = "/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Tcell_combined/"
# CD4_cellnames <- read.table(paste(savedir,"Table/CD4_cellnames.txt",sep = ""))[,1]
CD8_cellnames = read.table(paste(savedir,"Table/CD8_cellnames.txt",sep = ""))[,1]
CD8_combined = subset(combined, cells = CD8_cellnames)
CD8_combined = NormalizeData(CD8_combined)
CD8_combined = FindVariableFeatures(CD8_combined, selection.method = "vst", nfeatures = 4000)
CD8_combined = ScaleData(CD8_combined)
CD8_combined = RunPCA(CD8_combined)
CD8_combined = FindNeighbors(CD8_combined, dims = 1:30)
CD8_combined = RunUMAP(CD8_combined, dims = 1:30)

savedir = "/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Tcell_combined/CD8/"
dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE)
pdf(paste0(savedir,"UMAP/CD8_combined_UMAP.pdf"))
DimPlot(CD8_combined, group.by = "orig.ident") + NoLegend()
dev.off()

pdf(paste0(savedir,"UMAP/CD8_combined_UMAP_condition.pdf"))
DimPlot(CD8_combined, group.by = "condition", label = TRUE)
dev.off()

pdf(paste0(savedir,"UMAP/CD8_combined_UMAP_condition_splitted.pdf"))
DimPlot(CD8_combined, group.by = "condition", label = FALSE, split.by = "condition", ncol = 3)
dev.off()

CD8_combined <- FindClusters(CD8_combined, resolution = 0.8)

pdf(paste0(savedir,"UMAP/CD8_combined_UMAP_seurat_clus.pdf"))
DimPlot(CD8_combined, label = TRUE)
dev.off()

dir.create(paste0(savedir,"saveRDS_obj"), showWarnings = FALSE)
saveRDS(CD8_combined, paste0(savedir,"saveRDS_obj/CD8_combined.RDS"))

CD8_combined = JoinLayers(CD8_combined)
CD8_markers = FindAllMarkers(CD8_combined)

dir.create(paste0(savedir,"Table/"), showWarnings=FALSE)
write.table(CD8_markers, paste0(savedir,"Table/CD8_all_marker.txt"), sep = "\t", row.names = T, col.names = T, quote = F)

### Combined ADT
CD8_combined_ADT = subset(combined_ADT, cells = CD8_cellnames)
CD8_combined_ADT = NormalizeData(CD8_combined_ADT, normalization.method = 'CLR', margin = 2)
CD8_combined_ADT = FindVariableFeatures(CD8_combined_ADT, selection.method = "vst", nfeatures = 70)
CD8_combined_ADT = ScaleData(CD8_combined_ADT)
CD8_combined_ADT = RunPCA(CD8_combined_ADT)
savedir = "/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Tcell_combined/CD8/"
dir.create(paste(savedir, "elbow_plots/", sep = ""), showWarnings = FALSE)
pdf(paste(savedir, "elbow_plots/CD8_combined_ADT.pdf", sep = ""), width = 8, height = 8)
print(ElbowPlot(CD8_combined_ADT, ndims = 30, reduction = "pca"))
dev.off()
CD8_combined_ADT = FindNeighbors(CD8_combined_ADT, dims = 1:12)
CD8_combined_ADT = RunUMAP(CD8_combined_ADT, dims = 1:12)

savedir = "/mnt/data/projects/.immune/Mayo/Ines/Fibro_Tcells/analysis/downstream/seurat_R/Tcell_combined/CD8/"
dir.create(paste0(savedir,"UMAP"), showWarnings = FALSE)
pdf(paste0(savedir,"UMAP/CD8_combined_ADT_UMAP.pdf"))
DimPlot(CD8_combined_ADT, group.by = "orig.ident") + NoLegend()
dev.off()

CD8_combined_ADT@meta.data$condition= gsub(".*._","",CD8_combined_ADT@meta.data$orig.ident)
pdf(paste0(savedir,"UMAP/CD8_combined_ADT_UMAP_condition.pdf"))
DimPlot(CD8_combined_ADT, group.by = "condition", label = TRUE)
dev.off()

pdf(paste0(savedir,"UMAP/CD8_combined_ADT_UMAP_condition_splitted.pdf"))
DimPlot(CD8_combined_ADT, group.by = "condition", label = FALSE, split.by = "condition", ncol = 3)
dev.off()

CD8_combined_ADT <- FindClusters(CD8_combined_ADT, resolution = 0.8)

pdf(paste0(savedir,"UMAP/CD8_combined_ADT_UMAP_seurat_clus.pdf"))
DimPlot(CD8_combined_ADT, label = TRUE)
dev.off()

saveRDS(CD8_combined_ADT, paste0(savedir,"saveRDS_obj/CD8_combined_ADT.RDS"))

CD8_combined_ADT = JoinLayers(CD8_combined_ADT)
CD8_ADT_markers = FindAllMarkers(CD8_combined_ADT)

write.table(CD8_ADT_markers, paste0(savedir,"Table/CD8_all_marker_ADT.txt"), sep = "\t", row.names = T, col.names = T, quote = F)



# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_subsetting.R")
# RNA_obj_path = combined
# ADT_obj_path = combined_ADT
# ADT_main = combined_ADT
# CD4_cluster <- NULL
# objname = "CD4"
# Assay = "RNA"
# process = "subset"

# CD4_RNA <- RNA_subsetting(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path, ADT_main = ADT_main, subset_cluster = CD4_cluster,
#                           cellnames=CD4_cellnames, ngenes = 4000, saveDir = savedir, objname = objname, Assay = Assay, process = process)

# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_sctransform_V2.R")
# objname = "CD4"
# Assay = "RNA"
# process = "sctransform"
# CD4_RNA_sctransformed <- sctransform_V2_integration(obj = CD4_RNA, saveDir = savedir, ngenes = 4000,
#                                                     regress = c("mitoRatio","G2M.Score","S.Score","nCount_RNA"),
#                                                     dims = 30,
#                                                     Assay = Assay, process = process, objname = objname,
#                                                     split_by = "Run",
#                                                     reference = NULL,
#                                                     sample_tree = NULL)

# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
# objname = "CD4"
# process = "integration"
# Assay = "RNA"
# CD4_integrated_RNA <- RNA_integration(CD4_RNA_sctransformed, savedir, dims = 30, RNA_features = c("CD4","CD8A"),
#                                       Assay=Assay, process=process, objname=objname, ncol = 4, ndims = 50)

# CD4_integrated_RNA_NN <- FindNeighbors(CD4_integrated_RNA, dims = 1:30)

# source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
# res = c(0.2,0.4,0.6,0.8,1)
# for (i in 1:length(res)) {
#   process = paste("RNA_UMAP_QC",res[i],sep = "_")
#   Assay = "integrated"
#   CD4_integrated_RNA_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD4_integrated_RNA_NN, dims = 30, res = res[i], saveDir=savedir, Assay = Assay,
#                                                        QC_features = c("nCount_RNA","nFeature_RNA"), objname = "CD4_integrated_RNA_NN_cluster",
#                                                        process = process, col_sel = c("gender","Age","vaccine","condition","Run","condition","orig.ident"))
# }

# saveRDS(CD4_integrated_RNA_NN_cluster,paste(savedir,"saveRDS_obj/CD4_integrated_RNA_NN_cluster",sep = ""))


### CD4 ADT #######
CD4_ADT = readRDS(paste(savedir,"saveRDS_obj/CD4_subset_ADT.RDS",sep = ""))

CD4_ADT@meta.data$orig.ident <- gsub("_[A|T|G|C].*","",rownames(CD4_ADT@meta.data))
CD4_ADT@meta.data$condition = gsub("_.*","",CD4_ADT@meta.data$orig.ident)
CD4_ADT@meta.data$Run = gsub(".*_","",CD4_ADT@meta.data$orig.ident)
CD4_ADT@meta.data$vaccine <- gsub("DMSO_|gE_|_O.*.","",CD4_ADT@meta.data$orig.ident)
CD4_ADT@meta.data[grep("DMSO",CD4_ADT@meta.data$condition),"vaccine"] <- "DMSO"
# DMSO     S     Z
# 17642 37560 24146

CD4_ADT@meta.data$sample_id <- gsub(".*._B|_2022_.*|.*.I|_F.*.|_M.*.","",CD4_ADT@meta.data$orig.ident)
CD4_ADT@meta.data$gender <- gsub(".*._M_.*.","M",CD4_ADT@meta.data$orig.ident) %>% gsub(".*._F_.*.","F",.)
CD4_ADT@meta.data$Age <- gsub(".*._OO_.*.","OO",CD4_ADT@meta.data$orig.ident) %>% gsub(".*._O_.*.","O",.)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_ADT_merging.R")
objname = "CD4"
process = "integrated"
Assay = "ADT"

CD4_ADT_integrated <- ADT_merging(CD4_ADT, savedir, dims = 6, numfeatures=36,
                                  Assay=Assay, process=process, objname=objname,
                                  sample_tree = NULL, split_by = "Run",
                                  reference=NULL)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITESeq_RNA_intergation.R")
objname = "CD4"
process = "integration"
Assay = "ADT"
CD4_ADT_integrated <- RNA_integration(obj_path = CD4_ADT_integrated, saveDir = savedir, dims = 10,
                                      RNA_features = c("CD4-protein","CD8a-protein"),Assay = Assay,
                                      process = process, objname = objname, ncol = 4)

### Nearest Neighbour and Clustering
CD4_ADT_integrated_NN <- FindNeighbors(CD4_ADT_integrated, dims = 1:10)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/cluster_UMAP_QC.R")
res = c(0.2,0.4,0.6,0.8)
for (i in 1:length(res)) {
  process = paste("ADT_UMAP_QC",res[i],sep = "_")
  Assay = "integrated"
  CD4_ADT_integrated_NN_cluster <- cluster_UMAP_and_QC(obj_path = CD4_ADT_integrated_NN, dims = 10, res = res[i], saveDir=savedir, Assay = Assay,
                                                       QC_features = c("nCount_ADT","nFeature_ADT"), objname = "CD4_ADT_integrated_NN_cluster",
                                                       process = process,col_sel = c("gender","Age","vaccine","condition","Run","condition","orig.ident"))

}

saveRDS(CD4_ADT_integrated_NN_cluster,paste(savedir,"saveRDS_obj/CD4_ADT_integrated_NN_cluster.RDS",sep = ""))

#### Modality Integration ######
source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/scCITEseq_modality_integration.R")
# Find the last save object for RNA and ADT
# we have RNA already in the path
RNA_obj_path <- CD4_integrated_RNA_NN_cluster
ADT_obj_path <- CD4_ADT_integrated_NN_cluster

objname <- "CD4"
process <- "modality_integrate"
Assay <- "integrated"

CD4_modality_integrate <- modality_integration(RNA_obj_path = RNA_obj_path, ADT_obj_path = ADT_obj_path,
                                               RNA_dims = 30, ADT_dims = 10,
                                               saveDir = savedir, Assay = Assay, process = process, objname = objname)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/Modal_cluster_and_UMAP.R")
obj_path <- CD4_modality_integrate
RNA_features = c("CD4","CD8A")
ADT_features = c("CD4-protein","CD8a-protein")

res = c(0.2,0.4,0.6,0.8,1)
res=0.4
for (i in 1:length(res)) {
  objname = "CD4"
  Assay = "integrated"
  process = paste("modality_cluster_and_QC",res[i],sep = "_")
  CD4_modality_integrate_cluster <- Modal_Cluster_and_UMAP_QC(modal_integrate_obj_path = obj_path, saveDir = savedir,
                                                              res = res[i], RNA_features = RNA_features,
                                                              ADT_features = ADT_features,Assay = Assay, protein_assay = "ADT",
                                                              process = process, objname=objname, col_sel=c("gender","Age","vaccine","condition","Run","condition","orig.ident"))
}

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/scCITESeq/pipeline_functions/ADT_featureplot.R")
ADT_featureplot(obj = CD4_modality_integrate_cluster, savedir = savedir,objname = "CD4_modality_integrate_cluster",
                process = "featureplots",reduction = "wnn.umap",x = "wnnUMAP_1",y = "wnnUMAP_2")


CD4_samplewise_counts <- table(CD4_modality_integrate_cluster@meta.data$orig.ident) %>% as.data.frame()
write.table(CD4_samplewise_counts, paste(savedir,"Table/CD4_samplewise_counts.txt",sep = ""), sep = "\t",
            row.names = F, col.names = T, quote = F)

saveRDS(CD4_modality_integrate_cluster,paste(savedir,"saveRDS_obj/CD4_modality_integrate_cluster_0.6.RDS",sep = ""))

pdf(paste(savedir,"UMAP/CD4_modality_integrated_origident_splitted.pdf",sep = ""),width = 10, height = 8)
print(DimPlot(CD4_modality_integrate_cluster, split.by = "orig.ident", group.by="orig.ident", ncol = 3, reduction = "wnn.umap"))
dev.off()

pdf(paste(savedir,"UMAP/CD4_modality_integrated_0.8.pdf",sep = ""),width = 6, height = 6)
print(DimPlot(CD4_modality_integrate_cluster, reduction = "wnn.umap", label.size = 10, label = TRUE))
dev.off()
