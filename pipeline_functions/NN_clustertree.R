nearest_neigbour <- function(Obj_path,Dims,saveDir, Assay = "RNA", samplename, process){
  message(paste("\n\n######################################## Processsing",samplename,"#####################################################\n\n"))
  GEX_obj_3 <- Obj_path
  DefaultAssay(GEX_obj_3) <- Assay
  GEX_obj_3 <- FindNeighbors(GEX_obj_3, dims = 1:Dims)
  
  rm(forres)
  forres <- GEX_obj_3
  message(paste("\nRunning Cluster Tree for",samplename,"\n"))
  resolution <- c(0.001,0.01,0.03,0.05,0.07,0.09,0.1,0.125,0.15,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2)
  for(j in 1:length(resolution)){
    forres <- FindClusters(forres, resolution = resolution[j], verbose = FALSE)
  }
  dir.create(paste(saveDir,"clustertree/",sep = ""), showWarnings = FALSE)
  library(clustree)
  pdf(paste(saveDir,"clustertree/",samplename,"_",Assay,"_",process,".pdf",sep = ""), width = 7, height = 12)
  print(clustree(forres))
  dev.off()
  #message("Saving RDS file")
  #saveRDS(GEX_obj_3,paste(saveDir,"saveRDS_obj/",samplename,"_",process,"_",Assay,".RDS",sep = ""))
  #message("\n\nplease check for the resolution using clustertree is at this location: ",saveDir,"clustertree/\n\n")
  return(GEX_obj_3)
}
