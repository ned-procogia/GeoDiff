# Funtion extracts the variables for DE from the SeuratDsk

library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(HDF5Array)
library(hdf5r)

Seurat_GeoDiff<- function(SeuratDisk_file){
  # accepts the SeuratDisk_file as input 
  
  sem<- LoadH5Seurat("SeuratProject.h5seurat")
  # loading only the required slots from the Disk
  
  GeoDiff_var <- list(neg = Matrix(sem@assays$negprobes@counts, sparse=TRUE),
                      # negative probes
                      raw = Matrix(sem@assays$RNA@counts, sparse=TRUE),
                      # raw counts
                      annot = sem@misc$cell_metadata$rna
                      # cell meta data counts
  )
  # convert the matrices into a list of matrices 
  
  rm(sem)
  gc()
  return(GeoDiff_var)
  
}
