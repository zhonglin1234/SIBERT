library(Seurat)
library(SeuratDisk)
library(Matrix)
library(DropletUtils)

mat = ReadMtx(
  mtx = "/project/DPDS/YYang_lab/shared/sibert/SIBERT/data/visium_hd/pseudo_visium/filtered_feature_bc_matrix/matrix.mtx.gz",
  features = "/project/DPDS/YYang_lab/shared/sibert/SIBERT/data/visium_hd/pseudo_visium/filtered_feature_bc_matrix/features.tsv.gz",
  cells = "/project/DPDS/YYang_lab/shared/sibert/SIBERT/data/visium_hd/pseudo_visium/filtered_feature_bc_matrix/barcodes.tsv.gz"
)


write10xCounts(
  path = "/project/DPDS/YYang_lab/shared/sibert/SIBERT/data/visium_hd/pseudo_visium//filtered_feature_bc_matrix.h5",
  x = mat,
  type = "HDF5"
)

obj = Load10X_Spatial(data.dir = "/project/DPDS/YYang_lab/shared/sibert/SIBERT/data/visium_hd/pseudo_visium/")



umi <- Matrix::colSums(obj@assays$Spatial$counts)
sum(umi==0)

