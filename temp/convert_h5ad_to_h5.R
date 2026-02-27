# Use R 4.2
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(DropletUtils)
library(dplyr)
library(ggplot2)


data_folder = "/project/DPDS/YYang_lab/shared/sibert/SIBERT/data/visium_hd/"
hd_files = c()
for(f in list.files(data_folder))
{
  if(endsWith(f, "T_cell"))
  {
    hd_files = c(hd_files, f)
  }
}

for(f in hd_files[4:7])
{
  mat = ReadMtx(
    mtx = paste0(data_folder, f, "/filtered_feature_bc_matrix/matrix.mtx.gz"),
    features = paste0(data_folder, f, "/filtered_feature_bc_matrix/features.tsv.gz"),
    cells = paste0(data_folder, f, "/filtered_feature_bc_matrix/barcodes.tsv.gz")
  )
  
  
  write10xCounts(
    path = paste0(data_folder, f,"/filtered_feature_bc_matrix.h5"),
    x = mat,
    type = "HDF5"
  )
}


# obj = Load10X_Spatial(data.dir = paste0(data_folder, f))

for(f in hd_files)
{
  tissue_pos = read.csv(paste0(data_folder, f, "/spatial/tissue_positions.csv"))
  ct_comp = read.csv(paste0(data_folder, f, "/spot_cell_type_composition.csv"))
  df = merge(tissue_pos, ct_comp, by.x="barcode", by.y='spot_id')
  df$n_cells = df$Other + df$T + df$Tumor
  df$Tumor_frac = df$Tumor/df$n_cells
  df$T_frac = df$T/df$n_cells
  df$Other_frac = df$Other/df$n_cells
  df$Tumor_frac[which(is.na(df$Tumor_frac))] = 0
  df$T_frac[which(is.na(df$T_frac))] = 0
  df$Other_frac[which(is.na(df$Other_frac))] = 0
  
  df2 <- df %>%
    mutate(
      Tumor_frac = pmax(pmin(Tumor_frac, 1), 0),
      T_frac     = pmax(pmin(T_frac, 1), 0),
      Other_frac = pmax(pmin(Other_frac, 1), 0),
      s = Tumor_frac + T_frac + Other_frac,
      Tumor_frac = ifelse(s > 0, Tumor_frac / s, 0),
      T_frac     = ifelse(s > 0, T_frac / s, 0),
      Other_frac = ifelse(s > 0, Other_frac / s, 0),
      
      # base colors (0-1)
      r = 1.00 * Tumor_frac + 0.60 * Other_frac + 0.00 * T_frac,
      g = 0.00 * Tumor_frac + 0.60 * Other_frac + 0.00 * T_frac,
      b = 0.00 * Tumor_frac + 0.60 * Other_frac + 1.00 * T_frac,
      
      # transparency: empty spots invisible; optionally also fade low n_cells
      alpha = ifelse(n_cells <= 0, 0, 1),
      
      mix_col = rgb(r, g, b, alpha = alpha)
    )
  
  ggplot(df2, aes(x = pxl_row_in_fullres, y = pxl_col_in_fullres)) +
    geom_point(color = df2$mix_col, size = 1.8) +
    coord_fixed() +
    scale_y_reverse() +   
    theme_void()+
    ggtitle(f)
  
  ggsave(filename=paste0("/project/DPDS/YYang_lab/shared/sibert/misc/presentation/", f, ".png"),
         width = 9, height = 9)
}





legend_df <- tibble::tibble(
  label = c("Tumor", "T", "Other", "Tumor+T (50/50)", "All (1/3 each)"),
  col = c(
    rgb(1,0,0,1),
    rgb(0,0,1,1),
    rgb(0.6,0.6,0.6,1),
    rgb(0.5,0,0.5,1),
    rgb((1+0+0.6)/3, (0+0+0.6)/3, (0+1+0.6)/3, 1)
  )
)

ggplot(legend_df, aes(x=1, y=label)) +
  geom_point(aes(color = col), size=20) +
  scale_color_identity() +
  theme_minimal() +
  labs(x=NULL, y=NULL)







# umi <- Matrix::colSums(obj@assays$Spatial$counts)
# sum(umi==0)

