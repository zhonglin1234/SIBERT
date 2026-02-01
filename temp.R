library(hdf5r)
setwd("/project/DPDS/YYang_lab/shared/sibert/SIBERT")
func_dir <- file.path("functions")
func_files <- list.files(func_dir, pattern = "\\.R$", full.names = TRUE)
sapply(func_files, source)

data_dir <- file.path("data")
rdata_files <- list.files(data_dir, pattern = "\\.Rdata$", full.names = TRUE)
for (f in head(rdata_files, 2)) {
  message("Loading: ", f)
  load(f)
}


tcell_markers=c('CD3D','CD3E','CD3G','CD4','CD8A','CD8B2')
ex_markers=c('PDCD1','LAG3','HAVCR2','TIGIT','ENTPD1','CTLA4','TOX')

dat1=data_preparation_TRR(folder='data/visium_melanoma1',  
                          filenm = "CytAssist_FFPE_Human_Skin_Melanoma_filtered_feature_bc_matrix.h5", 
                          res=0.2,
                          loc_file='data/visium_melanoma1/spatial/tissue_positions.csv',
                          ex_markers=ex_markers,
                          tcell_markers=tcell_markers
)


sp.counts.r=dat1[[1]]
sp.counts.norm=dat1[[2]]
sp.counts.tpm=dat1[[3]]
vgg=dat1[[4]]
loc.raw=dat1[[5]]
meta_dat=dat1[[6]]
r.unit=dat1[[7]]


tmp.mean <- plot_seurat_clusters(
  loc = loc.raw,
  meta_dat = meta_dat,
  sp.counts.norm = sp.counts.norm,
  sp.counts.r = sp.counts.r,
  cancer.markers = cancer.markers$melanoma
)
tumor.clus <- c(0,4,5,6) #Define clusters that are tumors here
tumor.loc= loc.raw[rownames(meta_dat)[which(meta_dat$seurat_clusters %in% tumor.clus)],] # Get coordinates of your tumor regions here

## A better idea is to draw the tumor region mannually based on margins of the seurat clusters with high cancer marker expression.
tmp.spots <- SpatialPoints(loc.raw)
boundary_points <- locator(type = "l")  # Click on the plot to create a boundary
boundary_polygon <- Polygon(cbind(boundary_points$x, boundary_points$y))
boundary_sp <- SpatialPolygons(list(Polygons(list(boundary_polygon), ID = "boundary")))
tumor.loc<- data.frame(tmp.spots[boundary_sp, ]) 
plot(loc.raw)
points(tumor.loc,col='red')



j1.dat=estimate_J1_dat(genes='B2M',sp.dd=sp.counts.r,adj.umi=F) #omega is 0.84

sum_obs=j1.dat[[2]]
j1.dat=j1.dat[[1]]

loc.dat=j1.dat[[2]]
obs_ss=j1.dat[[3]]
nb.r=j1.dat[[4]]
nb.info=j1.dat[[5]]

first_dr<-find_optimal_mu_Q(mu_Q_candidates=seq(0.16, 0.44, by = 0.01),
                            init_x=sample(c(-1,1),size=length(obs_ss),prob=c(0.5,0.5),replace = T)  , #initial states
                            nb.list.cpp=nb.info, #neighbors
                            mu=0.3, sigma2=0.0049,sigma2_Q=0.0001, #prior of J
                            K=10, s_obs=sum_obs , n_flip=60000, M=10000) # can set mu_Q_candidates=seq(0.35, 0.41, by = 0.01) here to save time, as omega is 0.84

J1=as.numeric(names(first_dr)[which(first_dr==min(first_dr))]) 







Sys.setenv(HDF5_DIR="/project/apps/apps_rhel9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.4.1/hdf5-1.14.5-gy4hlf4r4igw23wwhiqqnyua4ennth33")
Sys.setenv(LD_LIBRARY_PATH = paste0("/project/apps/apps_rhel9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.4.1/hdf5-1.14.5-gy4hlf4r4igw23wwhiqqnyua4ennth33/lib:",
                                    Sys.getenv("LD_LIBRARY_PATH")))
install.packages("hdf5r")


hdf5_path <- "/project/apps/apps_rhel9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.4.1/hdf5-1.14.5-gy4hlf4r4igw23wwhiqqnyua4ennth33"

install.packages("hdf5r", configure.args = paste0("--with-hdf5=", hdf5_path, "/bin/h5cc"))

Sys.setenv(LD_LIBRARY_PATH = paste0(
  "/project/apps/apps_rhel9/spack/opt/spack/linux-rhel9-x86_64/gcc-11.4.1/hdf5-1.14.5-gy4hlf4r4igw23wwhiqqnyua4ennth33/lib:",
  Sys.getenv("LD_LIBRARY_PATH")
))

