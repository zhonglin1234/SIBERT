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


all.lines=get_lines(loc.raw,grid.type='hex',max_gap = 7,aa=0.2)
plot_lines(loc.raw,all.lines,grid.type='hex')


margin.spots=c()

for (i in 1:length(all.lines)) {
  tmp.lines=all.lines[[i]]
  margin.spots=c(margin.spots,unlist(lapply(tmp.lines,function(x) rownames(x)[c(1,nrow(x))])))
}
margin.spots=unique(margin.spots)

plot(loc.raw)
points(loc.raw[margin.spots,],col='red')

NN=nrow(sp.counts.r)
YY=rowSums(sp.counts.r[,ex_markers])
XX = cbind(rep(1, NN), rep(NA, NN),rowSums(sp.counts.r[,tcell_markers]))
rownames(XX)=names(YY)

all.output.tcell=get_cps_and_save(all.lines,
                                  YY=YY,
                                  XX=XX,
                                  K=0,
                                  alpha_vec=seq(0.1,1,length=10),
                                  method_cp='PG',
                                  dat_name='data/melanoma1_tcelladj',
                                  grid.type='hex')


load("data/melanoma1_tcelladj_exhst_cps_penalty_10_20_30_40_50_60_70_80_90_100.Rdata")

for (i in c(10)) {
  tmp.output=calculate_slopes(all_segs=all.segs.list[[i]],
                              all.lines=all.lines,
                              y_markers=ex_markers,
                              sp.counts.r,
                              method_beta='PG',
                              alpha=alpha_vec[i],
                              dat_name='data/melanoma1')
}



load("data/melanoma1_tcelladj_exhst_cps_penalty_10_20_30_40_50_60_70_80_90_100.Rdata")
load("data/melanoma1_exhst_slopes_penalty100.Rdata")

lambda.list=seq(1,10,length=5)
j2.dat=plot_J2(J1=J1,lambda.list=lambda.list)
mean.list=j2.dat[[1]]
n.add=j2.dat[[2]]
new.nb.list==j2.dat[[3]]

plot(mean.list, type = 'l', ylim = c(0, 20), xaxt = "n", xlab = "J2",ylab='% spots with Pr(TRR)>0.7')  # Suppress default x-axis
axis(1, at = seq_along(mean.list), labels = round(lambda.list,1))  # Add custom x-axis
#Find the transition point of the means and define it as J2


ss.mat=gibbs_sampling_new(ss=ss,
                          spot_names=names(obs_ss),
                          n_add=n.add*J2,
                          nb_list=nb.r,
                          new_nb_list=new.nb.list,
                          J=J1,
                          loc=loc.raw,
                          iterations = 5)
ss.mat[ss.mat==-1]=0
final.probs=calculate_final_prob(ss.mat, loc=loc.dat, burn_in = 2000,J=J1)



final.regions=get_TRR(final.probs)




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

