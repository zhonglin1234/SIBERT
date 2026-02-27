#' @title Data Preparation Functions for SIBERT
#' @description Functions for loading and preparing spatial transcriptomics data
#' @name data_preparation
NULL

#' Read Spatial Data from 10x Visium
#'
#' Loads 10x Visium spatial transcriptomics data into a Seurat object and
#' performs initial processing including normalization, PCA, and clustering.
#'
#' @param folder Path to the folder containing the Visium data
#' @param filenm Name of the HDF5 file (e.g., "filtered_feature_bc_matrix.h5")
#' @param resolut Resolution parameter for Seurat clustering (default: 0.3)
#' @param loc_file Path to the tissue positions CSV file
#' @param n.vgg Number of variable genes for SCTransform (default: 3000)
#'
#' @return A list containing:
#' \itemize{
#'   \item sp.counts.r - Raw count matrix (spots x genes)
#'   \item sp.counts.norm - Normalized count matrix
#'   \item sp.counts.scale - Scaled count matrix
#'   \item meta_dat - Seurat metadata
#'   \item loc.raw - XY coordinates of spots
#' }
#'
#' @importFrom Seurat Load10X_Spatial SCTransform RunPCA FindNeighbors FindClusters RunUMAP DimPlot
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- read.spatial.dat(
#'   folder = "data/visium_sample",
#'   filenm = "filtered_feature_bc_matrix.h5",
#'   resolut = 0.2,
#'   loc_file = "data/visium_sample/spatial/tissue_positions.csv"
#' )
#' }
read.spatial.dat <- function(folder, filenm, resolut = 0.3, loc_file, n.vgg = 3000) {
  spdat <- Seurat::Load10X_Spatial(folder,
                                    filename = filenm,
                                    assay = "Spatial",
                                    slice = "slice1",
                                    filter.matrix = TRUE,
                                    to.upper = FALSE,
                                    image = NULL)

  spdat <- Seurat::SCTransform(spdat, assay = 'Spatial', variable.features.n = n.vgg)
  spdat <- Seurat::RunPCA(spdat)
  spdat <- Seurat::FindNeighbors(spdat, dims = 1:30)
  spdat <- Seurat::FindClusters(spdat, resolution = resolut)
  spdat <- Seurat::RunUMAP(spdat, dims = 1:30)
  Seurat::DimPlot(spdat, reduction = 'umap', label = TRUE, repel = TRUE)

  sp.counts.r <- t(as.matrix(spdat@assays$Spatial$counts))
  sp.counts.norm <- t(as.matrix(as.matrix(spdat@assays$SCT@data)))
  sp.counts.r <- sp.counts.r[which(rownames(sp.counts.r) %in% rownames(sp.counts.norm)),
                              which(colnames(sp.counts.r) %in% colnames(sp.counts.norm))]
  meta_dat <- spdat@meta.data

  sp.counts.norm <- sp.counts.norm[rownames(meta_dat), ]
  sp.counts.r <- sp.counts.r[rownames(meta_dat), ]
  sp.counts.scale <- t(spdat@assays$SCT$scale.data)
  sp.counts.scale <- sp.counts.scale[rownames(meta_dat), ]

  row.sum <- rowSums(sp.counts.r)
  col.sum <- colSums(sp.counts.r)
  message(paste("Minimum number of transcripts per spot:", min(row.sum)))
  message(paste("Minimum number of transcripts per gene:", min(col.sum)))

  # Load location data
  tmp.loc <- utils::read.csv(loc_file, sep = ',', header = FALSE)
  tmp.loc <- tmp.loc[which(tmp.loc[, 1] %in% rownames(sp.counts.r)), ]
  rownames(tmp.loc) <- tmp.loc[, 1]
  tmp.loc <- tmp.loc[, 5:6]
  colnames(tmp.loc) <- c('x', 'y')
  loc.raw <- data.frame(tmp.loc)

  loc.raw$x <- as.numeric(loc.raw$x)
  loc.raw$y <- as.numeric(loc.raw$y)
  loc.raw <- loc.raw[rownames(meta_dat), ]

  out <- list(
    sp.counts.r = sp.counts.r,
    sp.counts.norm = sp.counts.norm,
    sp.counts.scale = sp.counts.scale,
    meta_dat = meta_dat,
    loc.raw = loc.raw
  )

  return(out)
}


#' Data Preparation for TRR Analysis
#'
#' Main data preparation function that loads Visium data and computes
#' exhaustion scores adjusted for T cell marker expression.
#'
#' @param folder Path to the folder containing the Visium data
#' @param filenm Name of the HDF5 file
#' @param loc_file Path to the tissue positions CSV file
#' @param res Resolution for Seurat clustering (default: 0.2)
#' @param tcell_markers Character vector of T cell marker genes
#' @param ex_markers Character vector of exhaustion marker genes
#' @param cutf Cutoff parameter (default: 100)
#'
#' @return A list containing:
#' \itemize{
#'   \item sp.count.r - Raw count matrix
#'   \item sp.counts.norm - Normalized count matrix
#'   \item sp.counts.tpm - TPM normalized counts
#'   \item vgg - Variable genes
#'   \item loc.raw - Spot coordinates
#'   \item meta_dat - Metadata with exhaustion scores
#'   \item r.unit - Unit distance between adjacent spots
#' }
#'
#' @importFrom MASS glm.nb
#' @importFrom stats residuals
#' @export
#'
#' @examples
#' \dontrun{
#' tcell_markers <- c('CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'CD8B2')
#' ex_markers <- c('PDCD1', 'LAG3', 'HAVCR2', 'TIGIT', 'ENTPD1', 'CTLA4', 'TOX')
#'
#' dat <- data_preparation_TRR(
#'   folder = "data/visium_sample",
#'   filenm = "filtered_feature_bc_matrix.h5",
#'   loc_file = "data/visium_sample/spatial/tissue_positions.csv",
#'   res = 0.2,
#'   tcell_markers = tcell_markers,
#'   ex_markers = ex_markers
#' )
#' }
data_preparation_TRR <- function(folder, filenm, loc_file, res = 0.2,
                                  tcell_markers, ex_markers, cutf = 100) {
  # Load data
  dat1 <- data_preparation(folder = folder,
                           filenm = filenm,
                           res = res,
                           loc_file = loc_file,
                           ex.list = ex_markers,
                           cell = tcell_markers,
                           cell_cut = 1,
                           cut.t = 1e-10,
                           n.vgg = 3000)

  sp.counts.r <- dat1[[1]]
  sp.counts.norm <- dat1[[2]]
  sp.counts.tpm <- dat1[[12]]
  vgg <- dat1[[3]]
  loc.raw <- dat1[[4]]
  meta_dat <- dat1[[5]]

  spot.sdist <- as.matrix(stats::dist(loc.raw))
  r.unit <- min(spot.sdist[which(spot.sdist != 0)]) * 1.02

  # Get exhaustion score adjusted for count of T cell markers from NB regression
  NN <- nrow(sp.counts.r)
  YY <- rowSums(sp.counts.r[, ex_markers])
  XX <- cbind(rep(1, NN), rep(NA, NN), rowSums(sp.counts.r[, tcell_markers]))
  rownames(XX) <- names(YY)
  dat <- data.frame(x = XX[, 3], y = YY)
  mm <- MASS::glm.nb(y ~ x, data = dat)
  yy_adjusted <- stats::residuals(mm, type = 'pearson')
  names(yy_adjusted) <- names(YY)
  meta_dat$ex.score.nb <- yy_adjusted[rownames(meta_dat)]

  return(list(
    sp.count.r = dat1[[1]],
    sp.counts.norm = dat1[[2]],
    sp.counts.tpm = sp.counts.tpm,
    vgg = dat1[[3]],
    loc.raw = dat1[[4]],
    meta_dat = meta_dat,
    r.unit = r.unit
  ))
}


#' Plot Seurat Clusters
#'
#' Visualizes Seurat clusters and computes mean expression of cancer markers
#' across clusters to help identify tumor regions.
#'
#' @param loc Spot coordinates data frame
#' @param meta_dat Seurat metadata
#' @param sp.counts.norm Normalized count matrix
#' @param sp.counts.r Raw count matrix
#' @param cancer.markers Character vector of cancer marker genes
#'
#' @return A list containing median and mean expression of cancer markers per cluster
#'
#' @export
plot_seurat_clusters <- function(loc, meta_dat, sp.counts.norm, sp.counts.r, cancer.markers) {
  loc <- loc[rownames(meta_dat), ]
  loc[, 'x'] <- as.numeric(loc[, 'x'])
  loc[, 'y'] <- as.numeric(loc[, 'y'])
  tmp.clusters <- meta_dat$seurat_clusters
  number.clusters <- length(table(tmp.clusters))
  col.list <- c('darkblue', 'blue', 'lightblue', 'cyan', 'green',
                'yellow', 'orange', 'chocolate', 'brown', 'red')

  graphics::plot(loc, col = 'white',
                 xlim = c(min(loc[, 'x']) * 0.9, max(loc[, 'x']) * 1.2),
                 main = 'Seurat clusters')
  for (i in 1:number.clusters) {
    graphics::points(loc[tmp.clusters == (i - 1), 'x'],
                     loc[tmp.clusters == (i - 1), 'y'],
                     col = col.list[i], pch = 19, cex = 0.5)
  }
  graphics::legend(x = max(loc[, 'x']) * 1.05, y = median(loc[, 'y']),
                   legend = 0:(number.clusters - 1),
                   col = col.list[1:number.clusters], pch = 19, cex = 0.7, bty = 'n')

  tmp <- clustermeans(count_dat = sp.counts.norm, meta_dat = meta_dat, cancer.markers = cancer.markers)
  message('Median expression of each cancer marker in normalized data:')
  print(tmp[[1]])
  message('Sum of median expression of each cancer marker in each cluster:')
  print(rowSums(tmp[[1]]))

  tmp <- clustermeans(count_dat = sp.counts.r, meta_dat = meta_dat, cancer.markers = cancer.markers)
  message('Median expression of each cancer marker in raw count data:')
  print(tmp[[1]])
  message('Sum of median expression per cluster (raw data):')
  print(rowSums(tmp[[1]]))

  return(tmp)
}


# Internal helper function for cluster means
clustermeans <- function(count_dat, meta_dat, cancer.markers) {
  tmp.markers <- cancer.markers[which(cancer.markers %in% colnames(count_dat))]
  tmp.dat <- count_dat[, tmp.markers]
  meta_dat <- meta_dat[rownames(tmp.dat), ]

  med <- apply(tmp.dat, 2, function(x)
    stats::aggregate(x ~ meta_dat$seurat_clusters, data = meta_dat, median)[, 2])
  mm <- apply(tmp.dat, 2, function(x)
    stats::aggregate(x ~ meta_dat$seurat_clusters, data = meta_dat, mean)[, 2])

  return(list(med, mm))
}


# Internal data preparation function
data_preparation <- function(folder, filenm, loc_file, res, ex.list, cell,
                             cell_cut, cut.t = 1e-10, n.vgg, clustering = FALSE) {
  # Load data
  seurat_items <- read.spatial.dat(folder = folder, filenm = filenm,
                                   resolut = res, loc_file = loc_file, n.vgg = n.vgg)
  sp.counts.r <- seurat_items[[1]]
  sp.counts.norm <- seurat_items[[2]]
  meta_dat <- seurat_items[[4]]
  loc.raw <- seurat_items[[5]]
  vgg <- colnames(seurat_items[[3]])

  sp.counts.tpm <- apply(sp.counts.r, 2, function(x)
    as.numeric(x) / meta_dat[rownames(sp.counts.r), 'nCount_Spatial'] * 1e6)
  rownames(sp.counts.tpm) <- rownames(sp.counts.r)

  # Add exhaustion score to metadata
  meta_dat <- add_sum_to_meta(ex.list = ex.list, cell = cell,
                              sp.counts.r = sp.counts.r, meta_dat = meta_dat)

  # Define T cell spots
  keep.ids <- rownames(meta_dat[which(meta_dat$cell.sum >= cell_cut), ])
  message(paste("Number of T cell spots:", length(keep.ids)))

  vgg.mean <- apply(sp.counts.r[keep.ids, vgg], 2, mean)
  vgg <- vgg[which(vgg.mean > 0.1)]

  # Add exhaustion score
  meta_dat <- add_pt_to_meta(keep.ids = keep.ids, meta_dat = meta_dat)
  meta_dat_cell <- meta_dat[keep.ids, ]
  pt <- meta_dat_cell$ex.score
  names(pt) <- rownames(meta_dat_cell)

  # Rescale coordinates
  loc.rsc <- rescale_loc(loc = loc.raw, pt = pt, aa = 1)
  loc1 <- loc.rsc[keep.ids, ]

  # Distance matrices
  spot.pt.dist <- as.matrix(stats::dist(pt))
  rownames(spot.pt.dist) <- rownames(loc1)
  colnames(spot.pt.dist) <- rownames(loc1)

  spot.sdist <- as.matrix(stats::dist(loc1, diag = TRUE))
  r.unit <- min(spot.sdist[which(spot.sdist != 0)]) * 1.015

  # Subclustering
  if (clustering) {
    sub.clusters <- assign_subcluster(cut.t = cut.t, nr = 1.1, spots = keep.ids,
                                      spot.sdist = spot.sdist, r.unit = r.unit,
                                      spot.pt.dist = spot.pt.dist)
  } else {
    sub.clusters <- as.list(keep.ids)
  }

  ll <- unlist(lapply(sub.clusters, length))
  message('Number of spots in each node:')
  print(table(ll))

  meta.cts <- meta_subclusters(sub.clusters = sub.clusters, loc = loc1,
                               pt = pt, meta_dat = meta_dat, keep.ids = keep.ids)
  names(sub.clusters) <- rownames(meta.cts)

  loc.cts <- meta.cts[, c('x', 'y')]
  pt.cts <- meta.cts$pt
  names(pt.cts) <- rownames(meta.cts)

  output <- list(
    sp.counts.r, sp.counts.norm, vgg, loc.raw, meta_dat,
    keep.ids, meta.cts, loc.cts, pt.cts, r.unit, sub.clusters, sp.counts.tpm
  )

  names(output) <- c('raw count data', 'normalized count data', 'Top variable genes',
                     'XY coordinate of the raw data', 'meta data of the spot data',
                     'Barcodes of T cell spots', 'meta data of the node data',
                     'Rescaled XY coordinate of the node data', 'Exhaustion scores of the nodes',
                     'Distance between adjacent nodes', 'Spots in nodes', 'TPM data')

  return(output)
}


# Internal helper functions
add_sum_to_meta <- function(ex.list, cell, sp.counts.r, meta_dat) {
  ex.sum <- rowSums(sp.counts.r[, ex.list])
  cell.sum <- rowSums(sp.counts.r[, cell])
  
  print(paste('Sum of ',paste(cell,collapse = ',')))
  print(quantile(cell.sum,probs=seq(0,1,length=21)))
  
  print(paste('Sum of ',paste(ex.list,collapse = ',')))
  print(quantile(round(ex.sum),probs=seq(0,1,length=201)))
  
  names(cell.sum) <- rownames(sp.counts.r)
  names(ex.sum) <- rownames(sp.counts.r)

  meta_dat$cell.sum <- cell.sum[rownames(meta_dat)]
  meta_dat$ex.sum <- ex.sum[rownames(meta_dat)]

  return(meta_dat)
}

add_pt_to_meta <- function(keep.ids, meta_dat) {
  meta_dat_cell <- meta_dat[keep.ids, ]
  tmp <- stats::lm(ex.sum ~ cell.sum, data = meta_dat_cell)
  ex.score <- tmp$residuals
  ex.score <- ex.score + abs(min(ex.score))
  meta_dat$ex.score <- NA
  meta_dat[keep.ids, 'ex.score'] <- ex.score
  return(meta_dat)
}

rescale_loc<-function(keep.ids=keep.ids,loc,pt,aa=1){ #aa act as a weight of physical distance
  loc=data.frame(loc[keep.ids,])
  loc$x=as.numeric(loc$x)
  loc$y=as.numeric(loc$y)
  
  loc$x=loc$x-min(loc$x)
  loc$y=loc$y-min(loc$y)
  loc$x=loc$x*max(pt)*aa/max(loc$x)
  loc$y=loc$y*max(pt)*aa/max(loc$y)
  plot(loc)
  return(loc)
}

assign_subcluster <- function(cut.t, spots, spot.sdist, nr, r.unit, spot.pt.dist) {
  tmp.list <- list()
  for (i in spots) {
    nbs <- names(which(spot.sdist[i, ] < nr * r.unit))
    if (length(nbs) == 0) {
      message('Error: no neighbors found')
    }
    tmp.clus <- nbs[which(spot.pt.dist[i, nbs] < cut.t)]
    tmp.list <- c(tmp.list, list(tmp.clus))
  }
  tmp.list <- unique(tmp.list)

  ll <- unlist(lapply(tmp.list, length))
  if (max(ll) <= 2) {
    return(tmp.list)
  }
  if (max(ll)>2) {
    final.list=c()
    for (j in 2:(max(ll)-1)){
      s.clus=tmp.list[ll==j]
      b.clus=tmp.list[ll>j]
      max.overlap=unlist(lapply(s.clus,function(x) max(n.overlap(v1=x,v.list=b.clus))))
      if (max(max.overlap)==j) {
        s.clus=s.clus[-which(max.overlap==j)]
      }
      final.list=c(final.list,s.clus)
    }
    final.list=c(tmp.list[ll==1],final.list,tmp.list[ll==max(ll)])
    print(sum(unlist(lapply(final.list,length))))
    print(table(unlist(lapply(final.list,length))))
    return(final.list)
  }
}

meta_subclusters<-function(sub.clusters,loc,pt,meta_dat,keep.ids){
  meta_dat_cell=meta_dat[keep.ids,]
  
  loc.cts=do.call(rbind,lapply(sub.clusters,function(x) colSums(loc[x,])/length(x)))
  rownames(loc.cts)=paste0('cts',1:nrow(loc.cts))
  colnames(loc.cts)=c('x','y')
  pt.cts=unlist(lapply(sub.clusters,function(x) mean(pt[x])))
  names(pt.cts)=rownames(loc.cts)
  
  tmp.cluster=lapply(sub.clusters,function(x) meta_dat_cell[x,'seurat_clusters'])
  tmp.cluster=lapply(tmp.cluster,function(x) prop.table(table(x)))
  tmp.cluster=unlist(lapply(tmp.cluster,function(x) names(x)[which(x==max(x))][1]))
  cts.meta=data.frame(loc.cts,pt.cts,tmp.cluster)
  colnames(cts.meta)=c('x','y','pt','seurat_clusters')
  return(cts.meta)
}
