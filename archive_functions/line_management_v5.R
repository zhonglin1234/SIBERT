
#Data preparation
require(Seurat)
require(cpop)
require(MASS)
require(BayesLogit)
require(coda)
require(changepoint)
require(changepoint.np)
require(ggplot2)
library(igraph)

colfunc<-colorRampPalette(c("royalblue","springgreen","yellow",'red'))
#plot(rep(1,50),col=(colfunc(50)), pch=19, cex=.8,cex=2)


data_preparation_TRR<-function(folder,filenm,loc_file,res=0.2, tcell_markers,ex_markers,cutf=100){
  #Load data in
  dat1=data_preparation(folder=folder, 
                        filenm = filenm,
                        res=res,
                        loc_file=loc_file,
                        ex.list=ex_markers,
                        cell=tcell_markers,
                        cell_cut=1,
                        cut.t=1e-10,
                        n.vgg=3000)
  
  
  sp.counts.r=dat1[[1]]
  sp.counts.norm=dat1[[2]]
  sp.counts.tpm=dat1[[12]]
  
  vgg=dat1[[3]]
  loc.raw=dat1[[4]] 
  meta_dat=dat1[[5]] 
  
  spot.sdist=as.matrix(dist(loc.raw))
  r.unit=min(spot.sdist[which(spot.sdist!=0)])*1.02
  
  #Get exhaustion score adjusted for count of T cell markers from NB regression
  NN=nrow(sp.counts.r)
  YY=rowSums(sp.counts.r[,ex_markers])
  XX = cbind(rep(1, NN), rep(NA, NN),rowSums(sp.counts.r[,tcell_markers]))
  rownames(XX)=names(YY)
  dat = data.frame(x = XX[,3], y = YY)
  mm = glm.nb(y ~ x, data = dat)
  yy_adjusted = residuals.glm(mm, type = 'pearson')  # Adjusted for the number of T cell markers
  names(yy_adjusted) = names(YY)
  meta_dat$ex.score.nb=yy_adjusted[rownames(meta_dat)]
  
  return(list(sp.count.r=dat1[[1]],
              sp.counts.norm=dat1[[2]],
              sp.counts.tpm=sp.counts.tpm,
              vgg=dat1[[3]],
              loc.raw=dat1[[4]],
              meta_dat=meta_dat,
              r.unit=r.unit))
}



################################
#Get all lines
################################

group_close_coords_square <- function(coords, cutoff) { #coords has two columns, first being x and second being y
  # Create a copy of the original indices to map back later
  original_order <- order(coords)
  
  # Sort the coordinates and initialize grouping
  sorted_coords <- sort(coords)
  group <- rep(1, length(sorted_coords))
  
  # Apply grouping based on the cutoff
  for (i in 2:length(sorted_coords)) {
    if (abs(sorted_coords[i] - sorted_coords[i - 1]) >= cutoff) {
      group[i] <- group[i - 1] + 1
    } else {
      group[i] <- group[i - 1]
    }
  }
  
  # Map the group assignments back to the original order
  group_in_original_order <- rep(NA, length(group))
  group_in_original_order[original_order] <- group
  
  return(group_in_original_order)
}

group_close_coords_hex <- function(coords, cutoff) {
  original_order <- order(coords)
  sorted_coords <- sort(coords)
  group <- rep(1, length(sorted_coords))
  
  for (i in 2:length(sorted_coords)) {
    if (abs(sorted_coords[i] - sorted_coords[i - 1]) >= cutoff) {
      group[i] <- group[i - 1] + 1
    } else {
      group[i] <- group[i - 1]
    }
  }
  
  # Map group back to original order
  group_in_original_order <- rep(NA, length(group))
  group_in_original_order[original_order] <- group
  
  return(group_in_original_order)
}

get_lines <- function(loc, grid.type, aa = 0.2, max_gap = 10) { 
  # grid.type can be 'square' or 'hexagon'
  
  loc <- data.frame(loc[,1:2])
  dist.mat <- as.matrix(dist(loc))
  min.dist <- min(dist.mat[upper.tri(dist.mat)])  # Smallest distance between points
  
  # Function to split lines by max_gap
  split_by_gap <- function(line_data, dist_cutoff) {
    if (nrow(line_data) < 2) return(list())  # Return empty list if fewer than 2 points
    
    segments <- list()
    current_segment <- line_data[1, , drop = FALSE]
    
    for (i in 2:nrow(line_data)) {
      point_dist <- sqrt(sum((line_data[i, ] - line_data[i - 1, ])^2))
      if (point_dist > dist_cutoff) {
        if (nrow(current_segment) >= 2) {  # Ensure segment has at least 2 points
          segments <- append(segments, list(current_segment))
        }
        current_segment <- line_data[i, , drop = FALSE]
      } else {
        current_segment <- rbind(current_segment, line_data[i, , drop = FALSE])
      }
    }
    
    if (nrow(current_segment) >= 2) {  # Add the last segment if it has at least 2 points
      segments <- append(segments, list(current_segment))
    }
    
    return(segments)
  }
  
  if (grid.type == 'square') {
    
    # Group for square grid
    loc$vertical_group <- group_close_coords_square(loc[, 1], cutoff = min.dist * aa)
    loc$horizontal_group <- group_close_coords_square(loc[, 2], cutoff = min.dist * aa)
    
    # Extract and order vertical lines by y coordinate
    vertical_lines <- split(loc[, 1:2], loc$vertical_group)
    vertical_lines <- lapply(vertical_lines, function(df) df[order(df[, 2]), ])  # Order by y-coordinate
    
    # Split vertical lines by max_gap, ensure at least two spots
    vertical_lines <- unlist(lapply(vertical_lines, function(line) split_by_gap(line, max_gap * min.dist)), recursive = FALSE)
    vertical_lines <- vertical_lines[sapply(vertical_lines, nrow) >= 2]  # Filter to ensure at least 2 points
    names(vertical_lines) <- paste0("Vertical Line ", seq_along(vertical_lines))
    
    # Extract and order horizontal lines by x coordinate
    horizontal_lines <- split(loc[, 1:2], loc$horizontal_group)
    horizontal_lines <- lapply(horizontal_lines, function(df) df[order(df[, 1]), ])  # Order by x-coordinate
    
    # Split horizontal lines by max_gap, ensure at least two spots
    horizontal_lines <- unlist(lapply(horizontal_lines, function(line) split_by_gap(line, max_gap * min.dist)), recursive = FALSE)
    horizontal_lines <- horizontal_lines[sapply(horizontal_lines, nrow) >= 2]  # Filter to ensure at least 2 points
    names(horizontal_lines) <- paste0("Horizontal Line ", seq_along(horizontal_lines))
    
    out.lines <- list(vertical_lines = vertical_lines, horizontal_lines = horizontal_lines)
    return(out.lines)
  }
  
  if (grid.type == 'hex') {
    # Rotation matrices for 30 degrees and 150 degrees
    rotate_30 <- matrix(c(cos(pi / 3), -sin(pi / 3), sin(pi / 3), cos(pi / 3)), ncol = 2)
    rotate_150 <- matrix(c(cos(2 * pi / 3), -sin(2 * pi / 3), sin(2 * pi / 3), cos(2 * pi / 3)), ncol = 2)
    
    # Project coordinates onto new axes (30-degree and 150-degree lines)
    loc_rotated_30 <- as.matrix(loc) %*% rotate_30
    loc_rotated_150 <- as.matrix(loc) %*% rotate_150
    
    # Group along the 30-degree and 150-degree directions
    loc$group_30 <- group_close_coords_hex(loc_rotated_30[, 1], cutoff = min.dist * aa)  # Projected x-axis (30 degrees)
    loc$group_150 <- group_close_coords_hex(loc_rotated_150[, 1], cutoff = min.dist * aa)  # Projected x-axis (150 degrees)
    
    # Group along the 90-degree (vertical) direction
    loc$group_90 <- group_close_coords_hex(loc[, 1], cutoff = min.dist * aa)
    
    # Extract and order lines by y coordinate
    lines_30 <- split(loc[, 1:2], loc$group_30)
    lines_30 <- lapply(lines_30, function(df) df[order(df[, 2]), ])  # Order by y-coordinate
    
    # Split 30-degree lines by max_gap, ensure at least two spots
    lines_30 <- unlist(lapply(lines_30, function(line) split_by_gap(line, max_gap * min.dist)), recursive = FALSE)
    lines_30 <- lines_30[sapply(lines_30, nrow) >= 2]  # Filter to ensure at least 2 points
    
    lines_150 <- split(loc[, 1:2], loc$group_150)
    lines_150 <- lapply(lines_150, function(df) df[order(df[, 2]), ])  # Order by y-coordinate
    
    # Split 150-degree lines by max_gap, ensure at least two spots
    lines_150 <- unlist(lapply(lines_150, function(line) split_by_gap(line, max_gap * min.dist)), recursive = FALSE)
    lines_150 <- lines_150[sapply(lines_150, nrow) >= 2]  # Filter to ensure at least 2 points
    
    lines_90 <- split(loc[, 1:2], loc$group_90)
    lines_90 <- lapply(lines_90, function(df) df[order(df[, 2]), ])  # Order by y-coordinate
    
    # Split 90-degree lines by max_gap, ensure at least two spots
    lines_90 <- unlist(lapply(lines_90, function(line) split_by_gap(line, max_gap * min.dist)), recursive = FALSE)
    lines_90 <- lines_90[sapply(lines_90, nrow) >= 2]  # Filter to ensure at least 2 points
    
    out.lines <- list(lines_30 = lines_30, lines_90 = lines_90, lines_150 = lines_150)
    return(out.lines)
  }
}

plot_lines<-function(loc,lines.list,grid.type){
  
  if (grid.type=='square') {
    par(mfrow=c(1,2))
    
    vertical_lines=lines.list$vertical_lines
    plot(loc[,1:2],xlab = "X Coordinate", ylab = "Y Coordinate",pch=19, cex=.3,col='grey',
         main='Vertical lines')
    for (i in seq_along(vertical_lines)) {
      line_data <- vertical_lines[[i]]
      if (nrow(line_data) > 1) {
        lines(line_data[,1], line_data[,2], col = "blue", lwd = .9)
      }
    }
    
    
    horizontal_lines=lines.list$horizontal_lines
    plot(loc[,1:2],xlab = "X Coordinate", ylab = "Y Coordinate",pch=19, cex=.3,col='grey',
         main='Horizontal lines')
    for (i in seq_along(horizontal_lines)) {
      line_data <- horizontal_lines[[i]]
      if (nrow(line_data) > 1) {
        lines(line_data[,1], line_data[,2], col = "blue", lwd = .9)
      }
    }
  }
  
  if (grid.type=='hex') {
    par(mfrow=c(1,3))
    
    lines_30=lines.list$lines_30
    plot(loc[,1:2],xlab = "X Coordinate", ylab = "Y Coordinate",pch=19, cex=.3,col='grey',
         main='30 degree lines (relative to x axis)')
    for (i in seq_along(lines_30)) {
      line_data <- lines_30[[i]]
      if (nrow(line_data) > 1) {
        lines(line_data[,1], line_data[,2], col = "blue", lwd = .9)
      }
    }
    
    lines_90=lines.list$lines_90
    plot(loc[,1:2],xlab = "X Coordinate", ylab = "Y Coordinate",pch=19, cex=.3,col='grey',
         main='90 degree lines (relative to x axis)')
    for (i in seq_along(lines_90)) {
      line_data <- lines_90[[i]]
      if (nrow(line_data) > 1) {
        lines(line_data[,1], line_data[,2], col = "blue", lwd = .9)
      }
    }
    
    lines_150=lines.list$lines_150
    plot(loc[,1:2],xlab = "X Coordinate", ylab = "Y Coordinate",pch=19, cex=.3,col='grey',
         main='150 degree lines (relative to x axis)')
    for (i in seq_along(lines_150)) {
      line_data <- lines_150[[i]]
      if (nrow(line_data) > 1) {
        lines(line_data[,1], line_data[,2], col = "blue", lwd = .9)
      }
    }
    
  }
  
}

plot_lines_with_cps<-function(loc,lines.list,grid.type,all.cps){
  
  
  if (grid.type=='square') {
    par(mfrow=c(1,2))
    
    vertical_lines=lines.list$vertical_lines
    plot(loc[,1:2],xlab = "X Coordinate", ylab = "Y Coordinate",pch=19, cex=.3,col='grey',
         main='Vertical lines')
    for (i in seq_along(vertical_lines)) {
      line_data <- vertical_lines[[i]]
      if (nrow(line_data) > 1) {
        lines(line_data[,1], line_data[,2], col = "blue", lwd = .9)
      }
    }
    
    
    horizontal_lines=lines.list$horizontal_lines
    plot(loc[,1:2],xlab = "X Coordinate", ylab = "Y Coordinate",pch=19, cex=.3,col='grey',
         main='Horizontal lines')
    for (i in seq_along(horizontal_lines)) {
      line_data <- horizontal_lines[[i]]
      if (nrow(line_data) > 1) {
        lines(line_data[,1], line_data[,2], col = "blue", lwd = .9)
      }
    }
  }
  
  if (grid.type=='hex') {
    
    changep.loc.l = lapply(all.cps, function(x) lapply(x, function(y) y$change_points))
    changep.spot = lapply(1:length(all.lines), function(x) lapply(1:length(all.lines[[x]]), function(y) {
      #print(x)
      #print(y)
      all.lines[[x]][[y]][changep.loc.l[[x]][[y]],]
    }
    ))
    
    cps1=do.call(rbind,changep.spot[[1]])
    cps2=do.call(rbind,changep.spot[[2]])
    cps3=do.call(rbind,changep.spot[[3]])
    
    par(mfrow=c(1,3))
    
    lines_30=lines.list$lines_30
    plot(loc[,1:2],xlab = "X Coordinate", ylab = "Y Coordinate",pch=19, cex=.3,col='grey',
         main='30 degree lines (relative to x axis)')
    for (i in seq_along(lines_30)) {
      line_data <- lines_30[[i]]
      if (nrow(line_data) > 1) {
        lines(line_data[,1], line_data[,2], col = "blue", lwd = .9)
      }
    }
    points(cps1,col='red',cex=.4)
    
    lines_90=lines.list$lines_90
    plot(loc[,1:2],xlab = "X Coordinate", ylab = "Y Coordinate",pch=19, cex=.3,col='grey',
         main='90 degree lines (relative to x axis)')
    for (i in seq_along(lines_90)) {
      line_data <- lines_90[[i]]
      if (nrow(line_data) > 1) {
        lines(line_data[,1], line_data[,2], col = "blue", lwd = .9)
      }
    }
    points(cps2,col='red',cex=.4)
    
    lines_150=lines.list$lines_150
    plot(loc[,1:2],xlab = "X Coordinate", ylab = "Y Coordinate",pch=19, cex=.3,col='grey',
         main='150 degree lines (relative to x axis)')
    for (i in seq_along(lines_150)) {
      line_data <- lines_150[[i]]
      if (nrow(line_data) > 1) {
        lines(line_data[,1], line_data[,2], col = "blue", lwd = .9)
      }
    }
    points(cps3,col='red',cex=.4)
  }
  
}




### When get_lines() does not work-use get_lines_visium()


#### Angle of a line segment
plot_lines_one_dir<-function(tmp.lines,spot.ids) {
  plot(loc.raw[spot.ids,],cex=0.2)
  for (i in 1:length(tmp.lines)) {
    lines(loc.raw[tmp.lines[[i]],],col='blue')
  }
}


angle <- function(p1,p2){
  return( (p2[2]-p1[2])/(p2[1]-p1[1]) )
}


#Get overlapped pairs (to form a line)
get_all_overlapped_pairs<-function(tmp.spots,ang.pps,all.spots) {
  tmp.loc=c()
  for (i in all.spots) {
    tmp.loc=c(tmp.loc,list(which(ang.pps[,1]==i|ang.pps[,2]==i)))
  }
  names(tmp.loc)=all.spots
  
  length2=1
  length1=0
  while(length2>length1) {
    length1=length(tmp.spots)
    tmp.comb=unique(do.call(c,tmp.loc[tmp.spots]))
    tmp.spots=names(table(ang.pps[tmp.comb,]))
    length2=length(tmp.spots)
  }
  return(tmp.spots)
} 


#### Get lines of the same direction
get_sep_lines<-function(ang.pps,loc){
  all.spots=names(table(ang.pps))
  line.list=c(list(NULL))
  for (i in all.spots) {
    if (i %in% do.call(c,line.list)) {next}
    tmp.line=get_all_overlapped_pairs(tmp.spots=i,ang.pps=ang.pps,all.spots=all.spots)
    tmp.loc=loc[tmp.line,]
    tmp.order=order(tmp.loc[,'y'],decreasing=T)
    tmp.line=tmp.line[tmp.order]
    line.list=c(line.list,list(tmp.line))
  }
  return(line.list)
}

#Order spots on a line
order_spots<-function(line,loc){
  tmp.order=order(loc[line,2])
  return(line[tmp.order])
}


#### Connect everything
get_lines_visium<-function(ids,dist.mat,nr,length.cut,loc,allspots){
  r.unit=min(dist.mat[which(dist.mat!=0)])*nr
  print(r.unit)
  #Get all pairs of adjecent spots
  all.pairs=c()
  for (i in ids){
    tmp.nbs=get_nbs(spot=i,dist.mat=dist.mat,nr=nr,r.unit=r.unit)
    tmp.pairs=cbind(rep(i,length(tmp.nbs)),tmp.nbs)
    all.pairs=rbind(all.pairs,tmp.pairs)
  }
  
  all.pairs=t(apply(all.pairs, 1, sort))
  all.pairs=unique(all.pairs)
  
  print(paste('Number of pairs of spots:',nrow(all.pairs)))
  
  all.ang=c()
  for (i in 1:nrow(all.pairs)) {
    ag=angle(p1=loc.raw[all.pairs[i,1],],p2=loc.raw[all.pairs[i,2],])
    all.ang=c(all.ang,ag)
  }
  
  angs=all.ang
  angs[which(all.ang>0 & all.ang<1)]=1
  angs[which(all.ang<0 & all.ang>-1)]=2
  angs[which(!angs %in% c(1,2))]=3
  table(as.numeric(angs))
  
  ang1.pps=all.pairs[which(angs==1),]
  ang2.pps=all.pairs[which(angs==2),]
  ang3.pps=all.pairs[which(angs==3),]
  
  tmp=get_sep_lines(ang.pps=ang1.pps,loc=loc)
  lines1=tmp[2:length(tmp)]
  tmp=get_sep_lines(ang.pps=ang2.pps,loc=loc)
  lines2=tmp[2:length(tmp)]
  tmp=get_sep_lines(ang.pps=ang3.pps,loc=loc)
  lines3=tmp[2:length(tmp)]
  
  ll1=unlist(lapply(lines1,length))
  ll2=unlist(lapply(lines2,length))
  ll3=unlist(lapply(lines3,length))
  
  lines1=lines1[which(ll1>length.cut)]
  lines2=lines2[which(ll2>length.cut)]
  lines3=lines3[which(ll3>length.cut)]
  
  #lines1=connect_left(lines=lines1,allspots=allspots,loc=loc,r.unit=r.unit)
  #lines2=connect_left(lines=lines2,allspots=allspots,loc=loc,r.unit=r.unit)
  #lines3=connect_left(lines=lines3,allspots=allspots,loc=loc,r.unit=r.unit)
  
  #Order spots on each line
  lines1=lapply(lines1,function(x) order_spots(x,loc)) #30 degree
  lines2=lapply(lines2,function(x) order_spots(x,loc)) #150 degree
  lines3=lapply(lines3,function(x) order_spots(x,loc)) #90 degree
  all.lines=c(list(lines1),list(lines3),list(lines2))
  names(all.lines)=c('lines_30','lines_90','lines_150')
  
  return(all.lines)
}




############################
#Get change points for lines
############################

#Get change points for one line
get_line_cps_vectorized<-function(tmp.line,YY,XX,method,alpha_vec,K=0){
  yy=YY[tmp.line]
  xx=XX[tmp.line,]
  print(yy)
  print(xx)
  xx[,2]=1:length(yy)
  tmp.cp=optimal_partitioning_pelt_vectorized(n_iter=20,Y=yy,X=xx,method=method,alpha_vec=alpha_vec,K=K)
  return(tmp.cp)
}


#Get segments
get_all_segs<-function(all.cps,all.lines){
  all.segs=c()
  for (i in 1:length(all.cps)) {
    tmp.outcome=all.cps[[i]]
    tmp.lines=all.lines[[i]]
    tmp.segs=lapply(1:length(tmp.outcome),function(x) lapply(tmp.outcome[[x]]$segmentation, function(seg) {
      #print(seg)
      #print(tmp.lines[[x]])
      tmp.lines[[x]][seg[1]:seg[2],]
    }))
    all.segs=c(all.segs,list(tmp.segs))
  }
  return(all.segs)
}


# Get change points and segments for all lines
get_cps_and_save<-function(all.lines,YY,XX,K=0,alpha_vec,method_cp='PG',dat_name,
                           grid.type){
  
  # Remove lines with less than 4 spots (each segment must have at least 2 data points)
  lines.ll = lapply(all.lines, function(x) unlist(lapply(x, nrow)))  # length of each line
  for (i in 1:length(lines.ll)) {
    tmp.ll = lines.ll[[i]]
    if (min(tmp.ll) < 4) {
      all.lines[[i]] = all.lines[[i]][-which(tmp.ll < 4)]
    }
  }
  
  all.cps = list()  # Initialize to store the results for all lines
  
  # Use vectorized get_line_cps
  for (i in 1:length(all.lines)) {
    tmp.lines = all.lines[[i]]
    tmp.cps = lapply(tmp.lines, function(x) 
      get_line_cps_vectorized(tmp.line = rownames(x), YY = YY, XX = XX, method = method_cp, alpha_vec = alpha_vec, K = K))
    all.cps = c(all.cps, list(tmp.cps))
  }
  
  # Initialize an empty list to hold the rearranged structure with 11 penalties
  
  n_alpha=length(alpha_vec)
  # Initialize the first output structure for the rearranged `all.cps` by penalty
  all.cps.rearranged <- vector("list", n_alpha)
  
  
  if (grid.type=='square') {
    # Loop over each penalty index (1 to 11)
    for (penalty_index in 1:n_alpha) {
      # Extract vertical lines for the current penalty
      vertical_lines_penalty <- lapply(all.cps[[1]], function(line) line[[1]][[penalty_index]])
      
      # Extract horizontal lines for the current penalty
      horizontal_lines_penalty <- lapply(all.cps[[2]], function(line) line[[1]][[penalty_index]])
      
      # Combine vertical and horizontal results for the current penalty
      all.cps.rearranged[[penalty_index]] <- list(
        vertical_lines = vertical_lines_penalty,
        horizontal_lines = horizontal_lines_penalty
      )
    }
    
    all.cost.mat <- list(
      vertical_lines = lapply(all.cps[[1]], function(line) line[[2]]),  # Extract `cost.mat` for each vertical line
      horizontal_lines = lapply(all.cps[[2]], function(line) line[[2]])  # Extract `cost.mat` for each horizontal line
    )
    
  }

  
  if (grid.type=='hex') {

    # Loop over each penalty index (1 to 11)
    for (penalty_index in 1:n_alpha) {
      # Extract 30-degree lines for the current penalty
      lines_30_penalty <- lapply(all.cps[[1]], function(line) line[[1]][[penalty_index]])
      
      # Extract 150-degree lines for the current penalty
      lines_150_penalty <- lapply(all.cps[[2]], function(line) line[[1]][[penalty_index]])
      
      # Extract 180-degree lines for the current penalty
      lines_180_penalty <- lapply(all.cps[[3]], function(line) line[[1]][[penalty_index]])
      
      # Combine results for each type of line for the current penalty
      all.cps.rearranged[[penalty_index]] <- list(
        lines_30 = lines_30_penalty,
        lines_150 = lines_150_penalty,
        lines_180 = lines_180_penalty
      )
    }
    
    # Initialize the second output structure for `cost.mat`
    all.cost.mat <- list(
      lines_30 = lapply(all.cps[[1]], function(line) line[[2]]),  # Extract `cost.mat` for each 30-degree line
      lines_150 = lapply(all.cps[[2]], function(line) line[[2]]), # Extract `cost.mat` for each 150-degree line
      lines_180 = lapply(all.cps[[3]], function(line) line[[2]])  # Extract `cost.mat` for each 180-degree line
    )
  }
  

  # Assign rearranged_all_cps back to all.cps if needed
  all.cps.list=all.cps.rearranged
  
  
  nn.cps.list =c()
  all.segs.list=c()
  all.segs.length.list=c()
  for (all.cps in all.cps.list) {
    nn.cps=lapply(all.cps, function(x) lapply(x, function(y) length(y$change_points)))
    all_segs = get_all_segs(all.cps, all.lines=all.lines)
    all_seg_length = lapply(all_segs, function(x) lapply(x, function(y) lapply(y, length)))
    
    nn.cps.list =c(nn.cps.list,list(nn.cps))
    all.segs.list=c(all.segs.list,list(all_segs))
    all.segs.length.list=c(all.segs.length.list,list(all_seg_length))
  }
  
  names(all.cps.list)=paste0('alpha_',sprintf("%02d", as.integer(alpha_vec * 100)))
  names(nn.cps.list)=paste0('alpha_',sprintf("%02d", as.integer(alpha_vec * 100)))
  names(all.segs.list)=paste0('alpha_',sprintf("%02d", as.integer(alpha_vec * 100)))
  names(all.segs.length.list)=paste0('alpha_',sprintf("%02d", as.integer(alpha_vec * 100)))
  
  # Save the specified variables before proceeding
  alpha_str <- paste0(sprintf("%02d", as.integer(alpha_vec * 100)), collapse = "_")
  save_file <- paste0("/Users/linzhong/Desktop/ST/segmentline/", dat_name, "_exhst_cps_penalty_", alpha_str, ".Rdata")
  
  save(all.lines, K, alpha_vec, all.cps.list, nn.cps.list, all.segs.list, all.segs.length.list,YY,XX,all.cost.mat, file = save_file)
  
  return(list(all_lines = all.lines, K = K, alpha_vec = alpha_vec, 
              all.cps.list=all.cps.list, nn.cps.list=nn.cps.list, 
              all.segs.list=all.segs.list, all.segs.length.list=all.segs.length.list,
              all.cost.mat=all.cost.mat,
              YY=YY, XX=XX))
}




#Get slope and means
calculate_slopes <- function(all_segs, all.lines, y_markers, sp.counts.r, method_beta='PG',alpha,dat_name) {
  
  if (length(y_markers) == 1) {
    YY = sp.counts.r[, y_markers]
    names(YY) = rownames(sp.counts.r)
  }
  if (length(y_markers) > 1) {
    YY = rowSums(sp.counts.r[, y_markers])
    names(YY) = rownames(sp.counts.r)
  }
  
  NN = length(YY)
  XX = cbind(rep(1, NN), rep(NA, NN), rowSums(sp.counts.r[, tcell_markers])) #adjusted for T cell marker sums
  rownames(XX) = names(YY)
  
  dat = data.frame(x = XX[,3], y = YY)
  mm = glm.nb(y ~ x, data = dat)
  
  # Get adjusted response values (residuals) using Pearson residuals
  yy_adjusted = residuals.glm(mm, type = 'pearson')  # Adjusted for the number of T cell markers
  names(yy_adjusted) = names(YY)
  
  # Get the slopes and adjusted means for all the segments
  all_betas = c()
  all_means = c()
  
  for (i in 1:length(all_segs)) {
    tmp.segs = all_segs[[i]]
    
    # Calculate slopes
    tmp.betas = lapply(tmp.segs, function(x) {
      lapply(x, function(seg) {
        if (length(unique(YY[rownames(seg)])) == 1) {
          return(0)  # Return 0 if there's only one unique element in 'seg'
        } else {
          # Modify the second column of XX[seg,] to be 1:length(seg)
          XX_modified = XX[rownames(seg), ]
          XX_modified[, 2] = 1:nrow(seg)
          
          # Apply get_NBR_ests with the modified XX
          return(get_NBR_ests(n_iter = 50, Y = YY[rownames(seg)], X = XX_modified, method = method_beta)$beta[2])
        }
      })
    })
    
    # Calculate adjusted means
    tmp.means = lapply(tmp.segs, function(x) {
      lapply(x, function(seg) {
        return(mean(yy_adjusted[rownames(seg)]))  # Adjusted means
      })
    })
    
    all_betas = c(all_betas, list(tmp.betas))
    all_means = c(all_means, list(tmp.means))
  }
  
  ### Assign slopes and means to each spot
  slope.mat = matrix(NA, nrow = nrow(loc.raw), ncol = length(all.lines))
  mean.mat = matrix(NA, nrow = nrow(loc.raw), ncol = length(all.lines))
  
  rownames(slope.mat) = rownames(loc.raw)
  rownames(mean.mat) = rownames(loc.raw)
    
  for (i in 1:length(all_segs)) {
    tmp_segs = all_segs[[i]]
    tmp_betas = all_betas[[i]]
    tmp_means = all_means[[i]]
    
    for (j in 1:length(tmp_segs)) {
      tmp.segs = tmp_segs[[j]]
      tmp.beta = tmp_betas[[j]]
      tmp.mean = tmp_means[[j]]
      
      l.segs = unlist(lapply(tmp.segs, nrow))
      
      # Set slope and mean of a segment with only 2 data points to NA
      if (min(l.segs) < 3) {
        tmp.beta[which(l.segs < 3)] = NA
        tmp.mean[which(l.segs < 3)] = NA
      }
      
      for (k in 1:length(tmp.segs)) {
        slope.mat[rownames(tmp.segs[[k]]), i] = tmp.beta[k][[1]]
        mean.mat[rownames(tmp.segs[[k]]), i] = tmp.mean[k][[1]]
      }
    }
  }
  
  slope.output=list(slope.mat=slope.mat,mean.mat=mean.mat)
  
  alpha_str <- sprintf("%02d", as.integer(alpha * 100))
  save_file <- paste0("/Users/linzhong/Desktop/ST/segmentline/", dat_name, "_exhst_slopes_penalty", alpha_str, ".Rdata")
  
  save(slope.output, file=save_file)
  return(slope.output)
}





#Get slope and means
calculate_slopes_sim <- function(all_segs, all.lines, YY,XX, method_beta='PG',alpha,dat_name,loc.raw) {
  
  NN = length(YY)
  
  dat = data.frame(x = XX[,3], y = YY)
  mm = glm.nb(y ~ x, data = dat)
  
  # Get adjusted response values (residuals) using Pearson residuals
  yy_adjusted = residuals.glm(mm, type = 'pearson')  # Adjusted for the number of T cell markers
  names(yy_adjusted) = names(YY)
  
  # Get the slopes and adjusted means for all the segments
  all_betas = c()
  all_means = c()
  
  for (i in 1:length(all_segs)) {
    tmp.segs = all_segs[[i]]
    
    # Calculate slopes
    tmp.betas = lapply(tmp.segs, function(x) {
      lapply(x, function(seg) {
        if (length(unique(YY[rownames(seg)])) == 1) {
          return(0)  # Return 0 if there's only one unique element in 'seg'
        } else if (nrow(seg)<3) { #if less than 3 elements in the segment, do not calcualte slope
          return(NA)
        }
        else
        {
          # Modify the second column of XX[seg,] to be 1:length(seg)
          XX_modified = XX[rownames(seg), ]
          XX_modified[, 2] = 1:nrow(seg)
          
          # Apply get_NBR_ests with the modified XX
          return(get_NBR_ests(n_iter = 50, Y = YY[rownames(seg)], X = XX_modified, method = method_beta)$beta[2])
        }
      })
    })
    
    # Calculate adjusted means
    tmp.means = lapply(tmp.segs, function(x) {
      lapply(x, function(seg) {
        return(mean(yy_adjusted[rownames(seg)]))  # Adjusted means
      })
    })
    
    all_betas = c(all_betas, list(tmp.betas))
    all_means = c(all_means, list(tmp.means))
  }
  
  ### Assign slopes and means to each spot
  slope.mat = matrix(NA, nrow = nrow(loc.raw), ncol = length(all.lines))
  mean.mat = matrix(NA, nrow = nrow(loc.raw), ncol = length(all.lines))
  
  rownames(slope.mat) = rownames(loc.raw)
  rownames(mean.mat) = rownames(loc.raw)
  
  for (i in 1:length(all_segs)) {
    tmp_segs = all_segs[[i]]
    tmp_betas = all_betas[[i]]
    tmp_means = all_means[[i]]
    
    for (j in 1:length(tmp_segs)) {
      tmp.segs = tmp_segs[[j]]
      tmp.beta = tmp_betas[[j]]
      tmp.mean = tmp_means[[j]]
      
      l.segs = unlist(lapply(tmp.segs, nrow))
      
      # Set slope and mean of a segment with only 2 data points to NA
      if (min(l.segs) < 3) {
        tmp.beta[which(l.segs < 3)] = NA
        tmp.mean[which(l.segs < 3)] = NA
      }
      
      for (k in 1:length(tmp.segs)) {
        slope.mat[rownames(tmp.segs[[k]]), i] = tmp.beta[k][[1]]
        mean.mat[rownames(tmp.segs[[k]]), i] = tmp.mean[k][[1]]
      }
    }
  }
  
  slope.output=list(slope.mat=slope.mat,mean.mat=mean.mat)
  
  alpha_str <- sprintf("%02d", as.integer(alpha * 100))
  save_file <- paste0("/Users/linzhong/Desktop/ST/segmentline/", dat_name, "_exhst_slopes_penalty", alpha_str, ".Rdata")
  
  save(slope.output, file=save_file)
  return(slope.output)
}



plot_contour<-function(){
  
  require(interp)
  
  YY = rowSums(sp.counts.r[, ex_markers])
  NN = length(YY)
  XX = cbind(rep(1, NN), rep(NA, NN), rowSums(sp.counts.r[, tcell_markers]))
  rownames(XX) = names(YY)
  dat = data.frame(x = XX[,3], y = YY)
  mm = glm.nb(y ~ x, data = dat)
  
  # Get adjusted response values (residuals) using Pearson residuals
  yy_new = residuals.glm(mm, type = 'pearson')  # Adjusted for the number of T cell markers
  names(yy_new) = names(YY)
  
  all.spots=lapply(all.lines,function(x) lapply(x, rownames))
  all.spots=unique(unlist(all.spots))
  
  dat=data.frame(loc.raw[all.spots,],yy_new[all.spots])
  colnames(dat)=c('x','y','z')
  
  interp_data <- with(dat, interp(x, y, z, duplicate = "mean"))
  # Plot the contour using the interpolated data
  contour(interp_data$x, interp_data$y, interp_data$z, 
          xlab = "X", ylab = "Y", main = "Contour Plot of Interpolated Data")
  
  filled.contour(interp_data$x, interp_data$y, interp_data$z,
                 color.palette = terrain.colors,
                 xlab = "X", ylab = "Y", main = "Colored Contour Plot of Interpolated Data")
  
}

plot_lines_with_colored_segments <- function(loc, lines.list, grid.type, all_segs) {
  colors <- c("purple", "green", "brown", "cyan", "blue")
  
  if (grid.type == 'square') {
    par(mfrow = c(1, 2))
    
    vertical_lines = lines.list$vertical_lines
    plot(loc[,1:2], xlab = "X Coordinate", ylab = "Y Coordinate", pch = 19, cex = .3, col = 'grey',
         main = 'Vertical lines with segments')
    for (i in seq_along(vertical_lines)) {
      line_segments <- all_segs[[1]][[i]]
      for (j in seq_along(line_segments)) {
        segment_data <- line_segments[[j]]
        if (nrow(segment_data) > 1) {
          lines(segment_data[,1], segment_data[,2], col = colors[(j - 1) %% length(colors) + 1], lwd = .9)
        }
      }
    }
    
    horizontal_lines = lines.list$horizontal_lines
    plot(loc[,1:2], xlab = "X Coordinate", ylab = "Y Coordinate", pch = 19, cex = .3, col = 'grey',
         main = 'Horizontal lines with segments')
    for (i in seq_along(horizontal_lines)) {
      line_segments <- all_segs[[2]][[i]]
      for (j in seq_along(line_segments)) {
        segment_data <- line_segments[[j]]
        if (nrow(segment_data) > 1) {
          lines(segment_data[,1], segment_data[,2], col = colors[(j - 1) %% length(colors) + 1], lwd = .9)
        }
      }
    }
  }
  
  if (grid.type == 'hex') {
    par(mfrow = c(1, 3))
    
    lines_30 = lines.list$lines_30
    plot(loc[,1:2], xlab = "X Coordinate", ylab = "Y Coordinate", pch = 19, cex = .3, col = 'grey',
         main = '30 degree lines with segments')
    for (i in seq_along(lines_30)) {
      line_segments <- all_segs[[1]][[i]]
      for (j in seq_along(line_segments)) {
        segment_data <- line_segments[[j]]
        if (nrow(segment_data) > 1) {
          lines(segment_data[,1], segment_data[,2], col = colors[(j - 1) %% length(colors) + 1], lwd = .9)
        }
      }
    }
    
    lines_90 = lines.list$lines_90
    plot(loc[,1:2], xlab = "X Coordinate", ylab = "Y Coordinate", pch = 19, cex = .3, col = 'grey',
         main = '90 degree lines with segments')
    for (i in seq_along(lines_90)) {
      line_segments <- all_segs[[2]][[i]]
      for (j in seq_along(line_segments)) {
        segment_data <- line_segments[[j]]
        if (nrow(segment_data) > 1) {
          lines(segment_data[,1], segment_data[,2], col = colors[(j - 1) %% length(colors) + 1], lwd = .9)
        }
      }
    }
    
    lines_150 = lines.list$lines_150
    plot(loc[,1:2], xlab = "X Coordinate", ylab = "Y Coordinate", pch = 19, cex = .3, col = 'grey',
         main = '150 degree lines with segments')
    for (i in seq_along(lines_150)) {
      line_segments <- all_segs[[3]][[i]]
      for (j in seq_along(line_segments)) {
        segment_data <- line_segments[[j]]
        if (nrow(segment_data) > 1) {
          lines(segment_data[,1], segment_data[,2], col = colors[(j - 1) %% length(colors) + 1], lwd = .9)
        }
      }
    }
  }
}


plot_spot_by_colors1 <- function(loc, values, leg_title) {
  # Ensure 'loc' has 'x' and 'y' columns
  if (!all(c("x", "y") %in% names(loc))) {
    stop("The 'loc' data frame must contain 'x' and 'y' columns.")
  }
  
  # Define the number of colors in the gradient
  n_colors <- 100
  
  # Create a color palette from blue to red
  color_palette <- colorRampPalette(c("blue", "red"))(n_colors)
  
  # Normalize 'values' to range between 1 and n_colors
  values_min <- min(values, na.rm = TRUE)
  values_max <- max(values, na.rm = TRUE)
  
  # Handle the case where all values are identical to prevent division by zero
  if (values_max == values_min) {
    normalized_values <- rep(1, length(values))
  } else {
    normalized_values <- (values - values_min) / (values_max - values_min)
  }
  
  # Assign colors based on the normalized values
  color_indices <- as.numeric(cut(normalized_values, breaks = n_colors, include.lowest = TRUE))
  colors_assigned <- color_palette[color_indices]
  
  # Plot the scatter plot
  plot(
    loc$x, loc$y,
    col = colors_assigned,
    pch = 16,               # Solid circles
    cex = 0.5,              # Point size
    xlab = "X Coordinate",
    ylab = "Y Coordinate",
    main=leg_title
  )
}
