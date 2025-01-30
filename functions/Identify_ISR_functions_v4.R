

###1. Get endpoints of 6 directions
get_endpoints<-function(all_segs,spots){
  endpoints=matrix(-1,nrow=length(spots),ncol=length(all_segs))
  rownames(endpoints)=spots
  colnames(endpoints)=c('degree_30','degree_90','degree_150')
  for (i in 1:length(all_segs)) {
    tmp.segs=all_segs[[i]]
    tmp.endpoints=unlist(lapply(tmp.segs,function(x) unlist(lapply(x, function(y) rownames(y)[c(1,nrow(y))]))))
    endpoints[tmp.endpoints,i]=1
  }
  endpoints=cbind(endpoints,endpoints)
  colnames(endpoints)[4:6]=c('degree_-150','degree_-90','degree_-30')
  return(endpoints) #endpoints of 3 directions
}


### 2. Assign direction 1-6 of the neighbor, direction is pointing the center spot to the neighbors
calculate_direction_center_to_nb <- function(spot_i, neighbors, loc) {
  spot_i_coords <- loc[spot_i, ]
  neighbor_coords <- loc[neighbors, ]
  
  # Calculate dx and dy from the center spot to the neighbors
  dx <- neighbor_coords[, 1] - as.numeric(spot_i_coords[1]) 
  dy <- neighbor_coords[, 2] - as.numeric(spot_i_coords[2]) 
  
  # Calculate the angle in degrees
  angle <- atan2(dy, dx) * 180 / pi # Convert radians to degrees
  angle <- round(angle, -1) # Round to the nearest 10 degrees
  angle <- paste0('degree_', angle)
  
  return(angle)
}


### 3. Get the list of neighbors of endpoints in each direction. The direction is pointing from endpoint to regular point

get_endpoints_neighbors<-function(endpoints.mat,slope.dat,nb.list,loc,margin.spots){
  
  nb.of.eps.list=c()
  endp.works.list=c() 
  angles=colnames(endpoints.mat)
  
  for (i in 1:6) {
    tmp.eps=rownames(endpoints.mat)[which(endpoints.mat[,i]==1)]
    tmp.eps=tmp.eps[which(!tmp.eps %in% margin.spots)] #Need to remove margin spots
    aa=angles[i] #angle
    nb_of_ends=c()
    endp_works=c()
    for (j in tmp.eps) {
      nb=nb.list[[j]]
      nb=nb[which(!nb %in% tmp.eps)] #exclude those who are the adjacent other endpoints
      direction=calculate_direction_center_to_nb(j,nb,loc)
      if (aa %in% direction) {
        if (slope.dat[nb[which(direction==aa)],aa]==1) {
          endp_works=c(endp_works,j)
          nb_of_ends=c(nb_of_ends,nb[which(direction==aa)])
        }
      }
    }
    names(nb_of_ends)=endp_works
    endp.works.list=c(endp.works.list,list(endp_works))
    nb.of.eps.list=c(nb.of.eps.list,list(nb_of_ends))
  }
  
  names(nb.of.eps.list)=colnames(endpoints.mat)
  names(endp.works.list)=colnames(endpoints.mat)
  
  #For spots.matter redo there neighbor list to include only those regular neighbors
  
  all_spots=unlist(nb.of.eps.list, use.names = T)
  adj_ep_spots=split(names(all_spots), all_spots)
  adj_ep_spots_noangle=lapply(adj_ep_spots, function(x) substr(x, nchar(x) - 17, nchar(x)))
  spots.add=unlist(lapply(adj_ep_spots_noangle,length))
  
  
  new.nb.list=c()
  for (spot.name in names(spots.add)) {
    nb=nb.list[[spot.name]]
    new.nb=nb[which(!nb %in% adj_ep_spots_noangle[[spot.name]])]
    #print(nb)
    #print(new.nb)
    #print('=======')
    new.nb.list=c(new.nb.list,list(new.nb))
  }
  return(list(n.add=spots.add,new.nb.list=new.nb.list))
}


###3. Assign Hamilton and probablity of being ss=1

get_cond_pp_new<-function(spot.name,n.add,nb.list,new.nb.list,ss,J){
  if (spot.name %in% names(n.add)) {
    nb=new.nb.list[[spot.name]]
    ss.nb=ss[nb]
    Ham1=-J*(length(which(ss.nb==1))+n.add[spot.name])
    Ham0=-J*(length(which(ss.nb==-1)))
    pp1=exp(-Ham1)/(exp(-Ham1)+exp(-Ham0))
  } else {
    nb=nb.list[[spot.name]]
    ss.nb=ss[nb]
    Ham1=-J*(length(which(ss.nb==1)))
    Ham0=-J*(length(which(ss.nb==-1)))
    pp1=exp(-Ham1)/(exp(-Ham1)+exp(-Ham0))
  }
  return(pp1)
}




### 4. Gibbs sampling to get updated ss
gibbs_sampling_new <- function(ss, spot_names, n_add,nb_list, new_nb_list,J, loc, iterations) {
  
  names(ss)=spot_names
  
  n <- length(ss)  # Number of elements in ss
  ss_samples <- matrix(0, nrow = iterations, ncol = n)  # Store all ss vectors across iterations
  colnames(ss_samples)=spot_names
  ss_samples[1, ] <- ss
  
  # Gibbs sampling iterations
  for (iter in 2:iterations) {
    if (iter %% 100==0) print(iter)
    current_ss <- ss_samples[iter - 1, ]  # Start with the previous ss vector
    # Update each element in ss
    for (i in 1:n) {
      # Calculate the conditional probability of ss[i] being 1
      prob_1=get_cond_pp_new(spot.name=spot_names[i],
                             n.add=n_add,
                             nb.list= nb_list,
                             new.nb.list=new_nb_list,
                             ss=current_ss, 
                             J=J)
      #print(prob_1)
      # Sample from a Bernoulli distribution based on prob_1
      current_ss[i] <- sample(c(1, -1), size = 1, prob = c(prob_1, 1 - prob_1))
    }
    # Store the updated ss for the current iteration
    ss_samples[iter, ] <- current_ss
  }
  return(ss_samples)
}






### 5. Calculate probability of ss=1
calculate_final_prob <- function(ss_samples, loc, burn_in,J) {
  # Remove burn-in samples
  post_burn_in_samples <- ss_samples[(burn_in + 1):nrow(ss_samples), ]
  
  # Calculate the final probability of being 1 for each element
  final_prob <- colMeans(post_burn_in_samples)
  
  # Plot the data
  plot(loc$x, loc$y, col = rgb(final_prob, 0, 1 - final_prob), pch = 16, cex = .5,xlim=c(min(loc[,1])-50,max(loc[,1])+10000),
       xlab = "X Coordinate", ylab = "Y Coordinate", 
       main = paste0("J=",J))
  
  # Determine legend position based on the range of x and y
  legend_x <- max(loc$x)+1000   # Shift legend to the right of the maximum x-coordinate
  legend_y <- max(loc$y)      # Place legend near the top of the y range
  
  # Add legend
  legend(legend_x, legend_y, legend = c("Low", "High"), 
         fill = c("blue", "red"), title = "Pr(ISR=1)", bty = "n",cex=.7,x.intersp = 0.5)
  
  return(final_prob)
}








get_expand_region<-function(region,tumor.loc,all.region.spots,tt=3){
  center.loc=colMeans(region)
  distances <- sqrt((tumor.loc$x - center.loc[1])^2 + (tumor.loc$y - center.loc[2])^2)
  exp.region=region
  exp.spots=rownames(exp.region)
  adj.spots=exp.spots[which(! exp.spots %in% rownames(region))]
  rr=1.4*r.unit
  while(length(adj.spots) <nrow(region)*tt){
    rr=rr+r.unit
    exp.region=tumor.loc[distances <= rr, ]
    exp.spots=rownames(exp.region)
    adj.spots=exp.spots[which(! exp.spots %in% rownames(region))]
    adj.spots=adj.spots[which(! adj.spots %in% all.region.spots)]
    #print(length(adj.spots))
  }
  return(adj.spots)
}




get_TRR<-function(final.probs){
  cutp=quantile(final.probs,probs = seq(0,1,length=11))[10]
  ltsa=final.probs
  ltsa[final.probs>cutp]=1
  ltsa[final.probs<=cutp]=0
  
  red_spots <- loc.raw[which(rownames(loc.raw) %in% names(ltsa)[which(ltsa == 1)]), ]
  
  # Create a graph where edges connect nearby points (e.g., within a certain distance threshold)
  threshold <- 1.5*r.unit  # Example distance threshold
  dist_matrix <- as.matrix(dist(red_spots))
  adj_matrix <- dist_matrix <= threshold & dist_matrix > 0  # Adjacency based on threshold
  graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  
  # Identify connected components (regions)
  components <- components(graph)
  
  # Filter out regions with fewer than 5 spots
  regions=lapply(1:length(components$csize), function(x) names(components$membership)[which(components$membership==x)])
  final.regions=regions[which(components$csize >= 10)]
  
  final.regions=lapply(final.regions,function(x) x[which(x %in% rownames(tumor.loc))])
  tmp.ll=unlist(lapply(final.regions,length))
  final.regions=final.regions[which(tmp.ll>=10)]
  
  par(mfrow=c(1,2))
  
  plot(loc.raw, main = "Spatial Plot", xlab = "X", ylab = "Y", pch = 16, col = "gray",ylim=c(0,max(loc.raw[,2]+2000)))
  points(tumor.loc,col='yellow',cex=0.6)
  points(loc.raw[which(rownames(loc.raw) %in% names(ltsa)[which(ltsa == 1)]), ], 
         col = "red", pch = 19, cex = 0.5)

  legend("topright", legend = c("LTSA", "Tumor Area"), 
         col = c("red", "yellow"), 
         pch = 19, 
         pt.cex = 1.2,  # Adjust point size in legend
         bty = "n")     # Remove legend box
  
  
  # Plot the updated red spots
  plot(loc.raw, main = "Filtered LTSA Regions", xlab = "X", ylab = "Y", 
       col = "gray", pch = 16,ylim=c(0,max(loc.raw[,2]+2000)))
  points(tumor.loc,col='yellow',cex=0.6)  # Blue with 50% transparency
  points(loc.raw[unlist(final.regions),], col = "red", pch = 19, cex = 0.7)

  legend("topright", legend = c("Entry region", "Tumor Area"), 
         col = c("red", "yellow"), 
         pch = 19, 
         pt.cex = 1.2,  # Adjust point size in legend
         bty = "n")     # Remove legend box
  
  par(mfrow=c(1,1))
  return(final.regions)
}



