if (!require("mclust")) install.packages("mclust")
library(mclust)

Rcpp::sourceCpp("functions/estimate_J_cpp_funcs_v1.cpp")
#Get neighbors for MC simulation


plot_spot_by_colors<-function(loc,values,leg){
  p1=ggplot(loc, aes(x = x, y = y, color = values))+
    geom_point(size = .5) +  # Adjust size to make dots smaller
    scale_color_gradient(low = "blue", high = "red") +  # Adjust color gradient
    theme_minimal() +  # Minimal theme
    labs(title = '',
         x = "X Coordinate",
         y = "Y Coordinate",
         color = leg) +
    theme(plot.title = element_text(hjust = 0.5))  # Center the plot title
  p1
}


get_nb<-function(xx,loc,threshold){
  dist_matrix=as.matrix(dist(loc))
  nb.list=lapply(1:nrow(loc), function(i) which(dist_matrix[i, ] <= threshold & dist_matrix[i, ] > 0))
  n.nb=unlist(lapply(nb.list,length))
  spots.keep=which(n.nb>0)
  xx.new=xx[spots.keep]
  loc.new=loc[spots.keep,]
  
  dist_matrix=as.matrix(dist(loc.new))
  nb.list=lapply(1:nrow(loc.new), function(i) which(dist_matrix[i, ] <= threshold & dist_matrix[i, ] > 0))
  nb.list.cpp=lapply(nb.list, function(x) as.numeric(x)-1) #this is to accomodate cpp
  
  nb.list=lapply(1:nrow(loc.new), function(i) rownames(loc.new)[which(dist_matrix[i, ] <= threshold & dist_matrix[i, ] > 0)])
  names(nb.list)=rownames(loc.new)
  names(nb.list.cpp)=rownames(loc.new)
    
  gmm <- Mclust(xx.new, G = 2)
  
  # View the summary of the fitted model
  print(summary(gmm))
  means <- gmm$parameters$mean
  variances <- gmm$parameters$variance$sigmasq
  print(list(means = means, variances = variances))
  
  # Cluster assignments (1 or 2)
  clusters <- gmm$classification
  ss=clusters
  ss[ss==2]=-1 #large mean
  print(table(ss))
  return(list(xx.new,loc.new,ss,nb.list,nb.list.cpp))
}


get_adjust_YY<-function(YY,XX=rowSums(sp.counts.r)){ #Adjust for umi if not specify
  NN = length(YY)
  dat = data.frame(x = XX, y = YY)
  mm = glm.nb(y ~ x, data = dat)
  # Get adjusted response values (residuals) using Pearson residuals
  yy_adjusted = residuals.glm(mm, type = 'pearson')  # Adjusted for the UMI
  names(yy_adjusted) = names(YY)
  return(yy_adjusted)
}



#use the Rcpp version
MC_sampling_x<-function(init_x,nb.list.cpp,J_k,n_flip,M){
  sum_old=calculate_sum_cpp(init_x,nb.list.cpp)
  n=length(init_x)
  sum.list=c()
  for (i in 1:n_flip) {
    new_x=init_x
    x_i=new_x[sample(1:n,1)]
    new_x[sample(1:n,1)]=-1*x_i
    sum_new=calculate_sum_cpp(new_x,nb.list.cpp)
    if (exp(J_k*(sum_new-sum_old))>runif(1)) {
      init_x=new_x #update configuration
      sum_old=sum_new #update sum
    }
    #print(exp(J_k*(sum_new-sum_old)))
    #print(sum_old)
    sum.list=c(sum.list,sum_old)
  }
  partial_der=mean(sum.list[c(n_flip-M:n_flip)])
  return(list(sum.list,partial_der))
}


test_J<-function(dat,n.sample,n.flip,mu_Q_candidates){ 
  #n.sample is number of samples draw
  #n.flip is number of flips
  
  unique.sums=unique(dat$obs_sum_list)
  p.list=table(dat$obs_sum_list)/length(dat$obs_sum_list)
  sample_sums=sample(unique.sums,prob=p.list,size=n.sample,replace=T)
  mu.list=c()
  for (i in 1:length(sample_sums)) {
    print(i)
    first_dr<-find_optimal_mu_Q(mu_Q_candidates=mu_Q_candidates,
                                init_x=sample(c(-1,1),size=nrow(dat$x_samples),replace=T), #initial states
                                nb.list.cpp=dat$nb.list.cpp, #neighbors
                                mu=0.3, sigma2=0.0049, #prior of J
                                sigma2_Q=0.0001,
                                K=10, s_obs=sample_sums[i] , n_flip=n.flip, M=10000) 
    mu_est=names(first_dr)[which(first_dr==min(first_dr))]
    mu.list=c(mu.list,mu_est)
  }
  return(as.numeric(mu.list))
}


#Randomly choose entry points
Ising_data_simulation<-function(n_rows,n_cols,J,n_flip){ #M is number of configurations wanted
  
  spot_diameter <- 1  # Distance between adjacent spots
  
  # Initialize lists to hold x and y coordinates
  x_coords <- numeric()
  y_coords <- numeric()
  
  # Generate coordinates in a triangular pattern
  for (row in 0:(n_rows - 1)) {
    for (col in 0:(n_cols - 1)) {
      # Calculate base x and y coordinates
      x <- col * spot_diameter
      y <- row * spot_diameter * sqrt(3) / 2  # sqrt(3)/2 factor for 60-degree angle
      
      # Offset every second row to create the triangular tiling
      if (row %% 2 != 0) {
        x <- x + (spot_diameter / 2)
      }
      
      # Append coordinates to lists
      x_coords <- c(x_coords, x)
      y_coords <- c(y_coords, y)
    }
  }
  
  # Combine into a data frame
  grid_data <- data.frame(x = x_coords, y = y_coords)
  
  grid_data=grid_data[,2:1]
  colnames(grid_data)=c('x','y')
  
  
  # View the first few coordinates
  rownames(grid_data)=paste0('spot_',1:nrow(grid_data))

  spots=rownames(grid_data)
  dist_matrix=as.matrix(dist(grid_data))
  r.unit=1
  
  threshold=1.2
  
  nb.list=lapply(1:nrow(grid_data), function(i) which(dist_matrix[i, ] <= threshold & dist_matrix[i, ] > 0))
  nb.list.cpp=lapply(nb.list, function(x) as.numeric(x)-1) #this is to accomodate cpp
  
  nb.list=lapply(1:nrow(grid_data), function(i) rownames(grid_data)[which(dist_matrix[i, ] <= threshold & dist_matrix[i, ] > 0)])
  
  names(nb.list)=spots
  
  names(nb.list.cpp)=spots
  
  #n.nb=unlist(lapply(nb.list,length))
  
  n=length(spots)
  init_x=sample(c(1,-1),prob=c(0.9,0.1),size=n,replace=T)
  sum_old=calculate_sum_cpp(init_x,nb.list.cpp)
  
  x.mat=c()
  sum.list=c()
  for (i in 1:n_flip) {
    new_x=init_x
    x_i=new_x[sample(1:n,1)]
    new_x[sample(1:n,1)]=-1*x_i
    sum_new=calculate_sum_cpp(new_x,nb.list.cpp)
    if (exp(J*(sum_new-sum_old))>runif(1)) {
      init_x=new_x #update configuration
      sum_old=sum_new #update sum
    }
    
    if (i>n_flip-5000) {
      x.mat=cbind(x.mat,init_x)
      sum.list=c(sum.list,sum_old)
    }

  }
    return(list(loc.dat=grid_data,nb.list.cpp=nb.list.cpp,J=J,obs_sum_list=sum.list,x_samples=x.mat))
}

find_optimal_mu_Q <- function(mu_Q_candidates=seq(0.16, 0.44, by = 0.01),init_x, nb.list.cpp, mu, sigma2,sigma2_Q, K, s_obs, n_flip, M) {
  
  # Initialize a vector to store the first derivative values for each candidate
  first_derivative_values <- numeric(length(mu_Q_candidates))
  names(first_derivative_values)=mu_Q_candidates
  first_derivative_values=c()
  
  # Loop through each candidate value
  for (i in seq_along(mu_Q_candidates)) {
    mu_Q <- mu_Q_candidates[i]
    
    # Sample epsilon_k ~ N(0, 1) and compute J_k
    J_k <- rnorm(K, mean = mu_Q, sd = sqrt(sigma2_Q)) # K samples from N(0, 1)
    
    d_logZ_dJk <- sapply(J_k, function(j) {
      MC_sampling_x_cpp(init_x = init_x, nb_list = nb.list.cpp, J_k = j, n_flip = n_flip, M = M)[[2]]
    })
    
    # Calculate the first derivative of the ELBO with respect to mu_Q
    first_derivative <- s_obs - (1 / K) * sum(d_logZ_dJk) - (mu_Q - mu) / sigma2
    
    # Store the first derivative value
    first_derivative_values=c(first_derivative_values,abs(first_derivative)) # Use absolute value to find closest to 0
    nn=length(first_derivative_values)
    if(nn>5) {
      if (min(diff(first_derivative_values[(nn-3):nn]))>0) {break}
    }
    print(first_derivative_values)
    #cat("mu_Q:",mu_Q)
    #cat("first derivative:",first_derivative)
    #print('==========')
  }
  names(first_derivative_values)=mu_Q_candidates[1:length(first_derivative_values)]
  return(first_derivative_values)
}




estimate_J1_dat<-function(genes,sp.dd=sp.counts.r,adj.umi=TRUE){
  par(mfrow=c(1,1))
  hla.c1=sp.dd[,genes]
  if (length(genes)>1) {hla.c1=rowSums(hla.c1)}
  
  if (adj.umi==TRUE)
  {hla.c1=get_adjust_YY(hla.c1)} #this function always adjust for UMI
  
  plot_spot_by_colors1(hla.c1,loc=loc.raw,leg='Level of antigen')
  
  dat=get_nb(xx=hla.c1,loc=loc.raw,threshold = r.unit*1.5)
  obs_ss=dat[[3]]
  loc.dat=dat[[2]]
  nb.r=dat[[4]]
  nb.info=dat[[5]]
  
  print("Number of numbers of each spot")
  print(table(unlist(lapply(nb.r,length))))
  
  n.pairs=sum(unlist(lapply(nb.info,length)))/2
  cat(paste0("number of pairs:",n.pairs,"; "))
  sum_obs=calculate_sum_cpp(obs_ss, nb.info)
  cat(paste0("number of pairs with same state:",sum_obs,'; '))
  omega=sum_obs/n.pairs
  cat(paste0("omega:",omega))
  return(list(dat,sum_obs))
}






plot_J2<-function(J1,lambda.list){
  #Locate impurities
  all.cps=all.cps.list[[10]]
  all_segs=all.segs.list[[10]]
  spots=rownames(loc.raw)
  endpoints.mat=get_endpoints(all_segs,spots=spots)
  
  slope.mat=slope.output[[1]]
  print("Matrix of slopes:")
  print(head(slope.mat))
  
  slope.dat=cbind(slope.mat,-slope.mat)
  colnames(slope.dat)=paste0('degree_',c(30,90,150,-150,-90,-30))
  slope.dat[is.na(slope.dat)]=0
  slope.dat[which(slope.dat> -0.01)]=-1
  slope.dat[which(slope.dat!= -1)]=1
  
  print("Number of slopes <-0.01:")
  print(table(slope.dat))
  
  tmp.dat=get_endpoints_neighbors(endpoints.mat,slope.dat,nb.r,loc=loc.raw,margin.spots=margin.spots)
  n.add=tmp.dat[[1]]
  new.nb.list=tmp.dat[[2]]
  
  #Get results for different J2's
  
  ss=sample(c(-1,1),size=length(obs_ss),prob=c(0.5,0.5),replace = T) 
  
  ss.mat.list=c()
  for (lambda in lambda.list) {
    print(lambda)
    tmp.ss.mat=gibbs_sampling_new(ss=ss, 
                                  spot_names=names(obs_ss), 
                                  n_add=n.add*lambda,
                                  nb_list=nb.r, 
                                  new_nb_list=new.nb.list,
                                  J=J1, 
                                  loc=loc.raw, 
                                  iterations = 1000) 
    ss.mat.list=c(ss.mat.list,list(tmp.ss.mat))
  }
  
  mean.list=c()
  for (i in 1:length(ss.mat.list)) {
    ss.mat=ss_samples=ss.mat.list[[i]]
    ss.mat[ss.mat==-1]=0
    tmp.final.probs=calculate_final_prob(ss.mat, loc=loc.dat, burn_in = 500,J=J1)
    mean.list=c(mean.list,length(which(tmp.final.probs>0.7))/length(ss))
  }
  
  mean.list=mean.list*100
  
  return(list(mean.list=mean.list,n.add=n.add, new.nb.list=new.nb.list))
}




























