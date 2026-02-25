require(rstan)
require(MASS)
require(BayesLogit)
require(Rcpp)
require(reshape2)

sourceCpp("functions/NBR_estimate_PG.cpp")

################
#Vectorized PELT
################


optimal_partitioning_pelt_vectorized <- function(n_iter=20, Y, X, method, alpha_vec, K=0) {
  
  #Get cost matrix
  cost.mat=optimal_partitioning_pelt(n_iter=n_iter,Y, X, method=method, alpha = .5 ,K=0)$cost_mat
  np=(ncol(X)+1) #number of parameters, beta and dispersion parameter
  
  results_list <- list()
  for (alpha in alpha_vec) {
    
    n <- length(Y)  # Length of the data
    
    # Initialize cost array and partition point array
    cost <- rep(Inf, n + 1)   # Cost for partitioning up to each point
    partition <- rep(0, n + 1)  # Best partition point for each endpoint
    
    # Base case: no cost for zero elements
    cost[1] <- 0
    
    #List of potential change points
    R=list()
    R[[1]]=0
    R[[2]]=1
    R[[3]]=1
    R[[4]]=c(1,3)
    R[[5]]=c(1,3,4)
    
    # Loop over all possible endpoints
    for (end in 2:n) {
      #1. Loop over all possible start points for each endpoint to get minimized cost for cost[end+1]
      for (start in R[[end]]) {  # Ensure at least 2 points in a segment
        if (start==2) {next}
        if (end - start + 1 >= 2) {  # Condition to ensure minimum segment length of 2
          if (!is.na(cost.mat[start,end])) {
            current_cost_llh=cost.mat[start,end]
            current_cost_aic=current_cost_llh+np*2*alpha
          } else {
            # Calculate the cost of the current segment
            current_cost <- segment_cost(n_iter=n_iter,Y=Y, X=X, start=start, end=end, method=method, alpha=alpha)
            current_cost_aic=current_cost[1]
            current_cost_llh=current_cost[2]
            cost.mat[start,end]=current_cost_llh
          }
          
          if (cost[start] + current_cost_aic < cost[end + 1]) {
            cost[end + 1] <- cost[start] + current_cost_aic
            partition[end + 1] <- start
            
          }
        }
      }
      #2. Pruning
      
      new_R=c() #change points that are kept after prunning for data point 'end'
      if (end>=5 & end < n)  { #At this time, we have cost[end+1] for F(end) 
        tmp.cp=R[[end]]
        for (start in tmp.cp[2:length(tmp.cp)]) {
          pruning_cost <- cost.mat[start,end] #cost with no penalty from of data point [start : end]
          f_t=cost[start] #Minimized cost (with penalty) of data point [1:(start-1)]
          
          # Pruning logic: keep start points that reduce the cost
          if (f_t + pruning_cost + K < cost[end+1]) {
            new_R <- c(new_R, start)
            #cat("Keeping start", start, "for endpoint", end, "\n")
          } else {
            #cat("Pruning start", start, "for endpoint", end, "\n")
          }
        }
        
        # Update the set of candidate changepoints for the next endpoint
        R[[end+1]] <- unique(c(1,new_R, end))
        #cat("New candidate changepoints for R[[", end + 1, "]]:", R[[end + 1]], "\n")
      }
    }
    
    # Reconstruct the partitioning solution and track change points
    segmentation <- list()
    change_points <- c()  # Store the locations of the change points
    end <- n
    while (end > 0) {
      start <- partition[end + 1]
      segmentation <- append(segmentation, list(c(start, end)))
      if (start != 1) {  # Only track the change points excluding 1
        change_points <- c(change_points, start)
      }
      end <- start - 1
    }
    
    results_list[[paste0("alpha_", alpha)]]=list(total_cost = cost[n + 1], segmentation = rev(segmentation), change_points = rev(change_points))
    
  }
  
  return(list(cp_output=results_list,cost_matrix=cost.mat))
  
}




#### OP 
optimal_partitioning <- function(n_iter=20,Y, X, method, alpha = 1) {
  n <- length(Y)  # Length of the data
  
  # Initialize cost array and partition point array
  cost <- rep(Inf, n + 1)   # Cost for partitioning up to each point
  partition <- rep(0, n + 1)  # Best partition point for each endpoint
  
  # Base case: no cost for zero elements
  cost[1] <- 0
  cost.mat<-matrix(NA,n,n) #Save cost (with no penalty, so just -2*llh)
  # Loop over all possible endpoints
  for (end in 2:n) {
    # Loop over all possible start points for each endpoint
    for (start in 1:(end - 1)) {  # Ensure at least 2 points in a segment
      if (start == 2) { next }
      if (end - start + 1 >= 2) {  # Condition to ensure minimum segment length of 2
        
        # Calculate the cost of the current segment
        current_cost <- segment_cost(n_iter=n_iter,Y=Y, X=X, start=start, end=end, method=method, alpha=alpha)[1]
        cost.mat[start,end]=current_cost_llh
        
        # Update the cost and partition arrays if we found a better partition
        if (cost[start] + current_cost < cost[end + 1]) {
          cost[end + 1] <- cost[start] + current_cost
          partition[end + 1] <- start
        }
      }
    }
  }
  
  # Reconstruct the partitioning solution
  segmentation <- list()
  change_points <- c()  # Store the change points (excluding 1)
  end <- n
  while (end > 0) {
    start <- partition[end + 1]
    segmentation <- append(segmentation, list(c(start, end)))
    
    if (start != 1) {  # Only track change points excluding 1
      change_points <- c(change_points, start)
    }
    
    # Print out the reconstructed segmentation
    print(paste("Segment from", start, "to", end))
    end <- start - 1
  }
  
  return(list(total_cost = cost[n + 1], segmentation = rev(segmentation), change_points = rev(change_points),
              cost_mat=))
}



####PELT
optimal_partitioning_pelt <- function(n_iter=20,Y, X, method, alpha = 1,K=0) {
  n <- length(Y)  # Length of the data
  
  # Initialize cost array and partition point array
  cost <- rep(Inf, n + 1)   # Cost for partitioning up to each point
  partition <- rep(0, n + 1)  # Best partition point for each endpoint
  
  # Base case: no cost for zero elements
  cost[1] <- 0
  
  #List of potential change points
  R=list()
  R[[1]]=0
  R[[2]]=1
  R[[3]]=1
  R[[4]]=c(1,3)
  R[[5]]=c(1,3,4)
  
  cost.mat<-matrix(NA,n,n) #Save cost (with no penalty, so just -2*llh)
  
  # Loop over all possible endpoints
  for (end in 2:n) {
    #1. Loop over all possible start points for each endpoint to get minimized cost for cost[end+1]
    for (start in R[[end]]) {  # Ensure at least 2 points in a segment
      if (start==2) {next}
      if (end - start + 1 >= 2) {  # Condition to ensure minimum segment length of 2
        
        # Calculate the cost of the current segment
        current_cost <- segment_cost(n_iter=n_iter,Y=Y, X=X, start=start, end=end, method=method, alpha=alpha)
        #print(start)
        #print(end)
        #print(current_cost)
        current_cost_aic=current_cost[1]
        current_cost_llh=current_cost[2]
        cost.mat[start,end]=current_cost_llh
        # Update the cost and partition arrays if we found a better partition
        
        #
        #print('debug')
        #print(start)
        #print(end)
        #print(current_cost_aic)
        #print(cost[start])
        #print(cost[end + 1])
        #print('==')
        if (cost[start] + current_cost_aic < cost[end + 1]) {
          cost[end + 1] <- cost[start] + current_cost_aic
          partition[end + 1] <- start
          
        }
      }
    }
    
    #2. Pruning
    
    new_R=c() #change points that are kept after prunning for data point 'end'
    if (end>=5 & end < n)  { #At this time, we have cost[end+1] for F(end) 
      tmp.cp=R[[end]]
      for (start in tmp.cp[2:length(tmp.cp)]) {
        pruning_cost <- cost.mat[start,end] #cost with no penalty from of data point [start : end]
        f_t=cost[start] #Minimized cost (with penalty) of data point [1:(start-1)]
        
        # Pruning logic: keep start points that reduce the cost
        if (f_t + pruning_cost + K < cost[end+1]) {
          new_R <- c(new_R, start)
          #cat("Keeping start", start, "for endpoint", end, "\n")
        } else {
          #cat("Pruning start", start, "for endpoint", end, "\n")
        }
      }
      
      # Update the set of candidate changepoints for the next endpoint
      R[[end+1]] <- unique(c(1,new_R, end))
      #cat("New candidate changepoints for R[[", end + 1, "]]:", R[[end + 1]], "\n")
    }
    
  }
  
  # Reconstruct the partitioning solution and track change points
  segmentation <- list()
  change_points <- c()  # Store the locations of the change points
  end <- n
  while (end > 0) {
    start <- partition[end + 1]
    segmentation <- append(segmentation, list(c(start, end)))
    if (start != 1) {  # Only track the change points excluding 1
      change_points <- c(change_points, start)
    }
    # Print out the reconstructed segmentation
    #print(paste("Segment from", start, "to", end))
    end <- start - 1
  }
  
  return(list(total_cost = cost[n + 1], segmentation = rev(segmentation), change_points = rev(change_points),
              cost_mat=cost.mat))
}





##### Calculate cost of a segment+penalty
segment_cost <- function(n_iter,Y,X,start,end,method=method,alpha=1) { #alpha is to adjust penalty
  yy=Y[start:end]
  xx=X[start:end,]
  #print(yy)
  #print(xx)
  out=get_NBR_ests(n_iter=n_iter,Y=yy,X=xx,method=method)
  #print(out)
  beta=out[[1]]
  phi=out[[2]]
  llh=out[[3]]
  #print(llh)
  aic=2*alpha*(length(beta)+1)-2*llh
  llh=-2*llh
  return(c(aic,llh))
}

segment_cost_vectorized <- function(n_iter, Y, X, start, end, method, alpha_vec) {
  # Extract the subset of data for the segment
  yy <- Y[start:end]
  xx <- X[start:end, ]
  
  # Compute the shared part (log-likelihood, beta, phi) just once
  out <- get_NBR_ests(n_iter=n_iter, Y=yy, X=xx, method=method)
  beta <- out[[1]]
  phi <- out[[2]]
  llh <- out[[3]]
  
  # Vectorize the AIC calculation over the provided alpha values
  aic_vec <- 2 * alpha_vec * (length(beta) + 1) - 2 * llh
  
  # Return a list with AIC and log-likelihood (llh) for all alpha values
  return(list(aic=aic_vec, llh=-2 * llh))
}


##### Get parameter estimates for negative binomial regression

get_NBR_ests<-function(n_iter,Y,X,method){ #method can be: glm.nb,MCMC,PG
  
  n_unique=apply(X,2,function(x) length(unique(x)))
  c.remove=which(n_unique==1)
  c.remove=c.remove[-1]
  c.keep=c(1,which(n_unique>1))
  ncol.x=ncol(X)
  
  if (length(c.remove)>0) {
    X=X[,-c.remove]
  }
  
  #print(X)
    
  if (ncol(X)==3) {
    X[,3]=X[,3]/mean(X[,3],na.rm=T)
  }
  if (ncol(X)>3) {
    x.rescale=apply(X[,3:ncol(X)],2,function(x) x/mean(x,na.rm=T))
    X[,3:ncol(X)]=x.rescale
  }
  
  if (method == 'glm.nb') {
    result <- tryCatch({
      # Prepare the data
      tmp.dat <- data.frame(Y = Y, X = X[,-1])  # Ensure X is correctly defined
      
      # Fit the negative binomial model
      fit <- glm.nb(Y ~ ., data = tmp.dat)  # Use "." to represent all variables
      
      # Calculate likelihood
      llh = nb_likelihood(Y, X, as.matrix(coef(fit)), 1 / fit$theta)
      
      # Return results
      list(beta = as.matrix(coef(fit)), r = 1 / fit$theta, llh = llh)
    }, error = function(e) {
      print('glm.nb fails')
      NULL  # Return NULL on failure, indicating glm.nb failed
    })
    
    # Check if result is NULL, indicating failure
    if (is.null(result)) {
      method = 'MCMC'  # Fall back to MCMC method
    }
  }
  
  if (method=='MCMC') {
    N=length(Y)
    K=ncol(X)
    stan_data <- list(
      N = N,      # Number of observations
      K = K,      # Number of predictors (including intercept)
      X = X,      # Predictor matrix
      Y = Y       # Count outcome
    )
    
    fit <- stan(model_code = stan_code, data = stan_data, iter = 2000, chains = 1, warmup = 500, thin = 1, seed = 1234,
                , verbose = TRUE, refresh = 0)
    # Get summary of parameter estimates
    summary_fit <- summary(fit)
    
    # Extract summary table for all parameters
    summary_table <- summary_fit$summary
    
    # Access mean, sd, and 95% credible intervals for beta
    beta_means <- as.matrix(summary_table[grep("beta", rownames(summary_table)), "mean"])
    beta_sd <- summary_table[grep("beta", rownames(summary_table)), "sd"]
    beta_credible_intervals <- summary_table[grep("beta", rownames(summary_table)), c("2.5%", "97.5%")]
    
    # Access mean and credible interval for r
    r_mean <- summary_table["r", "mean"]
    r_sd <- summary_table["r", "sd"]
    r_credible_interval <- summary_table["r", c("2.5%", "97.5%")]
    llh=nb_likelihood(Y,X,beta_means,r_mean)
    result=list(beta = beta_means, r = r_mean,llh=llh)
  }
  
  
  if (method=='PG') {
    phi=1
    beta=rep(0.001,ncol(X))
    for (i in 1:n_iter) {
      beta=cpp_update_beta(as.matrix(X),as.numeric(Y),
                           b=rep(0.001,ncol(X)),B=diag(c(0.5,rep(0.2,ncol(X)-1))),
                           beta=beta,phi=phi,niter_beta=20)
      phi=cpp_update_phi_MCMC(as.matrix(X),as.numeric(Y),beta=as.numeric(beta),
                              a_phi=0.01,b_phi=0.01,sd_phi=1,phi=phi,niter_phi=5000)
    }
    
    beta.list=cpp_update_beta_mat(as.matrix(X),as.numeric(Y),
                                  b=rep(0.001,ncol(X)),B=diag(c(0.5,rep(0.2,ncol(X)-1))),
                                  beta=beta,phi=phi,niter_beta=100)
    
    beta_means=apply(beta.list[31:100,],2,mean)
    
    phi.list=cpp_update_phi_MCMC_vec(as.matrix(X),as.numeric(Y),beta=as.numeric(beta_means),
                                     a_phi=0.01,b_phi=0.01,sd_phi=1,phi=phi,niter_phi=50000)
    
    
    r_mean=mean(phi.list[25000:50000])
    beta_credible_intervals=apply(beta.list,2,function(x) quantile(x[round(n_iter/2):n_iter], probs = c(0.025, 0.975)))
    r_credible_interval=quantile(phi.list[round(n_iter/2):n_iter], probs = c(0.025, 0.975))
    
    #Estimate likelihood
    llh=nb_likelihood_p(Y,X,beta_means,r_mean)

    if (length(c.remove)==0) {
      result=list(beta = as.matrix(beta_means), r = r_mean,llh=llh)
    }
    
    if (length(c.remove)>0) {
      beta.means.new=rep(0,ncol.x)
      beta.means.new[c.keep]=beta_means
      result=list(beta = as.matrix(beta.means.new), r = r_mean,llh=llh)
    }
  }
  
  return(result)
}



Bayes_NBR_ests<-function(n_iter=20,Y,X,method){ #method can be: MCMC,PG
  
  if (method=='MCMC') {
    N=length(Y)
    K=ncol(X)
    stan_data <- list(
      N = N,      # Number of observations
      K = K,      # Number of predictors (including intercept)
      X = X,      # Predictor matrix
      Y = Y       # Count outcome
    )
    
    fit <- stan(model_code = stan_code, data = stan_data, iter = 2000, chains = 1, warmup = 500, thin = 1, seed = 12,
                , verbose = FALSE, refresh = 0)
    # Get summary of parameter estimates
    summary_fit <- summary(fit)
    
    # Extract summary table for all parameters
    summary_table <- summary_fit$summary
    
    # Access mean, sd, and 95% credible intervals for beta
    beta_means <- as.matrix(summary_table[grep("beta", rownames(summary_table)), "mean"])
    beta_sd <- summary_table[grep("beta", rownames(summary_table)), "sd"]
    beta_credible_intervals <- summary_table[grep("beta", rownames(summary_table)), c("2.5%", "97.5%")]
    
    # Access mean and credible interval for r
    r_mean <- summary_table["r", "mean"]
    r_sd <- summary_table["r", "sd"]
    r_credible_interval <- summary_table["r", c("2.5%", "97.5%")]
    llh=nb_likelihood(Y,X,beta_means,r_mean)
    result=list(beta = as.matrix(beta_means), beta_credible_intervals=beta_credible_intervals,r = r_mean,r_credible_interval=r_credible_interval,llh=llh)
  }
  
  
  if (method=='PG') {
    phi=1
    beta=rep(0.001,ncol(X))
    for (i in 1:n_iter) {
      beta=cpp_update_beta(as.matrix(X),as.numeric(Y),
                           b=rep(0.001,ncol(X)),B=diag(c(0.5,rep(0.2,ncol(X)-1))),
                           beta=beta,phi=phi,niter_beta=20)
      phi=cpp_update_phi_MCMC(as.matrix(X),as.numeric(Y),beta=as.numeric(beta),
                              a_phi=0.01,b_phi=0.01,sd_phi=1,phi=phi,niter_phi=5000)
    }
    
    beta.list=cpp_update_beta_mat(as.matrix(X),as.numeric(Y),
                                  b=rep(0.001,ncol(X)),B=diag(c(0.5,rep(0.2,ncol(X)-1))),
                                  beta=beta,phi=phi,niter_beta=100)
    
    beta_means=apply(beta.list[round(n_iter/2):n_iter,],2,mean)
    
    phi.list=cpp_update_phi_MCMC_vec(as.matrix(X),as.numeric(Y),beta=as.numeric(beta_means),
                                     a_phi=0.01,b_phi=0.01,sd_phi=1,phi=phi,niter_phi=50000)
    
    
    r_mean=mean(phi.list[round(n_iter/2):n_iter])
    beta_credible_intervals=apply(beta.list,2,function(x) quantile(x[round(n_iter/2):n_iter], probs = c(0.025, 0.975)))
    r_credible_interval=quantile(phi.list[round(n_iter/2):n_iter], probs = c(0.025, 0.975))
    llh=nb_likelihood_p(Y,X,beta_means,r_mean)
    result=list(beta = as.matrix(beta_means), beta_credible_intervals=beta_credible_intervals,r = r_mean,r_credible_interval=r_credible_interval,llh=llh)
  }
  
  return(result)
}



nb_likelihood_p <- function(Y, X, beta, r) {
  X=as.matrix(X)
  # Step 1: Calculate logit(p) = X * beta
  logit_p <- X %*% as.matrix(beta)
  
  # Step 2: Convert logit(p) to probability p
  p <- 1 / (1 + exp(-logit_p))
  
  if (min(p)==0|max(p)==1) {
    print('probability is 0 or 1, return NaN for llh, check beta')
  }
  
  p[which(p<1e-10)]=1e-10
  p[which(p==1)]=0.999999999
  
  # Step 3: Calculate likelihood for each observation using dbinom for negative binomial
  likelihoods <- dnbinom(Y, size = r, prob = 1-p, log = TRUE)  # Log-likelihood for stability
  
  # Step 4: Sum the log-likelihoods to get total likelihood
  total_log_likelihood <- sum(likelihoods)
  
  return(total_log_likelihood)
}

nb_likelihood <- function(Y, X, beta, r) {
  # Step 1: Calculate the mean mu = exp(X * beta)
  X=as.matrix(X)
  log_mu <- X %*% beta
  mu <- exp(log_mu)
  
  # Step 2: Calculate the log-likelihood using the dnbinom function (log = TRUE for log-likelihood)
  log_likelihood <- sum(dnbinom(Y, size = r, mu = mu, log = TRUE))
  
  # Return the total log-likelihood
  return(log_likelihood)
}


stan_code <- "
data {
  int<lower=0> N;          // Number of observations
  int<lower=0> K;          // Number of predictors (including intercept)
  matrix[N, K] X;          // Predictor matrix
  int<lower=0> Y[N];       // Outcome variable (count data)
}

parameters {
  vector[K] beta;          // Regression coefficients
  real<lower=0> r;         // Dispersion parameter
}

model {
  vector[N] mu;

  // Priors
  beta ~ normal(0, 10);    // Prior for regression coefficients
  r ~ gamma(2, 0.1);       // Prior for dispersion parameter

  // Mean calculation
  mu = exp(X * beta);      // Link function: log(mu) = X * beta

  // Likelihood
  Y ~ neg_binomial_2(mu, r); // Negative binomial likelihood
}
"





#####Functions to check performance#############################################################################################################################################################


###1. Check NBR estimates

test.glm.nb<-function(yy,xx){
  tryCatch({glm.nb(yy~xx)},error=function(e) {
    cat("glm.nb() does not work on this data")
    return(NULL)})
}

get_NB_estimate_output<-function(dat.list){
  
  #Save the data
  Y.list=c()
  X.list=c()
  index.list=c()
  
  #Save the estimates
  intercept.mat=c()
  slope.mat=c()
  r.mat=c()
  llh.mat=c()
  
  for (i in 1:length(dat.list)) {
    print(i)
    tmp.dd=dat.list[[i]]
    
    Y=tmp.dd$adj_gene_exp
    X=tmp.dd[,c('intercept','id')]
    mm=test.glm.nb(yy=Y,xx=X[,2])
    if (is.null(mm) | length(unique(Y))<3) {next}
    
    tmp1=get_NBR_ests(n_iter=20,Y=Y,X=X,method='PG')
    tmp2=get_NBR_ests(n_iter=20,Y=Y,X=X,method='MCMC')
    tmp3=get_NBR_ests(n_iter=20,Y=Y,X=X,method='glm.nb')
    
    intercepts=c(tmp1[[1]][1],tmp2[[1]][1],tmp3[[1]][1])
    slopes=c(tmp1[[1]][2],tmp2[[1]][2],tmp3[[1]][2])
    phis=c(tmp1[[2]],tmp2[[2]],tmp3[[2]])
    llh=c(tmp1[[3]],tmp2[[3]],tmp3[[3]])
    
    Y.list=c(Y.list,list(Y))
    X.list=c(X.list,list(X))
    index.list=c(index.list,i)
    
    intercept.mat=rbind(intercept.mat,slopes)
    slope.mat=rbind(slope.mat,slopes)
    llh.mat=rbind(llh.mat,llh)
    r.mat=rbind(r.mat,phis)
  }
  
  colnames(intercept.mat)=c('PG','MCMC','glmnb')
  rownames(intercept.mat)=index.list
  
  colnames(slope.mat)=c('PG','MCMC','glmnb')
  rownames(slope.mat)=index.list
  
  colnames(llh.mat)=c('PG','MCMC','glmnb')
  rownames(llh.mat)=index.list
  
  colnames(r.mat)=c('PG','MCMC','glmnb')
  rownames(r.mat)=index.list
  
  output=list(Y.list,X.list,index.list,intercept.mat,slope.mat,r.mat,llh.mat)
  
  names(output)=c('Y','X','index','intercept','slope','phi','llh')
  
  return(output)
}


#(1) Get segments of NB distributed count data

get_NB_seg_data<-function() { #Get segmented data
  sim.dat=get_NB_simulated_data(beta.slope.list=c(-0.25,0.1,-0.1,0.25), #True slopes
                                beta.intercept.list=c(2,0.5,1,-10), #True intercepts
                                phi.list=c(5,5,5,5), #True dispersion
                                smooth=0)
  yy=sim.dat[[1]]
  nc.id=sim.dat[[3]]
  
  NN=length(yy)
  XX=cbind(rep(1,NN),1:NN)
  dat=data.frame(yy,XX,rep(1,NN))
  colnames(dat)=c('adj_gene_exp','intercept','id','SS')
  
  zeta=rep(0,NN)
  zeta[c(1,nc.id)]=1
  gg=get_gg_from_zeta(zeta)
  seg.dat=get_segmented_dat(dat=dat,gg=gg)
  return(seg.dat)
}

#(2) Get NB data with CP

get_NB_simulated_data<-function(beta.slope.list,
                                beta.intercept.list,
                                phi.list,
                                very_smooth=0,
                                smooth=1){
  nn.list=c(40,55,70)
  nc.list=c(list(c(1)),list(c(1:2)),list(c(1:3)))
  
  #Choose sample size
  index=sample(1:length(nn.list),1)
  nn=nn.list[index]
  
  #Choose number of change point
  nc=sample(nc.list[[index]],1)
  
  #get change points
  nc.id=sample(4:(nn-3),nc)
  if(length(nc.id)>1) {
    nc.id=nc.id[order(nc.id)]
    nc.lag=nc.id[2:length(nc.id)]-nc.id[1:(length(nc.id)-1)]
    nc.rm=which(nc.lag<3)
    while(length(nc.rm)>0) {
      nc.id=sample(4:(nn-3),nc)
      nc.id=nc.id[order(nc.id)]
      nc.lag=nc.id[2:length(nc.id)]-nc.id[1:(length(nc.id)-1)]
      nc.rm=which(nc.lag<5)
    }
  }
  
  #Get segmented dat
  
  segs=c(1,nc.id,nn)
  segs=cbind(segs[1:(length(segs)-1)],c((segs[2:(length(segs)-1)]-1),nn))
  n.segs=nrow(segs)
  gg=c()
  for (i in 1:n.segs) {
    gg=c(gg,rep(i,segs[i,2]-segs[i,1]+1))
  }
  
  if (smooth==1) {
    for (i in 2:(length(nc.id)+1)) {
      l.part=beta.intercept.list[(i-1)]+(nc.id[(i-1)]-1)*beta.slope.list[i-1]
      if (l.part<0) {l.part=0}
      left=l.part+log(phi.list[i-1])
      target.mean=left-(nc.id[(i-1)]-1)*beta.slope.list[i]-log(phi.list[i])
      
      #Whether add noise
      if (very_smooth==0) {
        beta.intercept.list[i]=rnorm(1,mean=target.mean,sd=abs(target.mean)/2)
      } else {beta.intercept.list[i]=target.mean}
    }
  }
  
  betas=cbind(beta.intercept.list,beta.slope.list)[1:(nc+1),]
  
  #Covariate matrix
  XX=cbind(rep(1,nn),1:nn)
  seg.XX=get_segmented_dat(gg=gg,dat=XX)
  
  yy1=c()
  for (i in 1:n.segs) {
    tmp=seg.XX[[i]] %*% betas[i,]
    tmp.nn=length(tmp)
    #print('eta')
    #print(tmp)
    pp=1/(1+exp(-tmp))
    #print('pp')
    #print(pp)
    mean.list=pp/(1-pp)*phi.list[i]
    #print('mean')
    #print(mean.list)
    #print('+++++++++++++++++++')
    tmp.yy=sapply(pp,function(x) rnbinom(1,size=phi.list[i],prob=1-x))
    yy1=c(yy1,tmp.yy)
  }
  
  
  plot(yy1,type='l')
  abline(v=nc.id,col='red')
  XX=cbind(rep(1,nn),1:nn)
  return(list(yy1,XX,nc.id,betas,phi.list[1:n.segs]))
} #intercept not fixed

get_segmented_dat<-function(gg,dat){
  n.seg=length(table(gg))
  ind.list=lapply(1:n.seg,function(x) which(gg==x))
  dat.list=lapply(ind.list,function(x) dat[x,])
  return(dat.list)
}

get_gg_from_zeta<-function(zeta){
  nn=length(zeta)
  nc.id=which(zeta==1)
  if (length(nc.id)==1) {gg=rep(1,nn)} else{
    segs=c(nc.id,nn)
    segs=cbind(segs[1:(length(segs)-1)],c((segs[2:(length(segs)-1)]-1),nn))
    n.segs=nrow(segs)
    gg=c()
    for (i in 1:n.segs) {
      gg=c(gg,rep(i,segs[i,2]-segs[i,1]+1))
    }
  }
  return(gg)
}



###2. Check performance of finding CPs

# Function to generate group IDs based on total number of data points and change points
generate_group_ids <- function(nn, change_points) {
  # Initialize group id vector
  group_ids <- numeric(nn)
  
  # Set the first segment
  start <- 1
  group <- 1
  
  # Loop through each change point to assign group ids
  for (cp in change_points) {
    group_ids[start:(cp - 1)] <- group
    group <- group + 1
    start <- cp
  }
  
  # Assign the last segment (after the last change point)
  group_ids[start:nn] <- group
  
  return(group_ids)
}


evaluate_CP<-function(nn,true.cp,estimated.cp){
  true.gg=generate_group_ids(nn,true.cp)
  if (is.null(estimated.cp)) {estimated.gg=rep(1,nn)} else {estimated.gg=generate_group_ids(nn,estimated.cp)}
  
  #aa: be in the same segment in both the true and estimated grouping
  #bb: be in differents segments in true but same segment in estimated grouping
  #cc: be in same segments in true but different segment in estimated grouping
  #dd: be in differents segments in both true and estimated grouping
  
  prs=combn(nn,m=2)
  tt=ncol(prs)
  seg.true=apply(prs,2,function(x) length(unique(true.gg[x])))
  seg.est=apply(prs,2,function(x) length(unique(estimated.gg[x])))
  aa=length(which(seg.true==1 & seg.est==1))
  bb=length(which(seg.true==2 & seg.est==1))
  cc=length(which(seg.true==1 & seg.est==2))
  dd=length(which(seg.true==2 & seg.est==2))
  
  ARI=(tt*(aa+dd)-((aa+bb)*(aa+cc)+(cc+dd)*(bb+dd)))/(tt^2-((aa+bb)*(aa+cc)+(cc+dd)*(bb+dd)))
  return(ARI)
}




check_CP<-function(dat,output) {
  ARI=evaluate_CP(nn=length(dat[[1]]),true.cp = dat[[3]],estimated.cp=output$change_points)
  return(ARI)
}

check_CP_performance<-function(dat.list,out.id,CP.output){
  ari.mat=c()
  for (i in 1:length(CP.output)) {
    tmp.dat=dat=dat.list[[i]]
    tmp.output=CP.output[[i]][out.id]
    tmp.ARI=unlist(lapply(tmp.output,function(x) check_CP(dat=dat.list[[i]],output=x)))
    ari.mat=rbind(ari.mat,tmp.ARI)
  }
  colnames(ari.mat)=c('PG','MCMC','glm.nb')
  return(ari.mat)
}



CP_other_methods<-function(dat) {
  exp=dat[[1]]
  nn=length(exp)
  true.cp=dat[[3]]
  id=1:nn
  
  #Bayesmile
  zeta.bs=detect_changepoint_poi_model(input.df=data.frame(id,exp), 
                                       h_vec=c(100000, 10),
                                       NITER = 20000,
                                       seed = 123)
  est.cp=which(zeta.bs==1)[-1]
  ARI.bs=evaluate_CP(nn,true.cp,est.cp)
  
  #cpm
  est.cp=processStream(exp, cpmType = "Exponential")$changePoints
  ARI.cpm=evaluate_CP(nn,true.cp,est.cp)
  
  
  #changepoint
  tmp=cpt.mean(exp)
  est.cp=cpts(tmp)
  ARI.cp=evaluate_CP(nn,true.cp,est.cp)
  
  #changepoint.np
  tmp=cpt.np(exp)
  est.cp=cpts(tmp)
  ARI.cp.np=evaluate_CP(nn,true.cp,est.cp)
  
  return(c(ARI.bs,ARI.cpm,ARI.cp,ARI.cp.np))
}





#() Get Poisson distributed count data
#Get Poisson data
get_poisson_simulated_data<-function(very_smooth=1){
  nn.list=c(40,55,70)
  nc.list=c(list(c(1)),list(c(1:2)),list(c(1:3)))
  
  #Choose sample size
  index=sample(1:length(nn.list),1)
  nn=nn.list[index]
  
  #Choose number of change point
  nc=sample(nc.list[[index]],1)
  
  #get change points
  nc.id=sample(4:(nn-3),nc)
  if(length(nc.id)>1) {
    nc.id=nc.id[order(nc.id)]
    nc.lag=nc.id[2:length(nc.id)]-nc.id[1:(length(nc.id)-1)]
    nc.rm=which(nc.lag<3)
    while(length(nc.rm)>0) {
      nc.id=sample(4:(nn-3),nc)
      nc.id=nc.id[order(nc.id)]
      nc.lag=nc.id[2:length(nc.id)]-nc.id[1:(length(nc.id)-1)]
      nc.rm=which(nc.lag<5)
    }
  }
  
  #Get segmented dat
  
  segs=c(1,nc.id,nn)
  segs=cbind(segs[1:(length(segs)-1)],c((segs[2:(length(segs)-1)]-1),nn))
  n.segs=nrow(segs)
  gg=c()
  for (i in 1:n.segs) {
    gg=c(gg,rep(i,segs[i,2]-segs[i,1]+1))
  }
  
  
  #Set beta-coefficients
  beta.slope.list=c(-0.2,0.15,-0.2,0.1)
  beta.intercept.list=c(4,0,0,0)
  for (i in 2:(length(nc.id)+1)) {
    l.part=beta.intercept.list[(i-1)]+(nc.id[(i-1)]-1)*beta.slope.list[i-1]
    if (l.part<0) {l.part=0}
    target.mean=l.part-(nc.id[(i-1)]-1)*beta.slope.list[i]
    if (very_smooth==0) {
      beta.intercept.list[i]=rnorm(1,mean=target.mean,sd=abs(target.mean)/3)
    } else {beta.intercept.list[i]=target.mean}
  }
  betas=cbind(beta.intercept.list,beta.slope.list)[1:(nc+1),]
  
  #Covariate matrix
  XX=cbind(rep(1,nn),1:nn)
  seg.XX=get_segmented_dat(gg=gg,dat=XX)
  
  yy1=c()
  for (i in 1:n.segs) {
    mu=exp(seg.XX[[i]] %*% betas[i,])
    #print(seg.XX[[i]] %*% betas[i,])
    #print(mu)
    yy1=c(yy1,rpois(n=length(mu),alpha=mu))
  }
  return(list(yy1,XX,nc.id,betas))
}









#####Functions with notations#####################################################################################################################################################################
# OP with notations
# OP with notations
optimal_partitioning_with_notes <- function(n_iter=20,Y, X, method, alpha = 1) {
  n <- length(Y)  # Length of the data
  
  # Initialize cost array and partition point array
  cost <- rep(Inf, n + 1)   # Cost for partitioning up to each point
  partition <- rep(0, n + 1)  # Best partition point for each endpoint
  
  # Base case: no cost for zero elements
  cost[1] <- 0
  
  # Loop over all possible endpoints
  for (end in 2:n) {
    # Loop over all possible start points for each endpoint
    for (start in 1:(end - 1)) {  # Ensure at least 2 points in a segment
      if (start==2) {next}
      if (end - start + 1 >= 2) {  # Condition to ensure minimum segment length of 2
        
        # Calculate the cost of the current segment
        current_cost <- segment_cost(n_iter=n_iter,Y=Y, X=X, start=start, end=end, method=method, alpha=alpha)[1]
        
        # Print out the intermediate values
        print(paste("Start:", start, "End:", end, "Current cost:", current_cost))
        
        # Update the cost and partition arrays if we found a better partition
        if (cost[start] + current_cost < cost[end + 1]) {
          cost[end + 1] <- cost[start] + current_cost
          partition[end + 1] <- start
          
          # Print out updated cost and partition information
          print(paste("Updated cost for endpoint", end, ": (cost[end+1])", cost[end + 1]))
          print(paste("Partition point at", end + 1, "set to", start))
          print('============================')
        }
      }
    }
  }
  
  # Reconstruct the partitioning solution
  segmentation <- list()
  change_points <- c()  # Store the change points (excluding 1)
  end <- n
  while (end > 0) {
    start <- partition[end + 1]
    segmentation <- append(segmentation, list(c(start, end)))
    if (start != 1) {  # Only track change points excluding 1
      change_points <- c(change_points, start)
    }
    # Print out the reconstructed segmentation
    print(paste("Segment from", start, "to", end))
    end <- start - 1
  }
  
  return(list(total_cost = cost[n + 1], segmentation = rev(segmentation), change_points = rev(change_points)))
}


# PELT with notations
optimal_partitioning_pelt_with_notes <- function(n_iter=20,Y, X, method, alpha = 1,K=0) {
  n <- length(Y)  # Length of the data
  
  # Initialize cost array and partition point array
  cost <- rep(Inf, n + 1)   # Cost for partitioning up to each point
  partition <- rep(0, n + 1)  # Best partition point for each endpoint
  
  # Base case: no cost for zero elements
  cost[1] <- 0
  
  # List of potential change points
  R=list()
  R[[1]]=0
  R[[2]]=1
  R[[3]]=1
  R[[4]]=c(1,3)
  R[[5]]=c(1,3,4)
  cost.mat=matrix(NA,n,n)
  
  # Loop over all possible endpoints
  for (end in 2:n) {
    #1. Loop over all possible start points for each endpoint to get minimized cost for cost[end+1]
    for (start in R[[end]]) {  # Ensure at least 2 points in a segment
      if (start==2) {next}
      if (end - start + 1 >= 2) {  # Condition to ensure minimum segment length of 2
        
        # Calculate the cost of the current segment
        current_cost <- segment_cost(n_iter=n_iter,Y=Y, X=X, start=start, end=end, method=method, alpha=alpha)
        current_cost_aic=current_cost[1]
        current_cost_llh=current_cost[2]
        cost.mat[start,end]=current_cost_llh
        
        # Print out the intermediate values
        print(paste("Start:", start, "End:", end, "Current cost with penalty:", current_cost_aic))
        
        # Update the cost and partition arrays if we found a better partition
        if (cost[start] + current_cost_aic < cost[end + 1]) {
          cost[end + 1] <- cost[start] + current_cost_aic
          partition[end + 1] <- start
          
          # Print out updated cost and partition information
          print(paste("Updated cost for endpoint", end, ": (cost[end+1])", cost[end + 1]))
          print(paste("Partition point at", end + 1, "set to", start))
          print('============================')
        }
      }
    }
    
    #2. Pruning
    new_R=c() #change points that are kept after pruning for data point 'end'
    if (end>=5 & end < n)  { #At this time, we have cost[end+1] for F(end)
      tmp.cp=R[[end]]
      for (start in tmp.cp[2:length(tmp.cp)]) {
        pruning_cost=cost.mat[start,end]
        f_t=cost[start]
        
        # Print pruning cost and condition details
        cat("Pruning cost (start to end):", pruning_cost, "\n")
        cat("F(t):", cost[start], "F(s):", cost[end+1], "\n")
        
        # Pruning logic: keep start points that reduce the cost
        if (f_t + pruning_cost + K < cost[end+1]) {
          new_R <- c(new_R, start)
          cat("Keeping start", start, "for endpoint", end, "\n")
        } else {
          cat("Pruning start", start, "for endpoint", end, "\n")
        }
      }
      
      # Update the set of candidate changepoints for the next endpoint
      R[[end+1]] <- unique(c(1,new_R, end))
      cat("New candidate changepoints for R[[", end + 1, "]]:", R[[end + 1]], "\n")
    }
    
  }
  
  # Reconstruct the partitioning solution
  segmentation <- list()
  change_points <- c()  # Store the change points (excluding 1)
  end <- n
  while (end > 0) {
    start <- partition[end + 1]
    segmentation <- append(segmentation, list(c(start, end)))
    if (start != 1) {  # Only track change points excluding 1
      change_points <- c(change_points, start)
    }
    # Print out the reconstructed segmentation
    print(paste("Segment from", start, "to", end))
    end <- start - 1
  }
  
  return(list(total_cost = cost[n + 1], segmentation = rev(segmentation), change_points = rev(change_points)))
}





stan_code_p <- "
data {
  int<lower=0> N;          // Number of observations
  int<lower=0> K;          // Number of predictors (including intercept)
  matrix[N, K] X;          // Predictor matrix
  int<lower=0> Y[N];       // Outcome variable (count data)
}

parameters {
  vector[K] beta;          // Regression coefficients
  real<lower=0> r;         // Dispersion parameter
}

transformed parameters {
  vector[N] p;             // Probability parameter for negative binomial
  vector[N] mu;            // Mean parameter derived from p and r
  
  // Calculate p using the logit link
  p = inv_logit(X * beta); 

  // Calculate mu based on r and p
  mu = r * (1.0 / (1.0 - p));
}

model {
  // Priors
  beta ~ normal(0, 10);    // Prior for regression coefficients
  r ~ gamma(2, 0.1);       // Prior for dispersion parameter

  // Likelihood
  Y ~ neg_binomial_2(mu, r); // Negative binomial likelihood using mu
}
"









