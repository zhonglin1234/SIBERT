#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector rpg_R(int n,double h,double z){
  Function f("rpg");   
  return f(n, h, z);
}

// [[Rcpp::export]]
arma::vec cpp_update_beta(arma::mat XX, 
                          arma::vec y_obs, 
                          arma::vec b,
                          arma::mat B, 
                          arma::vec beta, 
                          double phi,
                          int niter_beta) {
  
  int n = y_obs.size();
  int p = beta.size();
  
  arma::vec psi(n, arma::fill::zeros);
  arma::vec omega(n, arma::fill::zeros);
  
  arma::vec k(n, arma::fill::zeros);
  arma::mat Omega = arma::zeros<arma::mat>(n, n);
  arma::mat V = arma::zeros<arma::mat>(p, p);
  arma::vec mm(n, arma::fill::zeros);
  
  arma::mat beta_mat(niter_beta, p);
  
  
  // Generating multivariate normal random variates
  for (int i = 0; i < niter_beta; i++) {
    psi = XX * beta;
    
    for (int j = 0; j < n; j++) {
      omega[j] = rpg_R(1,phi + y_obs[j],psi[j])[0];
    }
    
    // if omega is all 0, then need to repeat the procedure until not
    while (arma::sum(omega) == 0) { 
      for (int j = 0; j < n; j++) {
        omega[j] = omega[j] = rpg_R(1,phi + y_obs[j],psi[j])[0];
      }
    }
    
    k = (y_obs - phi) / 2;
    Omega = arma::diagmat(omega);
    V = arma::inv(XX.t() * Omega * XX + arma::inv(B));
    mm = V * (XX.t() * k + arma::inv(B) * b);
    
    // Generating multivariate normal random variates
    Function mvrnorm("mvrnorm");
    beta = Rcpp::as<arma::vec>(mvrnorm(1, mm, V));
    beta_mat.row(i) = beta.t();
  }
  size_t start_row = beta_mat.n_rows / 2;
  
  //Extract the second half of the rows
  arma::mat second_half = beta_mat.rows(start_row, beta_mat.n_rows - 1);
  arma::vec post_mean = arma::mean(second_half, 0).t();
  return post_mean;
}


// [[Rcpp::export]]
arma::mat cpp_update_beta_mat(arma::mat XX, 
                              arma::vec y_obs, 
                              arma::vec b,
                              arma::mat B, 
                              arma::vec beta, 
                              double phi,
                              int niter_beta) {
  
  int n = y_obs.size();
  int p = beta.size();
  
  arma::vec psi(n, arma::fill::zeros);
  arma::vec omega(n, arma::fill::zeros);
  
  arma::vec k(n, arma::fill::zeros);
  arma::mat Omega = arma::zeros<arma::mat>(n, n);
  arma::mat V = arma::zeros<arma::mat>(p, p);
  arma::vec mm(n, arma::fill::zeros);
  
  arma::mat beta_mat(niter_beta, p);
  
  
  // Generating multivariate normal random variates
  for (int i = 0; i < niter_beta; i++) {
    psi = XX * beta;
    
    for (int j = 0; j < n; j++) {
      omega[j] = rpg_R(1,phi + y_obs[j],psi[j])[0];
    }
    
    // if omega is all 0, then need to repeat the procedure until not
    while (arma::sum(omega) == 0) { 
      for (int j = 0; j < n; j++) {
        omega[j] = omega[j] = rpg_R(1,phi + y_obs[j],psi[j])[0];
      }
    }
    
    k = (y_obs - phi) / 2;
    Omega = arma::diagmat(omega);
    V = arma::inv(XX.t() * Omega * XX + arma::inv(B));
    mm = V * (XX.t() * k + arma::inv(B) * b);
    
    // Generating multivariate normal random variates
    Function mvrnorm("mvrnorm");
    beta = Rcpp::as<arma::vec>(mvrnorm(1, mm, V));
    beta_mat.row(i) = beta.t();
  }
  return beta_mat;
}


// [[Rcpp::export]]
double cpp_update_phi_MCMC(arma::mat XX, 
                           arma::vec y_obs, 
                           arma::vec beta, 
                           double a_phi, // shape parameter, a small one
                           double b_phi, // rate parameter, a small one
                           double sd_phi,
                           double phi,
                           int niter_phi) {
  
  int n = y_obs.size();
  
  arma::vec psi = XX * beta;
  arma::vec pp = 1 / (1 + exp(-psi));
  
  // Initialize the storage for phi 
  arma::vec phi_list(niter_phi, arma::fill::zeros);
  
  // Initialize log likelihoods
  double llh_old = -10.0;
  double llh_new = -10.0;  
  
  // Initialize log_m_MH
  double log_m_MH = 0.0;
  
  // Initialize vectors for lambda
  arma::vec lambda(n, arma::fill::zeros);
  arma::vec lambda_new(n, arma::fill::zeros);
  
  // Initialize phi_new
  double phi_new = 1.0;
  
  for (int i = 0; i < niter_phi; i++) {
    
    // Calculate lambda and llh_old for current phi
    lambda = phi * pp / (1 - pp);
    llh_old = arma::accu(lgamma(y_obs + phi) - lgamma(y_obs + 1) - lgamma(phi) +
      phi * (log(phi) - log(lambda + phi)) +
      y_obs % (log(lambda) - log(lambda + phi)));
    
    // Propose new phi
    phi_new = phi + R::rnorm(0, sd_phi);
    
    // Ensure phi_new is non-negative
    while (phi_new < 0) {
      phi_new = phi + R::rnorm(0, sd_phi);
    }
    
    // Calculate lambda_new and llh_new for proposed phi_new
    lambda_new = phi_new * pp / (1 - pp);
    
    llh_new = arma::accu(lgamma(y_obs + phi_new) - lgamma(y_obs + 1) - lgamma(phi_new) +
      phi_new * (log(phi_new) - log(lambda_new + phi_new)) +
      y_obs % (log(lambda_new) - log(lambda_new + phi_new)));
    
    // Compute the log of the Metropolis-Hastings ratio
    log_m_MH = llh_new - llh_old + 
      R::dgamma(phi_new, a_phi, 1.0 / b_phi, true) - 
      R::dgamma(phi, a_phi, 1.0 / b_phi, true);
    
    // Acceptance step: accept with probability min(1, exp(log_m_MH))
    double acceptance_prob = exp(log_m_MH);
    if (R::runif(0, 1) < acceptance_prob) {
      phi = phi_new;  // Accept the new value
    }
    
    // Store the accepted phi value
    phi_list[i] = phi;
  }
  
  // Compute the mean of the second half of the samples
  arma::vec post_phis = phi_list.subvec(niter_phi / 2, niter_phi - 1);
  return mean(post_phis);
}


// [[Rcpp::export]]
arma::vec cpp_update_phi_MCMC_vec(arma::mat XX, 
                                  arma::vec y_obs, 
                                  arma::vec beta, 
                                  double a_phi, // shape parameter, a small one
                                  double b_phi, // rate parameter, a small one
                                  double sd_phi,
                                  double phi,
                                  int niter_phi) {
  
  int n = y_obs.size();
  
  // Compute psi and pp for all observations
  arma::vec psi = XX * beta;
  arma::vec pp = 1 / (1 + exp(-psi));
  
  // Initialize the storage for phi 
  arma::vec phi_list(niter_phi, arma::fill::zeros);
  
  // Initialize likelihoods
  double llh_old = -10.0;
  double llh_new = -10.0;
  
  // Initialize log_m_MH
  double log_m_MH = 0.0;
  
  // Initialize vectors lambda and lambda_new
  arma::vec lambda(n, arma::fill::zeros);
  arma::vec lambda_new(n, arma::fill::zeros);
  
  // Initialize phi_new
  double phi_new = 1.0;
  
  for (int i = 0; i < niter_phi; i++) {
    
    // Compute lambda and llh_old for current phi
    lambda = phi * pp / (1 - pp);
    llh_old = arma::accu(lgamma(y_obs + phi) - lgamma(y_obs + 1) - lgamma(phi) +
      phi * (log(phi) - log(lambda + phi)) +
      y_obs % (log(lambda) - log(lambda + phi)));
    
    // Propose new phi
    phi_new = phi + R::rnorm(0, sd_phi);
    
    // Ensure phi_new is non-negative
    while (phi_new < 0) {
      phi_new = phi + R::rnorm(0, sd_phi);
    }
    
    // Compute lambda_new and llh_new for proposed phi_new
    lambda_new = phi_new * pp / (1 - pp);
    llh_new = arma::accu(lgamma(y_obs + phi_new) - lgamma(y_obs + 1) - lgamma(phi_new) +
      phi_new * (log(phi_new) - log(lambda_new + phi_new)) +
      y_obs % (log(lambda_new) - log(lambda_new + phi_new)));
    
    // Compute the log of the Metropolis-Hastings ratio
    log_m_MH = llh_new - llh_old + 
      R::dgamma(phi_new, a_phi, 1.0 / b_phi, true) - 
      R::dgamma(phi, a_phi, 1.0 / b_phi, true);
    
    // Acceptance step: accept with probability min(1, exp(log_m_MH))
    double acceptance_prob = exp(log_m_MH);
    if (R::runif(0, 1) < acceptance_prob) {
      phi = phi_new;  // Accept the new value
    }
    
    // Store the accepted or current phi value in phi_list
    phi_list[i] = phi;
  }
  
  // Return the full vector of sampled phi values
  return phi_list;
}

