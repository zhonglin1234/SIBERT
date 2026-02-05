// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

//' Sample from Polya-Gamma Distribution
//'
//' Wrapper to call R's rpg function from the BayesLogit package.
//'
//' @param n Number of samples
//' @param h First parameter
//' @param z Second parameter
//' @return Numeric vector of Polya-Gamma samples
//' @keywords internal
// [[Rcpp::export]]
NumericVector rpg_R(int n, double h, double z) {
  Function f("rpg");
  return f(n, h, z);
}


//' Update Beta Coefficients using Polya-Gamma Augmentation
//'
//' Performs Gibbs sampling for beta coefficients in negative binomial
//' regression using Polya-Gamma data augmentation.
//'
//' @param XX Design matrix (n x p)
//' @param y_obs Response vector (counts)
//' @param b Prior mean for beta
//' @param B Prior covariance for beta
//' @param beta Initial beta values
//' @param phi Dispersion parameter
//' @param niter_beta Number of iterations
//'
//' @return Posterior mean of beta coefficients
//' @keywords internal
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

  for (int i = 0; i < niter_beta; i++) {
    psi = XX * beta;

    for (int j = 0; j < n; j++) {
      omega[j] = rpg_R(1, phi + y_obs[j], psi[j])[0];
    }

    while (arma::sum(omega) == 0) {
      for (int j = 0; j < n; j++) {
        omega[j] = rpg_R(1, phi + y_obs[j], psi[j])[0];
      }
    }

    k = (y_obs - phi) / 2;
    Omega = arma::diagmat(omega);
    V = arma::inv(XX.t() * Omega * XX + arma::inv(B));
    mm = V * (XX.t() * k + arma::inv(B) * b);

    Function mvrnorm("mvrnorm");
    beta = Rcpp::as<arma::vec>(mvrnorm(1, mm, V));
    beta_mat.row(i) = beta.t();
  }

  size_t start_row = beta_mat.n_rows / 2;
  arma::mat second_half = beta_mat.rows(start_row, beta_mat.n_rows - 1);
  arma::vec post_mean = arma::mean(second_half, 0).t();

  return post_mean;
}


//' Update Beta Coefficients - Return Full Chain
//'
//' Same as cpp_update_beta but returns the full MCMC chain.
//'
//' @param XX Design matrix (n x p)
//' @param y_obs Response vector (counts)
//' @param b Prior mean for beta
//' @param B Prior covariance for beta
//' @param beta Initial beta values
//' @param phi Dispersion parameter
//' @param niter_beta Number of iterations
//'
//' @return Matrix of beta samples (niter_beta x p)
//' @keywords internal
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

  for (int i = 0; i < niter_beta; i++) {
    psi = XX * beta;

    for (int j = 0; j < n; j++) {
      omega[j] = rpg_R(1, phi + y_obs[j], psi[j])[0];
    }

    while (arma::sum(omega) == 0) {
      for (int j = 0; j < n; j++) {
        omega[j] = rpg_R(1, phi + y_obs[j], psi[j])[0];
      }
    }

    k = (y_obs - phi) / 2;
    Omega = arma::diagmat(omega);
    V = arma::inv(XX.t() * Omega * XX + arma::inv(B));
    mm = V * (XX.t() * k + arma::inv(B) * b);

    Function mvrnorm("mvrnorm");
    beta = Rcpp::as<arma::vec>(mvrnorm(1, mm, V));
    beta_mat.row(i) = beta.t();
  }

  return beta_mat;
}


//' Update Dispersion Parameter using MCMC
//'
//' Performs Metropolis-Hastings sampling for the dispersion parameter
//' in negative binomial regression.
//'
//' @param XX Design matrix
//' @param y_obs Response vector
//' @param beta Regression coefficients
//' @param a_phi Shape parameter for gamma prior
//' @param b_phi Rate parameter for gamma prior
//' @param sd_phi Standard deviation for proposal distribution
//' @param phi Initial value
//' @param niter_phi Number of iterations
//'
//' @return Posterior mean of phi
//' @keywords internal
// [[Rcpp::export]]
double cpp_update_phi_MCMC(arma::mat XX,
                           arma::vec y_obs,
                           arma::vec beta,
                           double a_phi,
                           double b_phi,
                           double sd_phi,
                           double phi,
                           int niter_phi) {

  int n = y_obs.size();

  arma::vec psi = XX * beta;
  arma::vec pp = 1 / (1 + exp(-psi));

  arma::vec phi_list(niter_phi, arma::fill::zeros);

  double llh_old = -10.0;
  double llh_new = -10.0;
  double log_m_MH = 0.0;

  arma::vec lambda(n, arma::fill::zeros);
  arma::vec lambda_new(n, arma::fill::zeros);

  double phi_new = 1.0;

  for (int i = 0; i < niter_phi; i++) {
    lambda = phi * pp / (1 - pp);
    llh_old = arma::accu(lgamma(y_obs + phi) - lgamma(y_obs + 1) - lgamma(phi) +
      phi * (log(phi) - log(lambda + phi)) +
      y_obs % (log(lambda) - log(lambda + phi)));

    phi_new = phi + R::rnorm(0, sd_phi);

    while (phi_new < 0) {
      phi_new = phi + R::rnorm(0, sd_phi);
    }

    lambda_new = phi_new * pp / (1 - pp);
    llh_new = arma::accu(lgamma(y_obs + phi_new) - lgamma(y_obs + 1) - lgamma(phi_new) +
      phi_new * (log(phi_new) - log(lambda_new + phi_new)) +
      y_obs % (log(lambda_new) - log(lambda_new + phi_new)));

    log_m_MH = llh_new - llh_old +
      R::dgamma(phi_new, a_phi, 1.0 / b_phi, true) -
      R::dgamma(phi, a_phi, 1.0 / b_phi, true);

    double acceptance_prob = exp(log_m_MH);
    if (R::runif(0, 1) < acceptance_prob) {
      phi = phi_new;
    }

    phi_list[i] = phi;
  }

  arma::vec post_phis = phi_list.subvec(niter_phi / 2, niter_phi - 1);
  return mean(post_phis);
}


//' Update Dispersion Parameter - Return Full Chain
//'
//' Same as cpp_update_phi_MCMC but returns the full MCMC chain.
//'
//' @param XX Design matrix
//' @param y_obs Response vector
//' @param beta Regression coefficients
//' @param a_phi Shape parameter for gamma prior
//' @param b_phi Rate parameter for gamma prior
//' @param sd_phi Standard deviation for proposal
//' @param phi Initial value
//' @param niter_phi Number of iterations
//'
//' @return Vector of phi samples
//' @keywords internal
// [[Rcpp::export]]
arma::vec cpp_update_phi_MCMC_vec(arma::mat XX,
                                  arma::vec y_obs,
                                  arma::vec beta,
                                  double a_phi,
                                  double b_phi,
                                  double sd_phi,
                                  double phi,
                                  int niter_phi) {

  int n = y_obs.size();

  arma::vec psi = XX * beta;
  arma::vec pp = 1 / (1 + exp(-psi));

  arma::vec phi_list(niter_phi, arma::fill::zeros);

  double llh_old = -10.0;
  double llh_new = -10.0;
  double log_m_MH = 0.0;

  arma::vec lambda(n, arma::fill::zeros);
  arma::vec lambda_new(n, arma::fill::zeros);

  double phi_new = 1.0;

  for (int i = 0; i < niter_phi; i++) {
    lambda = phi * pp / (1 - pp);
    llh_old = arma::accu(lgamma(y_obs + phi) - lgamma(y_obs + 1) - lgamma(phi) +
      phi * (log(phi) - log(lambda + phi)) +
      y_obs % (log(lambda) - log(lambda + phi)));

    phi_new = phi + R::rnorm(0, sd_phi);

    while (phi_new < 0) {
      phi_new = phi + R::rnorm(0, sd_phi);
    }

    lambda_new = phi_new * pp / (1 - pp);
    llh_new = arma::accu(lgamma(y_obs + phi_new) - lgamma(y_obs + 1) - lgamma(phi_new) +
      phi_new * (log(phi_new) - log(lambda_new + phi_new)) +
      y_obs % (log(lambda_new) - log(lambda_new + phi_new)));

    log_m_MH = llh_new - llh_old +
      R::dgamma(phi_new, a_phi, 1.0 / b_phi, true) -
      R::dgamma(phi, a_phi, 1.0 / b_phi, true);

    double acceptance_prob = exp(log_m_MH);
    if (R::runif(0, 1) < acceptance_prob) {
      phi = phi_new;
    }

    phi_list[i] = phi;
  }

  return phi_list;
}
