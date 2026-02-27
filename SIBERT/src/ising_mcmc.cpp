// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

//' Compute Conditional Probability for Ising Model
//'
//' Computes the conditional probability for spin x_i = +1 given neighbors.
//'
//' @param J_k Coupling parameter
//' @param neighbor_sum Sum of neighboring spins
//' @return Probability of x_i = +1
//' @keywords internal
// [[Rcpp::export]]
double compute_prob_cpp(double J_k, double neighbor_sum) {
  double prob_1 = 1.0 / (1.0 + exp(-2.0 * J_k * neighbor_sum));
  return prob_1;
}


//' Calculate Sum of Spin Interactions
//'
//' Calculates the total sum of interactions between neighboring spins.
//'
//' @param x Spin configuration vector (-1 or 1)
//' @param nb_list List of neighbor indices (0-indexed)
//' @return Total sum of spin interactions (divided by 2 for double counting)
//' @keywords internal
// [[Rcpp::export]]
double calculate_sum_cpp(NumericVector x, List nb_list) {
  double total_sum = 0.0;

  for (int i = 0; i < x.size(); ++i) {
    NumericVector neighbors = nb_list[i];
    for (int j = 0; j < neighbors.size(); ++j) {
      total_sum += x[i] * x[neighbors[j]];
    }
  }
  total_sum /= 2.0;
  return total_sum;
}


//' Monte Carlo Sampling for Ising Model
//'
//' Performs Metropolis-Hastings sampling for the Ising model to estimate
//' the derivative of the log partition function.
//'
//' @param init_x Initial spin configuration
//' @param nb_list_cpp Neighbor list (0-indexed for C++)
//' @param J_k Coupling parameter
//' @param n_flip Number of proposed spin flips
//' @param M Number of samples for averaging
//'
//' @return List containing:
//'   - sum_list: Vector of running sums
//'   - partial_der: Estimated partial derivative
//'
//' @keywords internal
// [[Rcpp::export]]
List MC_sampling_x_cpp(NumericVector init_x, List nb_list_cpp, double J_k,
                        int n_flip, int M) {
  int n = init_x.size();
  double sum_old = calculate_sum_cpp(init_x, nb_list_cpp);
  NumericVector sum_list(n_flip);

  for (int i = 0; i < n_flip; i++) {
    NumericVector new_x = clone(init_x);
    int idx = floor(R::runif(0, n));
    new_x[idx] = -1 * new_x[idx];

    double sum_new = calculate_sum_cpp(new_x, nb_list_cpp);

    if (exp(J_k * (sum_new - sum_old)) > R::runif(0, 1)) {
      init_x = new_x;
      sum_old = sum_new;
    }

    sum_list[i] = sum_old;
  }

  double partial_der = 0;
  for (int i = n_flip - M; i < n_flip; i++) {
    partial_der += sum_list[i];
  }
  partial_der /= M;

  return List::create(Named("sum_list") = sum_list,
                      Named("partial_der") = partial_der);
}
