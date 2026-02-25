#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
double compute_prob_cpp(double J_k, double neighbor_sum) {
  // Compute the conditional probability for x_i = +1
  double prob_1 = 1.0 / (1.0 + exp(-2.0 * J_k * neighbor_sum));
  return prob_1;
}

// [[Rcpp::export]]
double calculate_sum_cpp(NumericVector x, List nb_list) {
  double total_sum = 0.0;
  
  for (int i = 0; i < x.size(); ++i) {
    // Get the neighbors of the current spin as a NumericVector
    NumericVector neighbors = nb_list[i];
    // Compute the contribution of x[i] with its neighbors
    for (int j = 0; j < neighbors.size(); ++j) {
      total_sum += x[i] * x[neighbors[j]];
    }
  }
  // Divide by 2 to account for double counting
  total_sum /= 2.0;
  return total_sum;
}


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List MC_sampling_x_cpp(NumericVector init_x, List nb_list_cpp, double J_k, int n_flip, int M) {
  // Initial variables
  int n = init_x.size();
  double sum_old = calculate_sum_cpp(init_x, nb_list_cpp);
  NumericVector sum_list(n_flip);
  
  // Main loop for flips
  for (int i = 0; i < n_flip; i++) {
    NumericVector new_x = clone(init_x); // Copy init_x
    int idx = floor(R::runif(0, n)); // Random index in [0, n-1]
    new_x[idx] = -1 * new_x[idx]; // Flip the sign of the selected element
    
    double sum_new = calculate_sum_cpp(new_x, nb_list_cpp);
    
    // Metropolis-Hastings acceptance step
    if (exp(J_k * (sum_new - sum_old)) > R::runif(0, 1)) {
      init_x = new_x;  // Update configuration
      sum_old = sum_new;  // Update sum
    }
    
    sum_list[i] = sum_old; // Store the sum
  }
  
  // Compute the partial derivative
  double partial_der = 0;
  for (int i = n_flip - M; i < n_flip; i++) {
    partial_der += sum_list[i];
  }
  partial_der /= M;
  
  // Return results
  return List::create(Named("sum_list") = sum_list,
                      Named("partial_der") = partial_der);
}


