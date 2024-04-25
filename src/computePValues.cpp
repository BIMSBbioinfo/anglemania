#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Compute p-values from a correlation matrix
// [[Rcpp::export]]
void computePValues(arma::mat &cor_mat, int df)
{
  // Compute the t-values directly in the cor_mat matrix to save space
  cor_mat = cor_mat % arma::sqrt(df / (1 - arma::square(cor_mat)));

  // Calculate p-values and store them directly in cor_mat
  cor_mat = 2 * (1 - arma::normcdf(arma::abs(cor_mat)));
}