#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::sp_mat matrixAddition(arma::sp_mat A, arma::sp_mat B) {
  // Perform matrix addition
  arma::sp_mat result = A + B;
  return result;
}

