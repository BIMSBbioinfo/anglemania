// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(bigstatsr, rmio)]]
#include <bigstatsr/BMAcc.h>
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using std::size_t;

// [[Rcpp::export]]
DataFrame select_genes_cpp(Environment BM_sn,
                           Environment BM_mean,
                           double zscore_sn_threshold,
                           double zscore_mean_threshold) {
  // Retrieve FBM pointers (using the "address_rw" slot)
  XPtr<FBM_RW> xpBM_sn = BM_sn["address_rw"];
  XPtr<FBM_RW> xpBM_mean = BM_mean["address_rw"];
  
  if (xpBM_sn->matrix_type() != 8 || xpBM_mean->matrix_type() != 8)
    stop("Only FBMs with double type are supported.");
  
  BMAcc_RW<double> sn_acc(xpBM_sn);
  BMAcc_RW<double> mean_acc(xpBM_mean);
  
  size_t n = sn_acc.nrow();
  if (sn_acc.ncol() != n)
    stop("The sn_zscore FBM must be square.");
  if (mean_acc.nrow() != n || mean_acc.ncol() != n)
    stop("The mean_zscore FBM must be square and match dimensions with sn_zscore.");
  
  // First pass: Count the number of pairs that pass the thresholds
  size_t count = 0;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      double sn_val = sn_acc(i, j);
      if (sn_val >= zscore_sn_threshold) {
        double mean_val = mean_acc(i, j);
        if (std::abs(mean_val) >= zscore_mean_threshold) {
          count++;
        }
      }
    }
  }
  
  // Preallocate vectors with the exact count
  IntegerVector rows(count);
  IntegerVector cols(count);
  NumericVector sn_vals(count);
  NumericVector mean_vals(count);
  
  // Second pass: Fill in the vectors
  size_t k = 0;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      double sn_val = sn_acc(i, j);
      if (sn_val >= zscore_sn_threshold) {
        double mean_val = mean_acc(i, j);
        if (std::abs(mean_val) >= zscore_mean_threshold) {
          // R uses 1-indexing
          rows[k] = i + 1;
          cols[k] = j + 1;
          sn_vals[k] = sn_val;
          mean_vals[k] = mean_val;
          k++;
        }
      }
    }
  }
  
  // Create the DataFrame
  DataFrame df = DataFrame::create(
    _["geneA"]  = rows,
    _["geneB"]  = cols,
    _["sn_zscore"]   = sn_vals,
    _["mean_zscore"] = mean_vals
  );
  
  return df;
}


// For a square FBM
// Adapted from bigstatsr::scaleK for export in this package
// Source: https://github.com/privefl/bigstatsr
// License: GPL-3
// [[Rcpp::export]]
void scaleK(Environment BM,
            const NumericVector& sums,
            const NumericVector& mu,
            const NumericVector& delta,
            int nrow) {

  XPtr<FBM_RW> xpBM = BM["address_rw"];
  BMAcc_RW<double> K(xpBM);

  size_t n = K.nrow();
  myassert_size(K.ncol(), n);

  for (size_t j = 0; j < n; j++) {
    for (size_t i = 0; i < n; i++) {
      K(i, j) -= sums[i] * mu[j] + mu[i] * sums[j];
      K(i, j) += nrow * mu[i] * mu[j];
      K(i, j) /= delta(i) * delta(j);
    }
  }
}