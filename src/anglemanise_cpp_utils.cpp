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

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(bigstatsr, rmio)]]
#include <bigstatsr/BMAcc.h>
#include <Rcpp.h>
using namespace Rcpp;
using std::size_t;

// [[Rcpp::export]]
NumericVector colmean_no_diag_FBM(Environment BM) {
  // Retrieve the FBM pointer from the "address_rw" slot
  XPtr<FBM_RW> xpBM = BM["address_rw"];
  
  // Check that the FBM is of double type (matrix_type == 8)
  if(xpBM->matrix_type() != 8)
    stop("Only FBMs with double type are supported.");
  
  // Create an accessor and get dimensions
  BMAcc_RW<double> acc(xpBM);
  size_t n = acc.nrow();
  size_t m = acc.ncol();
  
  // For a symmetric matrix with a diagonal, n should equal m
  if(n != m)
    stop("Matrix must be square when ignoring the diagonal.");
  
  // Prepare a NumericVector to store the column means
  NumericVector means(m);
  
  // Compute the mean for each column, ignoring the diagonal element
  for (size_t j = 0; j < m; j++) {
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) {
      if (i == j) continue;  // Skip the diagonal element
      sum += acc(i, j);
    }
    // Divide by (n-1) because we omitted one element
    means[j] = sum / (n - 1);
  }
  
  return means;
}


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(bigstatsr, rmio)]]
#include <bigstatsr/BMAcc.h>
#include <Rcpp.h>
#include <algorithm>  // for std::sort
#include <vector>
using namespace Rcpp;
using std::size_t;

// [[Rcpp::export]]
NumericVector colmedian_no_diag_FBM(Environment BM) {
  // Retrieve the FBM pointer from the "address_rw" slot.
  XPtr<FBM_RW> xpBM = BM["address_rw"];
  
  // Check that the FBM is of double type.
  if (xpBM->matrix_type() != 8)
    stop("Only FBMs with double type are supported.");
  
  // Create an accessor and get dimensions.
  BMAcc_RW<double> acc(xpBM);
  size_t n = acc.nrow();
  size_t m = acc.ncol();
  
  // The matrix should be square for a symmetric matrix.
  if (n != m)
    stop("Matrix must be square when ignoring the diagonal.");
  
  // Prepare a NumericVector to store the medians for each column.
  NumericVector medians(m);
  
  // Compute the median for each column, ignoring the diagonal element.
  for (size_t j = 0; j < m; j++) {
    std::vector<double> values;
    values.reserve(n - 1);  // Reserve space for off-diagonal elements.
    
    for (size_t i = 0; i < n; i++) {
      if (i == j) continue; // Skip the diagonal.
      values.push_back(acc(i, j));
    }
    
    // Sort the collected values.
    std::sort(values.begin(), values.end());
    
    // Calculate the median.
    size_t L = values.size();
    double median;
    if (L % 2 == 1) {
      median = values[L / 2];
    } else {
      median = (values[L / 2 - 1] + values[L / 2]) / 2.0;
    }
    
    medians[j] = median;
  }
  
  return medians;
}
