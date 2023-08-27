#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix dist2mat(NumericVector& x, int bf) {

  /* input validation */
  if (!x.inherits("dist")) stop("Input must be a 'dist' object");
  int n = x.attr("Size");
  if (n > 65536) stop("R cannot create a square matrix larger than 65536 x 65536");

  /* initialization */
  NumericMatrix A(n, n);

  /* use pointers */
  size_t j, i, jj, ni, nj; double *ptr_x = &x[0];
  double *A_jj, *A_ij, *A_ji, *col, *row, *end;

  /* fill in lower triangular part */
  for (j = 0; j < n; j++) {
    col = &A(j + 1, j); end = &A(n, j);
    while (col < end) *col++ = *ptr_x++;
    }

  /* cache blocking factor */
  size_t b = (size_t)bf;

  /* copy lower triangular to upper triangular; cache blocking applied */
  for (j = 0; j < n; j += b) {
    nj = n - j; if (nj > b) nj = b;
    /* diagonal block has size nj x nj */
    A_jj = &A(j, j);
    for (jj = nj - 1; jj > 0; jj--, A_jj += n + 1) {
      /* copy a column segment to a row segment */
      col = A_jj + 1; row = A_jj + n;
      for (end = col + jj; col < end; col++, row += n) *row = *col;
      }
    /* off-diagonal blocks */
    for (i = j + nj; i < n; i += b) {
      ni = n - i; if (ni > b) ni = b;
      /* off-diagonal block has size ni x nj */
      A_ij = &A(i, j); A_ji = &A(j, i);
      for (jj = 0; jj < nj; jj++) {
        /* copy a column segment to a row segment */
        col = A_ij + jj * n; row = A_ji + jj;
        for (end = col + ni; col < end; col++, row += n) *row = *col;
        }
      }
    }

  /* add row names and column names */
  A.attr("dimnames") = List::create(x.attr("Labels"), x.attr("Labels"));

  return A;
  }
