// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// select_genes_cpp
DataFrame select_genes_cpp(Environment BM_sn, Environment BM_mean, double zscore_sn_threshold, double zscore_mean_threshold);
RcppExport SEXP _anglemania_select_genes_cpp(SEXP BM_snSEXP, SEXP BM_meanSEXP, SEXP zscore_sn_thresholdSEXP, SEXP zscore_mean_thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM_sn(BM_snSEXP);
    Rcpp::traits::input_parameter< Environment >::type BM_mean(BM_meanSEXP);
    Rcpp::traits::input_parameter< double >::type zscore_sn_threshold(zscore_sn_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type zscore_mean_threshold(zscore_mean_thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(select_genes_cpp(BM_sn, BM_mean, zscore_sn_threshold, zscore_mean_threshold));
    return rcpp_result_gen;
END_RCPP
}
// colmean_no_diag_FBM
NumericVector colmean_no_diag_FBM(Environment BM);
RcppExport SEXP _anglemania_colmean_no_diag_FBM(SEXP BMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    rcpp_result_gen = Rcpp::wrap(colmean_no_diag_FBM(BM));
    return rcpp_result_gen;
END_RCPP
}
// colmedian_no_diag_FBM
NumericVector colmedian_no_diag_FBM(Environment BM);
RcppExport SEXP _anglemania_colmedian_no_diag_FBM(SEXP BMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Environment >::type BM(BMSEXP);
    rcpp_result_gen = Rcpp::wrap(colmedian_no_diag_FBM(BM));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_anglemania_select_genes_cpp", (DL_FUNC) &_anglemania_select_genes_cpp, 4},
    {"_anglemania_colmean_no_diag_FBM", (DL_FUNC) &_anglemania_colmean_no_diag_FBM, 1},
    {"_anglemania_colmedian_no_diag_FBM", (DL_FUNC) &_anglemania_colmedian_no_diag_FBM, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_anglemania(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
