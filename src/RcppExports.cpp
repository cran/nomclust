// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// anderberg_cpp
std::vector<double> anderberg_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const double freq_rel_r, const std::vector<double>& num_cat, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_anderberg_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP freq_relSEXP, SEXP freq_rel_rSEXP, SEXP num_catSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_rel(freq_relSEXP);
    Rcpp::traits::input_parameter< const double >::type freq_rel_r(freq_rel_rSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type num_cat(num_catSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(anderberg_cpp(r, s, dat_vec, freq_rel, freq_rel_r, num_cat, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// burnaby_cpp
std::vector<double> burnaby_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const double freq_rel_r, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_burnaby_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP freq_relSEXP, SEXP freq_rel_rSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_rel(freq_relSEXP);
    Rcpp::traits::input_parameter< const double >::type freq_rel_r(freq_rel_rSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(burnaby_cpp(r, s, dat_vec, freq_rel, freq_rel_r, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// eskin_cpp
std::vector<double> eskin_cpp(const double r, const double s, const std::vector<double>& num_cat, const std::vector<double>& dat_vec, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_eskin_cpp(SEXP rSEXP, SEXP sSEXP, SEXP num_catSEXP, SEXP dat_vecSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type num_cat(num_catSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(eskin_cpp(r, s, num_cat, dat_vec, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// gambaryan_cpp
std::vector<double> gambaryan_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const double freq_rel_r, const double num_cat_sum, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_gambaryan_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP freq_relSEXP, SEXP freq_rel_rSEXP, SEXP num_cat_sumSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_rel(freq_relSEXP);
    Rcpp::traits::input_parameter< const double >::type freq_rel_r(freq_rel_rSEXP);
    Rcpp::traits::input_parameter< const double >::type num_cat_sum(num_cat_sumSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(gambaryan_cpp(r, s, dat_vec, freq_rel, freq_rel_r, num_cat_sum, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// good1_cpp
std::vector<double> good1_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const std::vector<double>& freq_rel2, const double freq_rel_r, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_good1_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP freq_relSEXP, SEXP freq_rel2SEXP, SEXP freq_rel_rSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_rel(freq_relSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_rel2(freq_rel2SEXP);
    Rcpp::traits::input_parameter< const double >::type freq_rel_r(freq_rel_rSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(good1_cpp(r, s, dat_vec, freq_rel, freq_rel2, freq_rel_r, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// good2_cpp
std::vector<double> good2_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const std::vector<double>& freq_rel2, const double freq_rel_r, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_good2_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP freq_relSEXP, SEXP freq_rel2SEXP, SEXP freq_rel_rSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_rel(freq_relSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_rel2(freq_rel2SEXP);
    Rcpp::traits::input_parameter< const double >::type freq_rel_r(freq_rel_rSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(good2_cpp(r, s, dat_vec, freq_rel, freq_rel2, freq_rel_r, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// good3_cpp
std::vector<double> good3_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const double freq_rel_r, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_good3_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP freq_relSEXP, SEXP freq_rel_rSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_rel(freq_relSEXP);
    Rcpp::traits::input_parameter< const double >::type freq_rel_r(freq_rel_rSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(good3_cpp(r, s, dat_vec, freq_rel, freq_rel_r, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// good4_cpp
std::vector<double> good4_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const double freq_rel_r, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_good4_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP freq_relSEXP, SEXP freq_rel_rSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_rel(freq_relSEXP);
    Rcpp::traits::input_parameter< const double >::type freq_rel_r(freq_rel_rSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(good4_cpp(r, s, dat_vec, freq_rel, freq_rel_r, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// iof_cpp
std::vector<double> iof_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_abs, const double freq_abs_r, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_iof_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP freq_absSEXP, SEXP freq_abs_rSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_abs(freq_absSEXP);
    Rcpp::traits::input_parameter< const double >::type freq_abs_r(freq_abs_rSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(iof_cpp(r, s, dat_vec, freq_abs, freq_abs_r, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// lin1_cpp
std::vector<double> lin1_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const double freq_rel_r, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_lin1_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP freq_relSEXP, SEXP freq_rel_rSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_rel(freq_relSEXP);
    Rcpp::traits::input_parameter< const double >::type freq_rel_r(freq_rel_rSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(lin1_cpp(r, s, dat_vec, freq_rel, freq_rel_r, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// lin_cpp
std::vector<double> lin_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_rel, const double freq_rel_r, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_lin_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP freq_relSEXP, SEXP freq_rel_rSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_rel(freq_relSEXP);
    Rcpp::traits::input_parameter< const double >::type freq_rel_r(freq_rel_rSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(lin_cpp(r, s, dat_vec, freq_rel, freq_rel_r, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// of_cpp
std::vector<double> of_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_abs, const double freq_abs_r, const double n, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_of_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP freq_absSEXP, SEXP freq_abs_rSEXP, SEXP nSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_abs(freq_absSEXP);
    Rcpp::traits::input_parameter< const double >::type freq_abs_r(freq_abs_rSEXP);
    Rcpp::traits::input_parameter< const double >::type n(nSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(of_cpp(r, s, dat_vec, freq_abs, freq_abs_r, n, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// sm_cpp
std::vector<double> sm_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_sm_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(sm_cpp(r, s, dat_vec, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// smirnov_cpp
std::vector<double> smirnov_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& freq_abs, const double freq_abs_r, const double num_cat_sum, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_smirnov_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP freq_absSEXP, SEXP freq_abs_rSEXP, SEXP num_cat_sumSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type freq_abs(freq_absSEXP);
    Rcpp::traits::input_parameter< const double >::type freq_abs_r(freq_abs_rSEXP);
    Rcpp::traits::input_parameter< const double >::type num_cat_sum(num_cat_sumSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(smirnov_cpp(r, s, dat_vec, freq_abs, freq_abs_r, num_cat_sum, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// ve_cpp
std::vector<double> ve_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& norm_entropy, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_ve_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP norm_entropySEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type norm_entropy(norm_entropySEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(ve_cpp(r, s, dat_vec, norm_entropy, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}
// vm_cpp
std::vector<double> vm_cpp(const double r, const double s, const std::vector<double>& dat_vec, const std::vector<double>& norm_gini, const std::vector<double>& w, const double sum_w);
RcppExport SEXP _nomclust_vm_cpp(SEXP rSEXP, SEXP sSEXP, SEXP dat_vecSEXP, SEXP norm_giniSEXP, SEXP wSEXP, SEXP sum_wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type dat_vec(dat_vecSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type norm_gini(norm_giniSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const double >::type sum_w(sum_wSEXP);
    rcpp_result_gen = Rcpp::wrap(vm_cpp(r, s, dat_vec, norm_gini, w, sum_w));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nomclust_anderberg_cpp", (DL_FUNC) &_nomclust_anderberg_cpp, 8},
    {"_nomclust_burnaby_cpp", (DL_FUNC) &_nomclust_burnaby_cpp, 7},
    {"_nomclust_eskin_cpp", (DL_FUNC) &_nomclust_eskin_cpp, 6},
    {"_nomclust_gambaryan_cpp", (DL_FUNC) &_nomclust_gambaryan_cpp, 8},
    {"_nomclust_good1_cpp", (DL_FUNC) &_nomclust_good1_cpp, 8},
    {"_nomclust_good2_cpp", (DL_FUNC) &_nomclust_good2_cpp, 8},
    {"_nomclust_good3_cpp", (DL_FUNC) &_nomclust_good3_cpp, 7},
    {"_nomclust_good4_cpp", (DL_FUNC) &_nomclust_good4_cpp, 7},
    {"_nomclust_iof_cpp", (DL_FUNC) &_nomclust_iof_cpp, 7},
    {"_nomclust_lin1_cpp", (DL_FUNC) &_nomclust_lin1_cpp, 7},
    {"_nomclust_lin_cpp", (DL_FUNC) &_nomclust_lin_cpp, 7},
    {"_nomclust_of_cpp", (DL_FUNC) &_nomclust_of_cpp, 8},
    {"_nomclust_sm_cpp", (DL_FUNC) &_nomclust_sm_cpp, 5},
    {"_nomclust_smirnov_cpp", (DL_FUNC) &_nomclust_smirnov_cpp, 8},
    {"_nomclust_ve_cpp", (DL_FUNC) &_nomclust_ve_cpp, 6},
    {"_nomclust_vm_cpp", (DL_FUNC) &_nomclust_vm_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_nomclust(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
