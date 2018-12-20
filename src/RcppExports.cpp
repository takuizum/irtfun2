// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// Estep_girt
List Estep_girt(DataFrame x, NumericVector a0, NumericVector b0, NumericVector Xq, NumericVector AX, NumericVector Yr, NumericVector BY, double D);
RcppExport SEXP _irtfun2_Estep_girt(SEXP xSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP XqSEXP, SEXP AXSEXP, SEXP YrSEXP, SEXP BYSEXP, SEXP DSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Xq(XqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type AX(AXSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Yr(YrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type BY(BYSEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    rcpp_result_gen = Rcpp::wrap(Estep_girt(x, a0, b0, Xq, AX, Yr, BY, D));
    return rcpp_result_gen;
END_RCPP
}
// estip
List estip(DataFrame x, CharacterVector model0, const int N, const int bg0, int fc0, int ng, int gc0, const double eMLL, const double eEM, const double eM, const double emu, const double esd, const double D, const double ic, const double max, const double min, const double mu, const double sigma, const int Bayes, const String method, const double mu_a, const double sigma_a, const double mu_b, const double sigma_b, const double mu_c, const double w_c, const int fix, const int print, const double min_a, const double maxabs_b, const int maxiter_em, const int maxiter_j, const int maxskip_j, CharacterVector rm_list, const String thdist, const int EM_dist);
RcppExport SEXP _irtfun2_estip(SEXP xSEXP, SEXP model0SEXP, SEXP NSEXP, SEXP bg0SEXP, SEXP fc0SEXP, SEXP ngSEXP, SEXP gc0SEXP, SEXP eMLLSEXP, SEXP eEMSEXP, SEXP eMSEXP, SEXP emuSEXP, SEXP esdSEXP, SEXP DSEXP, SEXP icSEXP, SEXP maxSEXP, SEXP minSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP BayesSEXP, SEXP methodSEXP, SEXP mu_aSEXP, SEXP sigma_aSEXP, SEXP mu_bSEXP, SEXP sigma_bSEXP, SEXP mu_cSEXP, SEXP w_cSEXP, SEXP fixSEXP, SEXP printSEXP, SEXP min_aSEXP, SEXP maxabs_bSEXP, SEXP maxiter_emSEXP, SEXP maxiter_jSEXP, SEXP maxskip_jSEXP, SEXP rm_listSEXP, SEXP thdistSEXP, SEXP EM_distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type x(xSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type model0(model0SEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type bg0(bg0SEXP);
    Rcpp::traits::input_parameter< int >::type fc0(fc0SEXP);
    Rcpp::traits::input_parameter< int >::type ng(ngSEXP);
    Rcpp::traits::input_parameter< int >::type gc0(gc0SEXP);
    Rcpp::traits::input_parameter< const double >::type eMLL(eMLLSEXP);
    Rcpp::traits::input_parameter< const double >::type eEM(eEMSEXP);
    Rcpp::traits::input_parameter< const double >::type eM(eMSEXP);
    Rcpp::traits::input_parameter< const double >::type emu(emuSEXP);
    Rcpp::traits::input_parameter< const double >::type esd(esdSEXP);
    Rcpp::traits::input_parameter< const double >::type D(DSEXP);
    Rcpp::traits::input_parameter< const double >::type ic(icSEXP);
    Rcpp::traits::input_parameter< const double >::type max(maxSEXP);
    Rcpp::traits::input_parameter< const double >::type min(minSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const int >::type Bayes(BayesSEXP);
    Rcpp::traits::input_parameter< const String >::type method(methodSEXP);
    Rcpp::traits::input_parameter< const double >::type mu_a(mu_aSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma_a(sigma_aSEXP);
    Rcpp::traits::input_parameter< const double >::type mu_b(mu_bSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma_b(sigma_bSEXP);
    Rcpp::traits::input_parameter< const double >::type mu_c(mu_cSEXP);
    Rcpp::traits::input_parameter< const double >::type w_c(w_cSEXP);
    Rcpp::traits::input_parameter< const int >::type fix(fixSEXP);
    Rcpp::traits::input_parameter< const int >::type print(printSEXP);
    Rcpp::traits::input_parameter< const double >::type min_a(min_aSEXP);
    Rcpp::traits::input_parameter< const double >::type maxabs_b(maxabs_bSEXP);
    Rcpp::traits::input_parameter< const int >::type maxiter_em(maxiter_emSEXP);
    Rcpp::traits::input_parameter< const int >::type maxiter_j(maxiter_jSEXP);
    Rcpp::traits::input_parameter< const int >::type maxskip_j(maxskip_jSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type rm_list(rm_listSEXP);
    Rcpp::traits::input_parameter< const String >::type thdist(thdistSEXP);
    Rcpp::traits::input_parameter< const int >::type EM_dist(EM_distSEXP);
    rcpp_result_gen = Rcpp::wrap(estip(x, model0, N, bg0, fc0, ng, gc0, eMLL, eEM, eM, emu, esd, D, ic, max, min, mu, sigma, Bayes, method, mu_a, sigma_a, mu_b, sigma_b, mu_c, w_c, fix, print, min_a, maxabs_b, maxiter_em, maxiter_j, maxskip_j, rm_list, thdist, EM_dist));
    return rcpp_result_gen;
END_RCPP
}
// theta_pv
NumericMatrix theta_pv(DataFrame x, const int nofrands, NumericVector eap_apply, NumericVector const_apply, NumericVector map_apply, const int n, double maxtheta, double mintheta, NumericVector a, NumericVector b, NumericVector c, const double D, const double mu, const double sigma);
RcppExport SEXP _irtfun2_theta_pv(SEXP xSEXP, SEXP nofrandsSEXP, SEXP eap_applySEXP, SEXP const_applySEXP, SEXP map_applySEXP, SEXP nSEXP, SEXP maxthetaSEXP, SEXP minthetaSEXP, SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP DSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type nofrands(nofrandsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eap_apply(eap_applySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type const_apply(const_applySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type map_apply(map_applySEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type maxtheta(maxthetaSEXP);
    Rcpp::traits::input_parameter< double >::type mintheta(minthetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double >::type D(DSEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(theta_pv(x, nofrands, eap_apply, const_apply, map_apply, n, maxtheta, mintheta, a, b, c, D, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_irtfun2_Estep_girt", (DL_FUNC) &_irtfun2_Estep_girt, 8},
    {"_irtfun2_estip", (DL_FUNC) &_irtfun2_estip, 36},
    {"_irtfun2_theta_pv", (DL_FUNC) &_irtfun2_theta_pv, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_irtfun2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
