// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/bayeslm.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sharkfin_cpp_loop
List sharkfin_cpp_loop(arma::mat Y, arma::mat X, arma::vec prob_vec, arma::uvec penalize, arma::vec block_vec, arma::vec cc, int prior_type, double sigma, double s2, double kap2, int nsamps, int burn, int skip, double vglobal, bool sampling_vglobal, bool verb, bool icept, bool standardize, bool singular, bool scale_sigma_prior);
RcppExport SEXP _bayeslm_sharkfin_cpp_loop(SEXP YSEXP, SEXP XSEXP, SEXP prob_vecSEXP, SEXP penalizeSEXP, SEXP block_vecSEXP, SEXP ccSEXP, SEXP prior_typeSEXP, SEXP sigmaSEXP, SEXP s2SEXP, SEXP kap2SEXP, SEXP nsampsSEXP, SEXP burnSEXP, SEXP skipSEXP, SEXP vglobalSEXP, SEXP sampling_vglobalSEXP, SEXP verbSEXP, SEXP iceptSEXP, SEXP standardizeSEXP, SEXP singularSEXP, SEXP scale_sigma_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prob_vec(prob_vecSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type penalize(penalizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type block_vec(block_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cc(ccSEXP);
    Rcpp::traits::input_parameter< int >::type prior_type(prior_typeSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type kap2(kap2SEXP);
    Rcpp::traits::input_parameter< int >::type nsamps(nsampsSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< double >::type vglobal(vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type sampling_vglobal(sampling_vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< bool >::type icept(iceptSEXP);
    Rcpp::traits::input_parameter< bool >::type standardize(standardizeSEXP);
    Rcpp::traits::input_parameter< bool >::type singular(singularSEXP);
    Rcpp::traits::input_parameter< bool >::type scale_sigma_prior(scale_sigma_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(sharkfin_cpp_loop(Y, X, prob_vec, penalize, block_vec, cc, prior_type, sigma, s2, kap2, nsamps, burn, skip, vglobal, sampling_vglobal, verb, icept, standardize, singular, scale_sigma_prior));
    return rcpp_result_gen;
END_RCPP
}
// bayeslm_cpp_loop
List bayeslm_cpp_loop(arma::mat Y, arma::mat X, arma::uvec penalize, arma::vec block_vec, arma::vec cc, int prior_type, Rcpp::Nullable<Rcpp::Function> user_prior_function, double sigma, double s2, double kap2, int nsamps, int burn, int skip, double vglobal, bool sampling_vglobal, bool verb, bool icept, bool standardize, bool singular, bool scale_sigma_prior);
RcppExport SEXP _bayeslm_bayeslm_cpp_loop(SEXP YSEXP, SEXP XSEXP, SEXP penalizeSEXP, SEXP block_vecSEXP, SEXP ccSEXP, SEXP prior_typeSEXP, SEXP user_prior_functionSEXP, SEXP sigmaSEXP, SEXP s2SEXP, SEXP kap2SEXP, SEXP nsampsSEXP, SEXP burnSEXP, SEXP skipSEXP, SEXP vglobalSEXP, SEXP sampling_vglobalSEXP, SEXP verbSEXP, SEXP iceptSEXP, SEXP standardizeSEXP, SEXP singularSEXP, SEXP scale_sigma_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type penalize(penalizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type block_vec(block_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cc(ccSEXP);
    Rcpp::traits::input_parameter< int >::type prior_type(prior_typeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::Function> >::type user_prior_function(user_prior_functionSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type kap2(kap2SEXP);
    Rcpp::traits::input_parameter< int >::type nsamps(nsampsSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< double >::type vglobal(vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type sampling_vglobal(sampling_vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< bool >::type icept(iceptSEXP);
    Rcpp::traits::input_parameter< bool >::type standardize(standardizeSEXP);
    Rcpp::traits::input_parameter< bool >::type singular(singularSEXP);
    Rcpp::traits::input_parameter< bool >::type scale_sigma_prior(scale_sigma_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(bayeslm_cpp_loop(Y, X, penalize, block_vec, cc, prior_type, user_prior_function, sigma, s2, kap2, nsamps, burn, skip, vglobal, sampling_vglobal, verb, icept, standardize, singular, scale_sigma_prior));
    return rcpp_result_gen;
END_RCPP
}
// blasso_cpp_loop
List blasso_cpp_loop(arma::mat Y, arma::mat X, arma::uvec penalize, arma::vec block_vec, arma::vec cc, int prior_type, double sigma, double s2, double kap2, int nsamps, int burn, int skip, double vglobal, bool sampling_vglobal, bool verb, bool icept, bool standardize, bool singular, bool scale_sigma_prior);
RcppExport SEXP _bayeslm_blasso_cpp_loop(SEXP YSEXP, SEXP XSEXP, SEXP penalizeSEXP, SEXP block_vecSEXP, SEXP ccSEXP, SEXP prior_typeSEXP, SEXP sigmaSEXP, SEXP s2SEXP, SEXP kap2SEXP, SEXP nsampsSEXP, SEXP burnSEXP, SEXP skipSEXP, SEXP vglobalSEXP, SEXP sampling_vglobalSEXP, SEXP verbSEXP, SEXP iceptSEXP, SEXP standardizeSEXP, SEXP singularSEXP, SEXP scale_sigma_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type penalize(penalizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type block_vec(block_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cc(ccSEXP);
    Rcpp::traits::input_parameter< int >::type prior_type(prior_typeSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type kap2(kap2SEXP);
    Rcpp::traits::input_parameter< int >::type nsamps(nsampsSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< double >::type vglobal(vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type sampling_vglobal(sampling_vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< bool >::type icept(iceptSEXP);
    Rcpp::traits::input_parameter< bool >::type standardize(standardizeSEXP);
    Rcpp::traits::input_parameter< bool >::type singular(singularSEXP);
    Rcpp::traits::input_parameter< bool >::type scale_sigma_prior(scale_sigma_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(blasso_cpp_loop(Y, X, penalize, block_vec, cc, prior_type, sigma, s2, kap2, nsamps, burn, skip, vglobal, sampling_vglobal, verb, icept, standardize, singular, scale_sigma_prior));
    return rcpp_result_gen;
END_RCPP
}
// hs_gibbs
List hs_gibbs(arma::mat Y, arma::mat X, int nsamps, double a, double b, bool scale_sigma_prior);
RcppExport SEXP _bayeslm_hs_gibbs(SEXP YSEXP, SEXP XSEXP, SEXP nsampsSEXP, SEXP aSEXP, SEXP bSEXP, SEXP scale_sigma_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nsamps(nsampsSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type scale_sigma_prior(scale_sigma_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(hs_gibbs(Y, X, nsamps, a, b, scale_sigma_prior));
    return rcpp_result_gen;
END_RCPP
}
// hs_gibbs_2
List hs_gibbs_2(arma::mat Y, arma::mat X, int nsamps, double a, double b, bool scale_sigma_prior);
RcppExport SEXP _bayeslm_hs_gibbs_2(SEXP YSEXP, SEXP XSEXP, SEXP nsampsSEXP, SEXP aSEXP, SEXP bSEXP, SEXP scale_sigma_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nsamps(nsampsSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< bool >::type scale_sigma_prior(scale_sigma_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(hs_gibbs_2(Y, X, nsamps, a, b, scale_sigma_prior));
    return rcpp_result_gen;
END_RCPP
}
// horseshoe_cpp_loop
List horseshoe_cpp_loop(arma::mat Y, arma::mat X, arma::uvec penalize, arma::vec block_vec, arma::vec cc, int prior_type, double sigma, double s2, double kap2, int nsamps, int burn, int skip, double vglobal, bool sampling_vglobal, bool verb, bool icept, bool standardize, bool singular, bool scale_sigma_prior);
RcppExport SEXP _bayeslm_horseshoe_cpp_loop(SEXP YSEXP, SEXP XSEXP, SEXP penalizeSEXP, SEXP block_vecSEXP, SEXP ccSEXP, SEXP prior_typeSEXP, SEXP sigmaSEXP, SEXP s2SEXP, SEXP kap2SEXP, SEXP nsampsSEXP, SEXP burnSEXP, SEXP skipSEXP, SEXP vglobalSEXP, SEXP sampling_vglobalSEXP, SEXP verbSEXP, SEXP iceptSEXP, SEXP standardizeSEXP, SEXP singularSEXP, SEXP scale_sigma_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type penalize(penalizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type block_vec(block_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cc(ccSEXP);
    Rcpp::traits::input_parameter< int >::type prior_type(prior_typeSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type kap2(kap2SEXP);
    Rcpp::traits::input_parameter< int >::type nsamps(nsampsSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< double >::type vglobal(vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type sampling_vglobal(sampling_vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< bool >::type icept(iceptSEXP);
    Rcpp::traits::input_parameter< bool >::type standardize(standardizeSEXP);
    Rcpp::traits::input_parameter< bool >::type singular(singularSEXP);
    Rcpp::traits::input_parameter< bool >::type scale_sigma_prior(scale_sigma_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(horseshoe_cpp_loop(Y, X, penalize, block_vec, cc, prior_type, sigma, s2, kap2, nsamps, burn, skip, vglobal, sampling_vglobal, verb, icept, standardize, singular, scale_sigma_prior));
    return rcpp_result_gen;
END_RCPP
}
// inverseLaplace_cpp_loop
List inverseLaplace_cpp_loop(arma::mat Y, arma::mat X, double lambda, arma::uvec penalize, arma::vec block_vec, arma::vec cc, int prior_type, double sigma, double s2, double kap2, int nsamps, int burn, int skip, double vglobal, bool sampling_vglobal, bool verb, bool icept, bool standardize, bool singular, bool scale_sigma_prior);
RcppExport SEXP _bayeslm_inverseLaplace_cpp_loop(SEXP YSEXP, SEXP XSEXP, SEXP lambdaSEXP, SEXP penalizeSEXP, SEXP block_vecSEXP, SEXP ccSEXP, SEXP prior_typeSEXP, SEXP sigmaSEXP, SEXP s2SEXP, SEXP kap2SEXP, SEXP nsampsSEXP, SEXP burnSEXP, SEXP skipSEXP, SEXP vglobalSEXP, SEXP sampling_vglobalSEXP, SEXP verbSEXP, SEXP iceptSEXP, SEXP standardizeSEXP, SEXP singularSEXP, SEXP scale_sigma_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type penalize(penalizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type block_vec(block_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cc(ccSEXP);
    Rcpp::traits::input_parameter< int >::type prior_type(prior_typeSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type kap2(kap2SEXP);
    Rcpp::traits::input_parameter< int >::type nsamps(nsampsSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< double >::type vglobal(vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type sampling_vglobal(sampling_vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< bool >::type icept(iceptSEXP);
    Rcpp::traits::input_parameter< bool >::type standardize(standardizeSEXP);
    Rcpp::traits::input_parameter< bool >::type singular(singularSEXP);
    Rcpp::traits::input_parameter< bool >::type scale_sigma_prior(scale_sigma_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(inverseLaplace_cpp_loop(Y, X, lambda, penalize, block_vec, cc, prior_type, sigma, s2, kap2, nsamps, burn, skip, vglobal, sampling_vglobal, verb, icept, standardize, singular, scale_sigma_prior));
    return rcpp_result_gen;
END_RCPP
}
// nonlocal_cpp_loop
List nonlocal_cpp_loop(arma::mat Y, arma::mat X, arma::vec prior_mean, arma::uvec penalize, arma::vec block_vec, arma::vec cc, int prior_type, double sigma, double s2, double kap2, int nsamps, int burn, int skip, double vglobal, bool sampling_vglobal, bool verb, bool icept, bool standardize, bool singular, bool scale_sigma_prior);
RcppExport SEXP _bayeslm_nonlocal_cpp_loop(SEXP YSEXP, SEXP XSEXP, SEXP prior_meanSEXP, SEXP penalizeSEXP, SEXP block_vecSEXP, SEXP ccSEXP, SEXP prior_typeSEXP, SEXP sigmaSEXP, SEXP s2SEXP, SEXP kap2SEXP, SEXP nsampsSEXP, SEXP burnSEXP, SEXP skipSEXP, SEXP vglobalSEXP, SEXP sampling_vglobalSEXP, SEXP verbSEXP, SEXP iceptSEXP, SEXP standardizeSEXP, SEXP singularSEXP, SEXP scale_sigma_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type prior_mean(prior_meanSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type penalize(penalizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type block_vec(block_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cc(ccSEXP);
    Rcpp::traits::input_parameter< int >::type prior_type(prior_typeSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type kap2(kap2SEXP);
    Rcpp::traits::input_parameter< int >::type nsamps(nsampsSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< double >::type vglobal(vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type sampling_vglobal(sampling_vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< bool >::type icept(iceptSEXP);
    Rcpp::traits::input_parameter< bool >::type standardize(standardizeSEXP);
    Rcpp::traits::input_parameter< bool >::type singular(singularSEXP);
    Rcpp::traits::input_parameter< bool >::type scale_sigma_prior(scale_sigma_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(nonlocal_cpp_loop(Y, X, prior_mean, penalize, block_vec, cc, prior_type, sigma, s2, kap2, nsamps, burn, skip, vglobal, sampling_vglobal, verb, icept, standardize, singular, scale_sigma_prior));
    return rcpp_result_gen;
END_RCPP
}
// bridge_cpp_loop
List bridge_cpp_loop(arma::mat Y, arma::mat X, arma::uvec penalize, arma::vec block_vec, arma::vec cc, int prior_type, double sigma, double s2, double kap2, int nsamps, int burn, int skip, double vglobal, bool sampling_vglobal, bool verb, bool icept, bool standardize, bool singular, bool scale_sigma_prior);
RcppExport SEXP _bayeslm_bridge_cpp_loop(SEXP YSEXP, SEXP XSEXP, SEXP penalizeSEXP, SEXP block_vecSEXP, SEXP ccSEXP, SEXP prior_typeSEXP, SEXP sigmaSEXP, SEXP s2SEXP, SEXP kap2SEXP, SEXP nsampsSEXP, SEXP burnSEXP, SEXP skipSEXP, SEXP vglobalSEXP, SEXP sampling_vglobalSEXP, SEXP verbSEXP, SEXP iceptSEXP, SEXP standardizeSEXP, SEXP singularSEXP, SEXP scale_sigma_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type penalize(penalizeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type block_vec(block_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type cc(ccSEXP);
    Rcpp::traits::input_parameter< int >::type prior_type(prior_typeSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< double >::type kap2(kap2SEXP);
    Rcpp::traits::input_parameter< int >::type nsamps(nsampsSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type skip(skipSEXP);
    Rcpp::traits::input_parameter< double >::type vglobal(vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type sampling_vglobal(sampling_vglobalSEXP);
    Rcpp::traits::input_parameter< bool >::type verb(verbSEXP);
    Rcpp::traits::input_parameter< bool >::type icept(iceptSEXP);
    Rcpp::traits::input_parameter< bool >::type standardize(standardizeSEXP);
    Rcpp::traits::input_parameter< bool >::type singular(singularSEXP);
    Rcpp::traits::input_parameter< bool >::type scale_sigma_prior(scale_sigma_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(bridge_cpp_loop(Y, X, penalize, block_vec, cc, prior_type, sigma, s2, kap2, nsamps, burn, skip, vglobal, sampling_vglobal, verb, icept, standardize, singular, scale_sigma_prior));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bayeslm_sharkfin_cpp_loop", (DL_FUNC) &_bayeslm_sharkfin_cpp_loop, 20},
    {"_bayeslm_bayeslm_cpp_loop", (DL_FUNC) &_bayeslm_bayeslm_cpp_loop, 20},
    {"_bayeslm_blasso_cpp_loop", (DL_FUNC) &_bayeslm_blasso_cpp_loop, 19},
    {"_bayeslm_hs_gibbs", (DL_FUNC) &_bayeslm_hs_gibbs, 6},
    {"_bayeslm_hs_gibbs_2", (DL_FUNC) &_bayeslm_hs_gibbs_2, 6},
    {"_bayeslm_horseshoe_cpp_loop", (DL_FUNC) &_bayeslm_horseshoe_cpp_loop, 19},
    {"_bayeslm_inverseLaplace_cpp_loop", (DL_FUNC) &_bayeslm_inverseLaplace_cpp_loop, 20},
    {"_bayeslm_nonlocal_cpp_loop", (DL_FUNC) &_bayeslm_nonlocal_cpp_loop, 20},
    {"_bayeslm_bridge_cpp_loop", (DL_FUNC) &_bayeslm_bridge_cpp_loop, 19},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayeslm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
