#ifndef __BAYESLM_H__
#define __BAYESLM_H__


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <iostream>
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace R;



double log_normal_density_matrix(arma::mat x, arma::mat Sigma,bool singular);
double log_horseshoe_approx_prior(arma::mat beta, double v, double sigma, arma::uvec penalize, bool scale_sigma_prior);
arma::vec sample_exp(arma::vec lambda);
arma::mat scaling(arma::mat x);
double log_double_exp_prior(arma::mat beta, arma::vec v);
double log_cauchy_prior(arma::mat beta, arma::vec v);
double log_normal_prior(arma::mat beta, arma::vec v);
// double betaprior(arma::mat beta, arma::vec v, int prior, Rcpp::Nullable<Rcpp::Function> user_prior_function);
double user_prior_function_wrapper(arma::mat beta, arma::vec v, Rcpp::Function f);
arma::field<arma::mat> conditional_factors(arma::mat X, arma::vec V);
arma::field<arma::mat> conditional_factors_parallel(arma::mat X, arma::vec V);
double log_normal_density(double x, double mu, double sigma);
double log_cauchy_density(double x);
double log_nonlocal_prior(arma::mat beta, double vglobal, double sigma, arma::uvec penalize, arma::vec prob, bool scale_sigma_prior);
double log_ridge_prior(arma::mat beta, double lambda, double vglobal, double sigma, arma::uvec penalize, bool scale_sigma_prior);
double log_laplace_prior(arma::mat beta, double tau, double sigma, double vglobal, arma::uvec penalize);
double log_inverselaplace_prior(arma::mat beta, double lambda, double sigma, double vglobal, arma::uvec penalize);
double log_asymmetric_prior(arma::mat beta, double vglobal, double sigma, arma::vec prob, arma::uvec penalize, bool scale_sigma_prior);
arma::mat sampling_beta(arma::mat mu_n, arma::mat chol_Lambda_n_inv, double sigma, int p, bool scale_sigma_prior);
double sampling_sigma(double a_n, double b_0, arma::mat YY, arma::mat mu_n, arma::mat Lambda_n);
arma::vec sampling_lambda(arma::mat lambda, arma::mat beta, double sigma, double tau, int p, bool scale_sigma_prior);
double sampling_tau(arma::mat lambda, arma::mat beta, double sigma, double tau, bool scale_sigma_prior);

arma::mat sampling_beta_2(arma::mat old_beta, arma::mat mu_n, double sigma, int p, bool scale_sigma_prior, arma::mat lambda, double tau, arma::mat X, arma::mat Y);

//

#ifdef FALSE
   #undef FALSE
#endif

#endif
