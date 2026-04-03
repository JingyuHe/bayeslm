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



double log_normal_density_matrix(const arma::mat &x, const arma::mat &Sigma, bool singular);
double log_normal_density_matrix_2(const arma::mat &x, const arma::vec &Sigma_inv_diag, bool singular);
double log_normal_density_matrix_2(const arma::mat &x, const arma::vec &Sigma_inv_diag, double log_det_Sigma_inv, bool singular);
double log_horseshoe_approx_prior(const arma::mat &beta, double v, double sigma, const arma::uvec &penalize, bool scale_sigma_prior);
double log_horseshoe_approx_prior(const arma::mat &beta, double v, double sigma, const arma::uvec &penalize_index, bool all_penalized, bool scale_sigma_prior);
arma::vec sample_exp(const arma::vec &lambda);
arma::mat scaling(const arma::mat &x);
double log_double_exp_prior(const arma::mat &beta, const arma::vec &v);
double log_cauchy_prior(const arma::mat &beta, const arma::vec &v);
double log_normal_prior(const arma::mat &beta, const arma::vec &v);
// double betaprior(arma::mat beta, arma::vec v, int prior, Rcpp::Nullable<Rcpp::Function> user_prior_function);
double user_prior_function_wrapper(const arma::mat &beta, const arma::vec &v, Rcpp::Function f);
arma::field<arma::mat> conditional_factors(const arma::mat &X, const arma::vec &V);
arma::field<arma::mat> conditional_factors_parallel(const arma::mat &X, const arma::vec &V);
double scalar_conditional_mean(const arma::mat &mean_factor, const arma::vec &beta_dev, arma::uword block_start);
double log_normal_density(double x, double mu, double sigma);
double log_cauchy_density(double x);
double log_nonlocal_prior(const arma::mat &beta, double vglobal, double sigma, const arma::uvec &penalize, const arma::vec &prob, bool scale_sigma_prior);
double log_nonlocal_prior(const arma::mat &beta, double vglobal, double sigma, const arma::uvec &penalize_index, const arma::vec &prior_mean, bool all_penalized, bool scale_sigma_prior);
double log_ridge_prior(const arma::mat &beta, double lambda, double vglobal, double sigma, const arma::uvec &penalize, bool scale_sigma_prior);
double log_ridge_prior(const arma::mat &beta, double lambda, double vglobal, double sigma, const arma::uvec &penalize_index, bool all_penalized, bool scale_sigma_prior);
double log_laplace_prior(const arma::mat &beta, double tau, double sigma, double vglobal, const arma::uvec &penalize);
double log_laplace_prior(const arma::mat &beta, double tau, double sigma, double vglobal, const arma::uvec &penalize_index, bool all_penalized);
double log_inverselaplace_prior(const arma::mat &beta, double lambda, double sigma, double vglobal, const arma::uvec &penalize);
double log_inverselaplace_prior(const arma::mat &beta, double lambda, double sigma, double vglobal, const arma::uvec &penalize_index, bool all_penalized);
double log_asymmetric_prior(const arma::mat &beta, double vglobal, double sigma, const arma::vec &prob, const arma::uvec &penalize, bool scale_sigma_prior);
double log_asymmetric_prior(const arma::mat &beta, double vglobal, double sigma, const arma::vec &prob, const arma::uvec &penalize_index, bool all_penalized, bool scale_sigma_prior);
arma::mat sampling_beta(const arma::mat &mu_n, const arma::mat &chol_Lambda_n_inv, double sigma, int p, bool scale_sigma_prior);
double sampling_sigma(double a_n, double b_0, const arma::mat &YY, const arma::mat &mu_n, const arma::mat &Lambda_n);
arma::vec sampling_lambda(const arma::mat &lambda, const arma::mat &beta, double sigma, double tau, int p, bool scale_sigma_prior);
double sampling_tau(const arma::mat &lambda, const arma::mat &beta, double sigma, double tau, bool scale_sigma_prior);

arma::mat sampling_beta_2(const arma::mat &old_beta, const arma::mat &mu_n, double sigma, int p, bool scale_sigma_prior, const arma::mat &lambda, double tau, const arma::mat &X, const arma::mat &Y);

//

#ifdef FALSE
   #undef FALSE
#endif

#endif
