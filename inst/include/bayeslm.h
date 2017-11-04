#ifndef __BAYESLM_H__
#define __BAYESLM_H__


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <iostream>
using namespace Rcpp;
using namespace arma;
using namespace std;
using namespace R;



double log_normal_density_matrix(arma::mat x, arma::mat Sigma,bool singular);
double log_horseshoe_approx_prior(arma::mat beta, double v, arma::uvec penalize);
arma::vec sample_exp(arma::vec lambda);
arma::mat scaling(arma::mat x);
double log_double_exp_prior(arma::mat beta, arma::vec v);
double log_cauchy_prior(arma::mat beta, arma::vec v);
double log_normal_prior(arma::mat beta, arma::vec v);
// double betaprior(arma::mat beta, arma::vec v, int prior, Rcpp::Nullable<Rcpp::Function> user_prior_function);
double user_prior_function_wrapper(arma::mat beta, arma::vec v, Rcpp::Function f);
arma::field<arma::mat> conditional_factors(arma::mat X, arma::vec V);
double log_normal_density(double x, double mu, double sigma);
double log_asymmetric_prior(arma::mat beta, double vglobal, arma::vec prob);
double log_cauchy_density(double x);
double log_nonlocal_prior(arma::mat beta, double vglobal, arma::uvec penalize, arma::vec prob);
double log_ridge_prior(arma::mat beta, double lambda, double vglobal, arma::uvec penalize);
double log_laplace_prior(arma::mat beta, double tau, double sigma, double vglobal, arma::uvec penalize);
double log_asymmetric_prior(arma::mat beta, double vglobal, arma::vec prob, arma::uvec penalize);
//

#endif