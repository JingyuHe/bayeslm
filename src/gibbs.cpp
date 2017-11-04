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
inline arma::mat sampling_beta(arma::mat mu_n, arma::mat chol_Lambda_n_inv, double sigma, int p);
inline double sampling_sigma(double a_n, double b_0, arma::mat YY, arma::mat mu_n, arma::mat Lambda_n);
inline arma::vec sampling_lambda(arma::mat lambda, arma::mat beta, double sigma, double tau, int p);
inline double sampling_tau(arma::mat lambda, arma::mat beta, double sigma, double tau);
/*


Standard Gibbs Sampler for horseshoe regression


*/


// [[Rcpp::export]]
List hs_gibbs(arma::mat Y, arma::mat X, arma::mat beta_ini, int nsamps = 1000, double a = 1, double b = 1){
    clock_t t = clock();
    // data dimensions;
    int p = X.n_cols;
    int n = X.n_rows;
    
    // prior parameters
    double a_0 = a;
    double b_0 = b;

    // parameter storage
    arma::mat beta = zeros<mat>(nsamps, p);
    arma::vec sigma = zeros<vec>(nsamps);
    arma::mat lambda = zeros<mat>(nsamps, p);
    arma::mat gamma_l = zeros<mat>(nsamps, p);
    arma::vec tau = zeros<vec>(nsamps);
    arma::vec gamma_tt = zeros<vec>(nsamps);
    arma::mat Lambda0 = zeros<mat>(p,p);

    // initialize
    // beta.row(0) = zeros<mat>(1, p);
    beta.row(0) = beta_ini;
    lambda.row(0) = ones<mat>(1, p);
    Lambda0.diag() = lambda.row(0);
    sigma(0) = 1;
    tau(0) = 1;

    // compute sufficient statistics
    arma::mat XX = X.t() * X;
    arma::mat XY = X.t() * Y;
    arma::mat YY = Y.t() * Y;

    // initialize beta at OLS estimator
    // beta.row(0) = arma::trans(inv(XX) * XY);

    double a_n = a_0 + n / 2.0;

    arma::mat mu_n;
    arma::mat Lambda_n;
    arma::mat Lambda_n_inv;
    arma::mat chol_Lambda_n;
    arma::mat chol_Lambda_n_inv;
    // sampling

    t = clock() - t;
    cout << "fixed time " << t / CLOCKS_PER_SEC << endl;

    t = clock();
    for(int i = 1; i < nsamps; i ++){

        Lambda0.diag() = 1.0 / (pow(lambda.row(i-1), 2) * pow(tau(i-1), 2));
        
        Lambda_n = XX + Lambda0;

        Lambda_n_inv = inv(Lambda_n);

        chol_Lambda_n_inv = chol(Lambda_n_inv, "lower"); // LL' = inv(Lambda_n)

        mu_n = Lambda_n_inv * XY;

        beta.row(i) = sampling_beta(mu_n, chol_Lambda_n_inv, sigma(i-1), p).t();

        sigma(i) = sampling_sigma(a_n, b_0, YY, mu_n, Lambda_n);

        lambda.row(i) = sampling_lambda(lambda.row(i-1), beta.row(i), sigma(i), tau(i-1), p).t();

        tau(i) = sampling_tau(lambda.row(i), beta.row(i), sigma(i), tau(i-1));
        
    }

    t = clock() - t;
    cout << "float time" << t / CLOCKS_PER_SEC << endl;

    return List::create(Named("beta") = beta, Named("lambda") = lambda, Named("sigma") = sigma, Named("tau") = tau); 
}



inline arma::mat sampling_beta(arma::mat mu_n, arma::mat chol_Lambda_n_inv, double sigma, int p){

    arma::vec eps = Rcpp::rnorm(p);
    arma::mat output = sigma * chol_Lambda_n_inv * eps + mu_n;
    return output;

}


inline double sampling_sigma(double a_n, double b_0, arma::mat YY, arma::mat mu_n, arma::mat Lambda_n){
    double b_n = b_0 + 0.5 * as_scalar(YY - mu_n.t() * Lambda_n * mu_n);
    double output = 1.0 / sqrt(Rcpp::rgamma(1, a_n, 1.0 / b_n)[0]);
    return output;
}


inline arma::vec sampling_lambda(arma::mat lambda, arma::mat beta, double sigma, double tau, int p){
    // slice sampling for lambda
    // loop over all parameters
    arma::vec gamma_l(p);
    double u1;
    double trunc_limit;
    arma::mat mu2_j = pow(beta / (sigma * tau), 2);
    arma::mat rate_lambda = mu2_j / 2.0;
    arma::vec ub_lambda(p);
    double u2;
    for(int i = 0; i < p; i ++){
        gamma_l(i) = 1.0 / pow(lambda(i) , 2);
        u1 = Rcpp::runif(1, 0, 1.0 / (1.0 + gamma_l(i)))[0];
        trunc_limit = (1.0 - u1) / u1;
        ub_lambda(i) = R::pexp(trunc_limit, 1.0 / rate_lambda(i), 1, 0);
        u2 = Rcpp::runif(1, 0, ub_lambda(i))[0];
        gamma_l(i) = R::qexp(u2, 1.0 / rate_lambda(i), 1, 0);
    }
    gamma_l = 1.0 / arma::sqrt(gamma_l);
    return gamma_l;
}


inline double sampling_tau(arma::mat lambda, arma::mat beta, double sigma, double tau){
    // slice sampling for tau
    double shape_tau = 0.5 * (1.0 + lambda.n_elem);
    double gamma_tt = 1.0 / pow(tau, 2.0);
    double u1 = Rcpp::runif(1, 0, 1.0 / (1.0 + gamma_tt))[0];
    double trunc_limit_tau = (1.0 - u1) / u1;
    double mu2_tau = as_scalar(arma::sum(pow(beta / (sigma * lambda), 2), 1));
    double rate_tau = mu2_tau / 2.0;
    double ub_tau = R::pgamma(trunc_limit_tau, shape_tau, 1.0 / rate_tau, 1, 0);
    double u2 = Rcpp::runif(1, 0, ub_tau)[0];
    gamma_tt = R::qgamma(u2, shape_tau, 1.0 / rate_tau, 1, 0);
    double output = 1.0 / sqrt(gamma_tt);
    return output;
}
