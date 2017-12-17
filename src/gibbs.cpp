#include "../inst/include/bayeslm.h"
/*


Standard Gibbs Sampler for horseshoe regression


*/


// [[Rcpp::export]]
List hs_gibbs(arma::mat Y, arma::mat X, int nsamps = 1000, double a = 1, double b = 1, bool scale_sigma_prior = true){
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
    beta.row(0) = zeros<mat>(1, p);
    lambda.row(0) = ones<mat>(1, p);
    Lambda0.diag() = lambda.row(0);
    sigma(0) = 1;
    tau(0) = 1;

    // scaling
    double sdy = as_scalar(stddev(Y));
    arma::mat sdx = stddev(X, 0);
    X = scaling(X);
    Y = Y / sdy;
    X = X / double(sqrt(n - 1.0));         
    sdx = sdx * double(sqrt(n - 1.0));  

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
    Rcpp::Rcout << "fixed time " << t / CLOCKS_PER_SEC << endl;

    t = clock();
    for(int i = 1; i < nsamps; i ++){

        Lambda0.diag() = 1.0 / (pow(lambda.row(i-1), 2) * pow(tau(i-1), 2));
        
        Lambda_n = XX + Lambda0;

        Lambda_n_inv = inv(Lambda_n);

        chol_Lambda_n_inv = chol(Lambda_n_inv, "lower"); // LL' = inv(Lambda_n)

        mu_n = Lambda_n_inv * XY;

        beta.row(i) = sampling_beta(mu_n, chol_Lambda_n_inv, sigma(i-1), p, scale_sigma_prior).t();

        sigma(i) = sampling_sigma(a_n, b_0, YY, mu_n, Lambda_n);

        lambda.row(i) = sampling_lambda(lambda.row(i-1), beta.row(i), sigma(i), tau(i-1), p, scale_sigma_prior).t();

        tau(i) = sampling_tau(lambda.row(i), beta.row(i), sigma(i), tau(i-1), scale_sigma_prior);
        
    }

    // scale back
    for(int ll = 0; (unsigned) ll < beta.n_rows; ll ++){
        beta.row(ll) = beta.row(ll) / sdx * sdy;
    }


    t = clock() - t;
    Rcpp::Rcout << "float time" << t / CLOCKS_PER_SEC << endl;

    return List::create(Named("beta") = beta, Named("lambda") = lambda, Named("sigma") = sigma, Named("tau") = tau); 
}
