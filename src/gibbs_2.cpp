#include "../inst/include/bayeslm.h"
/*


Standard Gibbs Sampler for horseshoe regression


*/


// [[Rcpp::export]]
List hs_gibbs_2(arma::mat Y, arma::mat X, int nsamps = 1000, double a = 1, double b = 1, bool scale_sigma_prior = true){
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

        beta.row(i) = sampling_beta_2(beta.row(i-1), mu_n, (double) sigma(i-1), p, scale_sigma_prior, lambda.row(i-1), tau(i-1), X, Y);

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






arma::mat sampling_beta_2(arma::mat old_beta, arma::mat mu_n, double sigma, int p, bool scale_sigma_prior, arma::mat lambda, double tau, arma::mat X, arma::mat Y){
    // prior over beta_i is tau^2 * lambda_i^2 * sigma^2
    arma::mat new_beta;
    new_beta = old_beta;
    arma::mat res;
    arma::mat X_current;
    arma::mat X_others;
    arma::mat beta_others;
    arma::mat beta1;
    double M1;
    arma::vec eps;
    for(int i = 0; i < p; i ++ ){
        // loop over all betas 
        // compute residual Y - beta_{-i} * X_{-i}
        X_current = X.col(i);
        X_others = X;
        X_others.shed_col(i);
        beta_others = new_beta;
        beta_others.shed_col(i);
        res = Y - X_others * beta_others.t();
        M1 = as_scalar(1.0 / pow(tau, 2) / pow(lambda.col(i), 2)) + as_scalar(X_current.t() * X_current);
        beta1 = 1.0 / M1 * X_current.t() * res;
        // Rcpp::Rcout << beta1 << endl;
        // update
        eps = Rcpp::rnorm(1);
        new_beta.col(i) = beta1 + eps * sqrt(pow(sigma, 2) / M1);
    }
    return new_beta;
}

