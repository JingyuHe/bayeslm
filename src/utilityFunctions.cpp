#include "../inst/include/bayeslm.h"

double log_normal_density_matrix(const arma::mat &x, const arma::mat &Sigma_inv, bool singular)
{
    double output = 0.0;
    if (singular == true)
    {
        output = -0.5 * as_scalar(trans(x) * Sigma_inv * x);
        output = output + 0.5 * log(det(Sigma_inv));
    }
    else
    {
        output = 0.0;
    }
    return (output);
}

double log_normal_density_matrix_2(arma::mat x, arma::vec sig_diag, bool singular)
{
    // log normal density function, for diagonal covariance matrix only
    // input : sig_diag, diagonal elements of covariance matrix
    double output = 0.0;
    arma::mat deviation;
    if (singular == true)
    {
    }
    else
    {
        output = 0.0;
    }
    return (output);
}

struct block_factor : public Worker
{
    // input matrix, pass by reference
    const arma::mat &X;
    const arma::vec &V;
    const arma::uvec &block_size_vec;
    arma::field<arma::mat> &output;

    // constructor
    block_factor(const arma::mat &X, const arma::vec &V, const arma::uvec &block_size_vec, arma::field<arma::mat> &output) : X(X), V(V), block_size_vec(block_size_vec), output(output) {}

    // function call operator that work for specified index range
    void operator()(std::size_t begin, std::size_t end)
    {
        int n_blocks = V.n_elem;
        arma::mat temp;
        arma::mat temp11;
        arma::mat temp12;
        arma::mat temp22;
        int block_size_cumulate = 0;
        int block_size_cumulate_next = 0;
        for (std::size_t i = begin; i < end; i++)
        {
            block_size_cumulate = block_size_vec(i);
            block_size_cumulate_next = block_size_vec(i) + V(i);
            // other subsequent blocks, subscripts begin from block_size_cumulate to block_size_cumulate + V(i)
            temp11 = X.submat(block_size_cumulate, block_size_cumulate, block_size_cumulate_next - 1, block_size_cumulate_next - 1);
            temp12 = X.rows(block_size_cumulate, block_size_cumulate_next - 1);
            temp12.shed_cols(block_size_cumulate, block_size_cumulate_next - 1);
            temp22 = X;
            temp22.shed_rows(block_size_cumulate, block_size_cumulate_next - 1);
            temp22.shed_cols(block_size_cumulate, block_size_cumulate_next - 1);
            temp = solve(temp22, temp12.t()); // faster than inv(temp22) * temp21
            output(i + n_blocks) = chol(temp11 - temp12 * temp, "lower");
            output(i) = temp.t();
            // block_size_cumulate += V(i);
        }
    }
};

arma::field<arma::mat> conditional_factors_parallel(const arma::mat &X, const arma::vec &V)
{
    int n_blocks = V.n_elem;
    arma::vec V_cumsum = arma::cumsum(V);
    V_cumsum = V_cumsum - V(0);
    arma::uvec block_size_vec = conv_to<arma::uvec>::from(V_cumsum);
    arma::field<arma::mat> output(n_blocks * 2);
    block_factor block_factor(X, V, block_size_vec, output);
    parallelFor(0, n_blocks, block_factor);

    return output;
}

// double betaprior(arma::mat beta, arma::vec v, int prior, Rcpp::Nullable<Rcpp::Function> user_prior_function){

//     double output = 0;

//     if(user_prior_function.isNotNull()){

//         Function user_prior_function_2(user_prior_function);

//         output = user_prior_function_wrapper(beta, v, user_prior_function_2);
//     }
//     else{
//         switch (prior){
//             case 1:
//                 output = log_horseshoe_approx_prior(beta, v);
//                 break;
//             case 2:
//                 output = log_double_exp_prior(beta, v);
//                 break;
//             case 3:
//                 output = log_normal_prior(beta, v);
//                 break;
//             case 4:
//                 output = log_cauchy_prior(beta, v);
//                 break;
//             default:
//                 Rprintf("Wrong input of prior types.\n");
//             }
//     }
//     return output;
// }

double user_prior_function_wrapper(const arma::mat &beta, const arma::vec &v, Rcpp::Function f)
{
    SEXP result = f(beta, v);
    double output = Rcpp::as<double>(result);
    return output;
}

double log_horseshoe_approx_prior(const arma::mat &beta, double v, double sigma, const arma::uvec &penalize, bool scale_sigma_prior)
{
    if (scale_sigma_prior == true)
    {
        v = v * sigma;
    }
    arma::uvec penalize_index = find(penalize > 0);
    arma::vec beta2 = conv_to<vec>::from(beta.rows(penalize_index));
    double p = (double)beta2.n_elem;
    beta2 = beta2 / v;
    arma::vec temp = log(log(1.0 + 2.0 / (pow(beta2, 2.0))));
    double ll;
    ll = sum(temp) - log(v) * p;
    return ll;
}

arma::vec sample_exp(const arma::vec &lambda)
{
    arma::vec sample = randu<vec>(lambda.n_elem);
    sample = -log(1 - sample) / lambda;
    return (sample);
}

arma::mat scaling(const arma::mat &x)
{
    // This function normalize a matrix x by column
    arma::mat x_output = x;
    arma::rowvec mean_x = mean(x, 0);
    arma::rowvec sd_x = stddev(x, 0);
    x_output.each_row() -= mean_x;
    x_output.each_row() /= sd_x;
    return x_output;
}

double log_double_exp_prior(const arma::mat &beta, const arma::vec &v)
{
    // log density of double exponential prior
    arma::vec beta2 = conv_to<vec>::from(beta);
    beta2 = beta2 / v;
    arma::vec temp = (-1.0) * abs(beta2);
    double ll;
    ll = sum(temp) - sum(log(v));
    return ll;
}

double log_cauchy_prior(const arma::mat &beta, const arma::vec &v)
{
    // log density of Cauchy prior
    arma::vec beta2 = conv_to<vec>::from(beta);
    beta2 = beta2 / v;
    arma::vec temp = log(1.0 + pow(beta2, 2.0));
    double ll;
    ll = (-1.0) * sum(temp) - sum(log(v));
    return ll;
}

double log_normal_prior(const arma::mat &beta, const arma::vec &v)
{
    // log density of normal prior
    arma::vec beta2 = conv_to<vec>::from(beta);
    beta2 = beta2 / v;
    arma::vec temp = pow(beta2, 2.0);
    double ll;
    ll = (-1.0 / 2.0) * sum(temp) - sum(log(v));
    return ll;
}

arma::field<arma::mat> conditional_factors(const arma::mat &X, const arma::vec &V)
{
    /*
        This function computes cholesky factor (lower triangular) of conditional covariance for all blocks
        Arguments: X : the full covariance
                   V : vector indicates blocks V(i) is number of parameters in block i
        Return value:
            an arma::field objects with length 2 * n_blocks
            the first n_blocks objects are factor of conditional mean
            the next n_blocks objects are cholesky factors of conditional covariance
    */
    int n_blocks = V.n_elem;
    arma::field<arma::mat> output(n_blocks * 2);
    arma::mat temp;
    arma::mat temp11;
    arma::mat temp12;
    arma::mat temp22;
    int block_size_cumulate = 0;
    int block_size_cumulate_next = 0;
    for (int i = 0; i < n_blocks; i++)
    {
        block_size_cumulate_next = block_size_cumulate + V(i);
        // other subsequent blocks, subscripts begin from block_size_cumulate to block_size_cumulate + V(i)
        temp11 = X.submat(block_size_cumulate, block_size_cumulate, block_size_cumulate_next - 1, block_size_cumulate_next - 1);
        temp12 = X.rows(block_size_cumulate, block_size_cumulate_next - 1);
        temp12.shed_cols(block_size_cumulate, block_size_cumulate_next - 1);
        temp22 = X;
        temp22.shed_rows(block_size_cumulate, block_size_cumulate_next - 1);
        temp22.shed_cols(block_size_cumulate, block_size_cumulate_next - 1);
        temp = solve(temp22, temp12.t()); // faster than inv(temp22) * temp21
        output(i + n_blocks) = chol(temp11 - temp12 * temp, "lower");
        output(i) = temp.t();
        block_size_cumulate += V(i);
    }

    return output;
}

double log_normal_density(double x, double mu, double sigma)
{
    // returns log density of normal(mu, sigma)
    double output = -0.5 * log(2.0 * M_PI) - log(sigma) - pow((x - mu), 2) / 2.0 / pow(sigma, 2);
    return (output);
}

double log_asymmetric_prior(const arma::mat &beta, double vglobal, double sigma, const arma::vec &prob, const arma::uvec &penalize, bool scale_sigma_prior)
{
    if (scale_sigma_prior == true)
    {
        vglobal = vglobal * sigma;
    }
    arma::uvec penalize_index = find(penalize > 0);
    arma::vec beta_subset = conv_to<vec>::from(beta.rows(penalize_index));
    arma::vec prob_subset = prob.rows(penalize_index);
    // scale by vglobal
    beta_subset = beta_subset / vglobal;
    double p = (double)beta_subset.n_elem;
    arma::vec s = (1.0 - prob_subset) / prob_subset;
    arma::vec result(p);
    for (int i = 0; i < p; i++)
    {
        if (beta_subset(i) > 0)
        {
            result(i) = log(2) + log_cauchy_density(beta_subset(i) / s(i)) + log(1.0 - prob_subset(i)) - log(s(i));
        }
        else
        {
            result(i) = log(2) + log(prob_subset(i)) + log_cauchy_density(beta_subset(i));
        }
    }
    double output = as_scalar(sum(result));
    // Jacobian of scaling by vglobal;
    output = output - p * log(vglobal);
    return output;
}

double log_cauchy_density(double x)
{
    double output = -1.0 * log(M_PI) - log(1.0 + pow(x, 2));
    return output;
}

double log_nonlocal_prior(const arma::mat &beta, double vglobal, double sigma, const arma::uvec &penalize, const arma::vec &prior_mean, bool scale_sigma_prior)
{
    if (scale_sigma_prior == true)
    {
        vglobal = vglobal * sigma;
    }
    // scale by vglobal
    arma::uvec penalize_index = find(penalize > 0);
    arma::vec prior_mean_subset = prior_mean.rows(penalize_index);
    arma::vec beta_subset = conv_to<vec>::from(beta.rows(penalize_index)); // pick only penalized elements, otherwise flat (constant) prior
    double p = (double)beta_subset.n_elem;
    arma::vec result(p);
    for (int i = 0; i < p; i++)
    {
        result(i) = -log(2) + log_cauchy_density((beta_subset(i) / 0.25 - prior_mean_subset(i)) / vglobal) - log(2) + log_cauchy_density((beta_subset(i) / 0.25 + prior_mean_subset(i)) / vglobal);
    }
    double output = as_scalar(sum(result));
    output = output - p * log(vglobal);
    return output;
}

double log_ridge_prior(const arma::mat &beta, double lambda, double vglobal, double sigma, const arma::uvec &penalize, bool scale_sigma_prior)
{
    // lambda is in variance scale
    // beta ~ N(0, 1 / lambda)
    if (scale_sigma_prior == true)
    {
        vglobal = vglobal * sigma;
    }
    arma::uvec penalize_index = find(penalize > 0);
    arma::mat beta_subset = beta.rows(penalize_index); // only care about elements with penalization

    beta_subset = beta_subset / vglobal;
    double p = (double)beta_subset.n_elem;
    double output = arma::as_scalar(arma::sum(arma::pow(beta_subset, 2)));
    output = p * log(lambda) / 2.0 - lambda / 2.0 * output;
    output = output - p * log(vglobal);
    return output;
}

double log_laplace_prior(const arma::mat &beta, double tau, double sigma, double vglobal, const arma::uvec &penalize)
{
    arma::uvec penalize_index = find(penalize > 0);
    arma::vec beta_subset = conv_to<vec>::from(beta.rows(penalize_index));
    double p = (double)beta_subset.n_elem;
    beta_subset = beta_subset / vglobal;
    double out = p * log(tau / 2.0 / sigma);
    out = out - tau / sigma * arma::as_scalar(arma::sum(abs(beta_subset)));
    out = out - p * log(vglobal);
    return out;
}

double log_inverselaplace_prior(const arma::mat &beta, double lambda, double sigma, double vglobal, const arma::uvec &penalize)
{
    arma::uvec penalize_index = find(penalize > 0);
    arma::vec beta_subset = conv_to<vec>::from(beta.rows(penalize_index));
    double p = (double)beta_subset.n_elem;
    // beta = beta / vglobal;

    double out = 0.0;

    double temp = log(lambda / 2.0);

    for (size_t i = 0; i < p; i++)
    {
        if (beta_subset(i) != 0.0)
        {
            out = out + temp - 2 * log(abs(beta_subset(i))) - lambda / abs(beta_subset(i));
        }
    }

    return out;
}

arma::mat sampling_beta(const arma::mat &mu_n, const arma::mat &chol_Lambda_n_inv, double sigma, int p, bool scale_sigma_prior)
{

    arma::vec eps = Rcpp::rnorm(p);
    arma::mat output = sigma * chol_Lambda_n_inv * eps + mu_n;
    return output;
}

double sampling_sigma(double a_n, double b_0, const arma::mat &YY, const arma::mat &mu_n, const arma::mat &Lambda_n)
{
    double b_n = b_0 + 0.5 * as_scalar(YY - mu_n.t() * Lambda_n * mu_n);
    double output = 1.0 / sqrt(Rcpp::rgamma(1, a_n, 1.0 / b_n)[0]);
    return output;
}

arma::vec sampling_lambda(const arma::mat &lambda, const arma::mat &beta, double sigma, double tau, int p, bool scale_sigma_prior)
{
    // slice sampling for lambda
    // loop over all parameters
    arma::vec gamma_l(p);
    double u1;
    double trunc_limit;
    arma::mat mu2_j;
    if (scale_sigma_prior == true)
    {
        mu2_j = pow(beta / (sigma * tau), 2);
    }
    else
    {
        mu2_j = pow(beta / tau, 2);
    }

    arma::mat rate_lambda = mu2_j / 2.0;
    arma::vec ub_lambda(p);
    double u2;
    for (int i = 0; i < p; i++)
    {
        gamma_l(i) = 1.0 / pow(lambda(i), 2);
        u1 = Rcpp::runif(1, 0, 1.0 / (1.0 + gamma_l(i)))[0];
        trunc_limit = (1.0 - u1) / u1;
        ub_lambda(i) = R::pexp(trunc_limit, 1.0 / rate_lambda(i), 1, 0);
        u2 = Rcpp::runif(1, 0, ub_lambda(i))[0];
        gamma_l(i) = R::qexp(u2, 1.0 / rate_lambda(i), 1, 0);
    }
    gamma_l = 1.0 / arma::sqrt(gamma_l);
    return gamma_l;
}

double sampling_tau(const arma::mat &lambda, const arma::mat &beta, double sigma, double tau, bool scale_sigma_prior)
{
    // slice sampling for tau
    double shape_tau = 0.5 * (1.0 + lambda.n_elem);
    double gamma_tt = 1.0 / pow(tau, 2.0);
    double u1 = Rcpp::runif(1, 0, 1.0 / (1.0 + gamma_tt))[0];
    double trunc_limit_tau = (1.0 - u1) / u1;
    double mu2_tau;
    if (scale_sigma_prior == true)
    {
        mu2_tau = as_scalar(arma::sum(pow(beta / (sigma * lambda), 2), 1));
    }
    else
    {
        mu2_tau = as_scalar(arma::sum(pow(beta / lambda, 2), 1));
    }
    double rate_tau = mu2_tau / 2.0;
    double ub_tau = R::pgamma(trunc_limit_tau, shape_tau, 1.0 / rate_tau, 1, 0);
    double u2 = Rcpp::runif(1, 0, ub_tau)[0];
    gamma_tt = R::qgamma(u2, shape_tau, 1.0 / rate_tau, 1, 0);
    double output = 1.0 / sqrt(gamma_tt);
    return output;
}
