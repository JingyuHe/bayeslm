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

double log_normal_density_matrix_2(const arma::mat &x, const arma::vec &sig_diag, double log_det_sig_diag, bool singular)
{
    // log normal density function, for diagonal precision matrix only
    // input : sig_diag, diagonal elements of precision matrix
    if (singular == true)
    {
        arma::vec deviation = arma::vectorise(x);
        return -0.5 * arma::dot(sig_diag, arma::square(deviation)) + 0.5 * log_det_sig_diag;
    }
    return 0.0;
}

double log_normal_density_matrix_2(const arma::mat &x, const arma::vec &sig_diag, bool singular)
{
    // log normal density function, for diagonal precision matrix only
    // input : sig_diag, diagonal elements of precision matrix
    if (singular == true)
    {
        return log_normal_density_matrix_2(x, sig_diag, sum(log(sig_diag)), singular);
    }
    return 0.0;
}

double log_normal_density_scalar_2(double x, double sig_diag, double log_det_sig_diag, bool singular)
{
    // Scalar analogue of log_normal_density_matrix_2(), used by the
    // one-coefficient ESS fast path under singular designs.
    if (singular == true)
    {
        return -0.5 * sig_diag * x * x + 0.5 * log_det_sig_diag;
    }
    return 0.0;
}

namespace
{

int capped_parallel_threads()
{
    // CRAN policy requires packages to never use more than two threads on the
    // check farm. Respect a user request for 1 thread, but cap anything larger.
    const char *threads_env = std::getenv("RCPP_PARALLEL_NUM_THREADS");
    if (threads_env != nullptr)
    {
        const int requested = std::atoi(threads_env);
        if (requested > 0)
        {
            return std::min(requested, 2);
        }
    }
    return 2;
}

arma::vec select_penalized_beta(const arma::mat &beta, const arma::uvec &penalize_index, bool all_penalized)
{
    if (all_penalized)
    {
        return arma::vectorise(beta);
    }
    return conv_to<vec>::from(beta.rows(penalize_index));
}

arma::vec select_penalized_vec(const arma::vec &x, const arma::uvec &penalize_index, bool all_penalized)
{
    if (all_penalized)
    {
        return x;
    }
    return x.rows(penalize_index);
}

double abs_sum_impl(const arma::mat &beta, const arma::uvec &penalize_index, bool all_penalized)
{
    if (all_penalized)
    {
        return arma::accu(arma::abs(beta));
    }
    return arma::accu(arma::abs(beta.rows(penalize_index)));
}

double square_sum_impl(const arma::mat &beta, const arma::uvec &penalize_index, bool all_penalized)
{
    if (all_penalized)
    {
        return arma::accu(arma::square(beta));
    }
    return arma::accu(arma::square(beta.rows(penalize_index)));
}

arma::mat extract_off_block(const arma::mat &Q, arma::uword block_start, arma::uword block_end)
{
    const arma::uword block_size = block_end - block_start + 1;
    const arma::uword p = Q.n_cols;
    arma::mat Q12(block_size, p - block_size);

    if (block_start > 0)
    {
        Q12.head_cols(block_start) = Q.submat(block_start, 0, block_end, block_start - 1);
    }

    if (block_end + 1 < p)
    {
        Q12.tail_cols(p - block_end - 1) = Q.submat(block_start, block_end + 1, block_end, p - 1);
    }

    return Q12;
}

void compute_precision_block_factors(const arma::mat &Q, arma::uword block_start, arma::uword block_end, arma::mat &mean_factor, arma::mat &chol_factor)
{
    const arma::uword block_size = block_end - block_start + 1;
    const arma::uword p = Q.n_cols;

    // Work directly from the full precision matrix Q. This avoids first
    // materializing Sigma = Q^{-1} and then applying a Schur complement for
    // every block, which was the main fixed-cost bottleneck.
    if (block_size == 1)
    {
        const double q11 = Q(block_start, block_start);
        chol_factor.set_size(1, 1);
        chol_factor(0, 0) = 1.0 / std::sqrt(q11);

        if (p > 1)
        {
            mean_factor.set_size(1, p - 1);
            if (block_start > 0)
            {
                mean_factor.head_cols(block_start) = -Q.submat(block_start, 0, block_start, block_start - 1) / q11;
            }
            if (block_end + 1 < p)
            {
                mean_factor.tail_cols(p - block_end - 1) = -Q.submat(block_start, block_end + 1, block_start, p - 1) / q11;
            }
        }
        else
        {
            mean_factor.zeros(1, 0);
        }
        return;
    }

    arma::mat Q11 = Q.submat(block_start, block_start, block_end, block_end);
    arma::mat Q12 = extract_off_block(Q, block_start, block_end);
    arma::mat L = chol(Q11, "lower");

    // Q11^{-1/2} is the conditional covariance factor, and -Q11^{-1} Q12
    // maps the off-block deviation into the conditional mean.
    chol_factor.eye(block_size, block_size);
    chol_factor = solve(trimatl(L), chol_factor);

    if (Q12.n_cols > 0)
    {
        mean_factor = solve(trimatl(L), Q12);
        mean_factor = solve(trimatu(L.t()), mean_factor);
        mean_factor = -mean_factor;
    }
    else
    {
        mean_factor.zeros(block_size, 0);
    }
}

struct block_factor : public Worker
{
    // input precision matrix, pass by reference
    const arma::mat &Q;
    const arma::vec &V;
    const arma::uvec &block_size_vec;
    arma::field<arma::mat> &output;

    // constructor
    block_factor(const arma::mat &Q, const arma::vec &V, const arma::uvec &block_size_vec, arma::field<arma::mat> &output) : Q(Q), V(V), block_size_vec(block_size_vec), output(output) {}

    // function call operator that work for specified index range
    void operator()(std::size_t begin, std::size_t end)
    {
        int n_blocks = V.n_elem;
        for (std::size_t i = begin; i < end; i++)
        {
            const arma::uword block_start = block_size_vec(i);
            const arma::uword block_end = block_size_vec(i) + static_cast<arma::uword>(V(i)) - 1;
            compute_precision_block_factors(Q, block_start, block_end, output(i), output(i + n_blocks));
        }
    }
};

} // namespace

arma::field<arma::mat> conditional_factors_parallel(const arma::mat &X, const arma::vec &V)
{
    // Layout convention: the first n_blocks fields store conditional-mean
    // factors, and the next n_blocks fields store covariance Cholesky factors.
    int n_blocks = V.n_elem;
    arma::vec V_cumsum = arma::cumsum(V);
    V_cumsum = V_cumsum - V(0);
    arma::uvec block_size_vec = conv_to<arma::uvec>::from(V_cumsum);
    arma::field<arma::mat> output(n_blocks * 2);
    block_factor block_factor(X, V, block_size_vec, output);
    parallelFor(0, n_blocks, block_factor, 1, capped_parallel_threads());

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
    arma::uvec penalize_index = find(penalize > 0);
    return log_horseshoe_approx_prior(beta, v, sigma, penalize_index, penalize_index.n_elem == penalize.n_elem, scale_sigma_prior);
}

double log_horseshoe_approx_prior(const arma::mat &beta, double v, double sigma, const arma::uvec &penalize_index, bool all_penalized, bool scale_sigma_prior)
{
    if (scale_sigma_prior == true)
    {
        v = v * sigma;
    }
    arma::vec beta2 = select_penalized_beta(beta, penalize_index, all_penalized);
    double p = (double)beta2.n_elem;
    beta2 = beta2 / v;
    arma::vec temp = log(log(1.0 + 2.0 / (pow(beta2, 2.0))));
    double ll;
    ll = sum(temp) - log(v) * p;
    return ll;
}

double log_horseshoe_approx_prior_scalar(double beta, double v, double sigma, bool penalized, bool scale_sigma_prior)
{
    if (!penalized)
    {
        return 0.0;
    }
    if (scale_sigma_prior == true)
    {
        v = v * sigma;
    }
    const double beta_scaled = beta / v;
    return std::log(std::log(1.0 + 2.0 / (beta_scaled * beta_scaled))) - std::log(v);
}

arma::vec sample_exp(const arma::vec &lambda)
{
    int n = lambda.n_elem;
    arma::vec sample;
    sample = randu<vec>(n);
    sample = -log(1 - sample) / lambda;
    return (sample);
}

arma::mat scaling(const arma::mat &x)
{
    // This function normalize a matrix x by column
    int n = x.n_rows;
    // int p = x.n_cols;
    arma::mat x_output;
    arma::mat mean_x;
    arma::mat sd_x;
    // normalize each column
    x_output = x;
    mean_x = mean(x, 0);
    sd_x = stddev(x, 0);
    for (int i = 0; i < n; i++)
    {
        x_output.row(i) = (x.row(i) - mean_x) / sd_x;
    }
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
        This function computes block conditional factors from the full precision matrix
        Arguments: X : the full precision matrix
                   V : vector indicates blocks V(i) is number of parameters in block i
        Return value:
            an arma::field objects with length 2 * n_blocks
            the first n_blocks objects are factor of conditional mean
            the next n_blocks objects are cholesky factors of conditional covariance
    */
    int n_blocks = V.n_elem;
    arma::vec V_cumsum = arma::cumsum(V);
    V_cumsum = V_cumsum - V(0);
    arma::uvec block_size_vec = conv_to<arma::uvec>::from(V_cumsum);
    arma::field<arma::mat> output(n_blocks * 2);
    for (int i = 0; i < n_blocks; i++)
    {
        const arma::uword block_start = block_size_vec(i);
        const arma::uword block_end = block_size_vec(i) + static_cast<arma::uword>(V(i)) - 1;
        compute_precision_block_factors(X, block_start, block_end, output(i), output(i + n_blocks));
    }

    return output;
}

double scalar_conditional_mean(const arma::mat &mean_factor, const arma::vec &b, const arma::mat &beta_hat_fixed, arma::uword block_start)
{
    // mean_factor omits the current block, so we multiply it against the
    // concatenated left and right deviations from the unconditional mean.
    double output = 0.0;
    if (block_start > 0)
    {
        output += arma::as_scalar(mean_factor.cols(0, block_start - 1) * (b.head_rows(block_start) - beta_hat_fixed.head_rows(block_start)));
    }
    const arma::uword tail_size = b.n_elem - block_start - 1;
    if (tail_size > 0)
    {
        output += arma::as_scalar(mean_factor.cols(block_start, mean_factor.n_cols - 1) * (b.tail_rows(tail_size) - beta_hat_fixed.tail_rows(tail_size)));
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
    arma::uvec penalize_index = find(penalize > 0);
    return log_asymmetric_prior(beta, vglobal, sigma, prob, penalize_index, penalize_index.n_elem == penalize.n_elem, scale_sigma_prior);
}

double log_asymmetric_prior(const arma::mat &beta, double vglobal, double sigma, const arma::vec &prob, const arma::uvec &penalize_index, bool all_penalized, bool scale_sigma_prior)
{
    if (scale_sigma_prior == true)
    {
        vglobal = vglobal * sigma;
    }
    arma::vec beta_subset = select_penalized_beta(beta, penalize_index, all_penalized);
    arma::vec prob_subset = select_penalized_vec(prob, penalize_index, all_penalized);
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

double log_asymmetric_prior_scalar(double beta, double vglobal, double sigma, double prob, bool penalized, bool scale_sigma_prior)
{
    if (!penalized)
    {
        return 0.0;
    }
    if (scale_sigma_prior == true)
    {
        vglobal = vglobal * sigma;
    }
    const double beta_scaled = beta / vglobal;
    const double s = (1.0 - prob) / prob;
    double output;
    if (beta_scaled > 0.0)
    {
        output = std::log(2.0) + log_cauchy_density(beta_scaled / s) + std::log(1.0 - prob) - std::log(s);
    }
    else
    {
        output = std::log(2.0) + std::log(prob) + log_cauchy_density(beta_scaled);
    }
    return output - std::log(vglobal);
}

double log_cauchy_density(double x)
{
    double output = -1.0 * log(M_PI) - log(1.0 + pow(x, 2));
    return output;
}

double log_nonlocal_prior(const arma::mat &beta, double vglobal, double sigma, const arma::uvec &penalize, const arma::vec &prior_mean, bool scale_sigma_prior)
{
    arma::uvec penalize_index = find(penalize > 0);
    return log_nonlocal_prior(beta, vglobal, sigma, penalize_index, prior_mean, penalize_index.n_elem == penalize.n_elem, scale_sigma_prior);
}

double log_nonlocal_prior(const arma::mat &beta, double vglobal, double sigma, const arma::uvec &penalize_index, const arma::vec &prior_mean, bool all_penalized, bool scale_sigma_prior)
{
    if (scale_sigma_prior == true)
    {
        vglobal = vglobal * sigma;
    }
    // scale by vglobal
    arma::vec prior_mean_subset = select_penalized_vec(prior_mean, penalize_index, all_penalized);
    arma::vec beta_subset = select_penalized_beta(beta, penalize_index, all_penalized);
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

double log_nonlocal_prior_scalar(double beta, double vglobal, double sigma, double prior_mean, bool penalized, bool scale_sigma_prior)
{
    if (!penalized)
    {
        return 0.0;
    }
    if (scale_sigma_prior == true)
    {
        vglobal = vglobal * sigma;
    }
    double output = -std::log(2.0) + log_cauchy_density((beta / 0.25 - prior_mean) / vglobal) - std::log(2.0) + log_cauchy_density((beta / 0.25 + prior_mean) / vglobal);
    return output - std::log(vglobal);
}

double penalized_abs_sum(const arma::mat &beta, const arma::uvec &penalize_index, bool all_penalized)
{
    return abs_sum_impl(beta, penalize_index, all_penalized);
}

double penalized_square_sum(const arma::mat &beta, const arma::uvec &penalize_index, bool all_penalized)
{
    return square_sum_impl(beta, penalize_index, all_penalized);
}

double log_ridge_prior_from_stats(double beta_sq_sum, arma::uword penalized_count, double lambda, double vglobal, double sigma, bool scale_sigma_prior)
{
    if (scale_sigma_prior == true)
    {
        vglobal = vglobal * sigma;
    }
    double p = static_cast<double>(penalized_count);
    double output = beta_sq_sum / (vglobal * vglobal);
    output = p * log(lambda) / 2.0 - lambda / 2.0 * output;
    output = output - p * log(vglobal);
    return output;
}

double log_ridge_prior(const arma::mat &beta, double lambda, double vglobal, double sigma, const arma::uvec &penalize, bool scale_sigma_prior)
{
    arma::uvec penalize_index = find(penalize > 0);
    return log_ridge_prior(beta, lambda, vglobal, sigma, penalize_index, penalize_index.n_elem == penalize.n_elem, scale_sigma_prior);
}

double log_ridge_prior(const arma::mat &beta, double lambda, double vglobal, double sigma, const arma::uvec &penalize_index, bool all_penalized, bool scale_sigma_prior)
{
    const arma::uword penalized_count = all_penalized ? static_cast<arma::uword>(beta.n_elem) : static_cast<arma::uword>(penalize_index.n_elem);
    return log_ridge_prior_from_stats(square_sum_impl(beta, penalize_index, all_penalized), penalized_count, lambda, vglobal, sigma, scale_sigma_prior);
}

double log_ridge_prior_scalar(double beta, double lambda, double vglobal, double sigma, bool penalized, bool scale_sigma_prior)
{
    if (!penalized)
    {
        return 0.0;
    }
    return log_ridge_prior_from_stats(beta * beta, 1, lambda, vglobal, sigma, scale_sigma_prior);
}

double log_laplace_prior_from_stats(double beta_abs_sum, arma::uword penalized_count, double tau, double sigma, double vglobal)
{
    double p = static_cast<double>(penalized_count);
    double out = p * log(tau / 2.0 / sigma);
    out = out - tau / sigma / vglobal * beta_abs_sum;
    out = out - p * log(vglobal);
    return out;
}

double log_laplace_prior(const arma::mat &beta, double tau, double sigma, double vglobal, const arma::uvec &penalize)
{
    arma::uvec penalize_index = find(penalize > 0);
    return log_laplace_prior(beta, tau, sigma, vglobal, penalize_index, penalize_index.n_elem == penalize.n_elem);
}

double log_laplace_prior(const arma::mat &beta, double tau, double sigma, double vglobal, const arma::uvec &penalize_index, bool all_penalized)
{
    return log_laplace_prior_from_stats(abs_sum_impl(beta, penalize_index, all_penalized), all_penalized ? static_cast<arma::uword>(beta.n_elem) : static_cast<arma::uword>(penalize_index.n_elem), tau, sigma, vglobal);
}

double log_laplace_prior_scalar(double beta, double tau, double sigma, double vglobal, bool penalized)
{
    if (!penalized)
    {
        return 0.0;
    }
    return log_laplace_prior_from_stats(std::abs(beta), 1, tau, sigma, vglobal);
}

double log_inverselaplace_prior(const arma::mat &beta, double lambda, double sigma, double vglobal, const arma::uvec &penalize)
{
    arma::uvec penalize_index = find(penalize > 0);
    return log_inverselaplace_prior(beta, lambda, sigma, vglobal, penalize_index, penalize_index.n_elem == penalize.n_elem);
}

double log_inverselaplace_prior(const arma::mat &beta, double lambda, double sigma, double vglobal, const arma::uvec &penalize_index, bool all_penalized)
{
    arma::vec beta_subset = select_penalized_beta(beta, penalize_index, all_penalized);
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

double log_inverselaplace_prior_scalar(double beta, double lambda, double sigma, double vglobal, bool penalized)
{
    if (!penalized || beta == 0.0)
    {
        return 0.0;
    }
    return std::log(lambda / 2.0) - 2.0 * std::log(std::abs(beta)) - lambda / std::abs(beta);
}

arma::mat sampling_beta(const arma::mat &mu_n, const arma::mat &chol_Lambda_n_inv, double sigma, int p, bool scale_sigma_prior)
{

    arma::vec eps = Rcpp::rnorm(p);
    arma::mat output;
    if (scale_sigma_prior == true)
    {
        output = sigma * chol_Lambda_n_inv * eps + mu_n;
    }
    else
    {
        output = sigma * chol_Lambda_n_inv * eps + mu_n;
    }
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
