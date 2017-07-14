#include  "bayeslm.h"



double log_normal_density_matrix(arma::mat x, arma::mat Sigma_inv, bool singular){
    double output = 0.0;
    arma::mat deviation;
    if(singular == true){
        deviation = x;
        output = -0.5 * as_scalar(trans(deviation) * Sigma_inv * deviation);
        output = output + 0.5 * log(det(Sigma_inv));
    }else{
        output = 0.0;
    }
    return(output);
}


double betaprior(arma::mat beta, arma::vec v, int prior, Rcpp::Nullable<Rcpp::Function> user_prior_function){
    
    double output = 0;

    if(user_prior_function.isNotNull()){

        Function user_prior_function_2(user_prior_function);

        output = user_prior_function_wrapper(beta, v, user_prior_function_2);
    }
    else{
        switch (prior){
            case 1:
                output = log_horseshoe_approx_prior(beta, v);
                break;
            case 2:
                output = log_double_exp_prior(beta, v);
                break;
            case 3:
                output = log_normal_prior(beta, v);
                break;
            case 4:
                output = log_cauchy_prior(beta, v);
                break;
            default:
                Rprintf("Wrong input of prior types.\n");
            }
    }
    return output;
}

double user_prior_function_wrapper(arma::mat beta, arma::vec v, Rcpp::Function f)
{
    SEXP result = f(beta, v);
    double output = Rcpp::as<double>(result);
    return output;
}



double log_horseshoe_approx_prior(arma::mat beta, arma::vec v)
{
    arma::vec beta2 = conv_to<vec>::from(beta);
    beta2 = beta2 / v;
    arma::vec temp = log(log(1.0 + 2.0 / (pow(beta2, 2.0))));
    double ll;
    ll = sum(temp) - sum(log(v));
    return ll;
}

arma::vec sample_exp(arma::vec lambda)
{
    int n = lambda.n_elem;
    arma::vec sample;
    sample = randu<vec>(n);
    sample =  - log(1 - sample) / lambda;
    return (sample);
}

arma::mat scaling(arma::mat x){
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
    for(int i = 0; i < n; i++){
        x_output.row(i) = (x.row(i) - mean_x) / sd_x;
    }
    return x_output;
}



double log_double_exp_prior(arma::mat beta, arma::vec v){
    // log density of double exponential prior
    arma::vec beta2 = conv_to<vec>::from(beta);
    beta2 = beta2 / v;
    arma::vec temp = (-1.0) * abs(beta2);
    double ll;
    ll = sum(temp) - sum(log(v));
    return ll;
}

double log_cauchy_prior(arma::mat beta, arma::vec v){
    // log density of Cauchy prior
    arma::vec beta2 = conv_to<vec>::from(beta);
    beta2 = beta2 / v;
    arma::vec temp = log(1.0 + pow(beta2, 2.0));
    double ll;
    ll = (-1.0) * sum(temp) - sum(log(v));
    return ll;
}

double log_normal_prior(arma::mat beta, arma::vec v){
    // log density of normal prior
    arma::vec beta2 = conv_to<vec>::from(beta);
    beta2 = beta2 / v;
    arma::vec temp = pow(beta2, 2.0);
    double ll;
    ll = (- 1.0 / 2.0) * sum(temp) - sum(log(v));
    return ll;
}


arma::field<arma::mat> conditional_factors(arma::mat X, arma::vec V){
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
  for(int i = 0; i < n_blocks; i++){
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



double log_normal_density(double x, double mu, double sigma){
    // returns log density of normal(mu, sigma)
    double output = -0.5 * log(2.0 * M_PI) - log(sigma) - pow((x - mu), 2) / 2.0 / pow(sigma, 2);
    return(output);
}


double log_asymmetric_prior(arma::mat beta, double vglobal, arma::vec prob){
    // scale by vglobal
    beta = beta / vglobal;
    int p = beta.n_elem;
    arma::vec s = (1.0 - prob) / prob;
    arma::vec result = zeros<vec>(p);
    for(int i = 0; i < p; i ++){
        if(beta(i) > 0){
            result(i) = log(2) + log_cauchy_density(beta(i) / s(i)) + log(1.0 - prob(i)) - log(s(i));
        }else{
            result(i) = log(2) + log(prob(i)) + log_cauchy_density(beta(i));
        }
    }
    double output = as_scalar(sum(result));
    // Jacobian of scaling by vglobal;
    output = output - p * log(vglobal);
    return output;
}



double log_cauchy_density(double x){
    double output = -1.0 * log(M_PI) - log(1.0 + pow(x, 2));
    return output;
}


double log_nonlocal_prior(arma::mat beta, double vglobal, arma::vec prior_mean){
    // scale by vglobal
    int p = beta.n_elem;
    arma::vec result = zeros<vec>(p);
    for(int i = 0; i < p; i ++){
        result(i) = - log(2) + log_cauchy_density((beta(i) / 0.25 - prior_mean(i)) /vglobal) - log(2) + log_cauchy_density((beta(i) / 0.25  + prior_mean(i)) / vglobal);
    }
    double output = as_scalar(sum(result));
    output = output - p * log(vglobal);
    return output;
}






double log_ridge_prior(arma::mat beta, double lambda, double vglobal){
    // lambda is in variance scale
    // beta ~ N(0, 1 / lambda)
    beta = beta / vglobal;
    int p = beta.n_elem;
    double output = arma::as_scalar(arma::sum(arma::pow(beta, 2)));
    output = p * log(lambda) / 2.0 - lambda / 2.0 * output;
    output = output - p * log(vglobal);
    return output;
}





double log_laplace_prior(arma::mat beta, double tau, double sigma, double vglobal){
    int p = beta.n_elem;
    beta = beta / vglobal;
    double out = p * log(tau / 2.0 / sigma);
    out = out - tau / sigma * arma::as_scalar(arma::sum(abs(beta)));
    out = out - p * log(vglobal);
    return out;
}




