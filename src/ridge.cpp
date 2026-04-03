#include "../inst/include/bayeslm.h"
/*


blocked elliptical slice sampler, ridge prior


*/

// [[Rcpp::export]]
List bridge_cpp_loop(arma::mat Y, arma::mat X, arma::uvec penalize, arma::vec block_vec, arma::vec cc, int prior_type = 1, double sigma = 0.5, double s2 = 4, double kap2 = 16, int nsamps = 10000, int burn = 1000, int skip = 1, double vglobal = 1.0, bool sampling_vglobal = true, bool verb = false, bool icept = false, bool standardize = true, bool singular = false, bool scale_sigma_prior = true)
{

    auto t0 = std::chrono::high_resolution_clock::now();

    arma::vec beta_hat;
    arma::vec beta;

    // dimensions
    int n = X.n_rows;
    int p = X.n_cols;
    int N_blocks = block_vec.n_elem; // number of blocks

    // compute standard derivation
    double sdy;
    arma::mat sdx = stddev(X, 0);
    if (standardize == true)
    {
        sdy = as_scalar(stddev(Y));
    }
    else
    {
        // if do not standardize regressors, set SD to 1
        sdy = 1.0;
        sdx.fill(1.0);
    }

    arma::vec X_mean(X.n_cols);

    // intercepts
    if (icept == true)
    {
        // if add a column of ones for intercept
        if (standardize == true)
        {
            X_mean = arma::mean(X, 0).t();
            X = scaling(X);
            Y = Y / sdy;
        }
        X = arma::join_rows(arma::ones<mat>(n, 1), X);
        block_vec = join_cols(ones<vec>(1), block_vec); // put intercept in the first block
        p = p + 1;
        N_blocks = N_blocks + 1;
        sdx = join_rows(ones<mat>(1, 1), sdx);
        penalize = arma::join_cols(arma::zeros<uvec>(1), penalize); // add one indicator of penalization for intercept. Do not penalize intercept
    }
    else
    {
        if (standardize == true)
        {
            Y = scaling(Y);
            X = scaling(X);
        }
    }

    // compute sufficient statistics
    arma::mat YY = trans(Y) * Y;
    arma::mat YX = trans(Y) * X;
    arma::mat XX = trans(X) * X;
    double YY_scalar = as_scalar(YY);
    double sigma_shape = (n + kap2) / 2.0;
    arma::vec sdx_vec = arma::vectorise(sdx.t());
    burn = burn + 1;

    double s = sigma;
    double ssq;
    double ly;
    double thetaprop;
    double thetamin;
    double thetamax;
    arma::vec b;
    double vgprop;

    arma::vec eta = 1.0 / cc; // 1/c in the paper, precision of the prior
    arma::mat M0;
    arma::mat Precision;

    if (singular == true)
    {
        // if matrix X is singular, use the "conjugate regression" type adjustment
        M0 = arma::diagmat(eta);
        Precision = XX + M0;
        beta_hat = arma::solve(Precision, trans(YX), arma::solve_opts::likely_sympd);
    }
    else
    {
        Precision = XX;
        beta_hat = arma::solve(Precision, trans(YX), arma::solve_opts::likely_sympd);
    }

    // a initial value of the derivation from the mean
    beta = 0.1 * beta_hat;

    // initialize vectors to save posterior samples
    arma::mat bsamps(p, nsamps);
    arma::vec ssamps(nsamps);
    arma::vec vsamps(nsamps);
    int loopcount = 0;
    arma::vec loops(nsamps);
    double u;
    arma::mat nu;
    arma::vec eps(p);
    arma::vec betaprop;
    double priorcomp;
    int iter = 0;
    int h = 0;
    double ratio = 0.0;

    // initial value  deviation + mean
    b = beta + beta_hat;

    /*
        pre-loop computation for conditional mean and covariance matrix given other blocks
    */
    arma::field<arma::mat> output = conditional_factors_parallel(Precision, block_vec);
    arma::field<arma::mat> mean_factors = output.rows(0, N_blocks - 1);
    arma::field<arma::mat> chol_factors = output.rows(N_blocks, 2 * N_blocks - 1);

    arma::vec beta_condition;
    arma::mat beta_hat_block;
    arma::mat beta_hat_fixed = beta_hat;
    arma::vec beta_dev = b - beta_hat_fixed;
    arma::mat beta_block;

    // There are three beta related vectors
    // beta_hat_fixed, uncodintional mean of betas
    // beta deviace from the mean

    arma::uvec block_indexes(N_blocks + 1);
    arma::uvec block_tail_sizes(N_blocks);
    arma::vec block_cum_count = arma::cumsum(block_vec);
    block_indexes(0) = 0;
    for (int i = 0; i < N_blocks; i++)
    {
        block_indexes(i + 1) = block_cum_count(i);
        block_tail_sizes(i) = p - block_indexes(i + 1);
    }
    arma::field<arma::mat> beta_hat_fixed_blocks(N_blocks);
    arma::field<arma::uvec> penalize_blocks(N_blocks);
    arma::field<arma::vec> eta_blocks;
    arma::field<arma::mat> eta_diag_blocks;
    if (singular == true)
    {
        eta_blocks.set_size(N_blocks);
        eta_diag_blocks.set_size(N_blocks);
    }
    for (int i = 0; i < N_blocks; i++)
    {
        beta_hat_fixed_blocks(i) = beta_hat_fixed.rows(block_indexes(i), block_indexes(i + 1) - 1);
        penalize_blocks(i) = penalize.subvec(block_indexes(i), block_indexes(i + 1) - 1);
        if (singular == true)
        {
            eta_blocks(i) = eta.rows(block_indexes(i), block_indexes(i + 1) - 1);
            eta_diag_blocks(i) = arma::diagmat(eta_blocks(i));
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto time_fixed = 1.e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    if (verb)
    {
        Rcpp::Rcout << "fixed running time " << time_fixed << endl;
    }

    t0 = std::chrono::high_resolution_clock::now();
    double lambda_prop = 1.0;
    arma::vec lambdasamps(nsamps);
    double lambda = 2.0;

    while (h < nsamps)
    {

        if (verb == true && h % 1000 == 0 && h > 1 && iter % skip == 0)
        {
            Rprintf("%d\n", h);
        }

        // sampling beta
        // a gibbs sampler for all blocks
        for (int i = 0; i < N_blocks; i++)
        {
            // loop over all blocks, sample each block by ellipitical slice sampler
            // compute conditional mean
            const arma::uword block_start = block_indexes(i);
            const arma::uword block_end = block_indexes(i + 1) - 1;
            const arma::uword block_tail_size = block_tail_sizes(i);
            const bool skip_acceptance = (i == 0 && icept == true);
            const arma::uvec &penalize_block = penalize_blocks(i);
            const arma::mat b_block_current = b.rows(block_start, block_end);

            beta_condition.set_size(block_start + block_tail_size);
            if (block_start > 0)
            {
                beta_condition.head_rows(block_start) = beta_dev.head_rows(block_start);
            }
            if (block_tail_size > 0)
            {
                beta_condition.tail_rows(block_tail_size) = beta_dev.tail_rows(block_tail_size);
            }
            beta_hat_block = beta_hat_fixed_blocks(i) + mean_factors(i) * beta_condition;

            // define ellipse
            // eps = rnorm((uword) block_vec(i));
            eps = randn<vec>(block_vec(i), 1);
            nu = chol_factors(i) * eps;
            nu = s * nu;
            arma::mat Sigma_block_inv;
            if (singular == true && !skip_acceptance)
            {
                Sigma_block_inv = eta_diag_blocks(i) / pow(s, 2);
            }

            // acceptance threshold
            if (!skip_acceptance)
            {
                priorcomp = log_ridge_prior(b_block_current, lambda, vglobal, s, penalize_block, scale_sigma_prior);
                if (singular == true)
                {
                    priorcomp = priorcomp - log_normal_density_matrix(b_block_current, Sigma_block_inv, singular);
                }
            }

            u = arma::as_scalar(randu(1));

            if (!skip_acceptance)
            {
                ly = priorcomp + log(u);
            }

            // subtract the new conditional mean from last draw
            beta_block = b_block_current - beta_hat_block;

            thetaprop = arma::as_scalar(randu(1)) * 2 * M_PI;

            betaprop = beta_block * cos(thetaprop) + nu * sin(thetaprop);
            arma::mat proposed_beta = beta_hat_block + betaprop;

            thetamin = thetaprop - 2.0 * M_PI;

            thetamax = thetaprop;

            if (i == 0 && icept == true)
            {

                b.subvec(block_start, block_end) = proposed_beta;
                beta_dev.subvec(block_start, block_end) = arma::vectorise(proposed_beta - beta_hat_fixed_blocks(i));
            }
            else
            {
                while (true)
                {
                    priorcomp = log_ridge_prior(proposed_beta, lambda, vglobal, s, penalize_block, scale_sigma_prior);
                    if (singular == true)
                    {
                        priorcomp = priorcomp - log_normal_density_matrix(proposed_beta, Sigma_block_inv, singular);
                    }
                    if (priorcomp >= ly)
                    {
                        break;
                    }

                    loopcount += 1;

                    if (thetaprop < 0)
                    {

                        thetamin = thetaprop;
                    }
                    else
                    {

                        thetamax = thetaprop;
                    }

                    thetaprop = runif(1, thetamin, thetamax)[0];

                    betaprop = beta_block * cos(thetaprop) + nu * sin(thetaprop);
                    proposed_beta = beta_hat_block + betaprop;
                }

                b.subvec(block_start, block_end) = proposed_beta;
                beta_dev.subvec(block_start, block_end) = arma::vectorise(proposed_beta - beta_hat_fixed_blocks(i));
            }
        }

        // update lambda
        lambda_prop = exp(log(lambda) + arma::as_scalar(randn(1)) * 0.05);

        ratio = exp(log_ridge_prior(b, lambda_prop, vglobal, s, penalize, scale_sigma_prior) + log_normal_density(lambda_prop, 0.0, 100.0) - log_ridge_prior(b, lambda, vglobal, s, penalize, scale_sigma_prior) - log_normal_density(lambda, 0.0, 100.0) + log(lambda_prop) - log(lambda));

        if (as_scalar(randu(1)) < ratio)
        {
            lambda = lambda_prop;
        }

        if (sampling_vglobal)
        {
            // update the global shrinkage parameter
            vgprop = exp(log(vglobal) + arma::as_scalar(randn(1)) * 0.2);
            // if there is no intercept, pass the full vector
            ratio = exp(log_horseshoe_approx_prior(b, vgprop, s, penalize, scale_sigma_prior) - log_horseshoe_approx_prior(b, vglobal, s, penalize, scale_sigma_prior) + log(vgprop) - log(vglobal));

            if (as_scalar(randu(1)) < ratio)
            {
                vglobal = vgprop;
            }
        }

        // update sigma
        if (scale_sigma_prior == false)
        {
            ssq = YY_scalar - 2.0 * as_scalar(YX * (b)) + as_scalar(trans(b) * XX * (b));
        }
        else
        {
            ssq = YY_scalar - 2.0 * as_scalar(YX * (b)) + as_scalar(trans(b) * XX * (b)) + as_scalar(trans(b) * b);
        }

        s = 1.0 / sqrt(arma::as_scalar(arma::randg(1, distr_param(sigma_shape, 2.0 / (ssq + s2)))));

        iter = iter + 1;

        if (iter > burn)
        {
            if (iter % skip == 0)
            {
                lambdasamps(h) = lambda;
                bsamps.col(h) = b;
                ssamps(h) = s;
                vsamps(h) = vglobal;
                loops(h) = loopcount;
                h = h + 1;
            }
        }

        // re-count for the next round.
        loopcount = 0;
    }

    t1 = std::chrono::high_resolution_clock::now();
    auto time_sampling = 1.e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    if (verb)
    {
        Rcpp::Rcout << "sampling time " << time_sampling << endl;
    }

    // X and Y were scaled at the beginning, rescale estimations
    bsamps.each_col() /= sdx_vec;
    bsamps *= sdy;

    // adjust intercept if standardize all X variables
    if (icept && standardize)
    {
        for (size_t ll = 0; ll < bsamps.n_cols; ll++)
        {
            bsamps(0, ll) = bsamps(0, ll) - arma::sum(bsamps.submat(1, ll, p - 1, ll) % X_mean);
        }
    }

    ssamps = ssamps * sdy;

    bsamps = trans(bsamps);

    return List::create(Named("loops") = loops, Named("sigma") = ssamps, Named("vglobal") = vsamps, Named("beta") = bsamps, Named("lambda") = lambdasamps);
}
