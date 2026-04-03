#include "../inst/include/bayeslm.h"
/*


blocked elliptical slice sampler, inverseLaplace prior
https://arxiv.org/abs/2001.08327


*/

// [[Rcpp::export]]
List inverseLaplace_cpp_loop(arma::mat Y, arma::mat X, double lambda, arma::uvec penalize, arma::vec block_vec, arma::vec cc, int prior_type = 1, double sigma = 0.5, double s2 = 4, double kap2 = 16, int nsamps = 10000, int burn = 1000, int skip = 1, double vglobal = 1.0, bool sampling_vglobal = true, bool verb = false, bool icept = false, bool standardize = true, bool singular = false, bool scale_sigma_prior = true)
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

    /*
    the input of penalize is (0,1,1,0,1...)
    convert to indeces of 1
    */

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
        beta_hat = solve(Precision, trans(YX), arma::solve_opts::likely_sympd);
    }
    else
    {
        Precision = XX;
        beta_hat = solve(Precision, trans(YX), arma::solve_opts::likely_sympd);
    }

    // a initial value of the derivation from the mean
    beta = 0.1 * beta_hat;

    // initialize vectors to save posterior samples
    arma::mat bsamps(p, nsamps);
    bsamps.fill(0.0);
    arma::vec ssamps(nsamps);
    ssamps.fill(0.0);
    arma::vec vsamps(nsamps);
    vsamps.fill(0.0);
    int loopcount = 0;
    arma::vec loops(nsamps);
    loops.fill(0);
    double u;
    arma::mat nu;
    nu.fill(0.0);
    arma::vec eps(p);
    eps.fill(0.0);
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
    arma::mat beta_block;

    // There are three beta related vectors
    // beta_hat_fixed, uncodintional mean of betas
    // beta deviace from the mean

    arma::uvec block_indexes(N_blocks + 1);
    arma::vec block_cum_count = arma::cumsum(block_vec);
    block_indexes(0) = 0;
    for (int i = 0; i < N_blocks; i++)
    {
        block_indexes(i + 1) = block_cum_count(i);
    }
    arma::field<arma::uvec> penalize_index_blocks(N_blocks);
    arma::uvec penalize_all_blocks(N_blocks);
    arma::field<arma::vec> eta_blocks;
    arma::vec eta_logsum_blocks;
    if (singular == true)
    {
        eta_blocks.set_size(N_blocks);
        eta_logsum_blocks.set_size(N_blocks);
    }
    for (int i = 0; i < N_blocks; i++)
    {
        const arma::uword block_start = block_indexes(i);
        const arma::uword block_end = block_indexes(i + 1) - 1;
        penalize_index_blocks(i) = find(penalize.subvec(block_start, block_end) > 0);
        penalize_all_blocks(i) = penalize_index_blocks(i).n_elem == static_cast<arma::uword>(block_end - block_start + 1);
        if (singular == true)
        {
            eta_blocks(i) = eta.rows(block_start, block_end);
            eta_logsum_blocks(i) = arma::sum(arma::log(eta_blocks(i)));
        }
    }

    // double tau = 1.0;
    // double tau_prop = 0.0;
    // arma::vec tausamps(nsamps);
    // tausamps.fill(0.0);

    auto t1 = std::chrono::high_resolution_clock::now();
    auto time_fixed = 1.e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    if (verb == true)
    {
        Rcpp::Rcout << "fixed running time " << time_fixed << endl;
    }

    t0 = std::chrono::high_resolution_clock::now();
    while (h < nsamps)
    {

        if (verb == true && h % 1000 == 0 && h > 1 && iter % skip == 0)
        {
            Rprintf("%d\n", h);
        }

        double log_s = 0.0;
        double inv_s_sq = 0.0;
        if (singular == true)
        {
            log_s = std::log(s);
            inv_s_sq = 1.0 / (s * s);
        }

        // sampling beta
        // a gibbs sampler for all blocks
        for (int i = 0; i < N_blocks; i++)
        {
            // loop over all blocks, sample each block by ellipitical slice sampler
            // compute conditional mean
            const arma::uword block_start = block_indexes(i);
            const arma::uword block_end = block_indexes(i + 1) - 1;
            const bool scalar_block = (block_vec(i) == 1);
            const arma::uvec &penalize_index_block = penalize_index_blocks(i);
            const bool all_penalized_block = penalize_all_blocks(i);

            if (scalar_block)
            {
                // Scalar branch for the common one-coefficient block case.
                const double beta_hat_scalar = beta_hat_fixed(block_start) + scalar_conditional_mean(mean_factors(i), b, beta_hat_fixed, block_start);
                const double current_beta = b(block_start);
                eps = randn<vec>(block_vec(i), 1);
                const double nu_scalar = s * chol_factors(i)(0, 0) * eps(0);
                const bool penalized_scalar = all_penalized_block;
                const double Sigma_scalar_inv_diag = singular ? eta_blocks(i)(0) * inv_s_sq : 0.0;
                const double Sigma_scalar_log_det = singular ? (eta_logsum_blocks(i) - 2.0 * log_s) : 0.0;

                priorcomp = log_inverselaplace_prior_scalar(current_beta, lambda, s, vglobal, penalized_scalar);
                if (singular == true)
                {
                    priorcomp = priorcomp - log_normal_density_scalar_2(current_beta, Sigma_scalar_inv_diag, Sigma_scalar_log_det, singular);
                }

                u = arma::as_scalar(randu(1));
                ly = priorcomp + log(u);

                const double beta_block_scalar = current_beta - beta_hat_scalar;
                thetaprop = arma::as_scalar(randu(1)) * 2 * M_PI;
                double betaprop_scalar = beta_block_scalar * cos(thetaprop) + nu_scalar * sin(thetaprop);
                thetamin = thetaprop - 2.0 * M_PI;
                thetamax = thetaprop;

                if (i == 0 && icept == true)
                {
                    b(block_start) = betaprop_scalar + beta_hat_scalar;
                }
                else
                {
                    while (log_inverselaplace_prior_scalar(beta_hat_scalar + betaprop_scalar, lambda, s, vglobal, penalized_scalar) - log_normal_density_scalar_2(beta_hat_scalar + betaprop_scalar, Sigma_scalar_inv_diag, Sigma_scalar_log_det, singular) < ly)
                    {

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

                        betaprop_scalar = beta_block_scalar * cos(thetaprop) + nu_scalar * sin(thetaprop);
                    }

                    b(block_start) = betaprop_scalar + beta_hat_scalar;
                }
                continue;
            }

            // Grouped blocks still use the original vector-valued update.
            beta_condition = join_cols(b.head_rows(block_start) - beta_hat_fixed.head_rows(block_start), b.tail_rows(p - block_indexes(i + 1)) - beta_hat_fixed.tail_rows(p - block_indexes(i + 1)));
            beta_hat_block = beta_hat_fixed.rows(block_start, block_end) + mean_factors(i) * beta_condition;

            // define ellipse
            // eps = rnorm((uword) block_vec(i));
            eps = randn<vec>(block_vec(i), 1);
            nu = chol_factors(i) * eps;
            nu = s * nu;
            arma::vec Sigma_block_inv_diag;
            double Sigma_block_log_det = 0.0;
            if (singular == true)
            {
                Sigma_block_inv_diag = eta_blocks(i) * inv_s_sq;
                Sigma_block_log_det = eta_logsum_blocks(i) - 2.0 * static_cast<double>(eta_blocks(i).n_elem) * log_s;
            }

            // acceptance threshold
            priorcomp = log_inverselaplace_prior(b.rows(block_start, block_end), lambda, s, vglobal, penalize_index_block, all_penalized_block);
            if (singular == true)
            {
                priorcomp = priorcomp - log_normal_density_matrix_2(b.rows(block_start, block_end), Sigma_block_inv_diag, Sigma_block_log_det, singular);
            }

            u = arma::as_scalar(randu(1));

            ly = priorcomp + log(u);

            // subtract the new conditional mean from last draw
            beta_block = b.rows(block_start, block_end) - beta_hat_block;

            thetaprop = arma::as_scalar(randu(1)) * 2 * M_PI;

            betaprop = beta_block * cos(thetaprop) + nu * sin(thetaprop);

            thetamin = thetaprop - 2.0 * M_PI;

            thetamax = thetaprop;

            if (i == 0 && icept == true)
            {

                b.subvec(block_start, block_end) = betaprop + beta_hat_block;
            }
            else
            {
                while (log_inverselaplace_prior(beta_hat_block + betaprop, lambda, s, vglobal, penalize_index_block, all_penalized_block) - log_normal_density_matrix_2(beta_hat_block + betaprop, Sigma_block_inv_diag, Sigma_block_log_det, singular) < ly)
                {

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
                }

                b.subvec(block_start, block_end) = betaprop + beta_hat_block;
            }
        }

        // update tau
        // tau_prop = exp(log(tau) + arma::as_scalar(randn(1)) * 0.05 );

        // ratio = exp(log_inverselaplace_prior(b, tau_prop, s, vglobal, penalize) + log_normal_density(tau_prop, 0.0, 100.0)  - log_inverselaplace_prior(b, tau, s, vglobal, penalize) - log_normal_density(tau, 0.0, 100.0)  +log(tau_prop) - log(tau));

        // if(as_scalar(randu(1)) < ratio){
        //     tau = tau_prop;
        // }

        // if(sampling_vglobal){
        //     // update the global shrinkage parameter
        //     vgprop = exp(log(vglobal) + arma::as_scalar(randn(1)) * 0.2);
        //     // if there is no intercept, pass the full vector
        //     ratio = exp(log_horseshoe_approx_prior(b, vgprop, s, penalize, scale_sigma_prior) - log_horseshoe_approx_prior(b, vglobal, s, penalize, scale_sigma_prior) + log(vgprop) - log(vglobal));

        //     if(as_scalar(randu(1)) < ratio){
        //         vglobal = vgprop;
        //     }
        // }

        // update sigma
        if (scale_sigma_prior == false)
        {
            ssq = as_scalar(YY) - 2.0 * as_scalar(YX * (b)) + as_scalar(trans(b) * XX * (b));
        }
        else
        {
            ssq = as_scalar(YY) - 2.0 * as_scalar(YX * (b)) + as_scalar(trans(b) * XX * (b)) + as_scalar(trans(b) * b);
        }

        s = 1.0 / sqrt(arma::as_scalar(arma::randg(1, distr_param((n + kap2) / 2.0, 2.0 / (ssq + s2)))));

        iter = iter + 1;

        if (iter > burn)
        {
            if (iter % skip == 0)
            {
                // tausamps(h) = tau;
                bsamps.col(h) = b;
                ssamps(h) = s;
                vsamps(h) = vglobal;
                // ssq_out(h) = ssq;
                loops(h) = loopcount;
                h = h + 1;
            }
        }

        // re-count for the next round.
        loopcount = 0;
    }

    t1 = std::chrono::high_resolution_clock::now();
    auto time_sampling = 1.e-9 * std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
    if (verb == true)
    {
        Rcpp::Rcout << "sampling time " << time_sampling << endl;
    }

    // X and Y were scaled at the beginning, rescale estimations
    for (int ll = 0; (unsigned)ll < bsamps.n_cols; ll++)
    {
        bsamps.col(ll) = bsamps.col(ll) / trans(sdx) * sdy;
    }

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

    return List::create(Named("loops") = loops, Named("sigma") = ssamps, Named("vglobal") = vsamps, Named("beta") = bsamps);
}
