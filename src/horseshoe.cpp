#include "../inst/include/bayeslm.h"
/*


blocked elliptical slice sampler, horseshoe prior


*/

#include <chrono>

// [[Rcpp::export]]
List horseshoe(arma::mat Y, arma::mat X, arma::uvec penalize, arma::vec block_vec, int prior_type = 1, Rcpp::Nullable<Rcpp::Function> user_prior_function = R_NilValue, double sigma = 0.5, double s2 = 4, double kap2 = 16,  int nsamps = 10000, int burn = 1000, int skip = 1.0, double vglobal = 1.0, bool verb = false, bool icept = false, bool standardize = true, bool singular = false, bool scale_sigma_prior = true, arma::vec cc = NULL){

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
    if(standardize == true){
        sdy = as_scalar(stddev(Y));
    }else{
        // if do not standardize regressors, set SD to 1
        sdy = 1.0;
        sdx.fill(1.0);
    }
    
    // intercepts
    if(icept == true){
        // if add a column of ones for intercept
        if(standardize == true){
            X = scaling(X);
            Y = Y / sdy;
        }
        X = arma::join_rows(arma::ones<mat>(n, 1), X);
        block_vec = join_cols(ones<vec>(1), block_vec); // put intercept in the first block
        p = p + 1;
        N_blocks = N_blocks + 1;
        sdx = join_rows(ones<mat>(1,1), sdx);
        penalize = arma::join_cols(arma::zeros<uvec>(1), penalize); // add one indicator of penalization for intercept. Do not penalize intercept    
    }else{
        if(standardize == true){
            Y = scaling(Y);
            X = scaling(X);
        }
    }

 
    // compute sufficient statistics
    arma::mat YY = trans(Y) * Y;
    arma::mat YX = trans(Y) * X;
    arma::mat XX = trans(X) * X;
    arma::uvec penalize_index = find(penalize > 0);

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
    arma::vec v(p);
    v.fill(1.0);
    double vgprop;


    arma::mat Sigma_inv = XX;
    arma::vec eta = 1.0 / cc; // 1/c in the paper, precision of the prior
    arma::mat M0;
    arma::mat Sigma;
    
    if(singular == true){
        // if matrix X is singular, use the "conjugate regression" type adjustment
        M0 = arma::diagmat(eta);
        Sigma = inv(XX + M0);
        beta_hat = Sigma * (trans(YX));
    }else{
        Sigma = inv(XX);
        beta_hat = Sigma * (trans(YX));
    }

    // a initial value of the derivation from the mean
    beta = 0.1 * beta_hat;
    
    //initialize vectors to save posterior samples
    arma::mat bsamps(p, nsamps);
    bsamps.fill(0.0);
    arma::vec ssamps(nsamps);
    ssamps.fill(0.0);
    arma::vec vsamps(nsamps);
    vsamps.fill(0.0);
    arma::vec ssq_out(nsamps);
    ssq_out.fill(0.0);
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
    // arma::field<arma::mat> output = conditional_factors(Sigma, block_vec);

    // parallel versions
    arma::field<arma::mat> output = conditional_factors_parallel(Sigma, block_vec);

    arma::field<arma::mat> mean_factors = output.rows(0, N_blocks-1);
    arma::field<arma::mat> chol_factors = output.rows(N_blocks, 2 * N_blocks - 1);

    arma::vec beta_condition;
    arma::mat beta_hat_block;
    arma::mat beta_hat_fixed = beta_hat;
    arma::mat beta_block;

    // There are three beta related vectors
    // beta_hat_fixed, uncodintional mean of betas
    // beta_hat_lastround, conditional mean of previous round
    // beta deviace from the mean

    arma::uvec block_indexes(N_blocks + 1);
    arma::vec block_cum_count = arma::cumsum(block_vec);
    block_indexes(0) = 0;
    for(int i = 0; i < N_blocks; i ++){
        block_indexes(i+1) = block_cum_count(i);
    }

    
    arma::vec beta_hat_lastround = beta_hat;

    auto t1 = std::chrono::high_resolution_clock::now();
    auto time_fixed = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
    Rcpp::Rcout << "fixed running time " << time_fixed << endl;

    t0 = std::chrono::high_resolution_clock::now();
    while (h < nsamps)
    {

        if(verb == true && h % 1000 == 0 && h > 1 && iter % skip == 0){
            Rprintf("%d\n", h);
        }

        // sampling beta
        // a gibbs sampler for all blocks
        for(int i = 0; i < N_blocks; i ++){
            // loop over all blocks, sample each block by ellipitical slice sampler
            // compute conditional mean

            beta_condition = join_cols(b.head_rows(block_indexes(i)) - beta_hat_fixed.head_rows(block_indexes(i)), b.tail_rows(p - block_indexes(i+1)) - beta_hat_fixed.tail_rows(p - block_indexes(i+1)));
            beta_hat_block = beta_hat_fixed.rows(block_indexes(i), block_indexes(i+1)-1) + mean_factors(i) * beta_condition;
        
            // define ellipse
            //eps = rnorm((uword) block_vec(i));
            eps = randn<vec>(block_vec(i), 1);
            nu = chol_factors(i) * eps;
            nu = s * nu;

            // acceptance threshold
            priorcomp = log_horseshoe_approx_prior(b.rows(block_indexes(i), block_indexes(i+1)-1) , vglobal, s, penalize.rows(block_indexes(i), block_indexes(i+1)-1), scale_sigma_prior)- log_normal_density_matrix(b.rows(block_indexes(i), block_indexes(i+1)-1), arma::diagmat(eta.rows(block_indexes(i), block_indexes(i+1)-1)) / pow(s, 2), singular);

            u = arma::as_scalar(randu(1));

            ly = priorcomp + log(u);

            // subtract the new conditional mean from last draw
            beta_block = b.rows(block_indexes(i), block_indexes(i+1) - 1) - beta_hat_block;
            
            thetaprop = arma::as_scalar(randu(1)) * 2 * M_PI;

            betaprop = beta_block * cos(thetaprop) + nu * sin(thetaprop);

            thetamin = thetaprop - 2.0 * M_PI;

            thetamax = thetaprop;
            
            if(i == 0 && icept == true){

                b.subvec(block_indexes(i), block_indexes(i+1) - 1) = betaprop + beta_hat_block;

            }else{
                while (log_horseshoe_approx_prior(beta_hat_block + betaprop, vglobal, s, penalize.rows(block_indexes(i), block_indexes(i+1)-1), scale_sigma_prior) - log_normal_density_matrix(beta_hat_block + betaprop, arma::diagmat(eta.rows(block_indexes(i), block_indexes(i+1)-1)) / pow(s,2), singular) < ly){
                    
                    loopcount += 1;

                    if (thetaprop < 0){

                        thetamin = thetaprop;
                    
                    }else{
                        
                        thetamax = thetaprop;
                    
                    }

                    thetaprop = runif(1, thetamin, thetamax)[0];

                    betaprop = beta_block * cos(thetaprop) + nu * sin(thetaprop);

                }

                b.subvec(block_indexes(i), block_indexes(i+1) - 1) = betaprop + beta_hat_block;
            }
        }


        // update the global shrinkage parameter
        vgprop = exp(log(vglobal) + arma::as_scalar(randn(1)) * 0.2);


        // if there is no intercept, pass the full vector
        ratio = exp(log_horseshoe_approx_prior(b, vgprop, s, penalize, scale_sigma_prior) - log_horseshoe_approx_prior(b, vglobal, s, penalize, scale_sigma_prior) + log(vgprop) - log(vglobal));


        if(as_scalar(randu(1)) < ratio){
            vglobal = vgprop;
        }

        
        
        // update sigma
        if(scale_sigma_prior == false){
            ssq = as_scalar(YY) - 2.0 * as_scalar(YX * (b)) + as_scalar(trans(b) * XX * (b));
        }else{
            ssq = as_scalar(YY) - 2.0 * as_scalar(YX * (b)) + as_scalar(trans(b) * XX * (b)) + as_scalar(trans(b) * b);
        }
        
        s = 1.0 / sqrt(arma::as_scalar(arma::randg(1, distr_param((n + kap2) / 2.0, 2.0 / (ssq + s2)))));

        iter = iter + 1;

        if (iter > burn)
        {
            if (iter % skip == 0)
            {
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
    auto time_sampling = 1.e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(t1-t0).count();
    Rcpp::Rcout << "sampling time " << time_sampling << endl;

    // X and Y were scaled at the beginning, rescale estimations
    for(int ll = 0; (unsigned) ll < bsamps.n_cols; ll ++){
        bsamps.col(ll) = bsamps.col(ll) / trans(sdx) * sdy;
    }


    ssamps = ssamps * sdy;


    bsamps = trans(bsamps);


    return List::create(Named("loops") = loops, Named("sigma") = ssamps, Named("vglobal") = vsamps, Named("beta") = bsamps);
}
