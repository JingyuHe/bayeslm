bayeslm.default <- function(Y, X = FALSE, prior = "horseshoe", penalize = NULL, block_vec = NULL, sigma = NULL, s2 = 1, kap2 = 1, N = 20000L, burnin = 0L, thinning = 1L, vglobal = 1, verb = FALSE, icept = TRUE, standardize = TRUE, singular = FALSE, scale_sigma_prior = TRUE, prior_mean = NULL, prob_vec = NULL, cc = NULL, ...){

    Y = as.matrix(Y)

    X = as.matrix(X)

    if(is.null(cc)){
        cc = rep(1, dim(X)[2])
    }

    if(icept == TRUE){
        cc = c(1, cc)
    }

    try(if(dim(Y)[1] != dim(X)[1]) stop("Length of X and Y don't agree."))

    # scale the initial value
    if(is.null(sigma)){
        sigma = 0.5 * stats::sd(Y)
    }

    if(is.null(block_vec)){
        block_vec = rep(1, dim(X)[2])
    }

    if(is.null(penalize)){
        penalize = rep(1, dim(X)[2])
    }

    if(is.null(colnames(X))){
        Xnames = as.character(1:dim(X)[2])
        Xnames = paste("X", Xnames, sep = '')
    }else{
        Xnames = colnames(X)
    }

    if(icept == TRUE){
        Xnames = c("intercept", Xnames)
    }

    
    prior_type = 1L

    user_prior_function = NULL
    

    if(prior == "horseshoe"){
        cat("horseshoe prior \n")
        output = horseshoe_cpp_loop(Y, X, penalize, block_vec, prior_type, user_prior_function, sigma, s2, kap2, N, burnin, thinning, vglobal, verb, icept, standardize, singular, scale_sigma_prior, cc)
    }else if(prior == "laplace"){
        cat("laplace prior \n")
        output = blasso_cpp_loop(Y, X, penalize, block_vec, prior_type, sigma, s2, kap2, N, burnin, thinning, vglobal, verb, icept, standardize, singular, scale_sigma_prior, cc)
    }else if(prior == "ridge"){
        cat("ridge prior \n")
        output = bridge_cpp_loop(Y, X, penalize, block_vec, prior_type, sigma, s2, kap2, N, burnin, thinning, vglobal, verb, icept, standardize, singular, scale_sigma_prior, cc)
    }else if(prior == "nonlocal"){
        cat("nonlocal prior \n")
        if(is.null(prior_mean) == TRUE){
            cat("prior mean of nonlocal prior is not specified, use default value 1.5 for all coefficients \n")
            if(icept == TRUE){
                prior_mean = rep(0, dim(X)[2]+1)
            }else{
                prior_mean = rep(0, dim(X)[2])
            }
        }
        output = nonlocal_cpp_loop(Y, X, prior_mean, penalize, block_vec, prior_type, sigma, s2, kap2, N, burnin, thinning, vglobal, verb, icept, standardize, singular, scale_sigma_prior, cc)
    }else if(prior == "sharkfin"){
        cat("sharkfin prior \n")
        if(is.null(prob_vec) == TRUE){
            cat("prior parameter q of nonlocal prior is not specified, use default value 0.25 for all coefficients \n")
            if(icept == TRUE){
                prob_vec = rep(0, dim(X)[2] + 1)
            }else{
                prob_vec = rep(0, dim(X)[2])
            }
        }
        output = sharkfin_cpp_loop(Y, X, prob_vec, penalize, block_vec, prior_type, sigma, s2, kap2, N, burnin, thinning, vglobal, verb, icept, standardize, singular, scale_sigma_prior, cc)
    }else{
        cat("wrong prior \n")
    }

    attributes(output$sigma)$class = c("mcmc")
    attributes(output$vglobal)$class = c("mcmc")
    attributes(output$beta)$class = c("mcmc")
    attributes(output)$class = c("bayeslm.fit")
    colnames(output$beta) = Xnames
    colnames(output$sigma) = "sigma"    
    
    N = dim(output$beta)[1]

    if(icept==FALSE){
        fittedvalue = X %*% matrix(colMeans(output$beta[floor(0.2 * N):N, ]), ncol = 1)
    }else{
        fittedvalue = cbind(rep(1,dim(X)[1]), X) %*% matrix(colMeans(output$beta[floor(0.2 * N):N, ]), ncol = 1)
    }
    residuals = Y - fittedvalue

    output$call = NULL
    output$fitted.value = fittedvalue
    output$residuals = residuals


    return(output)
}