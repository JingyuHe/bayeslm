bayeslm.formula <- function(formula, data = list(), Y = FALSE, X = FALSE, prior = "horseshoe", penalize = NULL, block_vec = NULL, sigma = NULL, s2 = 1, kap2 = 1, N = 20000L, burnin = 0L, thinning = 1L, vglobal = 1, sampling_vglobal = TRUE, verb = FALSE, standardize = TRUE, singular = FALSE, scale_sigma_prior = TRUE, prior_mean = NULL, prob_vec = NULL, cc = NULL, lambda = NULL, ...){


    if (!inherits(formula, "formula"))
        stop("method is only for formula objects")

    this.call = match.call()

    # mf <- match.call(expand.dots = FALSE)

    # cat("mf \n")
    # print(mf)

    # m <- match(c("formula", "data", "weights", "offset"), names(mf), 0)
    # cat("m \n")
    # print(m)

    # mf <- mf[c(1, m)]
    # print(mf)

    # mf$drop.unused.levels <- TRUE
    # mf$na.action <- na.pass
    # mf[[1]] <- as.name("model.frame")
    # print(mf)

    # m <- mf
    # mf <- eval(mf, parent.frame())
    # Terms <- attr(mf, "terms")

    # Y <- model.response(mf)

    # W <- model.weights(mf)
    # offset <- model.offset(mf)

    # var.names <- attributes(Terms)$term.labels

    # X <- model.frame(terms(reformulate(var.names)), data, na.action=na.pass)

    # cat("dim ", dim(X), "\n")



    mf <- model.frame(formula = formula, data = data)

    # m <- match(c("formula", "data", "prior", "penalize", "block_vec", "sigma", "s2", "kap2", "N", "burnin", "thinning", "vglobal", "verb", "icept", "standardize", "singular", "prior_mean", "prob_vec"), names(mf), 0)

    # mf <- mf[m]

    # m <- mf

    # mf <- eval(mf, parent.frame())

    X <- model.matrix(attr(mf, "terms"), data = mf)

    Xnames = colnames(X)

    Y <- model.response(mf)

    if(sum(X[,1]) == dim(X)[1]){
        X = X[,-1]
        icept = TRUE
    }else{
        icept = FALSE
    }

    Y = as.matrix(Y)

    X = as.matrix(X)


    if(is.null(NULL)){
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
    
    prior_type = 1L
    

    if(prior == "horseshoe"){
        cat("horseshoe prior \n")
        output = horseshoe_cpp_loop(Y, X, penalize, block_vec, cc, prior_type, sigma, s2, kap2, N, burnin, thinning, vglobal, sampling_vglobal, verb, icept, standardize, singular, scale_sigma_prior)
    }else if(prior == "laplace"){
        cat("laplace prior \n")
        output = blasso_cpp_loop(Y, X, penalize, block_vec, cc, prior_type, sigma, s2, kap2, N, burnin, thinning, vglobal, sampling_vglobal, verb, icept, standardize, singular, scale_sigma_prior)
    }else if(prior == "ridge"){
        cat("ridge prior \n")
        output = bridge_cpp_loop(Y, X, penalize, block_vec, cc, prior_type, sigma, s2, kap2, N, burnin, thinning, vglobal, sampling_vglobal, verb, icept, standardize, singular, scale_sigma_prior)
    }else if(prior == "inverselaplace"){
        cat("inverse Laplace prior \n")
        if(is.null(lambda)){
            cat("lambda is missing for inverse Laplace prior. \n")
        }        
        output = inverseLaplace_cpp_loop(Y, X, lambda, penalize, block_vec, cc, prior_type, sigma, s2, kap2, N, burnin, thinning, vglobal, sampling_vglobal, verb, icept, standardize, singular, scale_sigma_prior) 
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
        output = nonlocal_cpp_loop(Y, X, prior_mean, penalize, block_vec, cc, prior_type, sigma, s2, kap2, N, burnin, thinning, vglobal, sampling_vglobal, verb, icept, standardize, singular, scale_sigma_prior)
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
        output = sharkfin_cpp_loop(Y, X, prob_vec, penalize, block_vec, cc, prior_type, sigma, s2, kap2, N, burnin, thinning, vglobal, sampling_vglobal, verb, icept, standardize, singular, scale_sigma_prior)
    }else{
        cat("wrong prior \n")
    }
    N = dim(output$beta)[1]

    if(icept==FALSE){
        fittedvalue = X %*% matrix(colMeans(output$beta[floor(0.2 * N):N, ]), ncol = 1)
    }else{
        fittedvalue = cbind(rep(1,dim(X)[1]), X) %*% matrix(colMeans(output$beta[floor(0.2 * N):N, ]), ncol = 1)
    }
    residuals = Y - fittedvalue

    attributes(output$sigma)$class = c("MCMC", "matrix")
    attributes(output$vglobal)$class = c("MCMC", "matrix")
    attributes(output$beta)$class = c("MCMC", "matrix")
    attributes(output)$class = c("bayeslm.fit")
    colnames(output$beta) = Xnames
    colnames(output$sigma) = "sigma"
    output$terms = terms(formula)
    
    output$call = this.call
    output$fitted.value = fittedvalue
    output$residuals = residuals
    output$icept = icept

    return(output)
}