predict.bayeslm.fit = function(object, data, burnin = NULL, X = NULL, ...){
    if(missing(data) && (!is.null(object$call$formula))){
        cat("You train the model with a formula.")
        cat("New data frame is missing, you need to provide one. \n")
        return()
    }

    if(missing(X) && (is.null(object$call$formula))){
        cat("You train the model with matrix input.")
        cat("New covariate matrix is missing, you need to provide one. \n")
    }

    # check input type
    if(!is(object,"bayeslm.fit")){
        cat("wrong input object, should be output of bayeslm.")
        return()
    }
    # compute posterior mean
    N = dim(object$beta)[1] # number of posterior draws
    if(is.null(burnin)){
        burnin = round(0.1 * N)
    }
    if(burnin >= N){
        cat("Error, burnin is larger than number of posterior draws.")
        return()
    }

    cat("Number of posterior draws is ", N, "\n")
    cat("Compute posterior mean with burn-in ", burnin, "\n")

    post_mean = colMeans(object$beta[(burnin+1):N, ])

    # detect type of the fit
    if(is.null(object$call$formula)){
        # when train the data, inputs are y, X matrices rather than formula
        if(object$icept == TRUE){
            if(ncol(X) != (length(post_mean) - 1)){
                cat("There is an intercept, but dimension of X and posterior draws do not match. \n")
                return()
            }
            X = cbind(1, X)
            pred = X %*% post_mean
        }else{
            if(ncol(X) != length(post_mean)){
                cat("No intercept, dimension of X and posterior draws do not match. \n")
                return()
            }
            pred = X %*% post_mean
        }
    }else{
        formula = object$call$formula
        tt = terms(as.formula(formula))
        TT = delete.response(tt)
        mf <- model.frame(TT, data = data)
        X <- model.matrix(attr(mf, "terms"), data = mf)
        pred = X %*% post_mean
    }

    return(pred)
}




