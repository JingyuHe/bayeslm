summary.bayeslm.fit = function(object, names, burnin=NULL,quantiles=FALSE,trailer=TRUE,...){
    # summary fitted object, returned by function bayeslm
    if(is.null(burnin)){
        burnin = 0.2 * dim(object$beta)[1]
    }
    cat("Average number of rejections before one acceptance : \n")
    cat(mean(object$loops), "\n")
    cat("Summary of beta draws \n")
    summary.MCMC(object$beta, names, burnin, quantiles, trailer, ...)
    cat("\n")
    cat("Summary of sigma draws \n")
    cat("Mean of standard deviation is ", mean(object$sigma[-c(1:burnin)]) ,"\n")
    cat("S.d. of standard deviation samples is ", sd(object$sigma[-c(1:burnin)]), "\n")
    cat("Effective sample size of s.d. is ", effectiveSize(object$sigma[-c(1:burnin)]), "\n")
}