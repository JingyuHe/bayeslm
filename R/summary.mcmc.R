summary.MCMC=function(object,names,burnin=trunc(.1*nrow(X)),quantiles=FALSE,trailer=TRUE,...){
#
# S3 method to compute and print posterior summaries for a matrix of draws
# Modified from package bayesm, author P. Rossi 2/07
#
    numEff= 
    function(x,m=as.integer(min(length(x),(100/sqrt(5000))*sqrt(length(x)))))
    {
    #
    # P. Rossi
    # revision history: 3/27/05
    #
    # purpose:
    #  compute N-W std error and relative numerical efficiency
    # 
    # Arguments:
    #  x is vector of draws
    #  m is number of lags to truncate acf
    #    def is such that m=100 if length(x)= 5000 and grows with sqrt(length)
    #
    # Output:
    #  list with numerical std error and variance multiple (f)
    #
    wgt=as.vector(seq(m,1,-1))/(m+1)
    z=acf(x,lag.max=m,plot=FALSE)
    f=1+2*wgt%*%as.vector(z$acf[-1])
    stderr=sqrt(var(x)*f/length(x))
    list(stderr=stderr,f=f,m=m)
    }
  X = object
  if(mode(X) == "list") stop("list entered \n Please extract sampling matrix from the list \n")
  if(mode(X) !="numeric") stop("Requires numeric argument \n")
  if(is.null(attributes(X)$dim)) X=as.matrix(X)
  nx=ncol(X)
  
  if(is.null(colnames(X))){
      if(missing(names)) names=as.character(1:nx)
  }else{
      names = colnames(X)
  }


  if(nrow(X) < 100) {cat("fewer than 100 draws submitted \n"); return(invisible())}
  X=X[(burnin+1):nrow(X),,drop=FALSE]
  mat=matrix(apply(X,2,mean),nrow=1)
  mat=rbind(mat,sqrt(matrix(apply(X,2,var),nrow=1)))
  num_se=double(nx); rel_eff=double(nx); eff_s_size=double(nx)
  for(i in 1:nx) 
     {out=numEff(X[,i])
      if(is.nan(out$stderr)) 
          {num_se[i]=-9999; rel_eff[i]=-9999; eff_s_size[i]=-9999} 
      else
          {#num_se[i]=out$stderr; rel_eff[i]=out$f; eff_s_size[i]=nrow(X)/ceiling(out$f)
            eff_s_size[i] = effectiveSize(X[,i])
          }
     }
  # mat=rbind(mat,num_se,rel_eff,eff_s_size)

  # construct credible intervals
  mat = t(mat)
  sig = rep(' ', dim(mat)[1])

  for(i in 1:dim(mat)[1]){
      if(abs(mat[i,1]) - qnorm(0.999) * mat[i,2] > 0){
          sig[i] = '***'
      }else if(abs(mat[i,1]) - qnorm(0.99) * mat[i,2] > 0){
          sig[i] = '**'
      }else if(abs(mat[i,1]) - qnorm(0.95) * mat[i,2] > 0){
          sig[i] = '*'
      }else if(abs(mat[i,1]) - qnorm(0.90) * mat[i,2] > 0){
          sig[i] = '.'
      }else{
          sig[i] = ' '
      }

  }
  
  mat = data.frame(mat, eff_s_size, sig)
  rownames(mat)=names
  colnames(mat)[1]="mean"
  colnames(mat)[2]="std dev"
  #rownames(mat)[3]="num se"
  #rownames(mat)[4]="rel eff"
  colnames(mat)[3]="eff sample size"
  colnames(mat)[4]=" "
  
  if(trailer) cat(paste("   based on ",nrow(X)," valid draws (number of burn in is ",burnin,") \n",sep=""))
  cat("Summary of Posterior draws ")
  cat("\n Moments \n")
  print((mat),digits=2)
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
  if(quantiles){
     qmat=apply(X,2,quantile,probs=c(.025,.05,.5,.95,.975))
     colnames(qmat)=names
     cat("\n Quantiles \n")
     print(t(qmat),digits=2)
   }
  return(invisible(1))
  # invisible(t(mat))
}