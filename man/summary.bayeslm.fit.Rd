\name{summary.bayeslm.fit}
\alias{summary.bayeslm.fit}

\title{Summarize fitted object of \code{bayeslm}}

\description{
\code{summary.bayeslm.fit} is an S3 method to summarize returned object of function \code{bayeslm}. The input should be \code{bayeslm.fit} object.
}

\usage{
\method{summary}{bayeslm.fit}(object,names,burnin=NULL,quantiles=FALSE,trailer=TRUE,...)
}

\arguments{
  \item{object    }{ \code{object} is a fitted object, returned by function \code{bayeslm}.}
  \item{names     }{ an optional character vector of names for all the coefficients.}
  \item{burnin    }{ number of draws to burn-in (if it is \code{NULL}, will set default value as \eqn{0.2*nrow(object\$beta)})}
  \item{quantiles }{ logical for should quantiles be displayed (def: \code{FALSE})}
  \item{trailer   }{ logical for should a trailer be displayed (def: \code{TRUE})}
  \item{...       }{ optional arguments for generic function}
}

\details{
This function summarize returned object of function \code{bayeslm}. It prints mean, std Dev, effective sample size (computed by function \code{effectiveSize} in package \code{coda}) coefficients posterior samples. If \code{quantiles=TRUE}, quantiles of marginal distirbutions in the columns of \eqn{X} are displayed.\cr

The function also returns significance level, defined by whether the symmetric posterior quantile-based credible interval excludes zero. For example, a regression coefficient with one * has 0.025 quantile and 0.975 quantile with the same sign. Similarly, '***' denotes 0.0005 and 0.9995, '**' denotes 0.005 and 0.995, '*' denotes 0.025 and 0.975, '.' denotes 0.05 and 0.95 quantiles with the same sign.
}

\author{Jingyu He}

\seealso{ \code{\link{summary.mcmc}}}

\examples{
x = matrix(rnorm(1000), 100, 10)
y = x \%*\% rnorm(10) + rnorm(100)
fit=bayeslm(y~x)
summary(fit)
}

\keyword{summary}