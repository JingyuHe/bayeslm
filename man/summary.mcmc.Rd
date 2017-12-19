\name{summary.mcmc}
\alias{summary.mcmc}

\title{Summarize posterior draws}

\description{
\code{summary.mcmc} is an S3 method to summarize posterior draws of the model. The input should be a matrix of draws.
}

\usage{
\method{summary}{mcmc}(object,names,burnin=trunc(.1*nrow(X)),quantiles=FALSE,trailer=TRUE,...)
}

\arguments{
  \item{object    }{ \code{object} is a matrix of draws, usually an object of class \code{mcmc}. It's same as \code{X}.}
  \item{names     }{ an optional character vector of names for the columns of \code{X}.}
  \item{burnin    }{ number of draws to burn-in (default value is \eqn{0.1*nrow(X)}).}
  \item{quantiles }{ logical for should quantiles be displayed (def: \code{FALSE}).}
  \item{trailer   }{ logical for should a trailer be displayed (def: \code{TRUE}).}
  \item{...       }{ optional arguments for generic function.}
}

\details{
This function is modified from package \code{bayesm} by Peter Rossi. It summarize object \code{mcmc}. Mean, Std Dev, effective sample size (computed by function \code{effectiveSize} in package \code{coda}) are displayed. If \code{quantiles=TRUE}, quantiles of marginal distirbutions in the columns of \eqn{X} are displayed.\cr

The function also returns significance level, defined by whether the symmetric posterior quantile-based credible interval excludes zero. For example, a regression coefficient with one * has 0.025 quantile and 0.975 quantile with the same sign. Similarly, '***' denotes 0.0005 and 0.9995, '**' denotes 0.005 and 0.995, '*' denotes 0.025 and 0.975, '.' denotes 0.05 and 0.95 quantiles with the same sign.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\seealso{ \code{\link{summary.bayeslm.fit}}}

\examples{
\dontrun{fit=bayeslm(y~x, data = data); summary(fit$beta)}
}

\keyword{univar}