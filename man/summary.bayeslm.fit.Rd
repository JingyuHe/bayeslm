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
  \item{quantiles }{ logical for should quantiles be displayed (def: \code{TRUE})}
  \item{trailer   }{ logical for should a trailer be displayed (def: \code{TRUE})}
  \item{...       }{ optional arguments for generic function}
}

\details{
This function summarize returned object of function \code{bayeslm}. It prints mean, std Dev, effective sample size (computed by function \code{effectiveSize} in package \code{coda}) coefficients posterior samples. If \code{quantiles=TRUE}, quantiles of marginal distirbutions in the columns of \eqn{X} are displayed.\cr
}

\author{Jingyu He}

\seealso{ \code{\link{summary.mcmc}}}

\examples{
\dontrun{fit=bayeslm(y~x, data = data); summary(fit)}
}

\keyword{univar}