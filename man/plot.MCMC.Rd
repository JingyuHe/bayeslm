\name{plot.MCMC}
\alias{plot.MCMC}

\title{Plot posterior summary}

\description{
\code{plot.MCMC} is an S3 method to plot empirical distribution of posterior draws. The input is a \code{MCMC} matrix
}

\usage{
\method{plot}{MCMC}(x,names,burnin=trunc(.1*nrow(X)),tvalues,TRACEPLOT=TRUE,DEN=TRUE,INT=TRUE,
      CHECK_NDRAWS=TRUE,... )
}

\arguments{
  \item{x         }{ A \code{MCMC} class matrix of posterior draws, such as \code{bayeslm\$beta}.}
  \item{names     }{ an optional character vector of names for the columns of \code{X}.}
  \item{burnin    }{ Number of draws to burn-in (default value is \eqn{0.1*nrow(X)}).}
  \item{tvalues   }{ vector of true values.}
  \item{TRACEPLOT }{ logical, TRUE provide sequence plots of draws and acfs (default: \code{TRUE})}
  \item{DEN       }{ logical, TRUE use density scale on histograms (default: \code{TRUE})}
  \item{INT       }{ logical, TRUE put various intervals and points on graph (default: \code{TRUE})}
  \item{CHECK_NDRAWS}{ logical, TRUE check that there are at least 100 draws (default: \code{TRUE})}
  \item{...       }{ optional arguments for generic function.}
}

\details{
This function is modified from package \code{bayesm} by Peter Rossi. It plots summary of posterior draws.
}

\author{Peter Rossi, Anderson School, UCLA, \email{perossichi@gmail.com}.}

\seealso{ \code{\link{summary.bayeslm.fit}}}

\examples{
x = matrix(rnorm(1000), 100, 10)
y = x \%*\% rnorm(10) + rnorm(100)
fit=bayeslm(y~x)
plot(fit$beta)
}

\keyword{summary}