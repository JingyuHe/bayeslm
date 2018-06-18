\name{hs_gibbs}
\alias{hs_gibbs}

\title{Gibbs sampler of horseshoe regression}

\description{
Standard Gibbs sampler of horseshoe regression.
}

\usage{
hs_gibbs(Y, X, nsamps, a, b, scale_sigma_prior)
}

\arguments{
  \item{Y}{Response of regression.}
  \item{X}{Matrix of regressors.}
  \item{nsamps}{Number of posterior samples.}
  \item{a}{Parameter of inverse Gamma prior on \eqn{\sigma}.}
  \item{b}{Parameter of inverse Gamma prior on \eqn{\sigma}.}
  \item{scale_sigma_prior}{Bool, if \code{TRUE}, use prior scaled by \eqn{\sigma}.}
}

\details{
This function implements standard Gibbs sampler of horseshoe regression. The prior is
\eqn{y \mid \beta, \sigma^2, X \sim MVN(X\beta, \sigma^2 I)}
\eqn{\beta_i \mid \tau, \lambda_i, \sigma \sim N(0, \lambda_i^2\tau^2\sigma^2)}
\eqn{\sigma^2\sim IG(a, b)}
\eqn{\tau \sim C^{+}(0,1)}
\eqn{\lambda_i \sim C^{+}(0,1)}
}

\author{Jingyu He}

\seealso{ \code{\link{summary.mcmc}}}

\examples{
x = matrix(rnorm(1000), 100, 10)
y = x \%*\% rnorm(10) + rnorm(100)
fit=hs_gibbs(y, x, 1000, 1, 1, TRUE)
summary(fit)
}

\keyword{univar}