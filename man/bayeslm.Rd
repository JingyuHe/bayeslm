\name{bayeslm}
\alias{bayeslm}
\alias{bayeslm.default}
\alias{bayeslm.formula}
\title{Efficient sampling for Gaussian linear model with arbitrary priors
}
\description{
  This package implements an efficient sampler for Gaussian Bayesian linear regression. The package uses elliptical slice sampler instead of regular Gibbs sampler. The function has several built-in priors and user can also provide their own prior function (written as a R function).
}
\usage{
\method{bayeslm}{default}(Y, X = FALSE, prior = "horseshoe", penalize = NULL, 
block_vec = NULL, sigma = NULL, s2 = 1, kap2 = 1, N = 20000L, burnin = 0L, 
thinning = 1L, vglobal = 1, sampling_vglobal = TRUE, verb = FALSE, icept = TRUE, 
standardize = TRUE, singular = FALSE, scale_sigma_prior = TRUE, prior_mean = NULL, 
prob_vec = NULL, cc = NULL, lambda = NULL, ...)

\method{bayeslm}{formula}(formula, data = list(), Y = FALSE, X = FALSE, 
prior = "horseshoe", penalize = NULL, block_vec = NULL, sigma = NULL, 
s2 = 1, kap2 = 1, N = 20000L, burnin = 0L, thinning = 1L, vglobal = 1, 
sampling_vglobal = TRUE, verb = FALSE, standardize = TRUE, singular = FALSE, 
scale_sigma_prior = TRUE, prior_mean = NULL, 
prob_vec = NULL, cc = NULL, lambda = NULL, ...)
}




\arguments{
  \item{formula}{\code{formula} of the model to fit.}
  \item{data}{an optional data frame containing the variables in the model.
    By default the variables are taken from the environment which
    \code{bayeslm} is called from.}
  \item{Y}{\code{data.frame}, \code{matrix}, or \code{vector} of inputs \code{Y}. Response variable. }
  \item{X}{\code{data.frame}, \code{matrix}, or \code{vector} of inputs \code{X}. Regressors. }
  \item{prior}{Indicating shrinkage prior to use. \code{"horseshoe"} for approximate horseshoe prior (default), \code{"laplace"} for laplace prior, \code{"ridge"} for ridge prior, \code{"sharkfin"} for "sharkfin" prior and \code{"nonlocal"} for nonlocal prior.}
  \item{block_vec}{A vector indicating number of regressors in each block. Sum of all entries should be the same as number of regressors. The default value is \code{block_vec = rep(1, p)}, put every regressor in its own block (slice-within-Gibbs sampler)}
  \item{penalize}{A vector indicating shrink regressors or not. It's length should be the same as number of regressors. \code{1} indicates shrink corresponding coefficient, \code{0} indicates no shrinkage. The default value is \code{rep(1, p)}, shrink all coefficients}
  \item{sigma}{Initial value of residual standard error. The default value is half of standard error of \code{Y}.}
  \item{s2, kap2}{Parameter of prior over sigma, an inverse gamma prior with rate s2 and shape s2.}
  \item{N}{Number of posterior samples (after burn-in).}
  \item{burnin}{Number of burn-in samples. If burnin > 0, the function will draw N + burnin samples and return the last N samples only.}
  \item{thinning}{Number of thinnings. \code{thinning = 1} means no thinning.}
  \item{vglobal}{Initial value of global shrinkage parameter. Default value is 1}
  \item{sampling_vglobal}{\code{Bool}, if \code{TRUE}, sampling the global shrinkage parameter by random walk Metropolis Hastings on log scale, otherwise always stay at the initial value \code{vglobal}.}
  \item{verb}{Bool, if \code{TRUE}, print out sampling progress.}
  \item{icept}{Bool, if the inputs are matrix \code{X} and \code{Y}, and \code{icept = TRUE}, the function will estimate intercept. Default value is \code{TRUE}. If the input is formula \code{Y ~ X}, option \code{icept} is useless, control intercept by formular \code{Y ~ X} or \code{Y ~ X - 1}.}
  \item{standardize}{Bool, if \code{TRUE}, standardize X and Y before sampling.}
  \item{singular}{Bool, if \code{TRUE}, take it as a rank-deficient case such as n < p or X'X is singular. See section 2.3.2 of the paper for details.}
  \item{scale_sigma_prior}{Bool, if \code{TRUE}, the prior of regression coefficient \eqn{\beta} is scaled by residual standard error \eqn{\sigma}.}
  \item{prior_mean}{\code{vector}, specify prior mean of nonlocal prior for each regressor. It should have length \code{p} (no intercept) or \code{p + 1} (intercept). The default value is 1.5 for all regressors.}
  \item{prob_vec}{\code{vector}, specify prior mean of sharkfin prior for each regressor. It should have length \code{p} (no intercept) or \code{p + 1} (intercept). The default value is 0.25 for all regressors.}
  \item{cc}{Only works when \code{singular == TRUE}, precision parameter of ridge adjustment. It should be a vector with length $p$. If it is \code{NULL}, it will be set as \code{rep(10, p)}.}
  \item{lambda}{The shrinkage parameter for Laplace prior only.}
  \item{...}{optional parameters to be passed to the low level function \code{bayeslm.default}.}
}
\details{
  For details of the approach, please see Hahn, He and Lopes (2017)
}
\value{
  \item{loops}{A \code{vector} of number of elliptical slice sampler loops for each posterior sample.}
  \item{sigma}{A \code{vector} of posterior samples of residual standard error.}
  \item{vglobal}{A \code{vector} of posterior samples of the global shrinkage parameter.}
  \item{beta}{A \code{matrix} of posterior samples of coefficients.}
  \item{fitted.values}{Fitted values of the regression model. Take posterior mean of coefficients with 20\% burnin samples.}
  \item{residuals}{Residuals of the regression model, equals \code{y - fitted.values}.}
}
\note{
  \code{horseshoe} is essentially call function \code{bayeslm} with \code{prior = "horseshoe"}. Same for \code{sharkfin}, \code{ridge}, \code{blasso}, \code{nonlocal}.
}

\references{
  Hahn, P. Richard, Jingyu He, and Hedibert Lopes. \emph{Efficient sampling for Gaussian linear regression with arbitrary priors.} (2017).
}

\author{ Jingyu He \email{jingyu.he@chicagobooth.edu} }


\examples{

p = 20
n = 100

kappa = 1.25
beta_true = c(c(1,2,3),rnorm(p-3,0,0.01))
sig_true = kappa*sqrt(sum(beta_true^2))


x = matrix(rnorm(p*n),n,p)
y = x \%*\% beta_true + sig_true * rnorm(n)


x = as.matrix(x)
y = as.matrix(y)
data = data.frame(x = x, y = y)

block_vec = rep(1, p) # slice-within-Gibbs sampler, put every coefficient in its own block

fitOLS = lm(y~x-1)

# call the function using formulas
fita = bayeslm(y ~ x, prior = 'horseshoe', 
        block_vec = block_vec, N = 10000, burnin = 2000)
# summary the results
summary(fita)
summary(fita$beta)

\donttest{
# put the first two coefficients in one elliptical sampling block
block_vec2 = c(2, rep(1, p-2))
fitb = bayeslm(y ~ x, data = data, prior = 'horseshoe', 
        block_vec = block_vec2, N = 10000, burnin = 2000)

# comparing several different priors

fit1 = bayeslm(y,x,prior = 'horseshoe', icept = FALSE, 
          block_vec = block_vec, N = 10000, burnin=2000)
beta_est1 = colMeans(fit1$beta)

fit2 = bayeslm(y,x,prior = 'laplace', icept = FALSE, 
          block_vec = block_vec, N = 10000, burnin=2000)
beta_est2 = colMeans(fit2$beta)

fit3 = bayeslm(y,x,prior = 'ridge', icept = FALSE, 
          block_vec = block_vec, N = 10000, burnin=2000)
beta_est3 = colMeans(fit3$beta)

fit4 = bayeslm(y,x,prior = 'sharkfin', icept = FALSE, 
          block_vec = block_vec, N = 10000, burnin=2000)
beta_est4 = colMeans(fit4$beta)

fit5 = bayeslm(y,x,prior = 'nonlocal', icept = FALSE, 
          block_vec = block_vec, N = 10000, burnin=2000)
beta_est5 = colMeans(fit5$beta)

plot(NULL,xlim=range(beta_true),ylim=range(beta_true), 
  xlab = "beta true", ylab = "estimation")
points(beta_true,beta_est1,pch=20)
points(beta_true,fitOLS$coef,col='red')
points(beta_true,beta_est2,pch=20,col='cyan')
points(beta_true,beta_est3,pch=20,col='orange')
points(beta_true,beta_est4,pch=20,col='pink')
points(beta_true,beta_est5,pch=20,col='lightgreen')

legend("topleft", c("OLS", "horseshoe", "laplace", "ridge", "sharkfin", 
  "nonlocal"), col = c("red", "black", "cyan", "orange", 
    "pink", "lightgreen"), pch = rep(1, 6))

abline(0,1,col='red')

rmseOLS = sqrt(sum((fitOLS$coef-beta_true)^2))
rmse1 = sqrt(sum((beta_est1-beta_true)^2))
rmse2 = sqrt(sum((beta_est2-beta_true)^2))
rmse3 = sqrt(sum((beta_est3-beta_true)^2))
rmse4 = sqrt(sum((beta_est4-beta_true)^2))
rmse5 = sqrt(sum((beta_est5-beta_true)^2))

print(cbind(ols = rmseOLS, hs = rmse1,laplace = rmse2,
ridge = rmse3,sharkfin = rmse4,nonlocal = rmse5))


}
}

\keyword{  linear regression }