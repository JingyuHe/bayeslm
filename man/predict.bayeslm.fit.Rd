\name{predict.bayeslm.fit}
\alias{predict.bayeslm.fit}

\title{Predict new data}

\description{
\code{predict.bayeslm.fit} is an S3 method to predict response on new data.
}

\usage{
\method{predict}{bayeslm.fit}(object,data,burnin,X,...)
}

\arguments{
  \item{object    }{ \code{object} is output of \code{bayeslm} function.}
  \item{data      }{ A data frame or list of new data to predict.}
  \item{burnin    }{ number of draws to burn-in (default value is \eqn{0.1*nrow(X)}).}
  \item{X         }{ If call \code{bayeslm} with matrices input \code{x} and {y} but not formula, pass new matrix to predict here. See example for details.}
  \item{...       }{ optional arguments for generic function.}
}

\details{
Make prediction on new data set, users are allowed to adjust number of burn-in samples.
}

\author{Jingyu He}

\examples{
x = matrix(rnorm(1000), 100, 10)
y = x \%*\% rnorm(10) + rnorm(100)
data = list(x = x, y = y)

# Train the model with formula input
fit1 = bayeslm(y ~ x, data = data)
# predict
pred1 = predict(fit1, data)

# Train the model with matrices input
fit2 = bayeslm(Y = y, X = x)
pred2 = predict(fit2, X = x)
}
\keyword{univar}