\name{summary.WEE.quantile}
\alias{summary.WEE.quantile}

\title{
Summary methods for WEE Quantile Regression
}
\description{
Returns a summary list for a WEE quantile regression fit.
}
\usage{
\method{summary}{WEE.quantile}(object, ...)
}
\arguments{
  \item{object}{
object produced by \code{\link{WEE.quantile}}.
}
  \item{...}{further arguments passed to or from other methods.}
}
\value{
a list is returned with the following components.
  \item{Coefficients }{a vector of coefficients}
  \item{StdErr }{bootstrap standard errors of the coefficients, only returned when boot > 0}
  \item{Wald }{Wald test statistics of the coefficients, only returned when boot > 0}
  \item{p.value }{p-values of the Wald test statistics, only returned when boot > 0}
  \item{Covariance}{the estimated covariance matrix for the coefficients in the model, provided that boot > 0 in the called sequence}
}

\seealso{
\code{\link{WEE.quantile}}
}

\examples{
\donttest{
## continued from WEE.quantile
summary(WEE.quantile(y ~ x, D, tau = 0.5, 
                     data = dat_quantile, pd_pop = pd))
summary(WEE.quantile(y ~ x + z, D, tau = c(0.25,0.5), 
                     data = dat_quantile, pd_pop = pd, boot=500))
}
}
