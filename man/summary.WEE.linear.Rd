\name{summary.WEE.linear}
\alias{summary.WEE.linear}

\title{
Summary methods for WEE linear Regression
}
\description{
Returns a summary list for a WEE linear regression fit.
}
\usage{
\method{summary}{WEE.linear}(object, ...)
}
\arguments{
  \item{object}{
object produced by \code{\link{WEE.linear}}. 
}
  \item{...}{further arguments passed to or from other methods.}
}
\value{
a list is returned with the following components.
  \item{Coefficients }{a vector of coefficients}
  \item{StdErr }{bootstrap standard errors of the coefficients, only returned when boot > 0}
  \item{Chisq }{Chi-squared test statistics of the coefficients, only returned when boot > 0}
  \item{p.value }{p-values of the chi-squared test statistics, only returned when boot > 0}
  \item{Covariance}{the estimated covariance matrix of the coefficients in the model, provided that boot > 0 in the called sequence.}
}

\seealso{
\code{\link{WEE.linear}}
}
\examples{
\donttest{
## continued from WEE.linear
## summary of WEE linear object
summary(WEE.linear(y ~ x + z, D, 
        data = dat_sample, pd_pop = pD))                                       
summary(WEE.linear(y ~ x + z, D, 
        data = dat_sample, pd_pop = pD, boot=500))
}
}

