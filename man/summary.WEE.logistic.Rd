\name{summary.WEE.logistic}
\alias{summary.WEE.logistic}

\title{
Summary methods for WEE logistic Regression
}
\description{
Returns a summary list for a WEE logistic regression fit.
}
\usage{
\method{summary}{WEE.logistic}(object, ...)
}
\arguments{
  \item{object}{
object produced by \code{\link{WEE.logistic}}.
}
  \item{...}{further arguments passed to or from other methods.}
}
\value{
a list is returned with the following components.
  \item{Coefficients }{a vector of coefficients}
  \item{Oddsratio }{the exponentiated coefficients, namely the odds ratio associated with the corresponding covariate}
  \item{StdErr }{bootstrap standard errors of the coefficients, only returned when boot > 0}
  \item{Wald }{Wald test statistics of the coefficients, only returned when boot > 0}
  \item{p.value }{p-values of the Wald test statistics, only returned when boot > 0}
  \item{Covariance}{the estimated covariance matrix for the coefficients in the model, provided that boot > 0 in the called sequence}
}

\seealso{
\code{\link{WEE.logistic}}
}

\examples{
\donttest{
## continued from WEE.logistic
## summary of WEE logistic object
summary(WEE.logistic(y ~ x + z, D, 
        data = dat_sample, pd_pop = pD)) 				      
summary(WEE.logistic(y ~ x + z, D, 
        data = dat_sample, pd_pop = pD, boot=500))
}
}


