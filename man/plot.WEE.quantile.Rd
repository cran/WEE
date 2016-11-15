\name{plot.WEE.quantile}
\alias{plot.WEE.quantile}
\title{
Plot coefficients estimated from WEE.quantile
}
\description{
Plot the estimated quantile coefficients and their pointwise confidence intervals from WEE.quantile regression
}
\usage{
\method{plot}{WEE.quantile}(x, CI = FALSE, level = 0.95, index = 1, \dots)
}

\arguments{
  \item{x}{
object produced by \code{\link{WEE.quantile}}.
}
  \item{CI}{
logical flag indicating whether to plot confidence interval: default is FALSE; if TRUE the pointwise confidence interval is plotted. Only set TRUE when boots  > 0 in the WEE.quantile fitting process.
}
  \item{level}{
confidence level.
}
  \item{index}{a vector to indicate the subset of coefficients to be plotted (e.g., 2 indicates the coefficient of the first covariate, 3 indicates the coefficient of the second covaraite). Default is 1, i.e. the intercept.}
  \item{\dots}{
further graphical parameters passed to \code{\link[graphics]{plot}}.
}
}
\seealso{
\code{\link{WEE.quantile}}
}
\examples{
\donttest{
## continued from WEE.quantile
## plot fitted model without pointwise confidence interval
plot(WEE.quantile(y ~ x, D, tau = 0.5, 
                  data = dat_quantile, pd_pop = pd),index = c(2,3)) 	

## plot fitted model with pointwise confidence interval
plot(WEE.quantile(y ~ x + z, D, tau = c(0.25,0.5), 
                  data = dat_quantile, pd_pop = pd, boot = 500), 
                  CI = TRUE)
}
}

