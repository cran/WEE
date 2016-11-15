\name{plot.predict.WEE.quantile}
\alias{plot.predict.WEE.quantile}
\title{
Plot predicted quantiles of WEE.quantile regression fit
}
\description{
Plot the predicted quantiles and their point-wise confidence intervals of a WEE-quantile fit on new dataset.
}
\usage{
\method{plot}{predict.WEE.quantile}(x, CI = FALSE, level = 0.95, index = 1, \dots)
}
\arguments{
  \item{x}{object produced by \code{\link{predict.WEE.quantile}}. }
  \item{CI}{logical flag indicating whether to plot confidence interval: default is FALSE; if TRUE the function not only plots point predictions for each of the 'newdata' points but also lower and upper confidence limits. Only set TRUE when boot > 0.}
  \item{level}{confidence level.}
  \item{index}{a vector to indicate the subset of newx to be plotted. Default is 1, i.e. the first combination of newx.}
  \item{\dots}{
further graphical parameters passed to \code{\link[graphics]{plot}}.
}
}

\seealso{
\code{\link{predict.WEE.quantile}}
}
\examples{
\donttest{
## continued from predict.WEE.quantile
## Plot prediction without confidence interval
plot(p1,index = c(2,3))

## Plot prediction with confidence interval
plot(p2, CI = TRUE)}
}