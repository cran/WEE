\name{predict.WEE.quantile}
\alias{predict.WEE.quantile}
\title{
WEE quantile Regression Prediction
}
\description{
Prediction on new dataset based on model fitted by WEE quantile regression
}
\usage{
\method{predict}{WEE.quantile}(object,newx, ...)
}
\arguments{
  \item{object}{
Object produced by \code{\link{WEE.quantile}}.
}
  \item{newx}{
A new data matrix in which to look for data with which to predict, newx cannot be omitted.
}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
Produces predicted values, obtained by evaluating the WEE quantile regression function on newx.
}
\value{
If in the WEE.quantile fitting procedure boot = 0, only point predictions are given.
If in the WEE.quantile fitting procedure boot > 0, standard errors of prediction are also given.
}

\seealso{
\code{\link{WEE.quantile}}
}
\examples{
\donttest{
## continued from WEE.quantile
## prediction based on newx
newx = dat[sample(500000,3, replace=F),][,c("x")]
p1 = predict(WEE.quantile(y ~ x, D, tau = 0.5, 
     data = dat_quantile, pd_pop = pd),newx)	
p1

newx = dat[sample(500000,3, replace=F),][,c("x","z")]
p2 = predict(WEE.quantile(y ~ x + z, D, tau = c(0.25,0.5), 
     data = dat_quantile, pd_pop = pd, boot = 500),newx)
p2
}
}

