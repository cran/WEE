\name{predict.WEE.linear}
\alias{predict.WEE.linear}
\title{
WEE Linear Regression Prediction
}
\description{
Prediction on new dataset based on model fitted by WEE linear regression
}
\usage{
\method{predict}{WEE.linear}(object,newx, ...)
}
\arguments{
  \item{object}{
Object produced by \code{\link{WEE.linear}}.
}
  \item{newx}{
A data matrix in which to look for variables with which to predict, newx cannot be omitted.
}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
Produces predicted values, obtained by evaluating the WEE linear regression function on newx.
}
\value{
If in the WEE.linear fitting procedure boot = 0, only point predictions are provided here.
If in the WEE.linear fitting procedurep boot > 0, standard errors of prediction are also provided.
}

\seealso{
\code{\link{WEE.linear}}
}
\examples{
\donttest{
## continued from WEE.linear
## predict outcome y based on newx
newx = dat[sample(500000,3, replace=F),][,c("x","z")]
predict(WEE.linear(y ~ x + z, D, 
data = dat_sample, pD),newx)                                        
predict(WEE.linear(y ~ x + z, D, 
data = dat_sample, pD, boot = 500),newx)
}
}
