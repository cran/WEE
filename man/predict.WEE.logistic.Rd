\name{predict.WEE.logistic}
\alias{predict.WEE.logistic}

\title{
WEE logistic Regression Prediction
}
\description{
Prediction on new dataset based on model fitted by WEE logistic regression
}
\usage{
\method{predict}{WEE.logistic}(object,newx, ...)
}
\arguments{
  \item{object}{
Object produced by \code{\link{WEE.logistic}}.
}
  \item{newx}{
A data matrix in which to look for variables with which to predict, newx cannot be omitted
}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
Produces predicted values, obtained by evaluating the WEE logistic regression function on newx.
}
\value{
If in the WEE.logistic fitting procedure boot = 0, linear predictor and predicted response of each newx are given.
If in the WEE.logistic fitting procedure boot > 0, standard errors of linear predictor and predicted response are given.
}
\seealso{
\code{\link{WEE.logistic}}
}
\examples{
\donttest{
## continued from WEE.logistic
## predict outcome y based on newx
newx = dat[sample(500000,3, replace=FALSE),][,c("x","z")]
predict(WEE.logistic(y ~ x + z, D, 
        data = dat_sample, pD),newx) 				      
predict(WEE.logistic(y ~ x + z, D, 
        data = dat_sample, pD, boot = 500),newx)
}
}

