\name{QRIndiv}
\alias{QRIndiv}
\title{
Individual quantile objective function
}
\description{
The estimation of conditional quantile curves using individual quantile objective function.
}
\usage{
QRIndiv(VecX, tau, times, subj, X, y, d, kn, degree, lambda, gam)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{VecX}{
The representative values for each covariate used to estimate the desired conditional quantile curves.
}
  \item{tau}{
The quantiles of interest.
}
  \item{times}{
The vector of the time variable.
}
  \item{subj}{
The vector of subjects/individuals.
}
  \item{X}{
The covariate matrix containing 1 as its first column (including intercept in the model).
}
  \item{y}{
The response vector.
}
  \item{d}{
The order of the differencing operator for each covariate.
}
  \item{kn}{
The number of knots for each covariate.
}
  \item{degree}{
The degree of the B-spline basis function for each covariate.
}
  \item{lambda}{
The grid for the smoothing parameter to control the trade of between fidelity and penalty term (use a fine grid of lambda).
}
  \item{gam}{
The power used in estimating the smooting parameter for each covariate (e.g. gam=1 or gam=0.5).
}
}
%\details{
%The grid of lambda contributes on the wiggliness of the conditional quantile curves estimators. Hence, we suggest to use a finer grid of lambda.
%}
\value{
\item{alpha}{The estimator of the coefficient vector of the basis B-splines.}
\item{hat_bt}{The varying coefficients estimators.}
\item{qhat}{The conditional quantile curves estimator.}
}
\references{
%Andriyana, Y. (2015). P-splines quantile regression in varying coefficient models. \emph{Doctoral Dissertation}, KU Leuven, Department of Mathematics, January 2015.

Andriyana, Y., Gijbels, I., and Verhasselt, A. P-splines quantile
regression estimation in varying coefficient models. \emph{Test} 23, 1 (2014a),153--194.
}
\author{
Yudhie Andriyana
}
\note{
Some warning messages are related to the function \code{rq.fit.sfn}.
%(See \url{http://www.inside-r.org/packages/cran/quantreg/docs/sfnMessage}).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[quantreg]{rq.fit.sfn}}
\code{\link[SparseM]{as.matrix.csr}}
\code{\link{truncSP}}
}
\examples{
%library(truncSP)
data(PM10)

PM10 = PM10[order(PM10$day,PM10$hour,decreasing=FALSE),]

y = PM10$PM10[1:200]
times = PM10$hour[1:200]
subj = PM10$day[1:200]
dim = length(y)
x0 = rep(1,200)
x1 = PM10$cars[1:200]
x2 = PM10$wind.speed[1:200]

X = cbind(x0, x1, x2)

VecX = c(1, max(x1), max(x2))

##########################
#### Input parameters ####
##########################
kn = c(10, 10, 10)
degree = c(3, 3, 3)
taus = seq(0.1,0.9,0.1)
lambdas = c(1,1.5,2)
d = c(1, 1, 1)
gam = 1/2
##########################


qhat = QRIndiv(VecX=VecX, tau=taus, times=times, subj=subj, X=X,
		y=y, d=d, kn=kn, degree=degree, lambda=lambdas, gam=gam)$qhat

qhat1 = qhat[,1]
qhat2 = qhat[,2]
qhat3 = qhat[,3]
qhat4 = qhat[,4]
qhat5 = qhat[,5]
qhat6 = qhat[,6]
qhat7 = qhat[,7]
qhat8 = qhat[,8]
qhat9 = qhat[,9]


i = order(times, y, qhat1, qhat2, qhat3, qhat4, qhat5, qhat6, qhat7,
    qhat8, qhat9);

times = times[i]; y = y[i]; qhat1 = qhat1[i]; qhat2=qhat2[i];
qhat3=qhat3[i]; qhat4=qhat4[i]; qhat5=qhat5[i]; qhat6=qhat6[i];
qhat7=qhat7[i]; qhat8=qhat8[i]; qhat9=qhat9[i];

ylim = range(qhat1, qhat9)
plot(qhat1~times, col="magenta", cex=0.2, lty=5, lwd=2, type="l",
    ylim=ylim, xlab="hour", ylab="PM10");
lines(qhat2~times, col="aquamarine4", cex=0.2, lty=4, lwd=2);
lines(qhat3~times, col="blue", cex=0.2, lty=3, lwd=3);
lines(qhat4~times, col="brown", cex=0.2, lty=2, lwd=2);
lines(qhat5~times, col="black", cex=0.2, lty=1, lwd=2);
lines(qhat6~times, col="orange", cex=0.2, lty=2, lwd=2)
lines(qhat7~times, col="darkcyan", cex=0.2, lty=3, lwd=3);
lines(qhat8~times, col="green", cex=0.2, lty=4, lwd=2);
lines(qhat9~times, col="red", cex=0.2, lty=5, lwd=3)

legend("bottom", c(expression(tau==0.9), expression(tau==0.8),
      expression(tau==0.7), expression(tau==0.6), expression(tau==0.5),
      expression(tau==0.4), expression(tau==0.3), expression(tau==0.2),
      expression(tau==0.1)), ncol=1, col=c("red","green","darkcyan",
      "orange","black","brown","blue","aquamarine4","magenta"),
      lwd=c(2,2,3,2,2,2,3,2,2), lty=c(5,4,3,2,1,2,3,4,5))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
