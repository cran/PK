\name{biexp}
\alias{biexp}
\title{Two-Phase Half-Life Estimation by Biexponential Model}
\description{Estimation of initial and terminal half-life by fitting a biexponential model.}

\usage{biexp(conc, time, log.scale=FALSE, baseline=0, tol=1E-9, maxit=500)}

\arguments{
  \item{conc}{ Levels of concentrations. }
  \item{time}{ Time points of concentration assessment. }
  \item{log.scale}{ Logical value indicating whether fitting is performed on the observed or log-scale (default=\code{FALSE}). }
  \item{baseline}{ Pre-dosing value (default=\code{0}). }
  \item{tol}{ Relative error tolerance (default=\code{1E-9}). } 
  \item{maxit}{ Maximum number of iterations (default=\code{500}).}
}

\details{
Estimation of initial and terminal half-life using the biexponential \code{y=a1*exp(-b1*x)+a2*exp(-b2*x)} model with a parameterization to ensure b1 > b2 > 0 fitted by the least squares criteria with function \code{optim} of package \code{base} with \code{method} "Nelder-Mead". Curve peeling (Foss, 1969) is used get start values for nonlinear model fitting. When no adequate starting values are determined by curve peeling, a single exponential model is fitted with starting values obtained from an OLS regression on log transformed values with a parameterization to ensure a slope > 0. \cr \cr

If the baseline value indicating that the intrinsic level is greater than zero, the pre-dosing value is subtracted from all concentration levels before calculation. Resulting values of zero or below zero are omitted.
}
 
\value{
A list of S3 class \code{"halflife"} containing the following components: \cr
  \item{parms}{ half-life and model estimates.}
  \item{time}{ time points of concentration assessments. }
  \item{conc}{ levels of concentrations. } 
  \item{method}{ "biexp". }
}

\seealso{\code{\link{lee}}}

\note{Records including missing values and values below or equal to zero are omitted. }

\references{
Foss S. D. (1969). A Method for Obtaining Initial Estimates of the Parameters in Exponential Curve Fitting. \emph{Biometrics}, 25:580-584. \cr

Pinheiro J. C. and Bates D. M. (2000). \emph{Mixed-Effects Models in S and S-PLUS}. Springer, New York. \cr
}

\author{Martin J. Wolfsegger and Thomas Jaki}

\examples{
## example from Pinheiro J.C. and Bates D.M. (2000, page 279)
## dataset Indometh of package datasets
require(datasets)
time <- Indometh$time
conc <- Indometh$conc

# fitting on observed and log-scale
res1 <- biexp(conc=conc, time=time, log.scale=FALSE)
res2 <- biexp(conc=conc, time=time, log.scale=TRUE)

print(res1$parms)
print(res2$parms)

plot(res1, xlim=c(0, max(time)))
plot(res2, xlim=c(0, max(time)), add=TRUE, lty=2)
legend(x=5, y=2.5, lty=c(1,2), legend=c("fitted on observed scale", "fitted on log-scale"))

# get residuals using function nls with maxiter=1 and tol=Inf
parms.obs <- list(a1=res1$parms[3,1], b1=res1$parms[2,1], a2=res1$parms[3,2], b2=res1$parms[2,2])
parms.log <- list(a1=res2$parms[3,1], b1=res2$parms[2,1], a2=res2$parms[3,2], b2=res2$parms[2,2])

resid.obs <- nls(conc ~ a1*exp(-b1*time) + a2*exp(-b2*time), start=parms.obs, 
             control=nls.control(maxiter=1, tol=Inf))
resid.log <- nls(log(conc) ~ log(a1*exp(-b1*time) + a2*exp(-b2*time)), start=parms.log, 
             control=nls.control(maxiter=1, tol=Inf))

plot.new()
plot(x=c(0, max(time)), y=c(-1,1), type="n", xlab="Time", ylab="Residuals")
abline(h=0)
points(y=summary(resid.obs)$res, x=time, col="red")
points(y=summary(resid.log)$res, x=time, col="blue")
legend(x=0, y=1, pch=c(21,21), col=c("red","blue"), legend=c("fitted on observed scale", 
       "fitted on log-scale")) 

# see also function SSbiexp of package stats 
parms.ssbiexp <- as.list(coef(nls(conc ~ SSbiexp(time, a1, b1, a2, b2))))
parms.ssbiexp <- list(a1=parms.ssbiexp$a1, b1=exp(parms.ssbiexp$b1), a2=parms.ssbiexp$a2, 
                 b2=exp(parms.ssbiexp$b2))
print(unlist(parms.ssbiexp))
print(unlist(parms.obs))
print(unlist(parms.log))
}
\keyword{misc}
