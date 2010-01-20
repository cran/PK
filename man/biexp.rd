\name{biexp}
\alias{biexp}
\title{Two-phase half-life estimation by biexponential model}
\description{Estimation of initial and terminal half-life by fitting a biexponential model.}

\usage{biexp(conc, time, log.scale=FALSE, baseline=0, tol=1E-9, maxit=500)}

\arguments{
  \item{conc}{ Levels of concentrations as a vector. }
  \item{time}{ Time points of concentration assessment as a vector. One time point for each concentration measured needs to be specified.}
  \item{log.scale}{ Logical value indicating whether fitting is performed on the observed or log-scale (default=\code{FALSE}). }
  \item{baseline}{ Pre-dosing value (default=\code{0}). }
  \item{tol}{ Relative error tolerance (default=\code{1E-9}). } 
  \item{maxit}{ Maximum number of iterations (default=\code{500}).}
}

\details{
Estimation of initial and terminal half-life using the biexponential \code{y=a1*exp(-b1*x)+a2*exp(-b2*x)} model with a parameterization to ensure b1 > b2 > 0 fitted by the least squares criteria with function \code{optim} of package \code{base} with \code{method} "Nelder-Mead". Curve peeling (Foss, 1969) is used get start values for nonlinear model fitting. When no adequate starting values are determined by curve peeling, a single exponential model is fitted with starting values obtained from an OLS regression on log transformed values with a parameterization to ensure a slope > 0. \cr\cr

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
Foss S. D. (1969). A Method for Obtaining Initial Estimates of the Parameters in Exponential Curve Fitting. \emph{Biometrics}, 25:580-584. \cr\cr

Pinheiro J. C. and Bates D. M. (2000). \emph{Mixed-Effects Models in S and S-PLUS}. Springer, New York. \cr\cr

Wolfsegger M. J. and Jaki T. (2009). Non-compartmental Estimation of Pharmacokinetic Parameters in Serial Sampling Designs. \emph{Journal of Pharmacokinetics and Pharmacodynamics}, 36(5):479-494. \cr	
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
parms.obs <- list(a1=res1$parms[3,1], b1=res1$parms[2,1], a2=res1$parms[3,2], 
                  b2=res1$parms[2,2])
parms.log <- list(a1=res2$parms[3,1], b1=res2$parms[2,1], a2=res2$parms[3,2], 
                  b2=res2$parms[2,2])

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


## comparing different non-linear fits
## example for a serial sampling data design from Wolfsegger and Jaki (2009)
conc <- c(0, 0, 0, 2.01, 2.85, 2.43, 0.85, 1.00, 0.91, 0.46, 0.35, 0.63, 0.39, 0.32, 
          0.45, 0.11, 0.18, 0.19, 0.08, 0.09, 0.06)
time <- c(rep(0,3), rep(5/60,3), rep(3,3), rep(6,3), rep(9,3), rep(16,3), rep(24,3))

res.biexp1 <- biexp(conc=conc, time=time, log=TRUE)
res.biexp2 <- biexp(conc=conc, time=time, log=FALSE)

print(res.biexp1$parms)
print(res.biexp2$parms)

split.screen(c(1,2)) 
screen(1)
plot(x=c(0,25), y=c(0,3), type='n', 
ylab='Plasma concentration (IU/mL)', xlab='Time (hours)')
points(x=time, y=conc, pch=21)
close.screen(1)
screen(2)
plot(x=c(0,25), y=c(0.01, 10), type='n', log='y', yaxt='n', 
ylab='Plasma concentration (IU/mL)', xlab='Time (hours)')
axis(side=2, at=c(0.01, 0.1, 1, 10), labels=c('0.01', '0.1', '1', '10'))
axis(side=2, at=seq(2,9,1), tcl=-0.25, labels=FALSE) 
axis(side=2, at=seq(0.2,0.9,0.1), tcl=-0.25, labels=FALSE) 
axis(side=2, at=seq(0.02,0.09,0.01), tcl=-0.25, labels=FALSE) 
points(x=time, y=conc, pch=21)
plot(res.biexp1, pch=NA, add=TRUE, lty=1)
plot(res.biexp2, pch=NA, add=TRUE, lty=2)
legend(x=25, y=10, xjust=1, col=c('black', 'black'), lty=c(1,2), 
        title='Nonlinear fitting:', 
        legend=c('logarithmic scale', 'observed scale'))
close.screen(2)

}
\keyword{misc}
