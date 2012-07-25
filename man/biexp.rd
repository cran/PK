\name{biexp}
\alias{biexp}
\title{Two-phase half-life estimation by biexponential model}
\description{Estimation of initial and terminal half-life by fitting a biexponential model.}

\usage{biexp(conc, time, log.scale=FALSE, tol=1E-9, maxit=500)}

\arguments{
  \item{conc}{ Levels of concentrations as a vector. }
  \item{time}{ Time points of concentration assessment as a vector. One time point for each concentration measured needs to be specified.}
  \item{log.scale}{ Logical value indicating whether fitting is performed on the observed or log-scale (default=\code{FALSE}). }
  \item{tol}{ Relative error tolerance (default=\code{1E-9}). } 
  \item{maxit}{ Maximum number of iterations (default=\code{500}).}
}

\details{
Estimation of initial and terminal half-life using the biexponential \code{y=a1*exp(-b1*x)+a2*exp(-b2*x)} model with a parameterization to ensure b1 > b2 > 0 fitted by the least squares criteria with function \code{optim} of package \code{base} with \code{method} "Nelder-Mead". Curve peeling (Foss, 1969) is used get start values for nonlinear model fitting. When no adequate starting values are 
determined by curve peeling, a single exponential model is fitted with starting values obtained from an OLS regression on log transformed values with a parameterization to ensure a slope > 0. \cr\cr

Fitting on the log-scale is based on the transform-both-sides approach described for example in chapter 4 of Bonate (2006) which is useful for some error distributions. An additional discussion regarding weighting schemes can be found in Gabrielsson and Weiner (2000, pages 368-374).
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
Bonate P. L. (2006). \emph{Pharmacokinetic-Pharmacodynamic Modeling and Simulation}. Springer, New York. \cr\cr

Gabrielsson J. and Weiner D. (2000). \emph{Pharmacokinetic and Pharmacodynamic Data Analysis: Concepts and Applications}. 4th Edition. Swedish Pharmaceutical Press, Stockholm. \cr\cr

Foss S. D. (1969). A Method for Obtaining Initial Estimates of the Parameters in Exponential Curve Fitting. \emph{Biometrics}, 25:580-584. \cr\cr

Pinheiro J. C. and Bates D. M. (2000). \emph{Mixed-Effects Models in S and S-PLUS}. Springer, New York. \cr\cr

Wolfsegger M. J. and Jaki T. (2009). Non-compartmental Estimation of Pharmacokinetic Parameters in Serial Sampling Designs. \emph{Journal of Pharmacokinetics and Pharmacodynamics}, 36(5):479-494. \cr	
}

\author{Martin J. Wolfsegger and Thomas Jaki}

\examples{
#### example from Pinheiro J.C. and Bates D.M. (2000, page 279)
#### dataset Indometh of package datasets
require(datasets)
data <- subset(Indometh, Subject==2) 
time <- data$time
conc <- data$conc 

## fitting on observed and log-scale
res.obs <- biexp(conc=conc, time=time, log.scale=FALSE)
res.log <- biexp(conc=conc, time=time, log.scale=TRUE)

print(res.obs$parms)
print(res.log$parms)

plot(res.obs, ylim=c(0,5), xlim=c(0, max(time)), las=1)
plot(res.log, ylim=c(0,5), xlim=c(0, max(time)), las=1, add=TRUE, lty=2)
legend(x=0, y=5, lty=c(1,2), legend=c("fitted on observed scale", "fitted on log-scale"))

## get residuals using function nls with tol=Inf
parms.obs <- list(a1=res.obs$parms[3,1], b1=res.obs$parms[2,1], a2=res.obs$parms[3,2], b2=res.obs$parms[2,2])
parms.log <- list(a1=res.log$parms[3,1], b1=res.log$parms[2,1], a2=res.log$parms[3,2], b2=res.log$parms[2,2])

mod.obs <- nls(conc ~ a1*exp(-b1*time) + a2*exp(-b2*time), start=parms.obs, control=nls.control(tol=Inf))
mod.log <- nls(conc ~ a1*exp(-b1*time) + a2*exp(-b2*time), start=parms.log, control=nls.control(tol=Inf))

## identical estimates to mod.log but different SEs  
summary(nls(log(conc)~log(a1*exp(-b1*time) + a2*exp(-b2*time)), start=parms.log, control=nls.control(tol=Inf)))

## different approach using weighted least squares (WLS) in nls 
mod.ols <- nls(conc ~ a1*exp(-b1*time) + a2*exp(-b2*time), start=parms.obs)
mod.wls1 <- nls(conc ~ a1*exp(-b1*time) + a2*exp(-b2*time), start=parms.obs, weight=1/predict(mod.ols)^1)
mod.wls2 <- nls(conc ~ a1*exp(-b1*time) + a2*exp(-b2*time), start=parms.obs, weight=1/predict(mod.ols)^2)

split.screen(c(2,2))
screen(1)
plot(ylim=c(-0.35,0.35), y=resid(mod.obs), x=predict(mod.obs), las=1, 
    main='Fitted using biexp on observed scale', xlab='Predicted', ylab='Residual')
abline(h=0)
screen(2)
plot(ylim=c(-0.35,0.35), y=resid(mod.log), x=predict(mod.log), las=1, 
     main='Fitted using biexp on log-scale', xlab='Predicted', ylab='Residual')
abline(h=0)
screen(3)
plot(ylim=c(-0.35,0.35), y=resid(mod.wls1), x=predict(mod.wls1), las=1, 
     main='Fitted using nls with weights 1/predict(mod.ols)^1', xlab='Predicted', ylab='Residual')
abline(h=0)
screen(4)
plot(ylim=c(-0.35,0.35), y=resid(mod.wls2), x=predict(mod.wls2), las=1, 
     main='Fitted using nls with weights 1/predict(mod.ols)^2', xlab='Predicted', ylab='Residual')
abline(h=0) 
close.screen(all.screens=TRUE)


#### example for a serial sampling data design from Wolfsegger and Jaki (2009)
conc <- c(2.01, 2.85, 2.43, 0.85, 1.00, 0.91, 0.46, 0.35, 0.63, 0.39, 0.32, 
          0.45, 0.11, 0.18, 0.19, 0.08, 0.09, 0.06)
time <- c(rep(5/60,3), rep(3,3), rep(6,3), rep(9,3), rep(16,3), rep(24,3))

res.biexp1 <- biexp(conc=conc, time=time, log=TRUE)
res.biexp2 <- biexp(conc=conc, time=time, log=FALSE)

print(res.biexp1$parms)
print(res.biexp2$parms)

split.screen(c(1,2)) 
screen(1)
plot(x=c(0,25), y=c(0,3), type='n', las=1, 
ylab='Plasma concentration (IU/mL)', xlab='Time (hours)')
points(x=time, y=conc, pch=21)
plot(res.biexp1, pch=NA, add=TRUE, lty=1)
plot(res.biexp2, pch=NA, add=TRUE, lty=2)
legend(x=25, y=3, xjust=1, col=c('black', 'black'), lty=c(1,2), 
        title='Nonlinear fitting with function biexp:', 
        legend=c('option: log=TRUE', 'option: log=FALSE'))
close.screen(1)
screen(2)
plot(x=c(0,25), y=c(0.01, 10), type='n', log='y', yaxt='n', 
ylab='Plasma concentration (IU/mL)', xlab='Time (hours)')
axis(side=2, at=c(0.01, 0.1, 1, 10), labels=c('0.01', '0.1', '1', '10'), las=1)
axis(side=2, at=seq(2,9,1), tcl=-0.25, labels=FALSE) 
axis(side=2, at=seq(0.2,0.9,0.1), tcl=-0.25, labels=FALSE) 
axis(side=2, at=seq(0.02,0.09,0.01), tcl=-0.25, labels=FALSE) 
points(x=time, y=conc, pch=21)
plot(res.biexp1, pch=NA, add=TRUE, lty=1)
plot(res.biexp2, pch=NA, add=TRUE, lty=2)
legend(x=25, y=10, xjust=1, col=c('black', 'black'), lty=c(1,2), 
        title='Nonlinear fitting with function biexp:', 
        legend=c('option: log=TRUE', 'option: log=FALSE'))
close.screen(all.screens=TRUE)


#### example from Gabrielsson and Weiner (2000, page 743) 
#### endogenous concentration is assumed to be constant over time  
dose <- 36630  
time <- c(-1, 0.167E-01, 0.1167, 0.1670, 0.25, 0.583, 0.8330, 1.083, 1.583, 2.083, 4.083, 8.083, 12, 23.5, 24.25, 26.75, 32)
conc <- c(20.34, 3683, 884.7, 481.1, 215.6, 114, 95.8, 87.89, 60.19, 60.17, 34.89, 20.99, 20.54, 19.28, 18.18, 19.39, 22.72)
data <- data.frame(conc,time)

## get starting values using function biexp using naive adjustment for endogenous concentration by subtraction of pre-value
data$concadj <- data$conc - data$conc[1]
data$concadj[min(which(data$concadj<0)):nrow(data)] <- NA
res.biexp <- biexp(conc=data$concadj[-1], time=data$time[-1])$parms 
start <- list(a1=res.biexp[3,1], k1=res.biexp[2,1], a2=res.biexp[3,2], k2=res.biexp[2,2])  

## specify indicator variable enabling inclusion of pre-dose concentration for fitting  
data$i1 <- ifelse(data$time <0, 1, 0)
data$i2 <- ifelse(data$time <0, 0, 1)

## assuming constant absolute error: ordinary least squares 
mod.ols <- nls(conc ~ i1*base + i2*(base + a1*exp(-k1*time) + a2*exp(-k2*time)), 
               start=c(base=20.34, start), data=data, trace=TRUE)
 
## assuming constant relative error (i.e. proportional error - weight of 2): weighted least squares  
mod.wls <- nls(conc ~ i1*base + i2*(base + a1*exp(-k1*time) + a2*exp(-k2*time)), 
               start=c(base=20.34, start), data=data, weight=1/predict(mod.ols)^2, trace=TRUE) 

## assuming constant relative error (i.e. proportional error - weight of 2): iteratively re-weighted least squares 
mod.irwls <- mod.wls 
for(i in 1:10){ 
   print(as.vector(coef(mod.irwls))) 
   mod.irwls <- nls(conc ~ i1*base + i2*(base + a1*exp(-k1*time) + a2*exp(-k2*time)), 
                    start=c(base=20.34, start), data=data, weight=1/predict(mod.irwls)^2)
}
 
summary(mod.ols)
summary(mod.wls)
summary(mod.irwls)

newdata <- data.frame(time=seq(0,32,0.01))
newdata$i1 <- ifelse(newdata$time <0, 1, 0)
newdata$i2 <- ifelse(newdata$time <0, 0, 1)
plot(conc ~ time, data=data, ylim=c(10,1E4), log='y', yaxt='n', 
     xlab='Time (hours)', ylab='Log of concentration (pmol/L)')
axis(side=2, at=c(10, 100, 1000, 10000), las=1)
axis(side=2, at=seq(1E1,1E2,1E1), tcl=-0.25, labels=FALSE) 
axis(side=2, at=seq(1E2,1E3,1E2), tcl=-0.25, labels=FALSE) 
axis(side=2, at=seq(1E3,1E4,1E3), tcl=-0.25, labels=FALSE) 
points(x=newdata$time, y=predict(mod.irwls, newdata), type='l')

## get total clearance (cls), inter-compartmental distribution (cld), and volume of distributions 
## from macro constant parametrization according to Gabrielsson and Weiner (2000) 
parm <- coef(mod.wls)[-1]
## get micro constants 
k21 <- as.real((parm[1]*parm[4]+parm[3]*parm[2])/(parm[1]+parm[3]))
k10 <- as.real(parm[2]*parm[4] / k21)
k12 <- as.real(parm[2]+parm[4] - k21 - k10)
## get cls, cld, vc, and vt 
cls <- as.real(dose / (parm[1]/parm[2] + parm[3]/parm[4])) 
vc <- as.real(dose / (parm[1] + parm[2]))
cld <- k12*vc 
vt <- cld / k21 
print(c(cls, cld, vc, vt))

## turnover model to account for endogenous baseline according to Gabrielsson and Weiner 
## using a biexponential (i.e. two-compartment) model parametrized in terms of clearance 

\dontrun{require(rgenoud)
require(deSolve)

k <- 2 # assuming proportional error - weighting in function objfun 
tinf <- 1/60 # duration of bolus in hours 
data <- subset(data, time>0) 

defun <- function(time, y, parms) { 
	rte1 <- ifelse(time <= tinf, dose/tinf, 0)
	dCptdt1 <- (rte1 + parms["synt"] - parms["cls"]*y[1] - parms["cld"]*y[1] + parms["cld"]*y[2]) / parms["vc"]
	dCptdt2 <- (parms["cld"]*y[1] - parms["cld"]*y[2])/parms["vt"] 
	list(c(dCptdt1, dCptdt2))
}

modfun <- function(time, synt, cls, cld, vc, vt) {
	out <- lsoda(y=c(synt/cls, synt/cls), times=c(0, data$time), defun, 
                parms=c(synt=synt, cls=cls, cld=cld, vc=vc, vt=vt), rtol=1e-5, atol=1e-5)[-1,2]
} 

objfun <- function(par) {
        out <- modfun(data$time, par[1], par[2], par[3], par[4], par[5])
        gift <- which(data$conc != 0 )
        sum((data$conc[gift]-out[gift])^2 / data$conc[gift]^k)		
}        

## grid search to get starting values for Nelder-Mead 
## increase values of pop.size and max.generation to get better starting values  
## values of 10 are used for illustration purpose only 
options(warn = -1) # omit warning when hard maximum limit is hit
gen <- genoud(objfun, nvars=5, max=FALSE, pop.size=10, max.generation=10,
             starting.value=c(1500, cls, cld, vc, vt), BFGS=FALSE,
             print.level=1, boundary.enforcement=2,
             Domains=matrix(c(0,0,0,0,0,1E4,1E3,1E3,1E3,1E3),5,2),
             MemoryMatrix=TRUE)
options(warn = 0) # set back to default  

opt <- optim(gen$par, objfun, method="Nelder-Mead") 

trn.wls <- nls(conc ~ modfun(time, synt, cls, cld, vc, vt), data=data, 
               start=list(synt=opt$par[1], cls=opt$par[2], cld=opt$par[3], vc=opt$par[4], vt=opt$par[5]), 
               trace=TRUE, nls.control(tol=Inf))

summary(trn.wls) 
}
}
\keyword{misc}
