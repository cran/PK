\name{lee}
\alias{lee}
\title{
Two-Phase Half-Life Estimation by Linear Fitting
}

\description{
Estimation of inital and terminal half-life by two-phase linear regression fitting.
}

\usage{
lee(time, conc, points=3, prev=0, method=c("lad", "ols", "hub", "npr")) 	     
}

\arguments{
  \item{time}{ time points of concentration assessments. }
  \item{conc}{ levels of concentrations. }
  \item{points}{ minimum data points in the terminal phase. }
  \item{prev}{ pre-dosing value. }
  \item{method}{ method of model fitting. } 
}

\details{
Estimation of inital and terminal half-life based on the method of Lee et al. (1990). This method uses a two-phase linear regression approach to break down the model into two straight lines based on the selection of the log10 transformed concentration values. For two-phase models the initial and terminal half-life were determined from the slopes of the regression lines. If a single-phase model is selected by this method, the half-life so determined is utilized as both initial and terminal phase half-life. Half-life is determined only for decreasing inital and terminal phases. \cr \cr

The method \code{ols} uses the ordinary least squares regression (OLS) to fit regression lines. \cr \cr

The method \code{lad} uses the absolute deviation regression (LAD) to fit regression lines by using the algorithm as described in Birkes and Dodge (chapter 4, 1993) for calculation of regression estimates.  \cr \cr

The method \code{hub} uses the Huber M regression to fit regression lines. Huber M-estimates are calculated by non-linear estimation by function \code{optim}, where OLS regression parameters are used as start values. The function that is minimized involved k = 1.5*1.483*MAD, where MAD is defined as the median of absolute deviation of residuals obtained by a least absolute deviation (LAD) regression based on the observed data. The initial value of MAD is used and not updated during iterations (Holland and Welsch, 1977). \cr \cr

The method \code{npr} uses the nonparametric regression to fit regression lines by using the algorithm as described in Birkes and Dodge (chapter 6, 1993) for calculation of regression estimates. \cr \cr

The selection criteria for the best tuple of regression lines is the sum of squared residuals for the \code{ols} method, the sum of Huber M residuals for the \code{hub} method, the sum of absolute residuals for the \code{lad} method and the sum of a function on ranked residuals suggest by Birkes and Dodge (page 115, 1993) for the \code{npr} method. \cr \cr 

If the pre-dosing value indicating the intrinsinc level is greater than 0, the pre-dosing value is subtracted from all concentration levels before calculation of inital and terminal half-life. 

}

\value{
  \item{parms}{ half-life and model estimates.}
  \item{chgpt}{ changepoint between inital and terminal phase. }
  \item{time}{ time points of concentration assessments. }
  \item{conc}{ levels of concentrations. } 
  \item{method}{ "lee". }
}

\note{Records including missing values and values below or equal to zero are omitted. }

\references{
Birkes D. and Dodge Y. (1993). Alternative Methods of Regression. Wiley, New York, Chichester, Brisbane, Toronto, Singapore.  \cr \cr
Holland P. W. and Welsch R. E. (1977). Robust regression using iteratively reweighted least-squares. Commun. Statist.-Theor. Meth. A6(9):813-827. \cr \cr
Lee M. L., Poon Wai-Yin, Kingdon H. S. (1990). A two-phase linear regression model for biologic half-life data. Journal of Laboratory and Clinical Medicine. 115(6):745-748. \cr \cr
}

\author{Martin Wolfsegger}

\examples{

## example for preparation 1 from Lee et. al (1990)
time <- c(0.5, 1.0, 4.0, 8.0, 12.0, 24.0)
conc <- c(75, 72, 61, 54, 36, 6)
result1 <- lee(conc=conc, time=time, method='ols', points=2)
print(result1$parms)
plot(result1)

## example for preparation 2 from Lee et. al (1990)
time <- c(0.5, 1.0, 2.0, 6.5, 8.0, 12.5, 24.0)
conc <- c(75, 55, 48, 51, 39, 9, 5)
result2 <- lee(conc=conc, time=time, method='ols', points=2)
print(result2$parms)
plot(result2)

## advanced plots 
xlim <- c(0,30)
ylim <- c(1,80)
ylab <- 'Log Concentration'
text1 <- paste('Initial half-life:', round(result1$parms[1,1],3), '   Terminal half-life:', round(result1$parms[1,2],3))
text2 <- paste('Initial half-life:', round(result2$parms[1,1],3), '   Terminal half-life:', round(result2$parms[1,2],3))
split.screen(figs=c(2,1)) 
screen(1)
plot(result1, ylab=ylab, main='Half-life: Preparation 1', xlim=xlim, ylim=ylim, log='y', sub=text1)
screen(2)
plot(result2, ylab=ylab, main='Half-life: Preparation 2', xlim=xlim, ylim=ylim, log='y', sub=text2)
close.screen(all=TRUE)
}

\keyword{misc}
