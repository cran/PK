\name{lee}
\alias{lee}
\title{Two-Phase Half-Life Estimation by Linear Fitting}
\description{Estimation of initial and terminal half-life by two-phase linear regression fitting.}
\usage{
lee(conc, time, points=3, baseline=0, method=c("lad", "ols", "hub", "npr"), lt=TRUE) 	     
}

\arguments{
  \item{conc}{ Levels of concentrations. }
  \item{time}{ Time points of concentration assessment. }
  \item{points}{ Minimum number of data points in the terminal phase (default=\code{3}). }
  \item{baseline}{ Pre-dosing value (default=\code{0}). }
  \item{method}{ Method of model fitting (default=\code{lad}). } 
  \item{lt}{ Logical value indicating whether requesting a longer terminal than initial half-life (default=\code{TRUE}).} 
}

\details{
Estimation of initial and terminal half-life based on the method of Lee \emph{et al.} (1990). This method uses a two-phase linear regression approach separate the model into two straight lines based on the selection of the log10 transformed concentration values. For two-phase models the initial and terminal half-life were determined from the slopes of the regression lines. If a single-phase model is selected by this method, the corresponding half-life is utilized as both initial and terminal phase half-life. Half-life is determined only for decreasing initial and terminal phases. \cr \cr

The method \code{ols} uses the ordinary least squares regression (OLS) to fit regression lines. \cr \cr

The method \code{lad} uses the absolute deviation regression (LAD) to fit regression lines by using the algorithm as described in Birkes and Dodge (chapter 4, 1993) for calculation of regression estimates.  \cr \cr

The method \code{hub} uses the Huber M regression to fit regression lines. Huber M-estimates are calculated by non-linear estimation using the function \code{optim}, where OLS regression parameters are used as starting values. The function that is minimized involved k = 1.5*1.483*MAD, where MAD is defined as the median of absolute deviation of residuals obtained by a least absolute deviation (LAD) regression based on the observed data. The initial value of MAD is used and not updated during iterations (Holland and Welsch, 1977). \cr \cr

The method \code{npr} uses the nonparametric regression to fit regression lines by using the algorithm as described in Birkes and Dodge (chapter 6, 1993) for calculation of regression estimates. \cr \cr

The selection criteria for the best tuple of regression lines is the sum of squared residuals for the \code{ols} method, the sum of Huber M residuals for the \code{hub} method, the sum of absolute residuals for the \code{lad} method and the sum of a function on ranked residuals for the \code{npr} method (see Birkes and Dodge (page 115, 1993)). Calculation details can be found in Wolfsegger (2006). \cr \cr

When \code{lt=TRUE}, the best two-phase model where terminal half-life >= initial half-life >= 0 is selected. When \code{lt=FALSE}, the best two-phase model among all possible tuples of regression is selected which can result in longer initial half-life than terminal half-life and/or in half-lifes < 0. \cr \cr

If the baseline value indicating that the intrinsic level is greater than zero, the pre-dosing value is subtracted from all concentration levels before calculation. Resulting values of zero or below zero are omitted.
}

\value{
A list of S3 class \code{"halflife"} containing the following components: \cr
  \item{parms}{ half-life and model estimates.}
  \item{chgpt}{ change point between initial and terminal phase. }
  \item{time}{ time points of concentration assessments. }
  \item{conc}{ levels of concentrations. } 
  \item{method}{ "lee". }
}

\note{Records including missing values and concentration values below or equal to zero are omitted. }

\references{
Birkes D. and Dodge Y. (1993). \emph{Alternative Methods of Regression}. Wiley, New York, Chichester, Brisbane, Toronto, Singapore.  \cr 

Holland P. W. and Welsch R. E. (1977). Robust regression using iteratively reweighted least-squares. \emph{Commun. Statist.-Theor. Meth.} A6(9):813-827. \cr

Lee M. L., Poon Wai-Yin, Kingdon H. S. (1990). A two-phase linear regression model for biologic half-life data. \emph{Journal of Laboratory and Clinical Medicine.} 115(6):745-748. \cr

Wolfsegger M. J. (2006). The R Package PK for Basic Pharmacokinetics. \emph{Biometrie und Medizin}, 5:61-68. \cr
}

\author{Martin J. Wolfsegger}

\examples{
## example for preparation 1 from Lee et al. (1990)
time <- c(0.5, 1.0, 4.0, 8.0, 12.0, 24.0)
conc <- c(75, 72, 61, 54, 36, 6)
res1 <- lee(conc=conc, time=time, method='ols', points=2, lt=TRUE)
res2 <- lee(conc=conc, time=time, method='ols', points=2, lt=FALSE)
plot(res1, log='y', ylim=c(1,100))
plot(res2, add=TRUE, lty=2)

## example for preparation 2 from Lee et al. (1990)
time <- c(0.5, 1.0, 2.0, 6.5, 8.0, 12.5, 24.0)
conc <- c(75, 55, 48, 51, 39, 9, 5)
res3 <- lee(conc=conc, time=time, method='ols', points=2, lt=FALSE)
print(res3$parms)
plot(res2, log='y', ylim=c(1,100), lty=1, pch=20)
plot(res3, add=TRUE, lty=2, pch=21)
legend(x=0, y=10, pch=c(20,21), lty=c(1,2), legend=c('Preperation 1','Preperation 2'))

## artificial examples
time <- c(5, 180,  360,  540,  720)/60
conc <- c(2.360, 0.061, 0.019, 0.012, 0.024)
res4 <- lee(conc=conc, time=time, method='lad', points=2, lt=FALSE)
print(res4$parms)
plot(res4)
plot(res4, log='y')

time <- seq(1:10)
conc <- c(1,2,3,4,5,5,4,3,2,1)
res5 <- lee(conc=conc, time=time, method='lad', points=2, lt=FALSE)
plot(res5, log='y', ylim=c(1,7), main='', xlab='', ylab='', pch=19)

## real life example
## dataset Indometh of package datasets
require(datasets)
res6 <- data.frame(matrix(ncol=3, nrow=length(unique(Indometh$Subject))))
colnames(res6) <- c('ID', 'initial', 'terminal')
row <- 1
for(i in unique(Indometh$Subject)){
   temp <- subset(Indometh, Subject==i)
   res6[row, 1] <- unique(temp$Subject)
   res6[row, c(2:3)] <- lee(conc=temp$conc, time=temp$time, method='lad')$parms[1,]
   row <- row + 1
}
print(res6)

# geometric means and corresponding two-sided CIs
exp(mean(log(res6$initial)))
exp(t.test(log(res6$initial), conf.level=0.95)$conf.int)

exp(mean(log(res6$terminal)))
exp(t.test(log(res6$terminal), conf.level=0.95)$conf.int)

}

\keyword{misc}
