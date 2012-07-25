\name{test}
\alias{test}
\alias{print.PKtest}
\alias{summary.PKtest}
\title{Function for hypothesis testing for objects of class PK}
\description{Generic function for hypothesis testing based on an object of class \code{PK}.}
\usage{
test(obj, theta=0, method = c("t", "fieller", "z", "resample"), nsample = 1000) 

\method{print}{PKtest}(x,hyp=FALSE,...)

\method{summary}{PKtest}(object,...)
}

\arguments{
  \item{obj}{ An output object of class PK.}
  \item{x}{ An output object of class PKtest.}
  \item{object}{ An output object of class PK test.}
  \item{theta}{ The reference value to be tested against. If multiple parameters are to be tested a vector can be supplied.}
  \item{method}{ A character string specifying the method for calculation of the test statistic. Possible values are \code{t} (the default) and \code{fieller} for a t-test based method, \code{z} for a z-test and \code{resample} for either a bootstrap or a permutation test.}
  \item{nsample}{ Number of resamples for the permutation/bootstrap test (default=\code{1000}). } 
  \item{hyp}{ Logical variable indicating if hypothesis tests should be printed explicitly (default=\code{FALSE}).}
  \item{...}{Arguments to be passed to methods, such as graphical.}
}

\details{
  Generic function to perform hypothesis test(s). \cr\cr

  The reference value for the test is to be specified in \code{theta}. If multiple tests are performed theta can be a vector.\cr\cr

  For method \code{"resample"} a permutation test is used for the difference of AUCs while a one-sample bootstrap test based on inverting a bootstrap-t statistic is implemented.

}

\seealso{\code{\link{auc}}, \code{\link{eqv}} and \code{\link{nca}}.}

\value{An object of the class PKtest containing the following components: \cr 
  \item{stat}{Test statistics.}
  \item{p.value}{p-values. }
  \item{theta}{Reference value(s) tested against.} 
  \item{conf.level}{Confidence level.}
  \item{alternative}{Type of alternative used.}
  \item{df}{Degrees of freedom of method \code{"t"}.}
  \item{design}{Sampling design used.}
  \item{method}{Type of test used.}}

\references{
Efron B and Tibshirani R. J. (1993). \emph{An introduction to the bootstrap}, Chapman and Hall, New York.\cr\cr

Holder D. J., Hsuan F., Dixit R. and Soper K. (1999). A method for estimating and testing area under the curve in serial sacrifice, batch, and complete data designs. \emph{Journal of Biopharmaceutical Statistics}, 9(3):451-464.\cr\cr

Wolfsegger M. J. and Jaki T. (2009) Assessing systemic drug exposure in repeated dose toxicity studies in the case of complete and incomplete sampling. \emph{Biometrical Journal}, 51(6):1017:1029.\cr
}


\author{Thomas Jaki}

\examples{
## example for a serial sampling data design from Wolfsegger and Jaki (2009)
conc <- c(0, 0, 0, 2.01, 2.85, 2.43, 0.85, 1.00, 0.91, 0.46, 0.35, 0.63, 0.39, 0.32, 
          0.45, 0.11, 0.18, 0.19, 0.08, 0.09, 0.06)
time <- c(rep(0,3), rep(5/60,3), rep(3,3), rep(6,3), rep(9,3), rep(16,3), rep(24,3))

obj <- nca(conc=conc, time=time, n.tail=4, dose=200, method="z", 
            conf.level=0.95, design="ssd")

## testing all parameters against different values using a z-test
res <- test(obj, theta=c(11, 12, 90, 7, 5, 16, 120), method="z")

print(res)

## a batch design example from Holder et al. (1999).
data(Rats)
data <- subset(Rats,Rats$dose==100)

obj <- auc(data=data,method=c('z','t'), design='batch')

## t-test
res <- test(obj, theta=100, method="t")

## making the hypothesis explicit
summary(res)

## bootstrap test for bioequivalence
# Note: This can take a few minutes
data(Glucose)
## one-sided permutation test
obj <- auc(conc=Glucose$conc, time=Glucose$time, group=Glucose$date, 
     method=c("t"), conf.level=0.90, alternative='less', 
     nsample=500, design="complete")

test(obj, theta=1, method="resample", nsample=500)

}

\keyword{htest}
