\name{eqv.ssd}
\alias{eqv.ssd}
\title{Establishing Bioequivalence in Serial Sampling Designs}
\description{Calculation of a confidence interval for the ratio of two AUCs assessed in a serial sampling design.}
\usage{
eqv.ssd(conc, time, group, method=c("fieller", "asymp", "boott"), 
        conf.level=0.90, strata=NULL, nsample=1000, data)	     
}

\arguments{
  \item{conc}{ vector of concentrations. }
  \item{time}{ vector of time points. }
  \item{group}{ vector of grouping variable. }
  \item{method}{ character string specifying the method for calculation of confidence intervals; default=\code{fieller}. }
  \item{conf.level}{ confidence level; default=\code{0.90}. }
  \item{strata}{ vector of one strata variable, only available for method \code{boott}. }
  \item{nsample}{ number of bootstrap iterations for method \code{boott}; default=\code{1000}. }
  \item{data}{ optional data frame containing variables named as \code{conc}, \code{time}, \code{group} and \code{strata}. } 
}

\details{
Calculation of a confidence interval for the ration of two AUCs (from 0 to the last time point) assessed in a serial sampling design.
 In a serial sampling design only one measurement is available per analysis subject at a specific time point. The AUC (from 0 to the last time point) is calculated using the linear trapezoidal rule on the arithmetic means at the different time points.
 Details for calculation can be found in Wolfsegger (2007) and in Jaki \emph{et al.} (in press). \cr \cr
The \code{fieller} method is based on Fieller's theorem (1954) using the critical value from a t-distribution with Satterthwaite's approximation (1946) to the degrees of freedom for calculation of confidence intervals. \cr \cr
The \code{asymp} method is based on the limit distribution for the ratio using the critical value from a normal distribution for calculation of confidence intervals. \cr \cr
The \code{boott} method uses the standard error of the limit distribution for the ratio where the critical value is obtained by the bootstrap-\emph{t} approach. Using \code{boott} an additional strata variable for bootstrapping can be specified. \cr
}

\value{
A data frame consisting of:
  \item{ratio}{ estimate for ratio of two AUCs.}
  \item{lower}{ lower limit of confidence interval. }
  \item{upper}{ upper limit of confidence interval. }
  \item{df}{ degrees of freedom when using method \code{fieller}. }
}

\seealso{\code{\link{ptest.ssd}}, \code{\link{auc.ssd}}.}

\note{Records including missing values are omitted. }

\references{
Fieller E. C. (1954). Some problems in interval estimation. \emph{Journal of the Royal Statistical Society, Series B}, 16:175-185. \cr 

Jaki T., Wolfsegger M. J. and Ploner M. Confidence intervals for ratios of AUC's in the case of serial sampling: A comparison of seven methods. \emph{Pharmaceutical statistics}, in press. \cr 

Nedelman J. R., Gibiansky E. and Lau D. T. W. (1995). Applying Bailer's method for AUC confidence intervals to sparse sampling. \emph{Pharmaceutical Research}, 12(1):124-128. \cr 

Satterthwaite F. E. (1946). An approximate distribution of estimates of variance components. \emph{Biometrics Bulletin}, 2:110-114. \cr

Wolfsegger M. J. (2007). Establishing bioequivalence in serial sacrifice designs. \emph{Journal of Pharmacokinetics and Pharmacodynamics}, 34(1):103-113. \cr 
}

\author{Martin J. Wolfsegger and Thomas Jaki}

\examples{
## example from Nedelman et al. (1995)
m.030 <- c(391, 396, 649, 1990, 3290, 3820, 844, 1650, 75.7, 288)
f.030 <- c(353, 384, 625, 1410, 1020, 1500, 933, 1030, 0, 80.5)
m.100 <- c(1910, 2550, 4230, 5110, 7490, 13500, 4380, 5380, 260, 326)
f.100 <- c(2790, 3280, 4980, 7550, 5500, 6650, 2250, 3220, 213, 636)
time <- c(1, 1, 2, 2, 4, 4, 8, 8, 24, 24)

data <- data.frame(conc=c(m.030, f.030, m.100, f.100), 
                   time=rep(time, 4), 
                   sex=c(rep('m', 10), rep('f', 10), rep('m', 10), rep('f', 10)),
                   dose=c(rep(30, 20), rep(100, 20)))

data$concadj <- data$conc / data$dose
eqv.ssd(conc=data$concadj, time=data$time, group=data$dose, method='asymp')
eqv.ssd(conc=data$concadj, time=data$time, group=data$dose, method='fieller')

set.seed(310578)
eqv.ssd(conc=data$concadj, time=data$time, group=data$dose, method='boott', nsample=1E4)

set.seed(310578)
eqv.ssd(conc=data$concadj, time=data$time, group=data$dose, method='boott', strata=data$sex, nsample=1E4)
}

\keyword{misc}