\name{auc.ssd}
\alias{auc.ssd}
\title{Estimation of AUC in Serial Sampling Designs}
\description{Calculation of a confidence interval for an AUC or for the difference between two AUCs assessed in a serial sampling design.}
\usage{
auc.ssd(conc, time, group=NULL, method=c("tang.burke", "bailer", "boott"), 
        alternative=c("two.sided", "less", "greater"), 
        conf.level=0.95, strata=NULL, nsample=1000, data)	     
}

\arguments{
  \item{conc}{ vector of concentrations. }
  \item{time}{ vector of time points. }
  \item{group}{ vector of grouping variable, if specified a confidence interval for the difference will be calculated; default=NULL. }
  \item{method}{ character string specifying the method for calculation of confidence intervals; default=\code{tang.burke}. }
  \item{alternative}{ character string specifying the alternative hypothesis; default=\code{two.sided}. } 
  \item{conf.level}{ confidence level; default=\code{0.95}. }
  \item{strata}{ vector of one strata variable, only available for method \code{boott}. }
  \item{nsample}{ number of bootstrap iterations for method \code{boott}; default=\code{1000}. } 
  \item{data}{ optional data frame containing variables named as \code{conc}, \code{time}, \code{group} and \code{strata}. }
}

\details{
Calculation of a confidence interval for an AUC (from 0 to the last time point) or for the difference between two AUCs assessed in a serial sampling design. \cr 
In a serial sampling design only one measurement is available per analysis subject at a specific time point. The AUC (from 0 to the last time point) is calculated using the linear trapezoidal rule on the arithmetic means at the different time points.
 If group=NULL a confidence interval for an AUC is calculated. If group specifies a factor variable (with two levels), a confidence interval for the difference between two AUCs is calculated. \cr \cr
The \code{tang.burke} method uses the critical value from a t-distribution with Satterthwaite's approximation (1946) to the degrees of freedom for calculation of confidence intervals as presented in Tang-Liu and Burke (1988) or in Nedelman \emph{et al.} (1995). \cr \cr
The \code{bailer} method uses the critical value from a normal distribution for calculation of confidence intervals as presented in Bailer (1988). \cr \cr
The \code{boott} method uses bootstrap-\emph{t} confidence intervals. Using \code{boott} an additional strata variable for bootstrapping can be specified. \cr
}

\value{
A data frame consisting of:
  \item{est}{ estimate for AUC or estimate for difference between two AUCs.}
  \item{stderr}{ standard error for estimate. }
  \item{lower}{ lower limit of confidence interval. }
  \item{upper}{ upper limit of confidence interval. }
  \item{df}{ degrees of freedom when using method \code{tang.burke}. }
}

\seealso{\code{\link{ptest.ssd}}, \code{\link{eqv.ssd}}.}

\note{Records including missing values are omitted. }

\references{
Bailer A. J. (1988). Testing for the equality of area under the curves when using destructive measurement techniques. \emph{Journal of Pharmacokinetics and Biopharmaceutics}, 16(3):303-309. \cr

Nedelman J. R., Gibiansky E. and Lau D. T. W. (1995). Applying Bailer's method for AUC confidence intervals to sparse sampling. \emph{Pharmaceutical Research}, 12(1):124-128. \cr 

Satterthwaite F. E. (1946). An approximate distribution of estimates of variance components. \emph{Biometrics Bulletin}, 2:110-114.  \cr 

Tang-Liu D. D.-S. and Burke J. P. (1988). The effect of azone on ocular levobunolol absorption: calculating the area under the curve and its standard error using tissue sampling compartments. \emph{Pharmaceutical Research}, 5(4):238-241. \cr 
}

\author{Martin J. Wolfsegger and Thomas Jaki}

\examples{
## example from Nedelman et al. (1995)
m.030 <- c(391, 396, 649, 1990, 3290, 3820, 844, 1650, 75.7, 288)
f.030 <- c(353, 384, 625, 1410, 1020, 1500, 933, 1030, 0, 80.5)
m.100 <- c(1910, 2550, 4230, 5110, 7490, 13500, 4380, 5380, 260, 326)
f.100 <- c(2790, 3280, 4980, 7550, 5500, 6650, 2250, 3220, 213, 636)
time <- c(1, 1, 2, 2, 4, 4, 8, 8, 24, 24)

auc.ssd(conc=m.030, time=time, method=c('bailer', 'tang.burke'))
auc.ssd(conc=f.030, time=time, method=c('bailer', 'tang.burke'))

auc.ssd(conc=m.100, time=time, method=c('bailer', 'tang.burke'))
auc.ssd(conc=f.100, time=time, method=c('bailer', 'tang.burke'))

data <- data.frame(conc=c(m.030, f.030, m.100, f.100), 
                   time=rep(time, 4), 
                   sex=c(rep('m', 10), rep('f', 10), rep('m', 10), rep('f', 10)),
                   dose=c(rep(30, 20), rep(100, 20)))

data$concadj <- data$conc / data$dose
auc.ssd(conc=data$concadj, time=data$time, group=data$dose, method=c('bailer', 'tang.burke'))

set.seed(260151)
auc.ssd(conc=data$concadj, time=data$time, group=data$dose, method='boott', strata=data$sex)

## example from Bailer (1988)
time <- c(rep(0, 4), rep(1.5, 4), rep(3, 4), rep(5, 4), rep(8, 4))
grp1 <- c(0.0658, 0.0320, 0.0338, 0.0438, 0.0059, 0.0030, 0.0084,
          0.0080, 0.0000, 0.0017, 0.0028, 0.0055, 0.0000, 0.0037,
          0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000)

grp2 <- c(0.2287, 0.3824, 0.2402, 0.2373, 0.1252, 0.0446, 0.0638,
          0.0511, 0.0182, 0.0000, 0.0117, 0.0126, 0.0000, 0.0440,
          0.0039, 0.0040, 0.0000, 0.0000, 0.0000, 0.0000)


grp3 <- c(0.4285, 0.5180, 0.3690, 0.5428, 0.0983, 0.0928, 0.1128,
          0.1157, 0.0234, 0.0311, 0.0344, 0.0349, 0.0032, 0.0052,
          0.0049, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000)


auc.ssd(conc=grp1, time=time, method='bailer')

auc.ssd(conc=grp2, time=time, method='bailer')

auc.ssd(conc=grp3, time=time, method='bailer')

data <- data.frame(conc=c(grp1, grp2, grp3), time=rep(time, 3),
                   group=c(rep(1, length(grp1)), rep(2, length(grp2)), rep(3, length(grp3))))

## function call with data frame using simultaneous confidence intervals based on bonferroni adjustment
auc.ssd(data=subset(data, group==1 | group==2), method=c('bailer', 'tang.burke'), conf.level=1-0.05/3)
auc.ssd(data=subset(data, group==1 | group==3), method=c('bailer', 'tang.burke'), conf.level=1-0.05/3)
auc.ssd(data=subset(data, group==2 | group==3), method=c('bailer', 'tang.burke'), conf.level=1-0.05/3)
}

\keyword{misc}