\name{auc.ci}
\alias{auc.ci}
\alias{auc.ssd.ci}
\alias{auc.batch.ci}
\title{Estimation of Confidence Intervals for Area Under the Concentration Time Curve in Complete and Incomplete Data Designs}
\description{Calculation of confidence intervals for an AUC or for the difference between two AUCs assessed in complete and incomplete data designs.}
\usage{
auc.ci(conc, time, group=NULL, method=c("t", "z", "boott"), 
        alternative=c("two.sided", "less", "greater"), 
        conf.level=0.95, strata=NULL, nsample=1000, 
        design=c('ssd','batch','complete'))	     
}

\arguments{
  \item{conc}{ Levels of concentrations. For batch designs a list is required, while a vector is expected otherwise. Multiple measurements per time point are expected to be adjacent to each other. }
  \item{time}{ Time points of concentration assessment. For batch designs a list with time points containing time points for each batch is required. Otherwise a vector is required. One time point for each concentration measured needs to be specified.}
  \item{group}{ A grouping variable (default=\code{NULL}). For batch designs a list is required, while a vector is expected otherwise. If specified, a confidence interval for the difference will be calculated. }
  \item{method}{ A character string specifying the method for calculation of confidence intervals (default=\code{c("t", "z", "boott")}). }
  \item{alternative}{ A character string specifying the alternative hypothesis (default=\code{"two.sided"}). } 
  \item{conf.level}{ Confidence level (default=\code{0.95}). }
  \item{strata}{ A vector of one strata variable (default=\code{NULL}). Only available for method \code{boott} in a serial sampling design. }
  \item{nsample}{ Number of bootstrap iterations for method \code{boott} (default=\code{1000}). } 
  \item{design}{ A character string indicating the type of design used. Possible values are \code{ssd} (the default) for a serial sampling design,
          \code{batch} for a batch design and \code{complete} for a complete data design. }
}

\details{
Calculation of confidence intervals for an AUC (from 0 to the last time point) or for the difference between two AUCs for serial sampling, batch and complete data designs. In a serial sampling design only one measurement is available per subject, while in a batch design multiple time points are measured for each subject. In a complete data design measurements are taken for all subjects at all time points. The AUC (from 0 to the last time point) is calculated using the linear trapezoidal rule on the arithmetic means at the different time points.\cr\cr

If group=NULL a confidence interval for an AUC is calculated. If group specifies a factor variable with exactly two levels, a confidence interval for the difference between two AUCs is calculated. \cr 

The \code{t} method uses the critical value from a t-distribution with Satterthwaite's approximation (Satterthwaite, 1946) to the degrees of freedom for calculation of confidence intervals as presented in Tang-Liu and Burke (1988), Nedelman et al (1995), Holder et al (1999) and in Jaki and Wolfsegger (in press). \cr 

The \code{z} method uses the critical value from a normal distribution for calculation of confidence intervals as presented in Bailer (1988) or in Jaki and Wolfsegger (in press). \cr 

The \code{boott} method uses bootstrap-\emph{t} confidence intervals as presented in Jaki and Wolfsegger (in press). Using \code{boott} an additional strata variable for bootstrapping can be specified in the case of serial sampling. \cr

For serial sampling designs missing data are omitted and unequal sample sizes per time point are allowed.\cr

For batch designs missing values are not permitted in batch designs and equal sample size per time point is required.
}

\seealso{\code{\link{auc}}, \code{\link{auc.test}}.}

\value{An object of the class PK including point estimates and confidence intervals.}

\note{This is a wrapper function for \code{\link{auc.complete.ci}, \link{auc.batch.ci}} and \code{\link{auc.ssd.ci}}. 
The function calculates point and interval estimates for AUC (from 0 to the last time point) using \code{\link{auc.complete}, \link{auc.batch}} and \code{\link{auc.ssd}} and corresponding confidence interval functions \code{\link{auc.complete.ci}, \link{auc.batch.ci}} and \code{\link{auc.ssd.ci}}.
}

\references{
Bailer A. J. (1988). Testing for the equality of area under the curves when using destructive measurement techniques. \emph{Journal of Pharmacokinetics and Biopharmaceutics}, 16(3):303-309. \cr

Holder D. J., Hsuan F., Dixit R. and Soper K. (1999). A method for estimating and testing area under the curve in serial sacrifice, batch, and complete data designs. \emph{Journal of Biopharmaceutical Statistics}, 9(3):451-464.\cr

Jaki T. and Wolfsegger M. J. (In press). A theoretical framework for estimation of AUCs in complete and incomplete sampling designs. \emph{Statistics in Biopharmaceutical Research}. \cr

Nedelman J. R., Gibiansky E. and Lau D. T. W. (1995). Applying Bailer's method for AUC confidence intervals to sparse sampling. \emph{Pharmaceutical Research}, 12(1):124-128. \cr

Satterthwaite F. E. (1946). An approximate distribution of estimates of variance components. \emph{Biometrics Bulletin}, 2:110-114.  \cr 

Tang-Liu D. D.-S. and Burke P. J. (1988). The effect of azone on ocular levobunolol absoprtion: Calculating the area under the curve and its standard error using tissue sampling compartments. \emph{Pharmaceutical Research}, 5(4):238-241. \cr
}



\author{Thomas Jaki and Martin J. Wolfsegger}

\examples{
## example from Nedelman et al. (1995)
m.030 <- c(391, 396, 649, 1990, 3290, 3820, 844, 1650, 75.7, 288)
f.030 <- c(353, 384, 625, 1410, 1020, 1500, 933, 1030, 0, 80.5)
m.100 <- c(1910, 2550, 4230, 5110, 7490, 13500, 4380, 5380, 260, 326)
f.100 <- c(2790, 3280, 4980, 7550, 5500, 6650, 2250, 3220, 213, 636)
time <- c(1, 1, 2, 2, 4, 4, 8, 8, 24, 24)

auc.ci(conc=m.030, time=time, method=c('z', 't'), design='ssd')
auc.ci(conc=f.030, time=time, method=c('z', 't'), design='ssd')

auc.ci(conc=m.100, time=time, method=c('z', 't'), design='ssd')
auc.ci(conc=f.100, time=time, method=c('z', 't'), design='ssd')

data <- data.frame(conc=c(m.030, f.030, m.100, f.100), 
                   time=rep(time, 4), 
                   sex=c(rep('m', 10), rep('f', 10), rep('m', 10), rep('f', 10)),
                   dose=c(rep(30, 20), rep(100, 20)))

# comparing the two dose levels
data$concadj <- data$conc / data$dose
auc.ci(conc=data$concadj, time=data$time, group=data$dose, method=c('z', 't'), design='ssd')

## caution this might take a few minutes
## comparing two dose level using stratified resampling per gender
set.seed(260151)
auc.ci(conc=data$concadj, time=data$time, group=data$dose, method='boott',
       strata=data$sex, design='ssd')

## example from Bailer (1988)
time <- c(rep(0,4), rep(1.5,4), rep(3,4), rep(5,4), rep(8,4))
grp1 <- c(0.0658, 0.0320, 0.0338, 0.0438, 0.0059, 0.0030, 0.0084,
          0.0080, 0.0000, 0.0017, 0.0028, 0.0055, 0.0000, 0.0037,
          0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000)

grp2 <- c(0.2287, 0.3824, 0.2402, 0.2373, 0.1252, 0.0446, 0.0638,
          0.0511, 0.0182, 0.0000, 0.0117, 0.0126, 0.0000, 0.0440,
          0.0039, 0.0040, 0.0000, 0.0000, 0.0000, 0.0000)

grp3 <- c(0.4285, 0.5180, 0.3690, 0.5428, 0.0983, 0.0928, 0.1128,
          0.1157, 0.0234, 0.0311, 0.0344, 0.0349, 0.0032, 0.0052,
          0.0049, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000)

auc.ci(conc=grp1, time=time, method='z', design='ssd')
auc.ci(conc=grp2, time=time, method='z', design='ssd')
auc.ci(conc=grp3, time=time, method='z', design='ssd')

data <- data.frame(conc=c(grp1, grp2, grp3), time=rep(time, 3),
                   group=c(rep(1, length(grp1)), rep(2, length(grp2)), rep(3, length(grp3))))

## function call with data frame using simultaneous confidence intervals based 
## on bonferroni adjustment
auc.ci(subset(data, group==1 | group==2)$conc, subset(data, group==1 | group==2)$time, 
       group=subset(data, group==1 | group==2)$group, method=c('z', 't'), 
       conf.level=1-0.05/3, design='ssd')
auc.ci(subset(data, group==1 | group==3)$conc, subset(data, group==1 | group==2)$time, 
       group=subset(data, group==1 | group==3)$group, method=c('z', 't'), 
       conf.level=1-0.05/3, design='ssd')
auc.ci(subset(data, group==2 | group==3)$conc, subset(data, group==1 | group==2)$time, 
       group=subset(data, group==2 | group==3)$group, method=c('z', 't'), 
       conf.level=1-0.05/3, design='ssd')

## a batch design example from Jaki and Wolfsegger (in press), originally in Holder et al. (1999).
conc<-list(batch1=c(0,0,0,1.75,2.2,1.58,4.63,2.99,1.52),batch2=c(3.03,1.98,2.22,3.34,1.3,1.22),
           batch3=c(3.54,2.84,2.55,0.3,0.0421,0.231))
time<-list(batch1=c(0,0,0,1,1,1,6,6,6),batch2=c(2,2,2,10,10,10),batch3=c(4,4,4,24,24,24))
auc.ci(conc,time,method=c('z','t'), design='batch')

## difference of two AUCs in batch design from Jaki and Wolfsegger (in press),
## originally in Holder et al. (1999).
conc<-list(batch1=c(0,0,0,0,0,0,1.75,2.2,1.58,3.21,2.01,2.73,4.63,2.99,1.52,3.25,3.31,2.22),
           batch2=c(3.03,1.98,2.22,3.76,3.91,4.37,3.34,1.3,1.22,1.66,6.83,2.98),
           batch3=c(3.54,2.84,2.55,2.14,4.03,2.78,0.3,0.0421,0.231,0.0294,1,0.29))
time<-list(batch1=c(rep(0,6),rep(1,6),rep(6,6)),batch2=c(rep(2,6),rep(10,6)),
           batch3=c(rep(4,6),rep(24,6)))
group<-list(batch1=c(1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2),
            batch2=c(1,1,1,2,2,2,1,1,1,2,2,2),
            batch3=c(1,1,1,2,2,2,1,1,1,2,2,2))
auc.ci(conc,time,group=group, design='batch')

}

\keyword{misc}
