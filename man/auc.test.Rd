\name{auc.test}
\alias{auc.test}
\alias{auc.batch.test}
\alias{auc.ssd.test}
\title{Hypothesis testing for Area Under the Concentration Time Curve in Complete and Incomplete Data Designs}
\description{Hypothesis testing for an AUC or for the difference between two AUCs assessed in complete and incomplete data designs.}
\usage{
auc.test(conc, time, theta=0, group=NULL, 
        alternative=c("two.sided", "less", "greater"), 
        conf.level=0.95, design=c("ssd","batch","complete"))	     
}

\arguments{
  \item{conc}{ Levels of concentrations. For batch designs a list is required, while a vector is expected otherwise. Multiple measurements per time point are expected to be adjacent to each other. }
  \item{time}{ Time points of concentration assessment. For batch designs a list with time points containing time points for each batch is required. Otherwise a vector is required. One time point for each concentration measured needs to be specified.}
  \item{theta}{ Value to test against (default=\code{0}).}
  \item{group}{ A grouping variable (default=\code{NULL}). For batch designs a list is required, while a vector is expected otherwise. If specified a confidence interval for the difference will be calculated. }
  \item{alternative}{ A character string specifying the alternative hypothesis (default=\code{two.sided}). } 
  \item{conf.level}{ Confidence level (default=\code{0.95}). }
  \item{design}{A character string indicating the type of design used. Possible values are \code{ssd} (the default) for a serial sampling design,
          \code{batch} for a batch design and \code{complete} for a complete data design. }
}

\details{
One- or two-sided hypothesis testing for an AUC (from 0 to the last time point) or for the difference between two AUCs based on the t-distribution for serial sampling, batch and complete data designs. In a serial sampling design only one measurement is available per subject, while in a batch design multiple time points are measured for each subject. In a complete data design measurements are taken for all subjects at all time points. The AUC (from 0 to the last time point) is calculated using the linear trapezoidal rule on the arithmetic means at the different time points.\cr\cr

If group=NULL a hypothesis test for an AUC is performed. If group specifies a factor variable with exactly two levels, the difference between two AUCs is tested. In addition the corresponding confidence intervals are also provided.\cr 

The test is performed using the t-distribution based method which uses the critical value from a t-distribution with Satterthwaite"s approximation (Satterthwaite, 1946) to the degrees of freedom for calculation of confidence intervals as presented in Tang-Liu and Burke (1988), Nedelman et al (1995), Holder et al (1999) and in Jaki and Wolfsegger (in press). \cr 

For serial sampling designs missing data are omitted, while missing values are not permitted in batch designs. \cr

Equal sample size per time point is required for batch designs. 
}

\value{An object of the class PK including confidence intervals and p-value.}

\seealso{\code{\link{auc}}, \code{\link{auc.ci}}.}

\references{
Bailer A. J. (1988). Testing for the equality of area under the curves when using destructive measurement techniques. \emph{Journal of Pharmacokinetics and Biopharmaceutics}, 16(3):303-309. \cr

Holder D. J., Hsuan F., Dixit R. and Soper K. (1999). A method for estimating and testing area under the curve in serial sacrifice, batch, and complete data designs. \emph{Journal of Biopharmaceutical Statistics}, 9(3):451-464.\cr

Jaki T. and Wolfsegger M. J. (In press). A theoretical framework for estimation of AUCs in complete and incomplete sampling designs. \emph{Statistics in Biopharmaceutical Research}. \cr

Nedelman J. R., Gibiansky E. and Lau D. T. W. (1995). Applying Bailer"s method for AUC confidence intervals to sparse sampling. \emph{Pharmaceutical Research}, 12(1):124-128. \cr

Satterthwaite F. E. (1946). An approximate distribution of estimates of variance components. \emph{Biometrics Bulletin}, 2:110-114.  \cr 

Tang-Liu D. D.-S. and Burke P. J. (1988). The effect of azone on ocular levobunolol absoprtion: Calculating the area under the curve and its standard error using tissue sampling compartments. \emph{Pharmaceutical Research}, 5(4):238-241. \cr
}

\author{Thomas Jaki}

\examples{
## example from Nedelman et al. (1995) for a serial sampling design
m.030 <- c(391, 396, 649, 1990, 3290, 3820, 844, 1650, 75.7, 288)
f.030 <- c(353, 384, 625, 1410, 1020, 1500, 933, 1030, 0, 80.5)
m.100 <- c(1910, 2550, 4230, 5110, 7490, 13500, 4380, 5380, 260, 326)
f.100 <- c(2790, 3280, 4980, 7550, 5500, 6650, 2250, 3220, 213, 636)
time <- c(1, 1, 2, 2, 4, 4, 8, 8, 24, 24)

## two sided test
auc.test(conc=m.030, time=time, theta=60000, 
         alternative="two.sided", design="ssd")

## one sided test
auc.test(conc=m.030, time=time, theta=60000, 
         alternative="greater", design="ssd")

data <- data.frame(conc=c(m.030, f.030, m.100, f.100), 
                   time=rep(time,4), 
                   sex=c(rep("m", 10), rep("f", 10), rep("m", 10), rep("f", 10)),
                   dose=c(rep(30, 20), rep(100, 20)))

# comparing the two dose levels
data$concadj <- data$conc / data$dose
auc.test(conc=data$concadj, time=data$time, theta=0, group=data$dose, design="ssd")

## a batch design example from Jaki and Wolfsegger (in press), originally in Holder et al. (1999).
conc<-list(batch1=c(0,0,0,1.75,2.2,1.58,4.63,2.99,1.52),batch2=c(3.03,1.98,2.22,3.34,1.3,1.22),
           batch3=c(3.54,2.84,2.55,0.3,0.0421,0.231))
time<-list(batch1=c(rep(0,3),rep(1,3),rep(6,3)), batch2=c(rep(2,3),rep(10,3)),
           batch3=c(rep(4,3),rep(24,3)))
auc.test(conc, time, theta=100, design="batch")

## difference of two AUCs in batch design from Jaki and Wolfsegger (in press),
## originally in Holder et al. (1999).
conc<-list(batch1=c(0,0,0,0,0,0,1.75,2.2,1.58,3.21,2.01,2.73,4.63,2.99,1.52,3.25,3.31,2.22),
           batch2=c(3.03,1.98,2.22,3.76,3.91,4.37,3.34,1.3,1.22,1.66,6.83,2.98),
           batch3=c(3.54,2.84,2.55,2.14,4.03,2.78,0.3,0.0421,0.231,0.0294,1,0.29))
time<-list(batch1=c(rep(0,6),rep(1,6),rep(6,6)), batch2=c(rep(2,6),rep(10,6)),
           batch3=c(rep(4,6),rep(24,6)))
group<-list(batch1=c(1,1,1,2,2,2,1,1,1,2,2,2,1,1,1,2,2,2),
            batch2=c(1,1,1,2,2,2,1,1,1,2,2,2),
            batch3=c(1,1,1,2,2,2,1,1,1,2,2,2))
auc.test(conc, time, theta=0, group=group, design="batch")

}

\keyword{misc}
