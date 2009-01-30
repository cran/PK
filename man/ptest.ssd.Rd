\name{ptest.ssd}
\alias{ptest.ssd}
\title{Comparing AUCs in Serial Sampling Designs }
\description{Comparing two AUCs assessed in a serial sampling design with a permutation test.}
\usage{
ptest.ssd(conc, time, group, alternative=c("two.sided", "less", "greater"), 
          nsample=1000, data)	     
}

\arguments{
  \item{conc}{ Levels of concentrations. }
  \item{time}{ Time points of concentration assessment. }
  \item{group}{ Grouping variable with two levels. }
  \item{alternative}{ Character string specifying the alternative hypothesis (default=\code{two.sided}). } 
  \item{nsample}{ Number of resampling iterations (default=\code{1000}). } 
  \item{data}{ Optional data frame containing variables named as \code{conc}, \code{time} and \code{group}. } 
}

\details{
Comparing two AUCs assessed in a serial sampling design with a permutation test using a Monte Carlo approximation for the permutation distribution. The difference between two AUCs is used as the test statistic and not the z-statistic as suggested in Bailer and Ruberg (1995). In a serial sampling design only one measurement is available per analysis subject at a specific time point. 
}

\value{
A data frame consisting of:
   \item{statistic}{ estimate for the difference.}
   \item{p.value}{ p-value. }
}

\seealso{\code{\link{auc.ci}}, \code{\link{eqv.ssd}}.}

\note{Records including missing values are omitted. }

\references{
Bailer A. J. (1988). Testing for the equality of area under the curves when using destructive measurement techniques. \emph{Journal of Pharmacokinetics and Biopharmaceutics}, 16(3):303-309. \cr

Bailer J. A. and Ruberg S. J. (1995). Randomization tests for assessing the equality of area under curves for studies using destructive sampling. \emph{Journal of Applied Toxicology}, 16(5):391-395. \cr 

Nedelman J. R., Gibiansky E. and Lau D. T. W. (1995). Applying Bailer"s method for AUC confidence intervals to sparse sampling. \emph{Pharmaceutical Research}, 12(1):124-128. \cr 
}

\author{Martin J. Wolfsegger}

\examples{
## example from Nedelman et al. (1995)
m.030 <- c(391, 396, 649, 1990, 3290, 3820, 844, 1650, 75.7, 288)
f.030 <- c(353, 384, 625, 1410, 1020, 1500, 933, 1030, 0, 80.5)
m.100 <- c(1910, 2550, 4230, 5110, 7490, 13500, 4380, 5380, 260, 326)
f.100 <- c(2790, 3280, 4980, 7550, 5500, 6650, 2250, 3220, 213, 636)
time <- c(1,1,2,2,4,4,8,8,24,24)

data <- data.frame(conc=c(m.030, f.030, m.100, f.100), 
                   time=rep(time, 4), 
                   sex=c(rep("m", 10), rep("f", 10), rep("m", 10), rep("f", 10)),
                   dose=c(rep(30, 20), rep(100, 20)))

data$concadj <- data$conc / data$dose

set.seed(523423)
ptest.ssd(conc=data$concadj, time=data$time, group=data$dose)

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

data <- data.frame(conc=c(grp1, grp2, grp3), time=rep(time, 3),
                   group=c(rep(1, length(grp1)), rep(2, length(grp2)), 
                   rep(3, length(grp3))))

## function call with data frame with subsequent multiple comparisons
set.seed(62432) 
pvalue <- rep(NA, 3)
pvalue[1] <- ptest.ssd(data=subset(data, group==1 | group==2), nsample=100)$p.value
pvalue[2] <- ptest.ssd(data=subset(data, group==1 | group==3), nsample=100)$p.value
pvalue[3] <- ptest.ssd(data=subset(data, group==2 | group==3), nsample=100)$p.value

pvalue
p.adjust(pvalue, method="holm")
}

\keyword{misc}
