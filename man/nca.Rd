\name{nca}
\alias{nca}
\alias{nca.ssd}
\alias{nca.batch}
\alias{nca.complete}
\title{Estimation of various PK parameters}
\description{Non-compartmental estimation of the area under the concentration versus time curve (AUC) to the last time point, AUC to infinity, area under the first moment curve (AUMC) to infinity, mean residence time (MRT), non-compartmental half-life, total clearance and volume of distribution at steady state.}
\usage{
nca(conc, time, n.tail=3, dose=0, method=c("z", "boott"), conf.level=0.95, 
     nsample=1000, design=c("ssd","batch","complete"), data)

nca.ssd(conc, time, n.tail=3, dose=0, method=c("z", "boott"), 
         conf.level=0.95, nsample=1000, data)

nca.batch(conc, time, n.tail=3, dose=0, method=c("z", "boott"), 
         conf.level=0.95, nsample=1000, data)

nca.complete(conc, time, n.tail=3, dose=0, method=c("z", "boott"), 
         conf.level=0.95, nsample=1000, data)
}

\arguments{
  \item{conc}{ Levels of concentrations. For batch designs a list is required, while a vector is expected otherwise.}
  \item{time}{ Time points of concentration assessment. For batch designs a list is required, while a vector is expected otherwise. One time point for each concentration measured needs to be specified.}
  \item{n.tail}{ Number of last data points used for tail area correction (default=\code{3}). }
  \item{dose}{ Dose administered as an IV bolus (default=\code{0}).}
  \item{method}{ A character string specifying the method for calculation of confidence intervals (default=\code{c("z", "boott")}). }
  \item{conf.level}{ Confidence level (default=\code{0.95}). }
  \item{nsample}{ Number of bootstrap iterations for bootstrap-\emph{t} interval (default=\code{1000}). } 
  \item{design}{ A character string indicating the type of design used. Possible values are \code{ssd} (the default) for a serial sampling design, \code{batch} for a batch design and \code{complete} for a complete data design. }
  \item{data}{Optional data frame containing variables named as \code{id}, \code{conc} and \code{time}.}
}

\details{
Estimation of the area under the concentration versus time curve from zero to the last time point (AUC 0-tlast), total area under the concentration versus time curve from zero to infinity (AUC 0-Inf), area under the first moment curve for zero to infinity (AUMC 0-Inf), mean residence time (MRT), non-compartmental half-life (HL), total clearance (CL) and volume of distribution at steady state (Vss). In a serial sampling design only one measurement is available per subject at a specific time point, while in a batch design multiple time points are measured for each subject. In a complete data design measurements are taken for all subjects at all time points.\cr\cr 

A constant coefficient of variation at the last \code{n.tail} time points is assumed.\cr\cr

The use of the standard errors and confidence intervals for the MRT and Vss from batch designs is depreciated due to very slow asymptotic that usually lead to severe undercoverage. For complete data designs only point estimates are provided. The parameters \code{method}, \code{conf.level} and \code{nsample=1000} are therefore not used. If data for only one subject is provided, the parameters are estimated for this subject while the geometric mean of the estimated parameters is found for multiple subjects (see Cawello, 2003, p. 114). \cr\cr

The AUC 0-tlast is calculated using the linear trapezoidal rule on the arithmetic means at the different time points while the extrapolation necessary for the AUC 0-Inf and AUMC 0-Inf is achieved assuming an exponential decay on the last \code{n.tail} time points. The other parameters are functions of these PK parameters and of the dosage and are defined as in Wolfsegger and Jaki (2009).\cr\cr

Two different confidence intervals are computed: an asymptotic confidence interval and a bootstrap-\emph{t} interval. The \code{z} method is based on the limit distribution of the parameter using the critical value from a normal distribution for calculation of confidence intervals together with asymptotic variances. The bootstrap-\emph{t} interval uses the same asymptotic variances, but while the critical value is obtained by the bootstrap-\emph{t} approach. If \code{nsample=0} only the asymptotic interval will be computed.\cr\cr

If \code{data} is specified the variable names \code{conc} and \code{time} are required and represent the corresponding variables.\cr\cr

Note that some estimators as provided assume IV bolus administration. If an oral administration is used\cr
- The clearance needs to be adjusted by the bioavailability, f. This can be achieved by either multiplying the obtained estimator by f or adjusting the dose parameter accordingly. \cr
- The MRT estimate produced corresponds to the mean transit time (MTT) which is the sum of MRT and mean absorption time (MAT).\cr
- HL and Vss are functions of the MRT and hence they will not be valid under oral administration.
}

\value{
An object of the class PK containing the following components: \cr 
  \item{est}{Point estimates.}
  \item{CIs}{Point estimates, standard errors and confidence intervals. }
  \item{conc}{Levels of concentrations. } 
  \item{conf.level}{Confidence level.}
  \item{design}{Sampling design used.}
  \item{time}{Time points measured.}}

\seealso{\code{\link{auc}}, \code{\link{estimator}}, \code{\link{ci}} and \code{\link{test}}.}

\note{At present only the option serial sampling design is available.} 

\references{
Cawello W. (2003). \emph{Parameters for compartment-free pharmacokinetics. Standardisation of study
design, data analysis and reporting}. Shaker Verlag, Aachen.\cr\cr

Gibaldi M. and Perrier D. (1982). \emph{Pharmacokinetics}. Marcel Dekker, New York and Basel.\cr\cr

Jaki T. and Wolfsegger M. J. (2012). Non-compartmental estimation of pharmacokinetic parameters for flexible sampling designs. \emph{Statistics in Medicine}, 31(11-12):1059-1073. \cr\cr

Wolfsegger M. J. and Jaki T. (2009). Non-compartmental Estimation of Pharmacokinetic Parameters in Serial Sampling Designs. \emph{Journal of Pharmacokinetics and Pharmacodynamics}, 36(5):479-494. \cr	
}

\author{Thomas Jaki and Martin J. Wolfsegger}

\examples{
#### serial sampling designs
## example for a serial sampling data design from Wolfsegger and Jaki (2009)
conc <- c(0, 0, 0, 2.01, 2.85, 2.43, 0.85, 1.00, 0.91, 0.46, 0.35, 0.63, 0.39, 0.32, 
          0.45, 0.11, 0.18, 0.19, 0.08, 0.09, 0.06)
time <- c(rep(0,3), rep(5/60,3), rep(3,3), rep(6,3), rep(9,3), rep(16,3), rep(24,3))

# Direct call of the function
# CAUTION: this might take a few minutes
# Note: 1E4 bootstrap replications were used in the example given 
# in Wolfsegger and Jaki (2009)
set.seed(34534)
nca.ssd(conc=conc, time=time, n.tail=4, dose=200, method=c("z","boott"), 
         conf.level=0.95, nsample=500)

# Call through the wrapper function using data
data <- data.frame(conc=conc, time=time)
nca(data=data, n.tail=4, dose=200, method="z", 
     conf.level=0.95, design="ssd")



#### batch design:
## a batch design example from Holder et al. (1999).
data(Rats)
data <- subset(Rats,Rats$dose==100)

# using the wrapper function
nca(data=data, n.tail=4, dose=100, method="z", 
     conf.level=0.95, design="batch")
# direct call
nca.batch(data=data, n.tail=4, dose=100, method="z", 
     conf.level=0.95)

## example with overlapping batches (Treatment A in Example of Jaki & Wolfsegger 2012)
conc <- list(batch1=c(0,0,0,0, 69.7,37.2,213,64.1,  167,306,799,406, 602,758,987,627,
                      1023,1124,1301,880, 1388,1374,1756,1120, 1481,1129,1665,1598,
                      1346,1043,1529,1481, 658,576,772,851, 336,325,461,492, 
                      84,75.9,82.6,116),
             batch2=c(0,0,0, 29.2,55.9,112.2, 145,153,169, 282,420,532, 727,1033,759, 
                      1360,1388,1425, 1939,1279,1318, 1614,1205,1542, 1238,1113,1386,
                      648,770,786, 392,438,511, 77.3,90.1,97.9))
time <- list(batch1=rep(c(0,0.5,0.75,1,1.5,2,3,4,8,12,24),each=4),
             batch2=rep(c(0,0.25,0.5,0.75,1,1.5,2,3,4,8,12,24),each=3))

nca.batch(conc,time,method="z",n.tail=4,dose=80)


#### complete data design
## example from Gibaldi and Perrier (1982, page 436) for individual PK parameters
time <- c(0, 0.165, 0.5, 1, 1.5, 3, 5, 7.5, 10)
conc <- c(0, 65.03, 28.69, 10.04, 4.93, 2.29, 1.36, 0.71, 0.38)
# using the wrapper function
nca(conc=conc, time=time, n.tail=3, dose=1E6, design="complete")
# direct call
nca.complete(conc=conc, time=time, n.tail=3, dose=1E6)
}

\keyword{htest}
