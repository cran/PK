\name{nca}
\alias{nca}
\alias{nca.ssd}
\title{Estimation of various PK parameters}
\description{Non-compartmental estimation of the area under the concentration versus time curve (AUC) to the last time point, AUC to infinity, area under the first moment curve (AUMC) to infinity, mean residence time (MRT), non-compartmental half-life, total clearance and volume of distribution at steady state.}
\usage{
nca(conc, time, n.tail=3, dose=0, method=c("z", "boott"), conf.level=0.95, 
     nsample=1000, design=c("ssd","batch","complete"), data)

nca.ssd(conc, time, n.tail=3, dose=0, method=c("z", "boott"), 
         conf.level=0.95, nsample=1000, data)
}

\arguments{
  \item{conc}{ Levels of concentrations as a vector.}
  \item{time}{ Time points of concentration assessment as a vector. One time point for each concentration measured needs to be specified.}
  \item{n.tail}{ Number of last data points used for tail area correction (default=\code{3}). }
  \item{dose}{ Dose administered as an IV bolus (default=\code{0}).}
  \item{method}{ A character string specifying the method for calculation of confidence intervals (default=\code{c("z", "boott")}). }
  \item{conf.level}{ Confidence level (default=\code{0.95}). }
  \item{nsample}{ Number of bootstrap iterations for bootstrap-\emph{t} interval (default=\code{1000}). } 
  \item{design}{ A character string indicating the type of design used. Possible values are \code{ssd} (the default) for a serial sampling design, \code{batch} for a batch design and \code{complete} for a complete data design. }
  \item{data}{Optional data frame containing variables named as \code{id}, \code{conc}, \code{time} and \code{group}.}
}

\details{
Estimation of the area under the concentration versus time curve from zero to the last time point (AUC 0-tlast), total area under the concentration versus time curve from zero to infinity (AUC 0-Inf), area under the first moment curve for zero to infinity (AUMC 0-Inf), mean residence time (MRT), non-compartmental half-life (HL), total clearance (CL) and volume of distribution at steady state (Vss). In a serial sampling design only one measurement is available per subject at a specific time point, while in a batch design multiple time points are measured for each subject. In a complete data design measurements are taken for all subjects at all time points. Note: Currently only the serial sampling design is implemented. Equal sample size per time point and a constant coefficient of variation at the last \code{n.tail} time points are assumed.\cr\cr

The AUC 0-tlast is calculated using the linear trapezoidal rule on the arithmetic means at the different time points while the extrapolation necessary for the AUC 0-Inf and AUMC 0-Inf is achieved assuming an exponential decay on the last \code{n.tail} time points. The other parameters are functions of these PK parameters and of the dosage and are defined as in Wolfsegger and Jaki (2009).\cr\cr

Two different confidence intervals are computed: an asymptotic confidence interval and a bootstrap-\emph{t} interval. The \code{z} method is based on the limit distribution of the parameter using the critical value from a normal distribution for calculation of confidence intervals together with asymptotic variances. The bootstrap-\emph{t} interval uses the same asymptotic variances, but while the critical value is obtained by the bootstrap-\emph{t} approach. If \code{nsample=0} only the asymptotic interval will be computed.\cr\cr

If \code{data} is specified the variable names \code{conc}, \code{time} and \code{group} are required and represent the corresponding variables. If \code{design} is \code{batch} an additional variable \code{id} is required to identify the subject.\cr\cr

Note that some estimators used are valid only under IV bolus. Straightforward modifications to the estimators are required for other routes of administrations.
}

\value{
An object of the class PK  containing the following components: \cr 
  \item{est}{Point estimates.}
  \item{CIs}{Point estimates, standard errors and confidence intervals. }
  \item{conc}{Levels of concentrations. } 
  \item{conf.level}{Confidence level.}
  \item{design}{Sampling design used.}
  \item{group}{Grouping variable.}
  \item{time}{Time points measured.}}

\seealso{\code{\link{auc}}, \code{\link{estimator}}, \code{\link{ci}} and \code{\link{test}}.}

\note{At present only the option serial sampling design is available.} 

\references{
Wolfsegger M. J. and Jaki T. (2009). Non-compartmental Estimation of Pharmacokinetic Parameters in Serial Sampling Designs. \emph{Journal of Pharmacokinetics and Pharmacodynamics}, 36(5):479-494. \cr	
}

\author{Thomas Jaki and Martin J. Wolfsegger}

\examples{
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
}

\keyword{htest}
