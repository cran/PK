\name{eqv.ssd}
\alias{eqv.ssd}
\title{Establishing Bioequivalence in Serial Sampling Designs}
\description{Confidence intervals for the ratio of independent or dependent area under the concentration versus time curves (AUCs) to the last time point in serial sampling designs.}
\usage{
eqv.ssd(conc, time, group, dependent=FALSE, method=c("fieller", "z", "boott"), 
        conf.level=0.90, strata=NULL, nsample=1000, data)	   
}

\arguments{
  \item{conc}{ Levels of concentrations as a vector. }
  \item{time}{ Time points of concentration assessment as a vector. One time point for each concentration measured needs to be specified.}
  \item{group}{A grouping variable as a vector.}
  \item{dependent}{ Logical variable indicating if concentrations are measured on the same subjects for both AUCs (default=\code{FALSE}). Not used.}
  \item{method}{A character string specifying the method for calculation of confidence intervals (default=\code{c("fieller", "z", "boott")}).}
  \item{conf.level}{Confidence level (default=\code{0.90}). }
  \item{strata}{A vector of one strata variable (default=\code{NULL}).}
  \item{nsample}{ Number of bootstrap iterations for method \code{boott} (default=\code{1000}). }
  \item{data}{Optional data frame containing variables named as \code{id}, \code{conc}, \code{time} and \code{group}.}
}

\details{
Calculation of confidence intervals for the ratio of (independent or dependent) AUCs (from 0 to the last time point) for serial sampling designs. In a serial sampling design only one measurement is available per subject. The AUC (from 0 to the last time point) is calculated using the linear trapezoidal rule on the arithmetic means at the different time points. \cr\cr

Confidence intervals are computed using the methods of Wolfsegger (2007) and Jaki \emph{et al.} (2009). Note that the option \code{dependent} is only included for similarity with similar calls but is not used in serial sampling designs as by definition only one sample is obtained per subject then.\cr\cr

The \code{fieller} method is based on Fieller's theorem (1954) which uses the asymptotic standard errors of the individual AUCs and a critical value from a t-distribution with Satterthwaite's approximation (1946) to the degrees of freedom for calculation of confidence intervals. The \code{z} method is based on the limit distribution for the ratio using the critical value from a normal distribution for calculation of confidence intervals. \cr\cr

The \code{boott} method uses the asymptotic standard errors of the ratio of two AUCs while the critical value is obtained by the bootstrap-\emph{t} approach and follows the idea discussed in the context of serial sampling designs in Jaki T. \emph{et al.} (2009). Using \code{boott} an additional strata variable for bootstrapping can be specified in serial sampling designs.\cr\cr

If \code{data} is specified the variable names \code{conc}, \code{time} and \code{group} are required and represent the corresponding variables.
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

\seealso{\code{\link{eqv}}, \code{\link{estimator}}, \code{\link{ci}} and \code{\link{test}}.}

\note{Records including missing values are omitted. }

\references{
Fieller E. C. (1954). Some problems in interval estimation. \emph{Journal of the Royal Statistical Society, Series B}, 16:175-185. \cr\cr

Jaki T., Wolfsegger M. J. and Ploner M. (2009). Confidence intervals for ratios of AUCs in the case of serial sampling: A comparison of seven methods. \emph{Pharmaceutical Statistics}, 8(1):12-24. \cr\cr

Nedelman J. R., Gibiansky E. and Lau D. T. W. (1995). Applying Bailer"s method for AUC confidence intervals to sparse sampling. \emph{Pharmaceutical Research}, 12(1):124-128. \cr\cr

Satterthwaite F. E. (1946). An approximate distribution of estimates of variance components. \emph{Biometrics Bulletin}, 2:110-114. \cr\cr

Wolfsegger M. J. (2007). Establishing bioequivalence in serial sacrifice designs. \emph{Journal of Pharmacokinetics and Pharmacodynamics}, 34(1):103-113. \cr 
}

\author{Martin J. Wolfsegger and Thomas Jaki}

\examples{
## example of a serial sampling design from Nedelman et al. (1995)
data(CPI975)
data <- subset(CPI975,dose>=30)

data$concadj <- data$conc / data$dose

# fieller and z-interval for ratio of dose-normalized AUCs
eqv.ssd(conc=data$concadj, time=data$time, group=data$dose, method=c("z","fieller"))

# bootstrap-t interval for ratio of dose-normalized AUCs stratified for sex
set.seed(310578)
eqv.ssd(conc=data$concadj, time=data$time, group=data$dose, method="boott",
        strata=data$sex, nsample=500)
}

\keyword{htest}
