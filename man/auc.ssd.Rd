\name{auc.ssd}
\alias{auc.ssd}
\title{Estimation of the Area Under the Concentration Time Curve in Serial Sampling Designs}
\description{Non-compartmental estimation of area under the concentration time curve (AUC) and area under the first moment curve (AUMC) in serial sampling designs.}
\usage{auc.ssd(conc, time, exact=NA, n.interpolate=2, n.tail=3) }

\arguments{
  \item{conc}{ A vector of concentration levels. }
  \item{time}{ A vector of time points of concentration assessment. One time point for each concentration measured needs to be specified.}
  \item{exact}{ Time point for linear interpolation/extrapolation (default=\code{NA}). Currently not implemented. }
  \item{n.interpolate}{ Number of last data points used for linear interpolation/extrapolation (default=\code{2}). Currently not implemented. }
  \item{n.tail}{ Number of last data points used for tail area estimation (default=\code{3}). }
}

\details{
Estimation of the area under the concentration time curve (AUC 0-tlast) and area under the first moment curve (AUMC 0-tlast) for serial sampling designs. In a serial sampling design only one measurement is available per subject at a specific time point, while in a batch design multiple time points are measured for each subject. 

Estimation of the linearly interpolated area under the concentration time curve (AUC 0-exact) and linearly interpolated area under the moment curve (AUMC 0-exact). The time point for linear interpolation/extrapolation must be after the last but one time point of concentration assessments. \cr

Estimation of the total area under the concentration curve (AUC 0-infinity) and total area under the first moment curve (AUMC 0-infinity) by using a tail area correction similar as suggested by Perrier and Gibaldi (Appendix D, 1982). \cr

}

\seealso{\code{\link{auc}}, \code{\link{auc.ci}}, \code{\link{auc.test}}.}

\value{An object of the class PK including AUC and AUMC estimates.}

\references{
Gibaldi M. and Perrier D. (1982). \emph{Pharmacokinetics. 2nd Edition}. Marcel Dekker, New York and Basel.  \cr 

Nedelman J. R., Gibiansky E. and Lau D. T. W. (1995). Applying Bailer's method for AUC confidence intervals to sparse sampling. \emph{Pharmaceutical Research}, 12(1):124-128. \cr

}

\author{Thomas Jaki}

\examples{ 
## example from Nedelman et al. (1995)
time <- c(1, 1, 2, 2, 4, 4, 8, 8, 24, 24)
m.030 <- c(391, 396, 649, 1990, 3290, 3820, 844, 1650, 75.7, 288)
f.030 <- c(353, 384, 625, 1410, 1020, 1500, 933, 1030, 0, 80.5)
auc.ssd(conc=m.030, time=time, n.tail=3)
auc.ssd(conc=f.030, time=time, n.tail=3)

}

\keyword{misc}
