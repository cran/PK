\name{AUC}
\alias{AUC}
\title{Area Under the Concentration Time Curve}
\description{Estimation of area under the concentration time curve (AUC) and area under the first moment curve (AUMC).}
\usage{AUC(conc, time, exact=NA, numintp=2, numtail=3, prev=0)	}

\arguments{
  \item{conc}{ levels of concentrations. }
  \item{time}{ time points of concentration assessment. }
  \item{exact}{ time point for linear interpolation/extrapolation. }
  \item{numintp}{ number of last data points used for linear interpolation/extrapolation. }
  \item{numtail}{ number of last data points used for tail area correction. }
  \item{prev}{ pre-dosing value. }
}

\details{
Estimation of the observed area under the concentration time curve (AUC 0-tlast) and observed area under the first moment curve (AUMC 0-tlast) using the linear trapezoidal method. \cr

Estimation of the linearly interpolated area under the concentration time curve (AUC 0-exact) and linearly interpolated area under the moment curve (AUMC 0-exact). The time point for linear interpolation/extrapolation must be after the last but one time point of concentration assessments. \cr

Estimation of the total area under the concentration curve (AUC 0-infinity) and total area under the first moment curve (AUMC 0-infinity) by using a tail area correction similar as suggested by Perrier and Gibaldi (Appendix D, 1982). \cr

If the pre-dosing value indicating that the intrinsic level is greater than 0, the pre-dosing value is subtracted from all concentration levels before calculation of AUC and AUMC.
}

\note{Records including missing values and values below zero are omitted. }
\value{Data frame including AUC and AUMC estimates.}

\references{
Cawello W. (2003). \emph{Parameters for Compartment-free Pharmacokinetics. Standardisation of Study Design, Data Analysis and Reporting}. Shaker Verlag, Aachen. \cr 

Gibaldi M. and Perrier D. (1982). \emph{Pharmacokinetics. 2nd Edition}. Marcel Dekker, New York and Basel.  \cr 
}

\author{Martin J. Wolfsegger and Thomas Jaki}

\examples{
## example from Cawello W. (2003, page 70 and 74)
time <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 6)
conc <- c(0, 5.67, 20.6, 28.7, 22.5, 17.4, 17.7, 13.4, 11.0, 8.23, 5.14, 2.84) 
AUC(conc=conc, time=time, exact=7, numtail=4)
}

\keyword{misc}