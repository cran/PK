\name{auc.complete}
\alias{auc.complete}
\title{Estimation of the Area Under the Concentration Time Curve in Complete Data Designs}
\description{Non-compartmental estimation of area under the concentration time curve (AUC) and area under the first moment curve (AUMC) in complete data designs.}
\usage{auc.complete(conc, time, exact=NA, n.interpolate=2, n.tail=3) }

\arguments{
  \item{conc}{ A vector of concentration levels. }
  \item{time}{ A vector of time points of concentration assessment. One time point for each concentration measured needs to be specified.}
  \item{exact}{ Time point for linear interpolation/extrapolation (default=\code{NA}). }
  \item{n.interpolate}{ Number of last data points used for linear interpolation/extrapolation (default=\code{2}). }
  \item{n.tail}{ Number of last data points used for tail area estimation (default=\code{3}). }
}

\details{
Estimation of the area under the concentration time curve (AUC 0-tlast) and area under the first moment curve (AUMC 0-tlast) for complete data designs. AUC and AUMC (from 0 to the last time point) are calculated using the linear trapezoidal rule at the different time points.\cr\cr

Estimation of the linearly interpolated area under the concentration time curve (AUC 0-exact) and linearly interpolated area under the moment curve (AUMC 0-exact). The time point for linear interpolation/extrapolation must be after the last but one time point of concentration assessments. \cr

Estimation of the total area under the concentration curve (AUC 0-infinity) and total area under the first moment curve (AUMC 0-infinity) by using a tail area correction similar as suggested by Perrier and Gibaldi (Appendix D, 1982). \cr\cr

}

\seealso{\code{\link{auc}}, \code{\link{auc.ci}}, \code{\link{auc.test}}.}


\value{An object of the class PK including AUC and AUMC estimates.}

\references{
Cawello W. (2003). \emph{Parameters for Compartment-free Pharmacokinetics. Standardisation of Study Design, Data Analysis and Reporting}. Shaker Verlag, Aachen. \cr 

Gibaldi M. and Perrier D. (1982). \emph{Pharmacokinetics. 2nd Edition}. Marcel Dekker, New York and Basel.  \cr 

}

\author{Martin J. Wolfsegger and Thomas Jaki}

\examples{
## example for a complete data design from Cawello W. (2003, page 70 and 74)
time <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 6)
conc <- c(0, 5.67, 20.6, 28.7, 22.5, 17.4, 17.7, 13.4, 11.0, 8.23, 5.14, 2.84) 

# using a vector for both concentration and time
auc.complete(conc=conc, time=time, exact=7, n.tail=4)

}

\keyword{misc}
