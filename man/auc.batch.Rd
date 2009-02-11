\name{auc.batch}
\alias{auc.batch}
\title{ Estimation of the Area Under the Concentration Time Curve in Batch Designs}
\description{ Non-compartmental estimation of area under the concentration time curve (AUC) and area under the first moment curve (AUMC) in batch designs.}
\usage{auc.batch(conc, time, exact=NA, n.interpolate=2, n.tail=3)}

\arguments{
  \item{conc}{ A list of concentration levels. One element in the list should correspond to one batch. Multiple measurements per time point are expected to be adjacent to each other. }
  \item{time}{ A list of time points of concentration assessment. One element in the list should correspond to one batch. One time point for each concentration measured needs to be specified.}
  \item{exact}{ Time point for linear interpolation/extrapolation (default=\code{NA}). Currently not implemented. }
  \item{n.interpolate}{ Number of last data points used for linear interpolation/extrapolation (default=\code{2}). Currently not implemented. }
  \item{n.tail}{ Number of last data points used for tail area estimation (default=\code{3}). }
}

\details{
Estimation of the area under the concentration time curve (AUC 0-tlast) and area under the first moment curve (AUMC 0-tlast) for batch designs. In a serial sampling design only one measurement is available per subject at a specific time point, while in a batch design multiple time points are measured for each subject. AUC and AUMC (from 0 to the last time point) are calculated using the linear trapezoidal rule at the different time points.\cr\cr

Estimation of the linearly interpolated area under the concentration time curve (AUC 0-exact) and linearly interpolated area under the moment curve (AUMC 0-exact). The time point for linear interpolation/extrapolation must be after the last but one time point of concentration assessments. \cr

Estimation of the total area under the concentration curve (AUC 0-infinity) and total area under the first moment curve (AUMC 0-infinity) by using a tail area correction similar as suggested by Perrier and Gibaldi (Appendix D, 1982). \cr

Equal sample size per time point is required. 
}

\seealso{\code{\link{auc}}, \code{\link{auc.ci}}, \code{\link{auc.test}}.}

\value{An object of the class PK including AUC and AUMC estimates.}

\references{
Gibaldi M. and Perrier D. (1982). \emph{Pharmacokinetics. 2nd Edition}. Marcel Dekker, New York and Basel.  \cr 

Holder D. J., Hsuan F., Dixit R. and Soper K. (1999). A method for estimating and testing area under the curve in serial sacrifice, batch, and complete data designs. \emph{Journal of Biopharmaceutical Statistics}, 9(3):451-464.\cr

Jaki T. and Wolfsegger M. J. (In press). A theoretical framework for estimation of AUCs in complete and incomplete sampling designs. \emph{Statistics in Biopharmaceutical Research}. \cr
}

\author{Thomas Jaki}

\examples{
## a batch design from Jaki and Wolfsegger (in press), originally in Holder et al. (1999).
conc <- list(batch1=c(0,0,0,1.75,2.2,1.58,4.63,2.99,1.52), 
             batch2=c(3.03,1.98,2.22,3.34,1.3,1.22),
             batch3=c(3.54,2.84,2.55,0.3,0.0421,0.231))
time <- list(batch1=c(0,0,0,1,1,1,6,6,6), batch2=c(2,2,2,10,10,10), batch3=c(4,4,4,24,24,24))
auc.batch(conc, time, n.tail=3)

}

\keyword{misc}
