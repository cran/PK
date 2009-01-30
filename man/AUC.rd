\name{auc}
\alias{auc}
\alias{print.PK}
\title{Estimation of the Area Under the Concentration Time Curve in Complete and Incomplete Data Designs}
\description{Non-compartmental estimation of area under the concentration time curve (AUC) and area under the first moment curve (AUMC).}
\usage{
auc(conc, time, exact=NA, n.interpolate=2, n.tail=3, 
     design=c('ssd','batch','complete'))}

\arguments{
  \item{conc}{ Levels of concentrations. For batch designs a list is required, while a vector is expected otherwise. Multiple measurements per time point are expected to be adjacent to each other.}
  \item{time}{ Time points of concentration assessment. For batch designs a list of batches containing time points for each batch is required. Otherwise a vector is required. One time point for each concentration measured needs to be specified.}
  \item{exact}{ Time point for linear interpolation/extrapolation (default=\code{NA}). Only implemented for a complete data design.}
  \item{n.interpolate}{ Number of last data points used for linear interpolation/extrapolation (default=\code{2}). Only implemented for a complete data design.}
  \item{n.tail}{ Number of last data points used for tail area correction (default=\code{3}). }
  \item{design}{ A character string indicating the type of design used. Possible values are \code{ssd} (the default) for a serial sampling design,
          \code{batch} for a batch design and \code{complete} for a complete data design. }
}

\details{
Estimation of the area under the concentration time curve (AUC 0-tlast) and area under the first moment curve (AUMC 0-tlast) for serial sampling, batch and complete data designs. In a serial sampling design only one measurement is available per subject at a specific time point, while in a batch design multiple time points are measured for each subject. In a complete data design all measurements are taken for all subjects at all time points. The AUC (from 0 to the last time point) is calculated using the linear trapezoidal rule on the arithmetic means at the different time points.\cr\cr

The total area under the concentration curve (AUC 0-infinity) and total area under the first moment curve (AUMC 0-infinity) is computed using a the tail area correction calculated similar as suggested by Perrier and Gibaldi (Appendix D, 1982). \cr

The linearly interpolated area under the concentration time curve (AUC 0-exact) and linearly interpolated area under the moment curve (AUMC 0-exact) is also estimated for complete data designs. The time point (\code{exact}) for linear interpolation/extrapolation must be after the second to last time point of concentration assessments. \code{n.interpolate} specifies the number of timepoints to be used for the linear interpolation.\cr

Equal sample size per time point is required for batch designs. 
}

\seealso{\code{\link{auc.ci}}, \code{\link{auc.test}}.}

\note{This is a wrapper function for \code{\link{auc.complete}}, \code{\link{auc.batch}} and \code{\link{auc.ssd}}. See the documentation of these functions for more details. }

\value{An object of the class PK including AUC and AUMC estimates.}

\references{
Cawello W. (2003). \emph{Parameters for Compartment-free Pharmacokinetics. Standardisation of Study Design, Data Analysis and Reporting}. Shaker Verlag, Aachen. \cr 

Gibaldi M. and Perrier D. (1982). \emph{Pharmacokinetics. 2nd Edition}. Marcel Dekker, New York and Basel.  \cr 

Holder D. J., Hsuan F., Dixit R. and Soper K. (1999). A method for estimating and testing area under the curve in serial sacrifice, batch, and complete data designs. \emph{Journal of Biopharmaceutical Statistics}, 9(3):451-464.\cr

Jaki T. and Wolfsegger M. J. (In press). A theoretical framework for estimation of AUCs in complete and incomplete sampling designs. \emph{Statistics in Biopharmaceutical Research}. \cr	

Nedelman J. R., Gibiansky E. and Lau D. T. W. (1995). Applying Bailer's method for AUC confidence intervals to sparse sampling. \emph{Pharmaceutical Research}, 12(1):124-128. \cr

}

\author{Thomas Jaki}

\examples{
## example for a complete data design from Cawello W. (2003, page 70 and 74)
time <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 6)
conc <- c(0, 5.67, 20.6, 28.7, 22.5, 17.4, 17.7, 13.4, 11.0, 8.23, 5.14, 2.84) 

# using a vector for both concentration and time
auc(conc=conc, time=time, exact=7, n.tail=4, design='complete')

## a batch design from Jaki and Wolfsegger (in press), originally in Holder et al. (1999).
conc <- list(batch1=c(0,0,0,1.75,2.2,1.58,4.63,2.99,1.52), 
             batch2=c(3.03,1.98,2.22,3.34,1.3,1.22),
             batch3=c(3.54,2.84,2.55,0.3,0.0421,0.231))
time <- list(batch1=c(0,0,0,1,1,1,6,6,6), batch2=c(2,2,2,10,10,10), 
             batch3=c(4,4,4,24,24,24))
auc(conc, time, n.tail=3, design='batch')

## example for a serial sampling design from Nedelman et al. (1995)
time <- c(1, 1, 2, 2, 4, 4, 8, 8, 24, 24)
m.030 <- c(391, 396, 649, 1990, 3290, 3820, 844, 1650, 75.7, 288)
f.030 <- c(353, 384, 625, 1410, 1020, 1500, 933, 1030, 0, 80.5)
auc(conc=m.030, time=time, n.tail=3, design='ssd')
auc(conc=f.030, time=time, n.tail=3, design='ssd')

}

\keyword{misc}
