\name{estimator}
\alias{estimator}
\title{Function to extract point estimate(s)}
\description{Generic function that extracts the point estimate(s) of an object of class PK.}
\usage{
estimator(obj,se=FALSE) }

\arguments{
  \item{obj}{ An output object of class PK.}
  \item{se}{ Logical variable indicating if the standard error should be provided as well (default=\code{FALSE}).}
}

\details{
  Generic function to allow easy extraction of point estimates.
}

\seealso{\code{\link{ci}} and \code{\link{test}}}

\value{A matrix containing the point estimate(s) and optionally the standard error(s).}

\author{Thomas Jaki}

\references{
Nedelman J. R., Gibiansky E. and Lau D. T. W. (1995). Applying Bailer's method for AUC confidence intervals to sparse sampling. \emph{Pharmaceutical Research}, 12(1):124-128. \cr
}

\examples{
# Example from Nedelman et al. (1995)
conc <- c(2790, 3280, 4980, 7550, 5500, 6650, 2250, 3220, 213, 636)
time <- c(1, 1, 2, 2, 4, 4, 8, 8, 24, 24)

obj <- auc(conc=conc, time=time, method=c('z', 't'), design='ssd')

estimator(obj,TRUE)
}

\keyword{manip}
