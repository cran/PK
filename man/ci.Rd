\name{ci}
\alias{ci}
\title{Function to extract confidence interval(s)}
\description{Generic function that extracts the confidence interval(s) of an object of class PK.}
\usage{
ci(obj, method=NULL) }

\arguments{
  \item{obj}{ An output object of class PK.}
  \item{method}{A character string specifying the method of the confidence interval. If \code{NULL} (default) all intervals are returned.}
}

\details{
  Generic function to allow easy extraction of confidence intervals.
}

\seealso{\code{\link{estimator}} and \code{\link{test}}}

\value{A matrix containing confidence interval bounds.}

\author{Thomas Jaki}

\references{
Nedelman J. R., Gibiansky E. and Lau D. T. W. (1995). Applying Bailer's method for AUC confidence intervals to sparse sampling. \emph{Pharmaceutical Research}, 12(1):124-128. \cr
}


\examples{
# Example from Nedelman et al. (1995)
conc <- c(2790, 3280, 4980, 7550, 5500, 6650, 2250, 3220, 213, 636)
time <- c(1, 1, 2, 2, 4, 4, 8, 8, 24, 24)

obj <- auc(conc=conc, time=time, method=c("z", "t"), design="ssd")

## all requested ci's
ci(obj)

## a specific ci
ci(obj, method="t")
}

\keyword{manip}
