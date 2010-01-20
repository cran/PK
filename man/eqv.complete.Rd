\name{eqv.complete}
\alias{eqv.complete}
\title{Establishing Bioequivalence in Complete Data Designs}
\description{Confidence intervals for the ratio of independent or dependent area under the concentration versus time curves (AUCs) to the last time point in complete data designs.}
\usage{
eqv.complete(conc, time, group, dependent=FALSE, 
              method=c("fieller", "z", "boott"),
              conf.level=0.90, nsample=1000, data)
}

\arguments{
  \item{conc}{ Levels of concentrations as a vector. }
  \item{time}{ Time points of concentration assessment as a vector. One time point for each concentration measured needs to be specified.}
  \item{group}{ A grouping variable as a vector. }
  \item{dependent}{ Logical variable indicating if concentrations are measured on the same subjects for both AUCs (default=\code{FALSE}).}
  \item{method}{A character string specifying the method for calculation of confidence intervals (default=\code{c("fieller", "z", "boott")}).}
  \item{conf.level}{Confidence level (default=\code{0.90}). }
  \item{nsample}{ Number of bootstrap iterations for method \code{boott} (default=\code{1000}). }
  \item{data}{Optional data frame containing variables named as \code{id}, \code{conc}, \code{time} and \code{group}.}
}

\details{
This function computes confidence intervals for the ratio of (independent or dependent) AUCs (from 0 to the last time point) in complete data designs. It does so by treating the complete data design as a batch design with a single batch. More information can therefore be found under \code{\link{eqv}} and \code{\link{eqv.batch}}.\cr\cr

The above approach, though correct, is often inefficient. A general implementation is not provided as the most efficient analysis strongly depends on the context. The interested reader is refered to chapter 8 of Cawello (2003) while some examples can be found under \code{\link{auc.complete}}.\cr\cr

If \code{data} is specified the variable names \code{id}, \code{conc}, \code{time} and \code{group} are required and represent the corresponding variables as well as the subject id.
}

\seealso{\code{\link{eqv}}, \code{\link{auc.complete}}, \code{\link{estimator}}, \code{\link{ci}} and \code{\link{test}}.}

\value{An object of the class PK  containing the following components: \cr 
  \item{est}{Point estimates.}
  \item{CIs}{Point estimates, standard errors and confidence intervals. }
  \item{conc}{Levels of concentrations. } 
  \item{conf.level}{Confidence level.}
  \item{design}{Sampling design used.}
  \item{group}{Grouping variable.}
  \item{time}{Time points measured.}}

\references{
Cawello W. (2003). \emph{Parameters for Compartment-free Pharmacokinetics. Standardisation of Study Design, Data Analysis and Reporting}. Shaker Verlag, Aachen. \cr\cr

Hand D. and Crowder M. (1996). \emph{Practical Longitudinal Data Analysis}, Chapman and Hall, London.\cr
}

\author{Thomas Jaki}

\examples{
## example of a complete data design from Hand and Crowder (1996)
data(Glucose)
set.seed(271184)
eqv.complete(conc=Glucose$conc, time=Glucose$time, group=Glucose$date, 
     dependent=TRUE, method=c("fieller", "z", "boott"), 
     nsample=500, conf.level=0.90)
}


\keyword{htest}
