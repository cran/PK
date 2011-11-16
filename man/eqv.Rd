\name{eqv}
\alias{eqv}
\title{Bioequivalence between AUCs}
\description{Confidence intervals for the ratio of independent or dependent area under the concentration versus time curves (AUCs) to the last time point.}
\usage{
eqv(conc, time, group, dependent=FALSE, method=c("fieller", "z", "boott"), 
     conf.level=0.90, strata=NULL, nsample=1000, 
     design=c("ssd","batch","complete"), data)

}

\arguments{
  \item{conc}{ Levels of concentrations. For batch designs a list is required, while a vector is expected otherwise.}
  \item{time}{ Time points of concentration assessment. For batch designs a list is required, while a vector is expected otherwise. One time point for each concentration measured needs to be specified.}
  \item{group}{ A grouping variable. For batch designs a list is required, while a vector is expected otherwise.}
  \item{dependent}{ Logical variable indicating if concentrations are measured on the same subjects for both AUCs (default=\code{FALSE}).}
  \item{method}{A character string specifying the method for calculation of confidence intervals (default=\code{c("fieller", "z", "boott")}). }
  \item{conf.level}{ Confidence level (default=\code{0.90}). }
  \item{strata}{ A vector of one strata variable (default=\code{NULL}). Only available for method \code{boott} in a serial sampling design. }
  \item{nsample}{ Number of bootstrap iterations for method \code{boott} (default=\code{1000}). }
  \item{design}{ A character string indicating the type of design used. Possible values are \code{ssd} (the default) for a serial sampling design, \code{batch} for a batch design and \code{complete} for a complete data design. }
  \item{data}{Optional data frame containing variables named as \code{id}, \code{conc}, \code{time} and \code{group}.}
}

\details{
Calculation of confidence intervals for the ratio of (independent or dependent) AUCs (from 0 to the last time point) for serial sampling, batch and complete data designs. In a serial sampling design only one measurement is available per subject, while in a batch design multiple time points are measured for each subject. In a complete data design measurements are taken for all subjects at all time points. The AUC (from 0 to the last time point) is calculated using the linear trapezoidal rule on the arithmetic means at the different time points. \cr\cr

\code{dependent} specifies if the AUCs are dependent, that is measured on the same subjects. If \code{FALSE}, the intervals are based on the method of Jaki \emph{et al.} (2009) for the serial sampling design and on Jaki \emph{et al.} (in press) for the batch design. For dependent AUCs the method of Wolfsegger and Jaki (in press), which assumes that animals, batches and time points are equal for both AUCs, is used. Note that the option \code{dependent} is not used in serial sampling designs as by definition only one sample is obtained per subject then.\cr\cr

The \code{fieller} method is based on Fieller's theorem (1954) which uses the asymptotic standard errors of the individual AUCs and a critical value from a t-distribution with Satterthwaite's approximation (1946) to the degrees of freedom for calculation of confidence intervals. The \code{z} method is based on the limit distribution for the ratio using the critical value from a normal distribution for calculation of confidence intervals. \cr\cr

The \code{boott} method uses the asymptotic standard errors of the ratio of two AUCs while the critical value is obtained by the bootstrap-\emph{t} approach and follows the idea discussed in the context of serial sampling designs in Jaki T. \emph{et al.} (2009).  An equivalent approach is used in batch designs as well. Using \code{boott} an additional strata variable for bootstrapping can be specified in serial sampling designs.\cr\cr

If \code{data} is specified the variable names \code{conc}, \code{time} and \code{group} are required and represent the corresponding variables. If \code{design} is \code{batch} an additional variable \code{id} is required to identify the subject.
}

\seealso{\code{\link{auc}}, \code{\link{estimator}}, \code{\link{ci}} and \code{\link{test}}.}

\value{An object of the class PK  containing the following components: \cr 
  \item{est}{Point estimates.}
  \item{CIs}{Point estimates, standard errors and confidence intervals. }
  \item{conc}{Levels of concentrations. } 
  \item{conf.level}{Confidence level.}
  \item{design}{Sampling design used.}
  \item{group}{Grouping variable.}
  \item{time}{Time points measured.}}

\note{This is a wrapper function for \code{\link{eqv.complete}, \link{eqv.batch}} and \code{\link{eqv.ssd}}. 
The function calculates point and interval estimates for the ratio of AUCs (from 0 to the last time point).}

\references{
Fieller E. C. (1954). Some problems in interval estimation. \emph{Journal of the Royal Statistical Society, Series B}, 16:175-185. \cr\cr

Hand D. and Crowder M. (1996). \emph{Practical Longitudinal Data Analysis}, Chapman and Hall, London.\cr\cr

Jaki T., Wolfsegger M. J. and Ploner M. (2009). Confidence intervals for ratios of AUCs in the case of serial sampling: A comparison of seven methods. \emph{Pharmaceutical Statistics}, 8(1):12-24. \cr\cr

Jaki T., Wolfsegger M. J. and Lawo J-P. (in press). Establishing bioequivalence in complete and incomplete data designs using AUCs. \emph{Journal of Biopharmaceutical Statistics}. \cr\cr

Nedelman J. R., Gibiansky E. and Lau D. T. W. (1995). Applying Bailer's method for AUC confidence intervals to sparse sampling. \emph{Pharmaceutical Research}, 12(1):124-128. \cr\cr

Satterthwaite F. E. (1946). An approximate distribution of estimates of variance components. \emph{Biometrics Bulletin}, 2:110-114.  \cr\cr

Wolfsegger M. J. and Jaki T. (2009) Assessing systemic drug exposure in repeated dose toxicity studies in the case of complete and incomplete sampling. \emph{Biometrical Journal}, 51(6):1017:1029.\cr\cr

Yeh, C. (1990). Estimation and significant tests of area under the curve derived from incomplete blood sampling. \emph{ASA Proceedings of the Biopharmaceutical Section}, 74-81.\cr

}

\author{Thomas Jaki}

\examples{
## example of a serial sampling design from Nedelman et al. (1995)
data(CPI975)
data <- subset(CPI975,dose>=30)

data$concadj <- data$conc / data$dose

# fieller and z-interval for ratio of dose-normalized AUCs
eqv(conc=data$concadj, time=data$time, group=data$dose, method=c("z","fieller"), 
     design="ssd")

# bootstrap-t interval for ratio of dose-normalized AUCs stratified for sex
set.seed(310578)
eqv(conc=data$concadj, time=data$time, group=data$dose, method="boott",
        strata=data$sex, nsample=500, design="ssd")


## Example of an independent batch design from Yeh (1990) 
conc <- list(batch1=c(0,0,0,0,0,0, 4.690,2.070,6.450,0.1,0.852,0.136, 
                      4.690,4.060,6.450,0.531,1.2,0.607),
           batch2=c(4,1.3,3.2,0.074,0.164,0.267, 6.68,3.83,6.08,0.669,1.21,0.878, 
                    8.13,9.54,6.29,0.923,1.65,1.04),
           batch3=c(9.360,13,5.48,1.090,1.370,1.430, 5.180,5.180,2.79,0.804,1.47,1.26, 
                    1.060,2.15,0.827,0.217,0.42,0.35))  
time <- list(batch1=c(rep(0,6),rep(1,6),rep(4,6)),
             batch2=c(rep(0.5,6),rep(2,6),rep(6,6)),
             batch3=c(rep(8,6),rep(12,6),rep(24,6)))
group <- list(batch1=rep(rep(c(1,2),each=3),3), batch2=rep(rep(c(1,2),each=3),3),
              batch3=rep(rep(c(1,2),each=3),3))

eqv(conc=conc, time=time, group=group, dependent=FALSE, method=c("fieller"), 
    conf.level=0.90, design="batch")


## example for a dependent batch data design from Wolfsegger and Jaki (2009)
conc <- list(batch1=c(0.46,0.2,0.1,0.1, 1.49,1.22,1.27,0.53, 0.51,0.36,0.44,0.28),
             batch2=c(1.51,1.80,2.52,1.91, 0.88,0.66,0.96,0.48),
             batch3=c(1.52,1.46,2.55,1.04, 0.54,0.61,0.55,0.27))
time <- list(batch1=c(0,0,0,0,1.5,1.5,1.5,1.5,10.5,10.5,10.5,10.5),
             batch2=c(5/60,5/60,5/60,5/60,4,4,4,4),
             batch3=c(0.5,0.5,0.5,0.5,7,7,7,7))
group <- list(batch1=c(1,1,2,2,1,1,2,2,1,1,2,2), batch2=c(1,1,2,2,1,1,2,2),
              batch3=c(1,1,2,2,1,1,2,2))

eqv(conc=conc, time=time, group=group, dependent=TRUE, method=c("fieller"), 
     conf.level=0.90, design="batch")


## example of a complete data design from Hand and Crowder (1996)
data(Glucose)
set.seed(271184)
eqv(conc=Glucose$conc, time=Glucose$time, group=Glucose$date, 
     dependent=TRUE, method=c("fieller", "z", "boott"), conf.level=0.90, 
     nsample=500, design="complete")
}

\keyword{htest}
