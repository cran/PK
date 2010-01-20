\name{auc}
\alias{auc}
\alias{auc.ssd}
\alias{auc.batch}
\title{Estimation of confidence intervals for the area under the concentration versus time curve in complete and incomplete data designs}
\description{Calculation of confidence intervals for an area under the concentration versus time curve (AUC) or for the difference between two AUCs assessed in complete and incomplete data designs.}
\usage{
auc(conc, time, group=NULL, method=c("t", "z", "boott"), 
     alternative=c("two.sided", "less", "greater"), 
     conf.level=0.95, strata=NULL, nsample=1000, 
     design=c("ssd","batch","complete"), data)	     

auc.ssd(conc, time, group=NULL, method=c("t", "z", "boott"), 
     alternative=c("two.sided", "less", "greater"), 
     conf.level=0.95, strata=NULL, nsample=1000, data)	     

auc.batch(conc, time, group=NULL, method=c("t", "z", "boott"), 
     alternative=c("two.sided", "less", "greater"), 
     conf.level=0.95, nsample=1000, data)	     

}

\arguments{
  \item{conc}{ Levels of concentrations. For batch designs a list is required, while a vector is expected otherwise.}
  \item{time}{ Time points of concentration assessment. For batch designs a list is required, while a vector is expected otherwise. One time point for each concentration measured needs to be specified.}
  \item{group}{ A grouping variable (default=\code{NULL}). For batch designs a list is required, while a vector is expected otherwise. If specified, a confidence interval for the difference of independent AUCs will be calculated. }
  \item{method}{ A character string specifying the method for calculation of confidence intervals (default=\code{c("t", "z", "boott")}). }
  \item{alternative}{ A character string specifying the alternative hypothesis. Possible values are \code{"less"}, \code{"greater"} and \code{"two.sided"} (the default).} 
  \item{conf.level}{ Confidence level (default=\code{0.95}). }
  \item{strata}{ A vector of one strata variable (default=\code{NULL}). Only available for method \code{boott} in a serial sampling design. }
  \item{nsample}{ Number of bootstrap iterations for method \code{boott} (default=\code{1000}). } 
  \item{design}{ A character string indicating the type of design used. Possible values are \code{"ssd"} for a serial sampling design, \code{"batch"} for a batch design and \code{"complete"} for a complete data design. }
  \item{data}{Optional data frame containing variables named as \code{id}, \code{conc}, \code{time} and \code{group}.}
}

\details{
Calculation of confidence intervals for an AUC (from 0 to the last time point) or for the difference between two AUCs for serial sampling, batch and complete data designs. In a serial sampling design only one measurement is available per subject, while in a batch design multiple time points are measured for each subject. In a complete data design measurements are taken for all subjects at all time points. The AUC (from 0 to the last time point) is calculated using the linear trapezoidal rule on the arithmetic means at the different time points.\cr\cr

If group=NULL a confidence interval for an AUC is calculated. If group specifies a factor variable with exactly two levels, a confidence interval for the difference between two independent AUCs is calculated. To obtain confidence intervals for dependent AUCs simply use the difference in concentrations for \code{conc}. See the example below.\cr\cr

The \code{t} method uses the critical value from a t-distribution with Satterthwaite's approximation (Satterthwaite, 1946) to the degrees of freedom for calculation of confidence intervals as presented in Tang-Liu and Burke (1988), Nedelman et al (1995), Holder et al (1999) and in Jaki and Wolfsegger (2009). The \code{z} method uses the critical value from a normal distribution for calculation of confidence intervals as presented in Bailer (1988) or in Jaki and Wolfsegger (2009). The \code{boott} method uses bootstrap-\emph{t} confidence intervals as presented in Jaki and Wolfsegger (2009). Using \code{boott} an additional strata variable for bootstrapping can be specified in the case of serial sampling. \cr\cr

For serial sampling designs missing data are omitted and unequal sample sizes per time point are allowed. For batch designs missing values are not permitted and equal sample size per time point is required.\cr\cr

If \code{data} is specified the variable names \code{conc}, \code{time} and \code{group} are required and represent the corresponding variables. If \code{design} is \code{batch} an additional variable \code{id} is required to identify the subject.\cr\cr

NOTE: Confidence intervals for AUCs assessed in complete data designs are found using a batch design with one batch based on the asymptotic normal distribution. Conventionally, AUCs are assumed to be log-normal distributed. See the help file \link{auc.complete} for some corresponding examples.
}

\seealso{\code{\link{auc.complete}}, \code{\link{nca}}, \code{\link{eqv}}, \code{\link{estimator}}, \code{\link{ci}} and \code{\link{test}}.}

\value{An object of the class PK  containing the following components: \cr 
  \item{est}{Point estimates.}
  \item{CIs}{Point estimates, standard errors and confidence intervals. }
  \item{conc}{Levels of concentrations. } 
  \item{conf.level}{Confidence level.}
  \item{design}{Sampling design used.}
  \item{group}{Grouping variable.}
  \item{time}{Time points measured.}
}

\note{This is a wrapper function for \code{\link{auc.complete}, \link{auc.batch}} and \code{\link{auc.ssd}}. 
The function calculates point and interval estimates for AUC (from 0 to the last time point).}

\references{
Bailer A. J. (1988). Testing for the equality of area under the curves when using destructive measurement techniques. \emph{Journal of Pharmacokinetics and Biopharmaceutics}, 16(3):303-309. \cr\cr

Holder D. J., Hsuan F., Dixit R. and Soper K. (1999). A method for estimating and testing area under the curve in serial sacrifice, batch, and complete data designs. \emph{Journal of Biopharmaceutical Statistics}, 9(3):451-464.\cr\cr

Jaki T. and Wolfsegger M. J. (2009). A theoretical framework for estimation of AUCs in complete and incomplete sampling designs. \emph{Statistics in Biopharmaceutical Research}, 1(2):176-184. \cr\cr

Nedelman J. R., Gibiansky E. and Lau D. T. W. (1995). Applying Bailer's method for AUC confidence intervals to sparse sampling. \emph{Pharmaceutical Research}, 12(1):124-128. \cr\cr

Satterthwaite F. E. (1946). An approximate distribution of estimates of variance components. \emph{Biometrics Bulletin}, 2:110-114.  \cr\cr

Tang-Liu D. D.-S. and Burke P. J. (1988). The effect of azone on ocular levobunolol absoprtion: Calculating the area under the curve and its standard error using tissue sampling compartments. \emph{Pharmaceutical Research}, 5(4):238-241. \cr\cr

Wolfsegger M. J. and Jaki T. (2009) Assessing systemic drug exposure in repeated dose toxicity studies in the case of complete and incomplete sampling. \emph{Biometrical Journal}, 51(6):1017:1029.\cr
}

\author{Thomas Jaki and Martin J. Wolfsegger}

\examples{
## example from Bailer (1988)
time <- c(rep(0,4), rep(1.5,4), rep(3,4), rep(5,4), rep(8,4))
grp1 <- c(0.0658, 0.0320, 0.0338, 0.0438, 0.0059, 0.0030, 0.0084,
          0.0080, 0.0000, 0.0017, 0.0028, 0.0055, 0.0000, 0.0037,
          0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000)

grp2 <- c(0.2287, 0.3824, 0.2402, 0.2373, 0.1252, 0.0446, 0.0638,
          0.0511, 0.0182, 0.0000, 0.0117, 0.0126, 0.0000, 0.0440,
          0.0039, 0.0040, 0.0000, 0.0000, 0.0000, 0.0000)

grp3 <- c(0.4285, 0.5180, 0.3690, 0.5428, 0.0983, 0.0928, 0.1128,
          0.1157, 0.0234, 0.0311, 0.0344, 0.0349, 0.0032, 0.0052,
          0.0049, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000)

auc(conc=grp1, time=time, method='z', design='ssd')
auc(conc=grp2, time=time, method='z', design='ssd')
auc(conc=grp3, time=time, method='z', design='ssd')

## function call with data frame using simultaneous confidence intervals based 
## on bonferroni adjustment
data <- data.frame(conc=c(grp1, grp2, grp3), time=rep(time, 3),
                   group=c(rep(1, length(grp1)), rep(2, length(grp2)), 
                   rep(3, length(grp3))))

auc(subset(data, group==1 | group==2)$conc, subset(data, group==1 | group==2)$time, 
    group=subset(data, group==1 | group==2)$group, method=c('z', 't'), 
    conf.level=1-0.05/3, design='ssd')

auc(subset(data, group==1 | group==3)$conc, subset(data, group==1 | group==2)$time, 
    group=subset(data, group==1 | group==3)$group, method=c('z', 't'), 
    conf.level=1-0.05/3, design='ssd')

auc(subset(data, group==2 | group==3)$conc, subset(data, group==1 | group==2)$time, 
    group=subset(data, group==2 | group==3)$group, method=c('z', 't'), 
    conf.level=1-0.05/3, design='ssd')

## example from Nedelman et al. (1995)
data(CPI975)
data <- CPI975[CPI975[,'dose']>=30 ,]

auc(data=subset(data,sex=='m' & dose==30), method=c('z', 't'), design='ssd')
auc(data=subset(data,sex=='f' & dose==30), method=c('z', 't'), design='ssd')

auc(data=subset(data,sex=='m' & dose==100), method=c('z', 't'), design='ssd')
auc(data=subset(data,sex=='f' & dose==100), method=c('z', 't'), design='ssd')

## comparing dose levels
data$concadj <- data$conc / data$dose
data.100 <- subset(data, dose==100)
data.030 <- subset(data, dose==30)
res.100 <- auc(conc=data.030$concadj, time=data.030$time, method='t', design='ssd')
res.030 <- auc(conc=data.100$concadj, time=data.100$time, method='t', design='ssd')
plot(res.030, ylim=c(0, 140), xlim=c(0,25), pch=19, ylab='Dose-normalized concentration', 
      main='Comparison of doses')
plot(res.100, col='red', pch=21, add=TRUE)
legend(x=25, y=140, xjust=1, lty=1, col=c('black','red'), 
       legend=c('Dose of 30', 'Dose of 100'))

res <- auc(conc=data$concadj, time=data$time, group=data$dose, method=c('t','z'), 
            design='ssd')
print(res)
summary(res)

## comparing two dose level using stratified resampling per gender
## caution this might take a few minutes
set.seed(260151)
auc(conc=data$concadj, time=data$time, group=data$dose, method='boott',
    strata=data$sex, design='ssd', nsample=500)

## a batch design example from Holder et al. (1999).
data(Rats)
data <- subset(Rats,Rats$dose==100)

# two-sided CI: data call
auc(data=data,method=c('z','t'), design='batch')
# one-sided CI: data call
auc(data=data,method=c('z','t'), alternative="less", design='batch')

## difference of two AUCs in batch design from Jaki and Wolfsegger (2009),
## originally in Holder et al. (1999).
data <- subset(Rats,Rats$dose==100 | Rats$dose==300 )
data$group <- data$dose
data$conc <- data$conc / data$dose

## data call
res1 <- auc(data=subset(data, dose==100), method='z', design='batch')
res2 <- auc(data=subset(data, dose==300), method='z', design='batch')
plot(res1, col='black', ylim=c(0,0.06), xlim=c(0,25), ylab='Dose-normalized concentration', 
      main='Comparison of doses')
plot(res2, col='red', add=TRUE)
legend(x=0, y=0.06, lty=1, col=c('black','red'), 
       legend=c('Dose of 100', 'Dose of 300'))

auc(data=data, method='z', design='batch')


## difference of two dependent AUCs in a batch design from Wolfsegger and Jaki (2009)
conc <- list(batch1=c(0.46,0.2,0.1,0.1, 1.49,1.22,1.27,0.53, 0.51,0.36,0.44,0.28),
             batch2=c(1.51,1.80,2.52,1.91, 0.88,0.66,0.96,0.48),
             batch3=c(1.52,1.46,2.55,1.04, 0.54,0.61,0.55,0.27))
time <- list(batch1=c(0,0,0,0,1.5,1.5,1.5,1.5,10.5,10.5,10.5,10.5),
             batch2=c(5/60,5/60,5/60,5/60,4,4,4,4),
             batch3=c(0.5,0.5,0.5,0.5,7,7,7,7))
group <- list(batch1=c(1,1,2,2,1,1,2,2,1,1,2,2),batch2=c(1,1,2,2,1,1,2,2),
              batch3=c(1,1,2,2,1,1,2,2))

# find difference in concentration and the corresponding times
dconc <- NULL
dtime <- NULL
grps <- unique(unlist(group))
B <- length(conc)
for(i in 1:B){
    dconc[[i]] <- conc[[i]][group[[i]]==grps[1]] - conc[[i]][group[[i]]==grps[2]]
    dtime[[i]] <- time[[i]][group[[i]]==grps[1]]
}
names(dconc) <- names(conc)

auc(conc=dconc, time=dtime, group=NULL, method="t", conf.level=0.90, design="batch")


## complete data design: example
## data Indometh
require(datasets)
Indometh$id <- as.character(Indometh$Subject)
Indometh <- Indometh[order(Indometh$id, Indometh$time),]
Indometh <- Indometh[order(Indometh$time),]
res <- auc.complete(conc=Indometh$conc, time=Indometh$time, method='t')
plot(res)

## more informative plot 
split.screen(c(1,2))
screen(1)
plot(x=c(0,8), y=c(0, 3), type='n', main='Observed concentration time-profiles', 
     xlab='Time', ylab='Concentration', las=1)
for(i in unique(Indometh$Subject)){
   temp <- subset(Indometh, Subject==i)
   points(x=temp$time, y=temp$conc, type='b')
}
screen(2)
plot(x=c(0,8), y=c(0.01, 9), type='n', main='Log-linear concentration time-profiles', 
     xlab='Time', ylab='Log of concentration', yaxt='n', log='y')
axis(side=2, at=c(0.01, 0.1, 1, 10), labels=c('0.01', '0.1', '1', '10'), las=1)
axis(side=2, at=seq(0.01, 0.1, 0.01), tcl=-0.2, labels=FALSE)
axis(side=2, at=seq(0.1, 1, 0.1), tcl=-0.2, labels=FALSE)
axis(side=2, at=seq(1, 10, 1), tcl=-0.2, labels=FALSE)
for(i in unique(Indometh$Subject)){
   temp <- subset(Indometh, Subject==i)
   points(x=temp$time, y=temp$conc, type='b')
}
close.screen(all = TRUE)
}

\keyword{htest}
