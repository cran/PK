\name{auc.complete}
\alias{auc.complete}
\title{Confidence intervals for the area under the concentration versus time curve in complete data designs}
\description{Examples to find confidence intervals for the area under the concentration versus time curve (AUC) in complete data designs.}
\usage{
auc.complete(conc, time, group=NULL, method=c("t", "z", "boott"), 
        alternative=c("two.sided", "less", "greater"), 
        conf.level=0.95, nsample=1000, data) }

\arguments{
  \item{conc}{ Levels of concentrations as a vector. }
  \item{time}{ Time points of concentration assessment as a vector. One time point for each concentration measured needs to be specified.}
  \item{group}{ A grouping variable as a vector (default=\code{NULL}). If specified, a confidence interval for the difference of independent AUCs will be calculated. }
  \item{method}{ A character string specifying the method for calculation of confidence intervals (default=\code{c("t", "z", "boott")}). }
  \item{alternative}{ A character string specifying the alternative hypothesis. Possible values are \code{"less"}, \code{"greater"} and \code{"two.sided"} (the default). } 
  \item{conf.level}{ Confidence level (default=\code{0.95}). }
  \item{nsample}{ Number of bootstrap iterations for method \code{boott} (default=\code{1000}). } 
  \item{data}{Optional data frame containing variables named as \code{conc}, \code{time} and \code{group}.}
}

\details{
This function computes confidence intervals for an AUC (from 0 to the last time point) or for the difference between two AUCs in complete data designs. \cr\cr

To compute confidence intervals in complete data designs the design is treated as a batch design with a single batch. More information can therefore be found under \code{\link{auc}}. A corresponding reminder message is produced if confidence intervals can be computed, ie when at least 2 measurements at each time point are available.\cr\cr

The above approach, though correct, is often inefficient and so we will illustrate alternative methods in this help file. A general implementation is not provided as the most efficient analysis strongly depends on the context. The interested reader is refered to chapter 8 of Cawello (2003).\cr\cr

If \code{data} is specified the variable names \code{conc}, \code{time} and \code{group} are required and represent the corresponding variables. 
}

\seealso{\code{\link{auc}}, \code{\link{estimator}}, \code{\link{ci}} and \code{\link{test}}}

\value{An object of the class PK  containing the following components: \cr 
  \item{est}{Point estimates.}
  \item{CIs}{Point estimates, standard errors and confidence intervals. }
  \item{conc}{Levels of concentrations. } 
  \item{conf.level}{Confidence level.}
  \item{design}{Sampling design used.}
  \item{group}{Grouping variable.}
  \item{time}{Time points measured.}
}

\references{
Cawello W. (2003). \emph{Parameters for Compartment-free Pharmacokinetics. Standardisation of Study Design, Data Analysis and Reporting}. Shaker Verlag, Aachen. \cr\cr

Gibaldi M. and Perrier D. (1982). \emph{Pharmacokinetics}. Marcel Dekker, New York and Basel.\cr
}

\author{Thomas Jaki and Martin Wolfsegger}

\examples{
## example from Gibaldi and Perrier (1982, page 436) for an individual AUC
time <- c(0, 0.165, 0.5, 1, 1.5, 3, 5, 7.5, 10)
conc <- c(0, 65.03, 28.69, 10.04, 4.93, 2.29, 1.36, 0.71, 0.38)
auc.complete(conc=conc, time=time)

## dataset Indometh of package datasets
## calculate individual AUCs
require(datasets)
row <- 1
res <- data.frame(matrix(nrow=length(unique(Indometh$Subject)), ncol=2))
colnames(res) <- c('id', 'auc')
for(i in unique(Indometh$Subject)){
   temp <- subset(Indometh, i==Subject)
   res[row, 1] <- i
   res[row, 2] <- auc.complete(data=temp[,c("conc","time")])$est[1,1]
   row <- row + 1
}
print(res)

# function to get geometric mean and corresponding CI
gm.ci <- function(x, conf.level=0.95){
   res <- t.test(x=log(x), conf.level=conf.level)
   out <- data.frame(gm=as.real(exp(res$estimate)), lower=exp(res$conf.int[1]), 
                      upper=exp(res$conf.int[2]))
   return(out)
}    

# geometric mean and corresponding CI: assuming log-normal distributed AUCs
gm.ci(res[,2], conf.level=0.95)

# arithmetic mean and corresponding CI: assuming normal distributed AUCs 
# or at least asymptotic normal distributed arithmetic mean 
t.test(x=res[,2], conf.level=0.95)
     
# alternatively: function auc.complete
set.seed(300874)
Indometh$id <- as.character(Indometh$Subject)
Indometh <- Indometh[order(Indometh$id, Indometh$time),]
Indometh <- Indometh[order(Indometh$time),]
auc.complete(conc=Indometh$conc, time=Indometh$time)


## example for comparing AUCs assessed in a repeated complete data design
## (dataset: Glucose)
## calculate individual AUCs
data(Glucose)
res <- data.frame(matrix(nrow=length(unique(Glucose$id))*2, ncol=3))
colnames(res) <- c('id', 'date', 'auc')
row <- 1
for(i in unique(Glucose$id)){
  for(j in unique(Glucose$date)){
     temp <- subset(Glucose, id==i & date==j)
     res[row, c(1,2)] <- c(i,j)
     res[row, 3] <- auc.complete(data=temp[,c("conc","time")])$est[1,1]
     row <- row + 1
  }
}
res <- res[order(res$id, res$date),]
print(res)

# assuming log-normally distributed AUCs
# geometric means and corresponding two-sided CIs per date           
tapply(res$auc, res$date, gm.ci)

# comparison of AUCs using ratio of geometric means and corresponding two-sided CI 
# repeated experiment
model <- t.test(log(auc)~date, data=res, paired=TRUE, conf.level=0.90)
exp(as.real(model$estimate))
exp(model$conf.int)
}

\keyword{htest}
