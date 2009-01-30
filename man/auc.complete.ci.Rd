\name{auc.complete.ci}
\alias{auc.complete.ci}
\title{Confidence intervals for the Area Under the Concentration Time Curve in Complete Data Designs.}
\description{Examples to find confidence intervals for the area under the concentration time curve (AUC) in complete data designs.}
\usage{
auc.complete.ci(conc, time) }

\arguments{
  \item{conc}{ Levels of concentrations. Not used. }
  \item{time}{ Time points of concentration assessment. Not used. }
}

\details{
This help file is for illustration of different methods to find confidence intervals for AUCs in complete data designs only. \cr\cr
For more information regarding analysis, presenting and reporting of PK parameters assessed in complete data designs refer for example to chapter 8 of Cawello (2003).
}

\seealso{\code{\link{auc}}, \code{\link{auc.test}} and \code{\link{auc.complete}}.}


\value{None.}

\references{
Cawello W. (2003). \emph{Parameters for Compartment-free Pharmacokinetics. Standardisation of Study Design, Data Analysis and Reporting}. Shaker Verlag, Aachen. \cr 
}

\author{Martin J. Wolfsegger}

\examples{
## example for comparing AUCs assessed in a repeated complete data design 
## dataset Glucose2 of package nlme
require(nlme)
Glucose2 <- Glucose2[order(Glucose2$Subject, Glucose2$Date, Glucose2$Time),]

# adjust for pre-infusion levels measured at time points -1 and 0
data <- NULL
for(i in unique(Glucose2$Subject)){
   for(j in unique(Glucose2$Date)){
      temp <- subset(Glucose2, Subject==i & Date==j)
      temp$Conc <- temp$glucose - mean(c(temp$glucose[1], temp$glucose[2]))
      temp$Conc <- ifelse(temp$Conc < 0 | temp$Time <= 0, 0, temp$Conc)
      # handle intermediate values > 0 
      index1 <- which.max(temp$Conc)
      index2 <- which.min(temp$Conc[-c(1:index1)]) + index1
      if(temp$Conc[index2]==0){temp$Conc[c(index2:nrow(temp))] <- 0}
      data <- rbind(data,temp)
   }
}     
data <- subset(data, Time >= 0, select=c('Subject', 'Date', 'Time', 'Conc'))
 
# calculate individual AUCs 
res <- data.frame(matrix(nrow=length(unique(data$Subject))*2, ncol=3))
colnames(res) <- c('Subject', 'Date', 'AUC')
row <- 1
for(i in unique(Glucose2$Subject)){
   for(j in unique(Glucose2$Date)){
      temp <- subset(data, Subject==i & Date==j)
      res[row, c(1,2)] <- c(i,j)
      # NOTE: warning messages are not relevant as the focus is on AUC from 0 to 30
      # and not from 0 to infinity      
      res[row, 3] <- suppressWarnings(auc.complete(conc=temp$Conc, time=temp$Time)$est[1,1])
      row <- row + 1      
   }
}     
res <- res[order(res$Subject, res$Date),]
print(res) 

# geometric means and corresponding two-sided CIs per date
gm.ci <- function(x, conf.level=0.95){
   gm <- function(x){exp(mean(log(x)))}
   ci  <- exp(t.test(x=log(x), conf.level=conf.level)$conf.int)
   res <- data.frame(gm=gm(x), lower=ci[1], upper=ci[2])
}               
tapply(res$AUC, res$Date, gm.ci)

# comparison of AUCs using ratio of geometric means and corresponding two-sided CI 
# repeated experiment
model <- t.test(log(AUC)~Date, data=res, paired=TRUE, conf.level=0.90)
exp(as.real(model$estimate))
exp(model$conf.int)


## example for estimation of mean residence time (MRT) for a drug given as IV bolus 
## dataset Indometh of package datasets
require(datasets)
Indometh$id <- as.character(Indometh$Subject)
Indometh <- Indometh[order(Indometh$id, Indometh$time),]
res <- data.frame(matrix(ncol=2, nrow=length(unique(Indometh$Subject))))
colnames(res) <- c('ID', 'MRT')
row <- 1
for(i in unique(Indometh$id)){
   temp <- subset(Indometh, id==i)
   res[row, 1] <- i
   # estimate number of time points for estimation of terminal elimination rate
   change <- lee(conc=temp$conc, time=temp$time, method='ols', points=2)$chgpt
   n.tail <- nrow(subset(temp, temp$time > change))
   # assuming baseline levels of zero 
   parm <- auc.complete(conc=c(0,temp$conc), time=c(0,temp$time), n.tail=n.tail)$est
   # calculate MRT: AUMC0-inf divided by AUC0-inf
   res[row, 2] <- parm[3,2] / parm[3,1]
   row <- row + 1
}
print(res)

# geometric mean and corresponding two-sided CI
print(gm.ci(x=res$MRT, conf.level=0.95))

}

\keyword{misc}
