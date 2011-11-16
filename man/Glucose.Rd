\name{Glucose}
\alias{Glucose}
\title{Baseline adjusted glucose levels following alcohol ingestion}
\description{
  The \code{Glucose} data frame has 196 rows and 4 columns. The dataset is originally in package nlme as Glucose2.
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{id}{
      a factor with levels
      \code{1} to \code{7} identifying the subject whose glucose
      level is measured. 
    }
    \item{date}{
      a factor with levels
      \code{1} 
      \code{2}
      indicating the occasion in which the experiment was conducted.
    }
    \item{time}{
      a numeric vector giving the time since alcohol ingestion (in min/10).
    }
    \item{conc}{
      a numeric vector giving the blood glucose level (in mg/dl) adjusted for baseline.
    }
  }
}
\details{
  Hand and Crowder (Table A.14, pp. 180-181, 1996) describe data on
  the blood glucose levels measured at 14 time points over 5 hours for 7
  volunteers who took alcohol at time 0. The same experiment was
  repeated on a second date with the same subjects but with a dietary
  additive used for all subjects.\cr\cr

  Dataset was corrected for baseline using the following code:

  \preformatted{## dataset Glucose2 of package nlme
require(nlme)
Glucose2 <- Glucose2[order(Glucose2$Subject, Glucose2$Date, Glucose2$Time),]
## adjust for pre-infusion levels measured at time points -1 and 0
data <- NULL
for(i in unique(Glucose2$Subject)){
  for(j in unique(Glucose2$Date)){
     temp <- subset(Glucose2, Subject==i & Date==j)
     temp$Conc <- temp$glucose - mean(c(temp$glucose[1], temp$glucose[2]))
     temp$Conc <- ifelse(temp$Conc < 0 | temp$Time <= 0, 0, temp$Conc)
     ## handle intermediate values > 0 
     index1 <- which.max(temp$Conc)
     index2 <- which.min(temp$Conc[-c(1:index1)]) + index1
     if(temp$Conc[index2]==0){temp$Conc[c(index2:nrow(temp))] <- 0}
     data <- rbind(data,temp)
   }
}    
Glucose <- subset(data, Time >= 0, 
                  select=c('Subject', 'Date', 'Time', 'Conc'))
names(Glucose) <- c("id","date","time","conc")
}
}
\source{
  Pinheiro, J. C. and Bates, D. M. (2000), \emph{Mixed-Effects Models in S
    and S-PLUS}, Springer, New York.  (Appendix A.10)

  Hand, D. and Crowder, M. (1996), \emph{Practical Longitudinal Data
    Analysis}, Chapman and Hall, London.
}
%\examples{}
\keyword{datasets}
