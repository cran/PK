auc.ssd <- function(conc, time, exact=NA, n.interpolate=2, n.tail=3) {		     
  ## check input parameters
  if (!is.vector(time) || !is.vector(conc)) {stop('Argument time and/or conc invalid')}
  if (n.tail < 2) {stop('Number of points for tail area estimation must be greater than 1')}

  if (length(time) != length(conc)) {stop('Length of conc and time differs')}

  weight <- function(time){
    if(length(time) == 0){return(0)}
    time <- unique(time)
    J <- length(time)
    i <- c(2:(J-1))
    w <- rep(NA,J)
    w[1] <- (time[2] - time[1])/2
    w[i] <- (time[i+1] - time[i-1])/2
    w[J] <- (time[J] - time[J-1])/2
    return(w)
  }
  data <- data.frame(conc=conc, time=time)
  data <- data[order(data$time), ]
  time <- data$time
  conc <- data$conc

  w<-weight(time)
  if(any(conc<0)){
    warning('concentration below zero were set to zero')
    conc[conc<0]<-0
  }

  meanconc<-tapply(conc,time,mean,na.rm=TRUE)

  auc <- sum(w*meanconc)
  aumc<- sum(unique(time)*w*meanconc)
  auc.interpol<-NA
  aumc.interpol<-NA

  temp<-auc.complete(as.vector(meanconc),unique(time),n.tail=n.tail)
  auc.infinity<-temp$est[3,1] 
  aumc.infinity<-temp$est[3,2]

  res <- list(est=data.frame(AUC=c(as.real(auc), as.real(auc.interpol), as.real(auc.infinity)), 
              AUMC=c(as.real(aumc), as.real(aumc.interpol), as.real(aumc.infinity))))
  res$design<-"ssd"
  res$CIs<-NULL
  res$test<-NULL
  rownames(res$est) <- c('observed', 'interpolated', 'infinity')
  class(res)<-"PK"
  return(res)     
}
