auc.batch <- function(conc, time, exact=NA, n.interpolate=2, n.tail=3) {		     
# check input parameters
  if (!is.list(time) || !is.list(conc)) stop('Both time and concentration need to be a list')
  if (length(time) != length(conc)) stop('Number of batches defined in time and conc are different')
  if (any(is.na(unlist(conc)))) stop('Missing values are not allowed for a batch design')
  if (sum(table(unlist(lapply(time,unique)))>1)>0) stop('Time points in each batch are not unique')
  if (sum(unlist(lapply(conc,length))!=unlist(lapply(time,length)))>0) stop('Length of conc and time differ')


  n<-as.numeric(table(time[[1]])[1])
  if(sum(as.vector(table(unlist(time)))!=rep(n,times=length(unique(unlist(time)))))>0) stop('Number of observations differs between batches')

  B<-length(time) ### number of batches
  if (any(unlist(conc) < 0)) {
    warning('concentration below zero are set to zero') 
    for(i in 1:B){
      conc[[i]][conc[[i]]<0]<-0    
    }
  }

  # sort concentrations by time order
  for(i in 1:B){
    data <- data.frame(conc=conc[[i]],time=time[[i]])
    data <- data[order(data$time), ]
    conc[[i]]<-data$conc
    time[[i]]<-data$time
  }

  # eliminates multiples in time points
  time <- lapply(time,unique)
  # finds the weights
  utimes<-sort(unlist(time))
  ntimes<-length(utimes)
  w<-rep(NA,ntimes)
  for(i in 2:(ntimes-1)){
    w[i]<-(utimes[i+1]-utimes[i-1])/2
  }
  w[1]<-(utimes[2]-utimes[1])/2
  w[ntimes]<-(utimes[ntimes]-utimes[ntimes-1])/2
  
  # finds the indices of time points
  tindex<-vector('list',B)
  rtimes<-rank(unlist(time))
  h<-1
  l<-unlist(lapply(time,length))
  for(i in 1:B){
    tindex[[i]]<-rtimes[seq(from=h,length=l[i])]
    h<-h+l[i]
  } 

  auc<-0
  aumc<-0
  for(i in 1:B){
    auc<-auc +  mean(colSums(w[tindex[[i]]]*matrix(conc[[i]],ncol=n,byrow=TRUE)))
    aumc<-aumc + mean(colSums(time[[i]]*w[tindex[[i]]]*matrix(conc[[i]],ncol=n,byrow=TRUE)))
  }
  auc.interpol<-NA
  aumc.interpol<-NA

  meanconc<-as.vector(tapply(unlist(conc),rep(unlist(time),each=n),mean,na.rm=TRUE))
  temp<-auc.complete(meanconc,unlist(time),n.tail=n.tail)
  auc.infinity<-temp$est[3,1] 
  aumc.infinity<-temp$est[3,2]
  res <- list(est=data.frame(AUC=c(as.real(auc), as.real(auc.interpol), as.real(auc.infinity)), 
              AUMC=c(as.real(aumc), as.real(aumc.interpol), as.real(aumc.infinity))))
  res$design<-"batch"
  res$CIs<-NULL
  res$test<-NULL
  rownames(res$est) <- c('observed', 'interpolated', 'infinity')
  class(res)<-"PK"
  return(res)
}

