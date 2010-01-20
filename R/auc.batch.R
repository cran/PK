auc.batch <- function(conc, time, group=NULL, method=c("t", "z", "boott"),  alternative=c("two.sided", "less", "greater"), conf.level=0.95, nsample=1000, data){
 
  # function for bootstrap resampling
  resample <- function(w, conc, time, group, alpha, nsample){
    B<-length(time)
    obsv.parm <- parms(w=w, conc=conc, time=time, group=group)
    boot.conc<-vector('list',B)
    boot.group<-vector('list',B)
    stat <- rep(NA, nsample) 
    for(j in 1:nsample){
       ## generate bootstrap data
       if(!is.null(group)){
         grps<-unique(unlist(group))
         n<-length(conc[[1]][group[[1]]==grps[1]])/length(time[[1]])
         m<-length(conc[[1]][group[[1]]==grps[2]])/length(time[[1]])
         for (i in 1:B){
           boot.conc[[i]]<-c(as.vector(matrix(conc[[i]][group[[i]]==grps[1]],nrow=n)[sample(1:n,replace=TRUE),]),as.vector(matrix(conc[[i]][group[[i]]==grps[2]],nrow=m)[sample(1:m,replace=TRUE),]))
           boot.group[[i]]<-c(rep(1,table(group[[i]])[1]),rep(2,table(group[[i]])[2]))
         }
       }else{
         n<-length(conc[[1]])/length(time[[1]])
         for (i in 1:B){
           boot.conc[[i]]<-as.vector(matrix(conc[[i]],nrow=n)[sample(1:n,replace=TRUE),])
         }
         boot.group<-NULL
       }
       boot.parm <- parms(w=w, conc=boot.conc, time=time, group=boot.group)
       stat[j] <- (boot.parm$est- obsv.parm$est) / sqrt(boot.parm$var)
    }
    lower <- quantile(stat, alpha/2, method=5)
    upper <- quantile(stat, 1-alpha/2, method=5)
    return(c(upper, lower*(-1)))
  }

  # function to calculate AUC and V(AUC)
  parms <- function(w, conc, time, group){
    B<-length(time) ### number of batches
    tindex<-vector('list',B)
    rtimes<-rank(unlist(time))
    h<-1
    l<-unlist(lapply(time,length))
    for(i in 1:B){
      tindex[[i]]<-rtimes[seq(from=h,length=l[i])]
      h<-h+l[i]
    } 
    ## difference of two AUC's
    if(!is.null(group)){
      conc1<-vector('list',B)
      conc2<-vector('list',B)
      grps<-unique(unlist(group))
      n<-length(conc[[1]][group[[1]]==grps[1]])/length(time[[1]])
      m<-length(conc[[1]][group[[1]]==grps[2]])/length(time[[1]])
  
      # seperating groups
      for(i in 1:B){
        conc1[[i]]<-conc[[i]][group[[i]]==grps[1]]
        conc2[[i]]<-conc[[i]][group[[i]]==grps[2]]
      }
      if(unlist(length(conc1))/n!=unlist(length(conc2))/m) stop('Number of observations does not match number of time points')
      auc1<-0
      auc2<-0
      aucvar1<-rep(NA,B)
      aucvar2<-rep(NA,B)
      for(i in 1:B){
        auc1<-auc1 +  mean(colSums(w[tindex[[i]]]*matrix(conc1[[i]],ncol=n,byrow=TRUE)))
        auc2<-auc2 +  mean(colSums(w[tindex[[i]]]*matrix(conc2[[i]],ncol=m,byrow=TRUE)))
        aucvar1[i]<-sum(w[tindex[[i]]]%o%w[tindex[[i]]]*cov(matrix(conc1[[i]],nrow=n)))/n
        aucvar2[i]<-sum(w[tindex[[i]]]%o%w[tindex[[i]]]*cov(matrix(conc2[[i]],nrow=m)))/m
      }
      res <- list(var=sum(aucvar1)+sum(aucvar2), est=auc1-auc2, df=(sum(aucvar1)/n)^2/sum(aucvar1^2/(n^2*(n-1)))+(sum(aucvar2)/m)^2/sum(aucvar2^2/(m^2*(m-1))))
      return(res)
    }else{      
      n<-length(conc[[1]])/length(time[[1]])
      auc<-0
      aucvar<-rep(NA,B)
      for(i in 1:B){
        auc<-auc +  mean(colSums(w[tindex[[i]]]*matrix(conc[[i]],ncol=n,byrow=TRUE)))
        aucvar[i]<-sum(w[tindex[[i]]]%o%w[tindex[[i]]]*cov(matrix(conc[[i]],nrow=n)))/n
      }
      res <- list(var=sum(aucvar), est=auc, df=(sum(aucvar)/n)^2/sum(aucvar^2/(n^2*(n-1))))
      return(res)
    }
  }

  # check input parameters  
  if(!missing(data)){
    cnames <- colnames(data)
    if(!any(cnames=='id')){stop("data does not contain a variable id")}
    if(!any(cnames=='conc')){stop("data does not contain a variable conc")}
    if(!any(cnames=='time')){stop("data does not contain a variable time")}
    temp <- .formattobatch(data)
    conc <- time <- group <- NULL
    for(i in 1:length(temp)){
      conc[[i]] <- temp[[i]]$conc 
      time[[i]] <- temp[[i]]$time
      if(any(names(temp[[1]])=='group')){
        group[[i]] <- temp[[i]]$group
      }
    }
  }  

  method <- match.arg(method,several.ok=TRUE)
  alternative <- match.arg(alternative)
  if(!is.null(group)){
    if(length(unlist(conc))!=length(unlist(group))){stop('different length of concentration and grouping lists')}
    if(nlevels(as.factor(unlist(group))) > 2){stop("limited for comparison of 2 groups")}
    if(nlevels(as.factor(unlist(group))) == 1) {group <- NULL}
  }
  
  # handle input parameters
  if(any(is.na(unlist(conc)))){stop('Missing values are not permitted in batch designs')}
  if (!is.list(time) || !is.list(conc)) stop('Both time and concentration need to be a list')
  if (!is.null(group) && !is.list(group)) stop('The parameter group needs to be a list or NULL')
  if (length(time) != length(conc)) stop('Number of batches defined in time and conc are different')
  if (length(time) != length(group) && !is.null(group)) stop('Number of batches defined in group does not match time and conc')

  n<-as.numeric(table(time[[1]])[1])
  if(sum(as.vector(table(unlist(time)))!=rep(n,times=length(unique(unlist(time)))))>0) stop('Number of observations differs between batches')
  if (sum(table(unlist(lapply(time,unique)))>1)>0) stop('Time points in each batch are not unique')

  grpfact <- levels(as.factor(unlist(group)))
  # specify alpha error
  alpha <- 1-conf.level

  if(is.null(group)){
    grp <- lapply(lapply(conc,'>',-Inf),as.numeric)
  }else{
    grp <- group
  }

  # sort concentrations by time order
  for(i in 1:length(time)){
    data <- data.frame(conc=conc[[i]], time=time[[i]], group=grp[[i]])
    data <- data[order(data$time), ]
    conc[[i]]<-data$conc
    time[[i]]<-data$time
    grp[[i]] <- data$group
  }

  if(!is.null(group)){
    group<-grp
  }  

  # eliminates multiples in time points
  time.unique <- lapply(time,unique)
  # calculate weights and observed parameters
  w <- .weight(time.unique)  

  # calculate confidence intervals
  obsv.parm <- parms(w=w, conc=conc, time=time.unique, group=group)
  if(alternative %in% c('less', 'greater')){alpha <- alpha*2}    
  z <- NULL
  
  method <- sort(method)  
  if (any(method=="boott")){
    z <- rbind(z,resample(w=w, conc=conc, time=time.unique, group=group, alpha=alpha, nsample=nsample))
  }

  if (any(method=="t")){
    z <- rbind(z,rep(qt(1-alpha/2, df=obsv.parm$df),2))
  }

  if (any(method=="z")) {
    z <- rbind(z,rep(qnorm(1-alpha/2),2))
  }
       
  est <- sum(obsv.parm$est[1], -obsv.parm$est[2], na.rm=TRUE)    
  lower <- est - sqrt(sum(obsv.parm$var))*z[,1]
  upper <- est + sqrt(sum(obsv.parm$var))*z[,2]
    
  switch(alternative,
    "less"={upper <- Inf},
    "greater"={lower <- -Inf},
    "two.sided"={},
  )

  df <- rep(NA, as.real(length(lower)))
  if(any(method=="t")){
    if(any(method=='z')){
      df=c(rep(NA, as.real(length(lower)-2)), obsv.parm$df, NA)
    }else{
      df=c(rep(NA, as.real(length(lower)-1)), obsv.parm$df)
    }
  } 
  res <- NULL
  res$est <- matrix(as.real(est),ncol=1)
  if(is.null(group)){
    rownames(res$est) <- 'AUC to tlast'
  }else{
    rownames(res$est) <- 'difference of AUCs to tlast'
  }
  colnames(res$est) <- 'est'
  res$design<-"batch"
  res$CIs<-data.frame(est=est, stderr=sqrt(sum(obsv.parm$var)), lower=lower, upper=upper, df=df,method=method)
  rownames(res$CIs) <- paste(conf.level*100,'% CI using a ', method,' distribution for AUC to tlast', sep='')
  res$conf.level <- conf.level
  res$conc <- conc
  res$time <- time
  res$group <- group
  class(res)<-"PK"
  return(res)
}


