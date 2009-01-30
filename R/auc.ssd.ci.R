auc.ssd.ci <- function(conc, time, group=NULL, method=c("t", "z", "boott"),  alternative=c("two.sided", "less", "greater"), conf.level=0.95, strata=NULL, nsample=1000){

  # function to calculate weights for linear trapezoidal rule
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
  # function for bootstrap resampling
  resample <- function(w, data, alpha, nsample, bystrata){
    obsv.parm <- parms(w=w, data=data)
    stat <- rep(NA, nsample) 
    for(i in 1:nsample){
      boot.conc <- unlist(by(data$conc, bystrata, sample, replace=TRUE))
      boot.data <- data.frame(time=unlist(split(data$time,bystrata[[3]])), conc=boot.conc, group=unlist(split(data$group,bystrata[[3]])))
      boot.parm <- parms(w=w, data=boot.data)
      stat[i] <- (sum(boot.parm$est[1], -boot.parm$est[2], na.rm=TRUE) - sum(obsv.parm$est[1], -obsv.parm$est[2], na.rm=TRUE)) /
                 sqrt(sum(boot.parm$var))
    }
    lower <- quantile(stat, alpha/2, method=5)
    upper <- quantile(stat, 1-alpha/2, method=5)
    return(c(upper, lower*(-1)))
  }

  # function to calculate AUC and V(AUC)
  parms <- function(w, data){
    n <- tapply(data$conc, list(data$time, data$group), length)
    var <- apply(w^2*tapply(data$conc, list(data$time, data$group), var)/n, 2, sum, na.rm=TRUE)
    est <- apply(w*tapply(data$conc, list(data$time, data$group), mean), 2, sum) 
    res <- list(n=n, var=var, est=est)
    return(res)
  }

  # check input parameters
  method <- match.arg(method,several.ok=TRUE)
  alternative <- match.arg(alternative)
  if (length(time) != length(conc) && is.null(group)) {stop('Length of conc and time differs')}
  if (!is.null(group) && length(time) != length(group)) {stop('Length of group does not match time and conc')}

  if(!is.null(strata)){
    if(length(conc)!=length(strata)){stop('different length of input vectors')}
    if(any(method!='boott')){warning("strata variable only applicable in method boott")}
  }
  if(!is.null(group)){
    if(length(conc)!=length(group)){stop('different length of input vectors')}
    if(nlevels(as.factor(group)) > 2){stop("limited for comparison of 2 groups")}
  }

  # handle input parameters
  if(!is.null(group)){
    if(nlevels(as.factor(group))!=2) stop('Exactly 2 groups need to be specified')
    if(length(conc)!=length(group)){stop('Different length of conc/time and group')}
  }

  data <- data.frame(conc=conc, time=time)
  if(is.null(strata)){strata <- rep(1, nrow(data))}
  if(is.null(group)){group <- rep(1, nrow(data))}
  data <- cbind(data, group=as.factor(group), strata=as.factor(strata))
  data <- na.omit(data[order(data$strata, data$group, data$time),])  
  data <- data[order(data$group,data$strata,data$time), ]
  grpfact <- levels(data$group)
  # check min number of observations
  if(any(na.omit(as.vector(tapply(data$conc, list(data$strata, data$group, data$time), length))) < 2)){
    stop('at least 2 observations per strata, group and time point required')
  }

  # handle strata variable for method boott
  if(nlevels(data$strata) == 1){bystrata <- list(data$time, data$group, data$strata)}
  if(nlevels(data$strata) >  1){bystrata <- list(data$time, data$group, data$strata)}
       
  # specify alpha error
  alpha <- 1-conf.level

  # calculate weights and observed parameters
  w1 <- weight(unlist(subset(data, data$group==grpfact[1], select='time')))
  w2 <- weight(unlist(subset(data, data$group==grpfact[2], select='time')))
  if(nlevels(data$group) == 2 & (length(w1) != length(w2) || !all(w1==w2))){stop('time points are not identical for both groups')}
  w <- weight(data$time)  

  # calculate confidence intervals
  z <- rep(NA,2)
  obsv.parm <- parms(w=w, data=data)
  if(alternative %in% c('less', 'greater')){alpha <- alpha*2}    
  z <- NULL
    
  method <- sort(method)  
  if (any(method=="boott")){
    z <- rbind(z,resample(w=w, data=data, alpha=alpha, nsample=nsample, bystrata=bystrata))
  }

  if (any(method=="t")){
    nom <- sum(obsv.parm$var)^2
    den <- sum(apply(w^4*tapply(data$conc,list(data$time,data$group),sd)^4/(obsv.parm$n^2*(obsv.parm$n-1)),2,sum,na.rm=TRUE))
    v <- nom / den
    z <- rbind(z,rep(qt(1-alpha/2, df=v),2))
  }

  if (any(method=="z")) {
    z <- rbind(z,rep(qnorm(1-alpha/2),2))
  }

  est <- sum(obsv.parm$est[1], -obsv.parm$est[2], na.rm=TRUE)
    
  lower <- est - sqrt(sum(obsv.parm$var))*z[,1]
  upper <- est + sqrt(sum(obsv.parm$var))*z[,2]
    
  switch(alternative,
    "less"={upper <- rep(Inf,length(lower))},
    "greater"={lower <- rep(-Inf,length(upper))},
    "two.sided"={},
  )

  df <- rep(NA, as.real(length(lower)))
  if(any(method=="t")){
    if(any(method=='z')){
      df=c(rep(NA, as.real(length(lower)-2)), v, NA)
    }else{
      df=c(rep(NA, as.real(length(lower)-1)), v)
    }
  }
  res <- list(est=data.frame(AUC=c(as.real(est), NA, NA), AUMC=c(NA, NA, NA)))
  rownames(res$est) <- c('observed', 'interpolated', 'infinity')
  res$design<-"ssd"
  res$CIs<-data.frame(est=est, stderr=sqrt(sum(obsv.parm$var)), lower=lower, upper=upper, df=df,method=method)
  rownames(res$CIs) <- paste(conf.level*100,'% CI using a ', method,' distribution', sep='')
  res$test<-NULL
  class(res)<-"PK"
  return(res)
}

