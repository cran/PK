eqv.ssd <- function(conc, time, group, method=c("fieller", "asymp", "boott"), conf.level=0.90, strata=NULL, nsample=1000, data=data){

    # function to calculate weights for AUC
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

    # function to calculate AUC from 0 to the last time point
    auc <- function(conc, time){
        xq <- tapply(conc, time, mean)
        sxq <- tapply(conc, time, var) / tapply(conc, time, length)
        est <- sum(weight(time)*xq)
        var <- sum(weight(time)^2*sxq)
        return(list(est=est, var=var))	
    }

    # function to calculate fieller confidence interval
    fieller <- function(data1, data2, w, alpha){
        # get estimates
        auc1 <- auc(data1$conc, data1$time)
        auc2 <- auc(data2$conc, data2$time)	
        est <- auc1$est / auc2$est

        # degress of freedom
        sd1 <- tapply(data1$conc, data1$time, sd)
        sd2 <- tapply(data2$conc, data2$time, sd)
        n1 <- tapply(data1$conc, data1$time, length)
        n2 <- tapply(data2$conc, data2$time, length)
        den <- (sum(w^2*sd1^2*n1^-1) + est^2*sum(w^2*sd2^2*n2^-1))^2
        nom <- sum(w^4*sd1^4*(n1^2*(n1-1))^-1) + est^4*sum(w^4*sd2^4*(n2^2*(n2-1))^-1)
        df <- den/nom
	
        # 1-alpha confidence interval
        t <- qt(p=1-alpha/2, df)
        lower <- auc1$est*auc2$est - sqrt((auc1$est*auc2$est)^2 - (auc1$est^2 - t^2*auc1$var)*(auc2$est^2 - t^2*auc2$var))
        upper <- auc1$est*auc2$est + sqrt((auc1$est*auc2$est)^2 - (auc1$est^2 - t^2*auc1$var)*(auc2$est^2 - t^2*auc2$var))
        lower <- lower / (auc2$est^2 - t^2*auc2$var)
        upper <- upper / (auc2$est^2 - t^2*auc2$var)

        # check condition
        res <- NA
        tmax <- (auc1$est^2*auc2$var + auc2$est^2*auc1$var)/(auc1$var*auc2$var)
        if(0 < t^2 & t^2 < tmax^2){res <- data.frame(ratio=est, lower=lower, upper=upper, df=df)}
        return(res)
    }

    # function to calculate asymptotic confidence interval
    asymp <- function(data1, data2, alpha){		
        # get estimates
        auc1 <- auc(data1$conc, data1$time)
        auc2 <- auc(data2$conc, data2$time)
        auc1$var <- auc1$var	
        auc2$var <- auc2$var

        # variance for ratio
        var.asymp <- auc1$var/(auc2$est^2) + (auc1$est^2)/(auc2$est^4)*auc2$var	

        # 1-alpha confidence interval
        z <- qnorm(p=1-alpha/2, mean=0, sd=sqrt(1))
        est <- auc1$est / auc2$est
        lower <- est - z*sqrt(var.asymp)
        upper <- est + z*sqrt(var.asymp)
        res <- data.frame(ratio=est, lower=lower, upper=upper)
        return(res)	
    }

    # function to calculate bootstrap-t confidence interval
    boott <- function(data1, data2, alpha, nsample, bystrata1, bystrata2){	
        # get estimates
        obsv.auc1 <- auc(data1$conc, data1$time)
        obsv.auc2 <- auc(data2$conc, data2$time)
        obsv.est <- obsv.auc1$est / obsv.auc2$est
        obsv.var <- obsv.auc1$var/(obsv.auc2$est^2) + (obsv.auc1$est^2)/(obsv.auc2$est^4)*obsv.auc2$var
			
        # bootstrap distribution
        boot.stat <- rep(NA,nsample)
        for(i in 1:nsample){
            boot.data1 <- data.frame(time=data1$time, conc=unlist(tapply(data1$conc, bystrata1, sample, replace=TRUE)))
            boot.data2 <- data.frame(time=data2$time, conc=unlist(tapply(data2$conc, bystrata2, sample, replace=TRUE)))
            boot.auc1 <- auc(boot.data1$conc, boot.data1$time)
            boot.auc2 <- auc(boot.data2$conc, boot.data2$time)
            boot.est <- boot.auc1$est / boot.auc2$est
            boot.var <- boot.auc1$var/(boot.auc2$est^2) + (boot.auc1$est^2)/(boot.auc2$est^4)*boot.auc2$var
            boot.stat[i] <- (boot.est - obsv.est) / sqrt(boot.var)
        }
        t.lb <- quantile(boot.stat, probs=c(alpha/2),   method=5, na.rm=TRUE)
        t.ub <- quantile(boot.stat, probs=c(1-alpha/2), method=5, na.rm=TRUE)
        base <- data.frame(est=obsv.est, t.lb=t.lb, t.ub=t.ub)	
        base$lower <- base$est - base$t.ub*sqrt(obsv.var)
        base$upper <- base$est - base$t.lb*sqrt(obsv.var)
        return(data.frame(ratio=base$est, lower=base$lower, upper=base$upper))
    }

    # check function call with data frame
    if(!missing(data)){
        conc <- data$conc
        time <- data$time
        group <- data$group
        strata <- data$strata
    }

    # check input parameters
    method <- match.arg(method)
    if(length(conc) != length(time)){stop('different length of input vectors')}
    if(!is.null(strata)){
        if(length(conc)!=length(strata)){stop('different length of input vectors')}
        if(method!='boott'){warning("strata variable only applicable in method boott")}
    }
    if(!is.null(group)){
        if(length(conc)!=length(group)){stop('different length of input vectors')}
        if(length(unique(group))!=2){stop("limited for comparison of 2 groups")}
    }

    # handle input parametersd
    data <- data.frame(conc=conc, time=time)
    if(is.null(strata)){strata <- rep(1, nrow(data))}
    data <- cbind(data, group=as.factor(group), strata=as.factor(strata))
    data <- na.omit(data[order(data$strata, data$group, data$time),])  
    grpfact <- levels(data$group)

    # check min number of observations
    if(any(na.omit(as.vector(tapply(data$conc, list(data$strata, data$group, data$time), length))) < 2)){
        stop('at least 2 observations per strata, group and time point required')
    }
	
    # handle strata variable for method boott
    if(method=='boott'){ 
        bystrata <- data.frame(time, group, strata)
        bystrata <- bystrata[order(bystrata$group, bystrata$strata, bystrata$time),]
        bystrata1 <- as.list(subset(bystrata, bystrata$group==grpfact[1]))
        bystrata2 <- as.list(subset(bystrata, bystrata$group==grpfact[2]))
    }
	
    # specify alpha error
    alpha <- 1-conf.level

    # calculate weights and observed parameters
    w1 <- weight(unlist(subset(data, data$group==grpfact[1], select='time')))
    w2 <- weight(unlist(subset(data, data$group==grpfact[2], select='time')))
    if(length(w1) != length(w2) || !all(w1==w2)){stop('time points are not identical for both groups')}
    w <- weight(data$time)

    data1 <- subset(data, data$group==grpfact[1])
    data2 <- subset(data, data$group==grpfact[2])
		 
    switch(method, 
        "asymp"={res <- asymp(data1=data1, data2=data2, alpha=alpha)}, 
        "fieller"={res <- fieller(data1=data1, data2=data2, w=w, alpha=alpha)},
        "boott"={res <- boott(data1=data1, data2=data2, alpha=alpha, nsample=nsample, bystrata1=bystrata1, bystrata2=bystrata2)
    },)
   
    rownames(res) <- paste(conf.level*100,'% CI using ', method, sep='')
    return(res)
}
