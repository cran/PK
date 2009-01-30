ptest.ssd <- function(conc, time, group, alternative=c("two.sided", "less", "greater"), nsample=1000, data){

    # function to calculate statistic 		
    statistic <- function(data, grpfact){	
        data1 <- subset(data, data$group == grpfact[1])
        data2 <- subset(data, data$group == grpfact[2])
        est1 <- auc.ci(conc=data1$conc, time=data1$time, method='z')$est
        est2 <- auc.ci(conc=data2$conc, time=data2$time, method='z')$est		
        return(est1[1,1]-est2[1,1])
    }

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
	
    # check function call with data frame
    if(!missing(data)){
        conc <- data$conc
        time <- data$time
        group <- data$group
    }

    # handle input data
    alternative <- match.arg(alternative)
    if(length(conc) != length(time) || length(conc) != length(group) || length(time) != length(group)){stop("different lengths of input vectors")}
    grpfact <- levels(as.factor(group))
    if(length(unique(group))!=2){stop("limited for comparison of 2 groups")}
    
    # check min number of observations
    if(any(na.omit(as.vector(tapply(conc, list(group, time), length))) < 2)){
        stop('at least 2 observations per group and time point required')
    }

    # check assumption of identical time points
    check <- data.frame(conc=conc, time=time, group=as.factor(group))
    check <- na.omit(check[order(check$group, check$time),])
    w1 <- weight(unlist(subset(check, check$group==grpfact[1], select='time')))
    w2 <- weight(unlist(subset(check, check$group==grpfact[2], select='time')))
    if(length(w1) != length(w2) || !all(w1==w2)){stop('time points are not identical for both groups')}
	
    data <- data.frame(conc=conc, time=time, group=factor(group))
    data <- data[order(data$time, data$group),]
	
    # calculate observed parameters
    t.obsv <- statistic(data=data, grpfact=grpfact)	
		
    intern <- data
    intern$group <- NULL
    t.perm <- rep(NA, nsample)
    for(i in 1:nsample){
        intern$ran <- runif(nrow(intern))
        intern <- intern[order(intern$time, intern$ran),]
        intern$group <- data$group
        t.perm[i] <- statistic(data=intern, grpfact=grpfact)
    }	
	
    switch(alternative, 
        "less"={p.value <- length(subset(t.perm, t.obsv <= t.perm))/nsample},	
        "greater"={p.value <- length(subset(t.perm, t.obsv >= t.perm))/nsample},
        "two.sided"={p.value <- length(subset(t.perm, abs(t.perm)>abs(t.obsv)))/nsample},
    )
    
    res <- data.frame(statistic=unlist(t.obsv), p.value=unlist(p.value))
    rownames(res) <- paste('Hypothesis:', alternative)
    return(res)
}
