auc.complete <- function(conc, time, exact=NA, n.interpolate=2, n.tail=3) {		     

	# function for linear interpolation/extrapolation
	linpol <- function(conc, time, exact){
		parms <- lm(conc~time)$coef
		return(parms[2] * exact + parms[1])
	}
         
	# function to add parts of auc and aumc
	add <- function(time, conc) {
		auc <- 0; aumc <- 0
		len <- length(conc)
		timediff <- diff(time)
		auc <- 0.5*sum(timediff*rowSums(cbind(conc[2:len],conc[1:(len-1)])))
		aumc <- 0.5*sum(timediff*rowSums(cbind(conc[2:len]*time[2:len],conc[1:(len-1)]*time[1:(len-1)])))
		return(list(auc=auc, aumc=aumc))
	}

	# check input parameters and exclude missing values
	if (!is.vector(time) || !is.vector(conc)) {stop('argument time and/or conc invalid')}
	if (length(time) != length(conc)) {stop('time and conc differ in length')}
	if (n.tail < 2) {stop('number of points for tail area correction must be greater than 1')}
	if (n.interpolate < 2) {stop('number of points for interpolation must be greater than 1')}
	data <- na.omit(data.frame(conc, time))
	
	if (any(data$conc < 0)) {
		data$conc[data$conc < 0] <- 0
		warning('concentration below zero were set to zero')
	}
	if (nrow(data) < 4) {stop('a minimum of 4 observations are required')}

	# use data as vectors
	data <- data[order(data$time),]
	n <- nrow(data)
	time <- data$time
	conc <- data$conc
	 	                
	# calculate observed auc and aumc
	auc.observed  <- add(time=time, conc=conc)$auc        
	aumc.observed <- add(time=time, conc=conc)$aumc

	# calculate auc from 0 to infinity and aumc from 0 to infinity by using last n.tail points above zero
	tail <- subset(data.frame(conc, time), conc > 0)
	tail <- tail[(nrow(tail)-n.tail+1) : nrow(tail), ]	

	lambda <- as.real(lm(log(tail$conc)~tail$time)$coef[2])*(-1)
	auc.infinity <- auc.observed + conc[n]/lambda 
	aumc.infinity <- aumc.observed + (conc[n]*time[n])/lambda + conc[n]/lambda**2
	if(lambda < 0){
		warning('tail area correction incorrect due to increasing concentration of last n.tail points')	
		auc.infinity <- NA
		aumc.infinity <- NA
	}

	# calculate auc and aumc from 0 to exact where exact must be greater than time[n-1]
	auc.interpol <- NA; aumc.interpol <- NA	
	if (!is.na(exact) & exact > time[n-1] & exact == time[n]) { # special case
		auc.interpol <- auc.observed; aumc.interpol <- aumc.observed
	}	
	if (!is.na(exact) & exact > time[n-1] & exact != time[n]) {
		conc[n] <- linpol(conc=conc[(n-n.interpolate+1):n], time=time[(n-n.interpolate+1):n], exact=exact)
		time[n] <- exact
		if(conc[n] < 0){warning('interpolated value below zero')}
		auc.interpol <- add(time=time, conc=conc)$auc
		aumc.interpol <- add(time=time, conc=conc)$aumc		
	} 

	# define output object
	res <- list(est=data.frame(AUC=c(as.real(auc.observed), as.real(auc.interpol), as.real(auc.infinity)), 
		AUMC=c(as.real(aumc.observed), as.real(aumc.interpol), as.real(aumc.infinity))))
        res$design<-"complete"
        res$CIs<-NULL
        res$test<-NULL
	rownames(res$est) <- c('observed', 'interpolated', 'infinity')
        class(res)<-"PK"
	return(res)      
}
