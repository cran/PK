

AUC <- function(time, conc, exact=NA, numintp=2, numtail=3, prev=0) {		     

	# function for linear interpolation/extrapolation
	linpol <- function(conc, time, exact){
		parms <- lm(conc~time)
		return(parms$coef[2] * exact + parms$coef[1])
	}
         
	# function to add parts of auc and aumc
	add <- function(time, conc) {
		auc <- 0; aumc <- 0
		for (i in 2:length(conc)) {
			auc <- auc  + 1/2 * (time[i]-time[i-1]) * (conc[i]+conc[i-1])
			aumc <- aumc + 1/2 * (time[i]-time[i-1]) * (conc[i]*time[i] + conc[i-1]*time[i-1])
		}
		return(list(auc=auc, aumc=aumc))
	}
	
	# remove missing values
	data <- na.omit(data.frame(conc, time))
	data <- data[order(data$time),]
	time <- data$time
	conc <- data$conc
	n <- nrow(data) 	           

	# subtraction of pre dosing concentration for single dose studies
	if (prev > 0) {conc <- conc - prev}

	# remove values below zero
	if (any(conc < 0)) {
		for (i in 1:length(conc)) {
			if (conc[i] < 0) {conc[i] <- NA}
		}
		warning('concentration below zero were omitted')
		data <- na.omit(data.frame(conc, time))
		time <- data$time
		conc <- data$conc	
		n <- nrow(data) 
	}

	# check input parameters
	if (numtail < 2) {stop('number of points for tail area correction must be greater than 1')}
	if (numintp < 2) {stop('number of points for interpolation must be greater than 1')}
	         
	# calculate observed auc and aumc
	auc.observed  <- add(time=time, conc=conc)$auc        
	aumc.observed <- add(time=time, conc=conc)$aumc
    
	# calculate auc from 0 to infinity and aumc from 0 to infinity  by using last numtail points above zero
	tail <- subset(data.frame(conc, time), conc > 0)
	tail <- tail[(nrow(tail)-numtail+1) : nrow(tail), ]	
	lamda <- as.real(lm(log(tail$conc)~tail$time)$coef[2])*(-1)
	auc.infinity <- auc.observed + conc[n]/lamda 
	aumc.infinity <- aumc.observed + (conc[n]*time[n])/lamda + conc[n]/lamda**2
	if(lamda < 0){
		warning('tail area correction incorrect due to increasing concentration of last numtail points')	
		auc.infinity <- NA
		aumc.infinity <- NA
	}

	# calculate auc and aumc from 0 to exact where exact must be greater than time[n-1]
	auc.interpol <- NA; aumc.interpol <- NA	
	if (!is.na(exact) & exact > time[n-1] & exact == time[n]) { # special case
		auc.interpol <- auc.observed; aumc.interpol <- aumc.observed
	}	
	if (!is.na(exact) & exact > time[n-1] & exact != time[n]) {
		conc[n] <- linpol(conc=conc[(n-numintp+1):n], time=time[(n-numintp+1):n], exact=exact)
		time[n] <- exact
		if(conc[n] < 0){warning('interpolated value below zero')}
		auc.interpol <- add(time=time, conc=conc)$auc
		aumc.interpol <- add(time=time, conc=conc)$aumc		
	} 
	

	# define output object
	res <- data.frame(AUC=c(auc.observed, auc.interpol, auc.infinity), 
		AUMC=c(aumc.observed, aumc.interpol, aumc.infinity))
	rownames(res) <- c('observed', 'interpolated', 'infinity')
	return(res)      
}      
  