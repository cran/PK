
biexp <- function(conc, time, prev=0, tol=1E-9){
	
	# get start values for optim by curve peeling
	curve.peeling <- function(x, y){
	
		n <- length(y)	
		res <- NA
		Fmin <- Inf
	
		for(i in (n-3):3) { 
			parms <- lm(log(y[(i+1):n])~x[(i+1):n])
			b2 <- parms$coef[2]*(-1)
			a2 <- exp(parms$coef[1])
			ynew <- abs(y - a2*exp(-b2*x))
			if(!any(is.na(ynew)) && !any(ynew==0) && !any(ynew==Inf)){
				parms <- lm(log(ynew[1:i])~x[1:i])
				b1 <- parms$coef[2]*(-1)
				a1 <- exp(parms$coef[1])
				F <- sum((y-(a1*exp(-b1*x)+a2*exp(-b2*x)))*(y-(a1*exp(-b1*x)+a2*exp(-b2*x))))
				if (!is.na(F) && F < Fmin && !is.na(all(a1,a2,b1,b2)) && all(b1>0,b2>0,b1>b2)) {
					res <- as.real(c(a1=a1, dl=log(b1)-log(b2), a2=a2, b2=log(b2)))
					Fmin <- F
				}

			}
		}


		if (is.na(any(res))){
		parms <- lm(log(y)~x)
			b <- parms$coef[2]*(-1)
			a <- exp(parms$coef[1])
			F <- sum(parms$resid*parms$resid)
			if (!is.na(F) && F < Fmin && !is.na(b) && b > 0){
				res <- as.real(c(a=a, b=log(b)))
				Fmin <- F
			}
		}
		return(res)
	}	

	# remove missing values
	data <- na.omit(data.frame(conc, time))
			
	# check input 
	if(!is.real(prev)){stop('argument prev invalid')}
	if(prev<0){stop('pre-dosing value must be greater 0')}
	if (prev > 0) {data$conc <- data$conc - prev}

	# check input parameters and remove values below or equal to zero
	if(!is.vector(data$time)){stop('argument time invalid')}
	if(!is.vector(data$conc)){stop('argument conc invalid')}
	if (any(data$time < 0)) {stop('timepoint below zero')}
	
	time <- data$time
	conc <- data$conc

	# remove values below or equal to zero
	if (any(conc <= 0)) {
		for (i in 1:length(conc)) {
			if (conc[i] <= 0) {conc[i] <- NA}
		}
		warning('concentration below or equal to zero were omitted')
		data <- na.omit(data.frame(conc, time))
		conc <- data$conc
		time <- data$time		
	}
	if (nrow(data) < 4) {stop('a minimum of 4 observations are required')}

	biexploss <- function(par){
		a1 <- par[1]
		dl <- par[2]
		a2 <- par[3]
		b2 <- par[4]
		sum((conc - a1*exp(-(exp(b2) + exp(dl))*time) - a2*exp(-exp(b2)*time))^2)
	}

	singleloss <- function(par){
		a <- par[1]
		b <- par[2]
		sum((conc - a*exp(-exp(b)*time))^2)
	}

	start <- curve.peeling(y=conc, x=time)
	type <- as.character(length(start))	

	if(type == "4" && sum((conc - start[1]*exp(-(exp(start[4]) + exp(start[2]))*time) - 
			start[3]*exp(-exp(start[4])*time))^2 == Inf)) {type <- "1"}

	if(type == "2" && sum((conc - start[1]*exp(-exp(start[2])*time))^2 == Inf)){type <- "1"}

	switch(type, 
		"4" = {sol <- optim(par=start, fn=biexploss,  method=c("Nelder-Mead"), control=list(reltol=tol))	
			b1 <- (exp(sol$par[4]) + exp(sol$par[2]))
			a1 <- sol$par[1]
			b2 <- exp(sol$par[4])
			a2 <- sol$par[3]},
		"2" = {sol <- optim(par=start, fn=singleloss,  method=c("Nelder-Mead"), control=list(reltol=tol))
			b1 <- exp(sol$par[2]) 
			a1 <- sol$par[1]
			b2 <- exp(sol$par[2])
			a2 <- sol$par[1]},
		"1" = {a1 <- NA; b1 <- NA; a2<- NA; b2 <- NA},
	)

	# calculate halflife
	init.hl <- log(2) / b1
	term.hl <- log(2) / b2

	# format output object
	parms <- data.frame(initial=as.real(c(init.hl, b1, a1)),
			terminal=as.real(c(term.hl, b2, a2)))
	rownames(parms) <- c('halflife', 'slope', 'intercept')
	res <- list(parms=parms, conc=conc, time=time, method="biexp")
	class(res) <- 'halflife'
	return(res)
}
