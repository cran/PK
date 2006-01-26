biexp <- function(conc, time, prev=0, tol=1E-9, maxit=500){

	# get start values for optim by curve peeling
	curve.peeling <- function(x, y){
	
		n <- length(y)	
		res <- NA
		Fmin <- Inf
	
		for(i in (n-3):3) { 
			parms1 <- tryCatch(lm(log(y[(i+1):n])~x[(i+1):n])$coef, error=function(e) rep(NA,2))
			if(!is.na(all(parms1))){
				b2 <- parms1[2]*(-1)
				a2 <- exp(parms1[1])
				ynew <- abs(y - a2*exp(-b2*x))
				parms2 <- tryCatch(lm(log(ynew[1:i])~x[1:i])$coef, error=function(e) rep(NA,2))
				if(!is.na(all(parms2))){
					b1 <- parms2[2]*(-1)
					a1 <- exp(parms2[1])
					F <- sum((y-(a1*exp(-b1*x)+a2*exp(-b2*x)))*(y-(a1*exp(-b1*x)+a2*exp(-b2*x))))
					if (F < Fmin && all(b1>0,b2>0,b1>b2)) {
						res <- as.real(c(a1=a1, dl=log(b1)-log(b2), a2=a2, b2=log(b2)))
						Fmin <- F
					}

				}
			}

		}

		if (is.na(all(res))){
			parms <- tryCatch(lm(log(y)~x), error=function(e) rep(NA,2))
			if(!is.na(all(parms))){
				b <- parms$coef[2]*(-1)
				a <- exp(parms$coef[1])
				F <- sum(parms$resid*parms$resid)
				if (b > 0){
					res <- as.real(c(a=a, b=log(b)))
					Fmin <- F
				}
			}
		}
		return(res)
	}	



	# check input parameters and exclude missing values
	if (!is.vector(time) || !is.vector(conc)) {stop('argument time and/or conc invalid')}
	if (length(time) != length(conc)) {stop('time and conc differ in length')}
	if (any(time < 0)) {stop('at least one timepoint below zero')}
	data <- na.omit(data.frame(conc, time))
	
	# check input parameters and remove values below or equal to zero
	if (prev < 0) {stop('pre-dosing value must be greater 0')}
	if (prev > 0) {data$conc <- data$conc - prev}
	if (any(data$conc <= 0)) {
		data$conc[data$conc <= 0] <- NA
		warning('concentration below or equal to zero were omitted')
		data <- na.omit(data)	
	}
	if (nrow(data) < 4) {stop('a minimum of 4 observations are required')}

	# use data as vectors
	n <- nrow(data)
	time <- data$time
	conc <- data$conc

	# loss function for biexp
	biexploss <- function(par){
		a1 <- par[1]
		dl <- par[2]
		a2 <- par[3]
		b2 <- par[4]
		sum((conc - a1*exp(-(exp(b2) + exp(dl))*time) - a2*exp(-exp(b2)*time))^2)
	}

	# loss function for single exp
	singleloss <- function(par){
		a <- par[1]
		b <- par[2]
		sum((conc - a*exp(-exp(b)*time))^2)
	}

	# get starting values
	start <- curve.peeling(y=conc, x=time)
	type <- as.character(length(start))

	# check sum of squared residuals using estimates obtained by curve peeling
	switch(type,
		"4" = {sum.resid <- biexploss(par=start)},
		"2" = {sum.resid <- singleloss(par=start)},
	)
	if(sum.resid == Inf){type <- "1"}

	switch(type, 
		"4" = {sol <- optim(par=start, fn=biexploss,  method=c("Nelder-Mead"), control=list(reltol=tol, maxit=maxit))$par	
			b1 <- (exp(sol[4]) + exp(sol[2]))
			a1 <- sol[1]
			b2 <- exp(sol[4])
			a2 <- sol[3]},
		"2" = {sol <- optim(par=start, fn=singleloss,  method=c("Nelder-Mead"), control=list(reltol=tol, maxit=maxit))$par
			b1 <- exp(sol[2]) 
			a1 <- sol[1]
			b2 <- exp(sol[2])
			a2 <- sol[1]},
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







