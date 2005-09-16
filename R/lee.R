
lee <- function(time, conc, points=3, prev=0, method=c("lad", "ols", "hub", "npr"), longer.terminal=TRUE) {

	# function for lad regression
	lad <- function(y, x) {
		resid <- Inf
		for (i in 1:length(y)) {
			for (j in 1:length(y)) {
				if ((i != j) & (x[j] != x[i])) {
					slope <- (y[j] - y[i]) / (x[j] - x[i])
					intct <- y[j] - slope*x[j] 
					absresid <- abs(y - (intct + slope*x))
					if (sum(absresid) < resid) {
						mad <- median(absresid)
						resid <- sum(absresid)
						k <- slope
						d <- intct
					}
				}
			}
		}
		return(list(k=as.real(k), d=as.real(d), resid=as.real(resid), mad=as.real(mad)))
	}
	
	# function for huber m regression 
	# acknowledgment to werner engl
	hub <- function(y, x, mad=lad(y=y,x=x)$mad, sigmafactor=1.483, kfactor=1.5) { 
		hubloss <- function(kd) { # Huber loss for k=kd[1], d=kd[2]
			absresid <- abs(y-kd[[1]]*x-kd[[2]])
			khuber <- kfactor*sigmafactor*mad
			sum(ifelse(absresid < khuber, absresid*absresid, khuber*(2*absresid-khuber))) 
		}
		start <- as.vector(c(lm(y~x)$coef[2], lm(y~x)$coef[1]))
		res <- optim(start, hubloss, method="Nelder-Mead", control=c(reltol=1e-9))		
		return(list(k = as.real(res$par[1]), d = as.real(res$par[2]), resid = as.real(res$value)))		
	}

	# function for nonparametric regression
	npr <- function(y, x) { 

		weighted.median <- function(w, x) { 
			data <- data.frame(x, w)
			data <- data[order(data$x),]
			i <- 1; while(sum(data$w[1:i]) <= 0.5) {i <- i + 1}
			ifelse (sum(data$w[1:i-1]) == 0.5, return((data$x[i-1]+data$x[i])/2), return(data$x[i]))
		}

		total <- 0
		for(i in 1:(length(y)-1)) {
			for(j in (i+1):length(y)){total <- total + abs(x[i]-x[j])}
		}

		l <- 1
		b <- array(1:(length(y)*(length(y)-1)/2))
		w <- array(1:(length(y)*(length(y)-1)/2))
		for(i in 1:(length(y)-1)) {
			for(j in (i+1):length(y)){	
				b[l] <- ((y[i]-y[j])/(x[i]-x[j]))
				w[l] <- abs(x[i]-x[j]) / total			 
				l <- l + 1									
			}
		}

		data <- subset(data.frame(w=as.vector(w), b=as.vector(b)), b != Inf & b != -Inf)
		k <- weighted.median(w=data$w, x=data$b)
		d <- median(y-k*x)
		e <- y-k*x-d
		resid <- sum((rank(e) - 1/2*(length(y)+1))*e)
		return(list(k=as.real(k), d=as.real(d), resid=as.real(resid)))
	}
	
	# function for internal ols regression
	ols <- function(y, x){
		res <- lm(y~x)
		return(list(k = as.real(res$coef[2]), 
	    		d = as.real(res$coef[1]), 
	    		resid = as.real(sum(res$resid*res$resid))))

	}

	# exclude missing values
	data <- na.omit(data.frame(conc, time))
	
	# subtraction of pre administration concentration for single dose studies
	if(!is.real(prev)){stop('argument prev invalid')}
	if(prev<0){stop('pre-dosing value must be greater 0')}
	if (prev > 0) {data$conc <- data$conc - prev}

	# check input parameters and remove values below or equal to zero
        method = match.arg(method)
	if(!is.vector(data$time)){stop('argument time invalid')}
	if(!is.vector(data$conc)){stop('argument conc invalid')}
	if (any(data$time < 0)) {stop('timepoint below zero')}
	if (points < 2) {stop('not enough points in terminal phase')}
	if(!is.logical(longer.terminal)){stop('argument longer.terminal invalid')}
	if(!is.real(points) || points%%1!=0){stop('argument points invalid')}
	if(length(unique(data$time))!=length(data$time)){stop('limited for one observation per time point')}
	
	
	# remove values below or equal to zero
	if (any(data$conc <= 0)) {
		for (i in 1:nrow(data)) {
			if (data$conc[i] <= 0) {data$conc[i] <- NA}
		}
		warning('concentration below or equal to zero were omitted')
		data <- na.omit(data)	
	}
	if (nrow(data) < 4) {stop('a minimum of 4 observations are required')}
  
	# transform data by logarithm at base 10
	n <- nrow(data)
 	data$conc <- log10(data$conc)
	conc <- data$conc
	time <- data$time
	
	# calculate parameters of one-phase model
	switch(method, 
		"lad"={model <- lad(y=conc, x=time)}, 
		"ols"={model <- ols(y=conc, x=time)}, 
		"hub"={model <- hub(y=conc, x=time, mad=lad(y=conc, x=time)$mad)},
		"npr"={model <- npr(y=conc, x=time)
	},)

	# inital halflife = terminal halflife for one-phase model
	final.term.model <- model
	final.init.model <- model
	resid <- model$resid
	final.chgpt <- NA

	# check special cases
	if(model$k >= 0) {
		resid <- Inf
		final.term.model$k <- NA
		final.init.model$k <- NA
	}
		
	# calculate parameters of two-phase models
	if (points > n-2) {stop('not enough points for inital phase')}
	for (i in 2:(n-points)) {

		init.conc <- conc[1:i]
		init.time <- time[1:i]
		term.conc <- conc[(i+1):n]
		term.time <- time[(i+1):n]

		# calculate parameters of two-phase model
		switch(method, 
			"lad"={
				init.model <- lad(y=init.conc, x=init.time)
				term.model <- lad(y=term.conc, x=term.time)
		},	"ols"={
				init.model <- ols(y=init.conc, x=init.time)
				term.model <- ols(y=term.conc, x=term.time)
		}, 	"hub"={
				init.model <- hub(y=init.conc, x=init.time, mad=lad(y=init.conc, x=init.time)$mad)
				term.model <- hub(y=term.conc, x=term.time, mad=lad(y=term.conc, x=term.time)$mad)
		}, 	"npr"={
				init.model <- npr(y=init.conc, x=init.time)
				term.model <- npr(y=term.conc, x=term.time)
		},)


		# check changeover criteria and for negative slopes 
		lower <- data$time[i]
		upper <- data$time[i+1]
		chgpt <- (init.model$d - term.model$d) / (term.model$k - init.model$k)
		 if (!(chgpt <= lower | chgpt >= upper) &
			(term.model$k < 0) & (init.model$k < 0)) {

			if(!longer.terminal){
  				if (sum(term.model$resid, init.model$resid) < resid) {
					final.init.model <- init.model
					final.term.model <- term.model
					final.chgpt <- as.real(chgpt)
					resid <- sum(term.model$resid, init.model$resid)
				}  
			}

			if(longer.terminal){
  				if (init.model$k <= term.model$k & sum(term.model$resid, init.model$resid) < resid) {
					final.init.model <- init.model
					final.term.model <- term.model
					final.chgpt <- as.real(chgpt)
					resid <- sum(term.model$resid, init.model$resid)
				}  
			}


		}
	}	

	init.hl <- -log10(2)/final.init.model$k
	term.hl <- -log10(2)/final.term.model$k
	if(is.na(init.hl) | is.na(term.hl)){warning('No model evaluated')}

	# format output objects
	parms <- data.frame(initial=as.real(c(init.hl, final.init.model$k, final.init.model$d)),
			terminal=as.real(c(term.hl, final.term.model$k, final.term.model$d)))
	rownames(parms) <- c('halflife', 'slope', 'intercept')
	res <- list(parms=parms, chgpt=as.real(final.chgpt), conc=10**conc, time=time, method='lee')
	class(res) <- 'halflife'
	return(res)
	
}

