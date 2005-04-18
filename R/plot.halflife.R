

# plot function for halflife objects 
plot.halflife <- function(x, xlab='Time', ylab='Concentration', main='Half-life Estimation', xlim=NULL, ylim=NULL, ...) {

	init.k <- x$parms[2,1] 
	init.d <- x$parms[3,1] 
	term.k <- x$parms[2,2] 
	term.d <- x$parms[3,2] 

	if(is.null(xlim)){xlim <- c(min(x$time), max(x$time))}
	if(is.null(ylim)){ylim <- c(min(x$conc), max(x$conc))}

	plot(x$time, x$conc, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, ...)

	switch(x$method,
		"lee"={
			if (is.na(x$chgpt)) {x$chgpt <- min(x$time)}
			plot(function(x) 10**(init.k*x+init.d), xlim[1], x$chgpt, add=TRUE)
			plot(function(x) 10**(term.k*x+term.d), x$chgpt, xlim[2], add=TRUE) 	
	},	"biexp"={
			if (init.k == term.k){
				plot(function(x) init.d*exp(-init.k*x), xlim[1], xlim[2], add=TRUE)
			}
			if (init.k != term.k){
				plot(function(x) init.d*exp(-init.k*x)+term.d*exp(-term.k*x), xlim[1], xlim[2], add=TRUE)
			}
	},)

		
}	

