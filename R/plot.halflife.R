plot.halflife <- function (x, xlab="Time", ylab="Concentration", main="Half-life Estimation", xlim=NULL, ylim=NULL, add=FALSE, ...) {

    b1 <- x$parms[2,1]
    a1 <- x$parms[3,1]
    b2 <- x$parms[2,2]
    a2 <- x$parms[3,2]

    if (is.null(xlim)) {xlim <- c(min(x$time), max(x$time))}
    if (is.null(ylim)) {ylim <- c(min(x$conc), max(x$conc))}

    if(add==FALSE){plot(x=x$time, y=x$conc, xlim=c(xlim[1], xlim[2]), ylim=c(ylim[1], ylim[2]), xlab=xlab, ylab=ylab, main=main, ...)}
    if(add==TRUE){points(x=x$time, y=x$conc, ...)}

    switch(x$method, lee = {
        if (is.na(x$chgpt)) {x$chgpt <- min(x$time)}
        curve(10^(b1 * x + a1), from=xlim[1], to=x$chgpt, add=TRUE, ...)
        curve(10^(b2 * x + a2), from=x$chgpt, to=xlim[2], add=TRUE, ...)
      }, biexp = {
        if (b1 == b2) {curve(a1 * exp(-b1 * x), from=xlim[1], to=xlim[2], add=TRUE, ...)}
        if (b1 != b2) {
           curve(a1 * exp(-b1 * x) + a2 * exp(-b2 * x), from=xlim[1], to=xlim[2], add=TRUE, ...)
        }
    }, )
}