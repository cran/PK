auc <- function(conc, time, exact=NA, n.interpolate=2, n.tail=3, design=c("ssd","batch","complete")) {
  design <- match.arg(design)
  if(design=='complete'){
    if(is.list(conc) || is.list(time)) stop('Both time and concentration need to be a vector')
    auc.complete(conc=conc, time=time, exact=exact, n.interpolate=n.interpolate, n.tail=n.tail)
  }else{
    if(design=='batch'){
      auc.batch(conc=conc, time=time, exact=exact, n.interpolate=n.interpolate, n.tail=n.tail)
    }else{
      if(is.list(conc) || is.list(time)) stop('Both time and concentration need to be a vector')
      auc.ssd(conc=conc, time=time, exact=exact, n.interpolate=n.interpolate, n.tail=n.tail)
    }
  }
}
