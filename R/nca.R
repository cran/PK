nca <- function(conc, time, n.tail=3, dose=0, method=c("z", "boott"), conf.level=0.95, nsample=1000, design=c("ssd","batch","complete"), data) {
  if(missing(design)) {stop("A design needs to be specified")}
  design <- match.arg(design)
  if(design=='complete'){
    print('Complete data design currently not implemented')
  }else{
    if(design=='batch'){
      print('Batch design currently not implemented')
    }else{
      if(!missing(conc) && (is.list(conc) || is.list(time))) stop('Both time and concentration need to be a vector')
      nca.ssd(conc, time, n.tail=n.tail, dose=dose, method=method, conf.level=conf.level, nsample=nsample, data)
    }
  }
}
