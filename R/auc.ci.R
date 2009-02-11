auc.ci <- function(conc, time, group=NULL, method=c("t", "z", "boott"), alternative=c("two.sided", "less", "greater"), conf.level=0.95, strata=NULL, nsample=1000, design=c('ssd','batch','complete')){

  design <- match.arg(design)
  if(design=='complete'){
    auc.complete.ci(conc=conc, time=time)
  }else{
    if(design=='batch'){
      if(!is.null(strata)) warning('Stratification will be ignored in a batch design')
      auc.batch.ci(conc=conc, time=time, group=group, method=method, alternative=alternative, conf.level=conf.level, nsample=nsample)
    }else{
      if(is.list(conc) || is.list(time)) stop('Both time and concentration need to be a vector')
      auc.ssd.ci(conc=conc, time=time, group=group, method=method, alternative=alternative, conf.level=conf.level, strata=strata, nsample=nsample)
    }
  }
}

