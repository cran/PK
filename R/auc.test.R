auc.test <- function(conc, time, theta=0, group=NULL, alternative=c("two.sided", "less", "greater"), conf.level=0.95, design=c('ssd','batch','complete')){

  if(conf.level<0.5) {
    conf.level<-1-conf.level
    warning(paste('It appears that you provided a significance level rather than a confidence level. The confidence level is therefore changed to ',conf.level,sep=''))
  }

  alternative <- match.arg(alternative)
  design <- match.arg(design)
  if(design=='complete'){
    print('Hypothesis testing is not implemented for complete data designs. Please see the help file auc.complete.ci for examples on how to find confidence intervals for this case.')
  }else{
    if(design=='batch'){
      cis <- auc.batch.ci(conc=conc, time=time, group=group, method='t', alternative=alternative, conf.level=conf.level, nsample=1)
    }else{
      if(is.list(conc) || is.list(time)) stop('Both time and concentration need to be a vector')
      cis <- auc.ssd.ci(conc=conc, time=time, group=group, method='t', alternative=alternative, conf.level=conf.level, nsample=1)
    }
  } 
  df<-min(cis$CIs$df,na.rm=TRUE)
  t.stat <- (cis$CIs$est[1]-theta)/cis$CIs$stderr[1]
  if(alternative=='two.sided'){
    pval <- 2*(1-pt(abs(t.stat),df=df))
  }else{
    if(alternative=='less'){
      pval <- 1-pt(t.stat,df=df)
    }else{
      pval <- pt(t.stat,df=df)
    }
  }
  res <- list(est=cis$est)
  res$design<-design
  res$CIs<-cis$CIs
  res$test<-data.frame(p.value=pval,theta=theta,alternative=alternative,conf.level=conf.level)
  class(res)<-"PK"
  return(res)
}
