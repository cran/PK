print.PK<-function(x,...){

  if(is.null(x$CIs) && is.null(x$test)){
    if(x$design=="complete"){
      cat("Estimation for a complete data design\n\n")
      print(x$est)
    }else{
      if(x$design=="batch"){
        cat("Estimation for a batch design\n\n")
        print(x$est)
      }else{
        cat("Estimation for a serial sampling design\n\n")
        print(x$est)
      }
    }
  }
  if (!is.null(x$CIs) && is.null(x$test)){
    if(x$design=="complete"){
      cat("Confidence intervals for a complete data design\n\n")
      print(x$CIs[,1:5])
    }else{
      if(x$design=="batch"){
        cat("Confidence intervals for a batch design\n\n")
        print(x$CIs[,1:5])
      }else{
        cat("Confidence intervals for a serial sampling design\n\n")
        print(x$CIs[,1:5])
      }
    }
  }
  if (!is.null(x$CIs) && !is.null(x$test)){
    if(x$test$alternative=='less') {
      hyp0<-paste('     H0: theta >=',x$test$theta,'\n')
      hyp1<-paste('     H1: theta < ',x$test$theta,'\n\n')
    }else{
      if(x$test$alternative=='greater') {
        hyp0<-paste('     H0: theta <=',x$test$theta,'\n')
        hyp1<-paste('     H1: theta > ',x$test$theta,'\n\n')
      }else{
        hyp0<-paste('     H0: theta = ',x$test$theta,'\n')
        hyp1<-paste('     H1: theta <>',x$test$theta,'\n\n')
      }
    }
    if(x$design=="complete"){
      cat("Hypothesis testing for a complete data design\n\n")
      cat(hyp0)
      cat(hyp1)
      print(x$CIs)
    }else{
      if(x$design=="batch"){
        cat("Hypothesis testing for a batch design\n\n")
      }else{
        cat("Hypothesis testing for a serial sampling design\n\n")
      }
    }
    cat(hyp0)
    cat(hyp1)
    cat('p-value:            ',x$test$p.value,'\n')
    cat('point estimate:     ',x$est[1,1],'\n')
    cat('degrees of freedom: ',x$CIs$df[which(!is.na(x$CIs$df))],'\n')
    cat('\n',x$test$conf.level*100,'% confidence interval:\n\n',sep='')
    for(i in 1:nrow(x$CIs)){
      cat(as.character(x$CIs[i,6]),'-interval: (',round(x$CIs[i,3],2),' ; ',round(x$CIs[i,4],2),')\n',sep='')
    }
    ifelse(x$test$p.value<=1-x$test$conf.level,reject<-'reject',reject<-'fail to reject')
    pval<-x$test$p.value
    ifelse(pval<0.0001,pval<-'<0.0001',pval<-format(pval,digits=1,nsmall=4))
    cat('\n\nThe p-value for this test is ',pval,' and so we ',reject,'\nthe null hypothesis at a significance level of ',1-x$test$conf.level,'.\n\n',sep='')
  }
}


