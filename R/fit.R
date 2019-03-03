#Likelihood rato test
#' @noRd
#' @export
lrTest<-function(LLRestricted,LLFull,DoF,alpha=0.05,boundaryF=F){
  Ratio   =  2 * (LLFull - LLRestricted)
  pValue  =  1 - stats:::pchisq(Ratio , DoF)
  criticalValue  =  stats:::qchisq(1 - alpha , DoF)
  if (boundaryF==T){
    pValue = pValue/2
    criticalValue  = NA
  }
  return(list(Ratio=Ratio,pValue=pValue,criticalValue=criticalValue))
}

#information criteria and their weights
#' @noRd
#' @export
AIC.BIC<-function(vLL,vnPar, Observ=NULL){
  res=NULL
  if (is.null(Observ)){
    AICVec=-2*vLL + 2*vnPar
    BICVec=NULL; BICWeights=NULL;
    DeltaVec=AICVec-min(AICVec);
    AICWeights=exp(-0.5*DeltaVec)/sum(exp(-0.5*DeltaVec));
  } else {    
    AICVec=-2*vLL + 2*vnPar + (2*vnPar*(vnPar+1))/(Observ-vnPar-1);
    DeltaVec=AICVec-min(AICVec);
    AICWeights=exp(-0.5*DeltaVec)/sum(exp(-0.5*DeltaVec));
    BICVec=-2*vLL+vnPar*log(Observ);
    DeltaVec=BICVec-min(BICVec);
    BICWeights=exp(-0.5*DeltaVec)/sum(exp(-0.5*DeltaVec));
  }
  res$AIC.vec=AICVec
  res$AIC.weights=AICWeights
  res$BIC.vec=BICVec
  res$BIC.weights=BICWeights
  res=as.data.frame(res)
  return(res)
}

#keep your estimates in bounds
#' @noRd
CorrectEstm<-function(Estm,lbvec,ubvec){
  if (any(is.na(Estm))) { #NA or NaN
    cat(paste('Wrong MLEs:',Estm,'\n'))
    stop
  }
  Estm[Estm<lbvec]<-lbvec[Estm<lbvec]
  Estm[Estm>ubvec]<-ubvec[Estm>ubvec]
  return(Estm)
}

#negative poisson likelihood function
#' @noRd
#' @export
neg.Poisson.LLike <- function(MortFunc,par,x,I,E,lbvec,ubvec) { 
  par=CorrectEstm(par,lbvec,ubvec)
  miun <- MortFunc(par,x)
  out <- -sum(I*log(miun)-E*miun) #Giancarlo used this: sum(I*log(miun*E)-E*miun), but this should give similar results     
  if (is.na(out)){
    out=1e200
  }
  if (abs(out)==Inf){
    out=1e200
  } 
  return(out)
}

#' @noRd
#' @export
optLL<-function(inipar.list,MortFunc,x,I,E,fn=neg.Poisson.LLike,lbvec,ubvec,doprecmax=25,doEA=T){
  #remove empty elements from the list
  inipar.list[sapply(1:length(inipar.list),function(X) is.null(inipar.list[[X]]))]=NULL
  if (sum((sapply(1:length(inipar.list),function(X) length(inipar.list[[X]])))-length(inipar.list[[1]]))>0){
    cat('Unequal length of initial parameters\n')
    stop
  }
  if (doEA){
    tmp=DEoptim(fn=fn,lower=lbvec,upper=ubvec,MortFunc=MortFunc, x=x, I=I, E=E, lbvec=lbvec, ubvec=ubvec,
                control = DEoptim.control(itermax=max(150,length(lbvec)*75),trace=F))
    inipar.list=c(list(tmp$optim$bestmem),inipar.list)
    EAvalue=-tmp$optim$bestval
    EApar=tmp$optim$bestmem
  } else{
    EAvalue=NA
    EApar=NA
  }   
  MaxLL=-1e200
  for (k in 1:length(inipar.list)){
    mypar=inipar.list[[k]]
    cat(k)
    LastLL=-1e200
    for (moreprec in 1:doprecmax){
      opt<-optim(par=mypar, MortFunc=MortFunc, x=x, I=I, E=E, fn=fn, lbvec=lbvec, ubvec=ubvec, hessian=F)
      mypar=CorrectEstm(opt$par,lbvec=lbvec, ubvec=ubvec)
      LL=-opt$value
      if (LL==LastLL) break
      cat(paste('(',LL-LastLL,')',sep=''))
      LastLL=-opt$value
    }
    LL=-opt$value
    if (MaxLL<LL){
      Max.opt=mypar
      MaxLL=LL
    }
  }
  cat('\n')
  return(list(value=MaxLL,par=Max.opt,EAvalue=EAvalue,EApar=EApar))
}
