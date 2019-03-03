#' @noRd
#' @export
CalcLT<-function(Times,Events=rep(1,length(Times)),First.Bin=0,Bin.Step=1,Last.Bin=NULL,Lost.Events.Status=1){
  if (length(Times)!=length(Events)) stop
  K=sort(Times,index.return=T)
  Times=K$x
  Events=Events[K$ix]
  if (!is.null(Last.Bin)) {
    RemovedTimes=sum((Times[Times>=Last.Bin[1]]-Last.Bin[1]))
    RemovedTimesR=sum(floor((Times[Times>=Last.Bin[1]]-Last.Bin[1]))+Bin.Step/2)
    RemovedTimesNb=length(Times[Times>=Last.Bin[1]])
    Events[Times>=Last.Bin[1]]<-(Lost.Events.Status)*Events[Times>=Last.Bin[1]]
    Times[Times>=Last.Bin[1]]<-Last.Bin[1]
  } else {
    RemovedTimes=0
    RemovedTimesR=0
    RemovedTimesNb=0
  }
  if (First.Bin>min(Times)) First.Bin=floor(min(Times))
  x=seq(from=First.Bin,to=max(Times),by=Bin.Step)
  xn=paste(x)
  #xnh=paste(x+Bin.Step/2)
  if (!is.null(Last.Bin)) {
    if (x[length(x)]==Last.Bin) {
      xn[length(xn)]<-paste(xn[length(xn)],'+',collapse='',sep='')
      censored=T
    } else censored=F
  }
  #if (!is.null(Last.Bin)) if (x[length(x)]==Last.Bin) xnh[length(xnh)]<-paste(xnh[length(xnh)],'+',collapse='',sep='')
  if (x[length(x)]<max(Times)) x=c(x,x[length(x)]+Bin.Step)
  L.Times=lapply(X=1:(length(x)),FUN<-function(X) Times[(Times-x[X]-Bin.Step<0) & (Times-x[X]>=0)])
  L.Events=lapply(X=1:(length(x)),FUN<-function(X) Events[(Times-x[X]-Bin.Step<0) & (Times-x[X]>=0)])
  
  DX=sapply(X=1:length(L.Events),FUN<-function(X) sum(L.Events[[X]]))
  N=sapply(X=1:length(L.Events),FUN<-function(X) length(L.Events[[X]]))
  
  NewDeaths=unlist(lapply(X=1:length(x),FUN<-function(X) rep(x[X]+Bin.Step/2,DX[X])))
  NewCens=unlist(lapply(X=1:length(x),FUN<-function(X) rep(x[X]+Bin.Step/2,N[X]-DX[X])))
  NewTimes=c(NewDeaths,NewCens)
  NewEvents=c(rep(1,length(NewDeaths)),rep(1,length(NewCens)))
  
  NX=(rep(sum(N),length(N)+1)-c(0,cumsum(N)))[1:length(N)]
  
  ZX=sapply(X=1:length(L.Times),FUN<-function(X) sum(L.Times[[X]]-x[X]))
  
  EX=NX-N+ZX
  EX.R=NX-N*Bin.Step/2
  EX[length(EX)]=RemovedTimes*Lost.Events.Status
  EX.R[length(EX.R)]=RemovedTimesR*Lost.Events.Status
  MU=DX/EX
  LT<-data.frame(X=x,NX=NX,DX=DX,EX=EX,EX.R=EX.R,MU=MU,MU.X=x+Bin.Step/2,N.Events=N)
  LT<-LT[NX>0,]
  rownames(LT)<-xn
  
  names(L.Times)<-paste('x_',xn,sep='')
  names(L.Events)<-paste('x_',xn,sep='')
  return(list(Test=list(L.T=L.Times,L.E=L.Events),xn=xn,newdata=list(Times=NewTimes,Events=Events),
              LT=LT,RemovedEX=RemovedTimes,RemovedEvents=RemovedTimesNb,stacked=censored))
}

