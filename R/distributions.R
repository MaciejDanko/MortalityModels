#Basic distributions to be rewritten in c++

#Gompertz
gomp<-function(a,b,x){
  MU=a*exp(b*x)
  SX=exp(-a*(exp(b*x)-1)/b)  
  return(list(MuBar=MU,SBar=SX))
}

#Gompertz-Makeham
gomp.mak<-function(a,b,c,x){
  MU=a*exp(b*x)+c
  SX=exp(-a*(exp(b*x)-1)/b)*exp(-c*c)  
  return(list(MuBar=MU,SBar=SX))
}

#PH-Gamma-Gompertz  
PH.gam.gomp<-function(a,b,sa,x) {
  if (sa>0){
    MU=a*exp(b*x)/(1+a*sa*(exp(b*x)-1)/b)
    SX=(1/(1+a*sa*(exp(b*x)-1)/b))^(1/sa)
  } else {
    MU=a*exp(b*x)
    SX=exp(-a*(exp(b*x)-1)/b)
  }  
  return(list(MuBar=MU,SBar=SX))
}

#PH-Gamma-Gompertz-Makeham  
PH.gam.gomp.mak<-function(a,b,c,sa,x) {
  if (sa>0){
    MU=(a*exp(b*x)/(1+a*sa*(exp(b*x)-1)/b))+c
    SX=exp(-c*x)*((1/(1+a*sa*(exp(b*x)-1)/b))^(1/sa))
  } else {
    MU=a*exp(b*x)+c
    SX=exp(-c*x)*exp(-a*(exp(b*x)-1)/b)
  }  
  return(list(MuBar=MU,SBar=SX))
}

#AFT-Gamma-Gompertz  , no heterogeneity in a
AFT.gam.gomp<-function(a,b,sb,x,mean_zb=1,steps=500,lo.zb=10^-10,hi.zb=5){
  gomp.mu<-function(A,B,x) A*exp(B*x)
  gomp.sx<-function(A,B,x) exp(-A*(exp(B*x)-1)/B)
  if (sb<=0){
    MU=gomp.mu(a,b*mean_zb,x)
    SX=gomp.sx(a,b*mean_zb,x)
    Correction=1
  } else {
    dgam.a=(mean_zb*mean_zb)/(sb^2) 
    dgam.b=(sb^2)/mean_zb
    zb=seq(from=lo.zb,to=hi.zb,length.out=steps)
    db=diff(zb)[1]
    gpdf=dgamma(x=zb,shape=dgam.a,scale=dgam.b)*db
    Correction=sum(gpdf)  #should be 1 anyway
    gpdf=gpdf/Correction
    zb=zb[gpdf>0]
    gpdf=gpdf[gpdf>0]
    A=zb*a
    B=zb*b
    A=repmat(X=A,m=length(x),n=1)
    B=repmat(X=B,m=length(x),n=1)
    gpdf=repmat(X=gpdf,m=length(x),n=1)
    X=t(repmat(X=x,m=length(zb),n=1))
    H=gomp.mu(A,B,x=X)
    S=gomp.sx(A,B,x=X)
    SX=rowSums(S*gpdf)
    Fx=rowSums(S*H*gpdf)
    MU=Fx/SX
  }
  return(list(MuBar=MU,SBar=SX,DistrMissed=1-Correction))
}

#AFT-Gamma-Gompertz-Makeham  , no heterogeneity in a
AFT.gam.gomp.mak<-function(a,b,c,sb,x,mean_zb=1,steps=500,lo.zb=10^-10,hi.zb=5){
  gomp.mu<-function(A,B,x) A*exp(B*x)
  gomp.sx<-function(A,B,x) exp(-A*(exp(B*x)-1)/B)
  if (sb<=0){
    MU.M=gomp.mu(a,b*mean_zb,x)+c
    SX.M=gomp.sx(a,b*mean_zb,x)*exp(-x*c)
    Correction=1
  } else {
    dgam.a=(mean_zb*mean_zb)/(sb^2) 
    dgam.b=(sb^2)/mean_zb
    zb=seq(from=lo.zb,to=hi.zb,length.out=steps)
    db=diff(zb)[1]
    gpdf=dgamma(x=zb,shape=dgam.a,scale=dgam.b)*db
    Correction=sum(gpdf)
    gpdf=gpdf/Correction
    zb=zb[gpdf>0]
    gpdf=gpdf[gpdf>0]
    A=zb*a
    B=zb*b
    A=repmat(X=A,m=length(x),n=1)
    B=repmat(X=B,m=length(x),n=1)
    gpdf=repmat(X=gpdf,m=length(x),n=1)
    X=t(repmat(X=x,m=length(zb),n=1))
    H=gomp.mu(A,B,x=X)
    S=gomp.sx(A,B,x=X)
    SX=rowSums(S*gpdf)
    Fx=rowSums(S*H*gpdf)
    MU.M=Fx/SX+c
    SX.M=exp(-x*c)*SX
  }
  return(list(MuBar=MU.M,SBar=SX.M,z=zb,pdf=gpdf,DistrMissed=1-Correction))
}
