#' @noRd
Inv.C.surv<-function(Sx,a) -log(Sx)/a

#' @noRd
Inv.G.surv<-function(Sx,a,b) return(log((1-log(Sx)*(b/a)))/b)

#' @noRd
Inv.GG.surv<-function(Sx,a,b,s) if (s<1e-12) {
  return(Inv.G.surv(Sx,a,b))
} else return(log((Sx^(-s)-1)*b/a/s+1)/b)

#' @noRd
Inv.GM.surv<-function(Sx,a,b,c) {
  z=(a*exp(-(b*log(Sx) - a)/c))/c #detect numerical problems, they are close to Gompertz
  Index=((c<=0) | (z>1e280))
  InvSur=Sx*NA
  InvSur[Index]=log((1-log(Sx[Index])*(b/a)))/b #Gompertz
  InvSur[!Index]= -(log(Sx[!Index]) + (c*gsl::lambert_W0( z[!Index]))/b - (a/b))/c
  return(InvSur)
}

#' @noRd
#' @export
invGomp<-function(a,b,N,Sx=runif(N)) log((1-log(Sx)*(b/a)))/b

#' @noRd
#' @export
invGompMak<-function(a,b,c,N,Sx=runif(N)) Inv.GM.surv(Sx,a,b,c)

#' @noRd
#' @export
invPHGamGomp<-function(a,b,sa,N,Sx=runif(N)) log((a*sa - b + b/Sx^sa)/(a*sa))/b

#' @noRd
#' @export
invPHGamGompMak<-function(a,b,c,sa,N,Sx=runif(N),mean_zb=1){
  dgam.a=(mean_zb*mean_zb)/(sa) 
  dgam.b=(sa)/mean_zb
  rangam=rgamma(n=N,shape=dgam.a,scale=dgam.b)
  times=invGompMak(a=a*rangam,b=b,c=c,N=N)
  return(times)
}

#' @noRd
#' @export
invAFTGamGomp<-function(a,b,sb,N,mean_zb=1){
  dgam.a=(mean_zb*mean_zb)/(sb^2) 
  dgam.b=(sb^2)/mean_zb
  rangam=rgamma(n=N,shape=dgam.a,scale=dgam.b)
  times=invGomp(a=a*rangam,b=b*rangam,N=N)
  return(times)
}

#' @noRd
#' @export
invAFTGamGompMak<-function(a,b,c,sb,N,mean_zb=1){
  dgam.a=(mean_zb*mean_zb)/(sb^2) 
  dgam.b=(sb^2)/mean_zb
  rangam=rgamma(n=N,shape=dgam.a,scale=dgam.b)
  times=invGompMak(a=a*rangam,b=b*rangam,c=c,N=N)
  return(times)
}
