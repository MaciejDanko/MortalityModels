#' @noRd
#' @export
r_gompertz<-function(a,b,N,Sx=runif(N)) log((1-log(Sx)*(b/a)))/b

#' @noRd
#' @export
#' @importFrom gsl lambert_W0
r_gompertz_makeham <- function(a,b,c,N,Sx=runif(N)) {
  z <- (a*exp(-(b*log(Sx) - a)/c))/c 
  #detect numerical problems, they are close to Gompertz
  Index <- ((c<=0) | (z>1e280))
  InvSur <- rep(NA, length(Sx))
  InvSur[Index] <- log((1-log(Sx[Index])*(b/a)))/b #Gompertz
  InvSur[!Index] <- -(log(Sx[!Index]) + 
                        (c*gsl::lambert_W0( z[!Index]))/b - (a/b))/c
  InvSur
}

#' @noRd
#' @export
r_phg_gompertz <- function(a,b,sa,N,Sx=runif(N)) 
  log((a*sa - b + b/Sx^sa)/(a*sa))/b

#' @noRd
#' @export
#' @importFrom stats rgamma
r_phg_gompertz_makeham <- function(a,b,c,sa,N,Sx=runif(N),mean_zb=1){
  dgam.a <- (mean_zb*mean_zb)/(sa) 
  dgam.b <- sa/mean_zb
  rangam <- stats::rgamma(n=N,shape=dgam.a,scale=dgam.b)
  r_gompertz_makeham(a=a*rangam,b=b,c=c,N=N)
}

#' @noRd
#' @export
#' @importFrom stats rgamma
r_aftg_gompertz<-function(a,b,sb,N,mean_zb=1){
  dgam.a <- (mean_zb*mean_zb)/(sb^2) 
  dgam.b <- (sb^2)/mean_zb
  rangam <- stats::rgamma(n=N,shape=dgam.a,scale=dgam.b)
  r_gompertz(a=a*rangam,b=b*rangam,N=N)
}

#' @noRd
#' @export
#' @importFrom stats rgamma
r_aftg_gompertz_makeham<-function(a,b,c,sb,N,mean_zb=1){
  dgam.a <- (mean_zb*mean_zb)/(sb^2) 
  dgam.b <- (sb^2)/mean_zb
  rangam <- stats::rgamma(n=N,shape=dgam.a,scale=dgam.b)
  r_gompertz_makeham(a=a*rangam,b=b*rangam,c=c,N=N)
}
