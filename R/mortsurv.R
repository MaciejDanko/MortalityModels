#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
aftg_gompertz_makeham_mus <-function(a,b,c,g,x,mzb=1,lzb=1e-10,hzb=5,steps=500)
  C_aftg_gompertz_makeham_mus(a,b,c,g,x,
                              mzb,lzb,hzb,steps)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
aftg_gompertz_makeham_mu <-function(a,b,c,g,x,mzb=1,lzb=1e-10,hzb=5,steps=500)
  C_aftg_gompertz_makeham_mus(a,b,c,g,x,mzb,lzb,hzb,steps)$mu

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
aftg_gompertz_makeham_s <-function(a,b,c,g,x,mzb=1,lzb=1e-10,hzb=5,steps=500)
  C_aftg_gompertz_makeham_mus(a,b,c,g,x,mzb,lzb,hzb,steps)$s

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
gompertz_mu<-function(a,b,x) C_gompertz_mu(a,b,x)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
gompertz_s<-function(a,b,x) C_gompertz_s(a,b,x)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
gompertz_makeham_mu<-function(a,b,c,x) C_gompertz_makeham_mu(a,b,c,x)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
gompertz_makeham_s<-function(a,b,c,x) C_gompertz_makeham_s(a,b,c,x)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
phg_gompertz_mu<-function(a,b,g,x) C_phg_gompertz_mu(a,b,c,g,x)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
phg_gompertz_s<-function(a,b,g,x) C_phg_gompertz_s(a,b,g,x)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
phg_gompertz_makeham_mu<-function(a,b,c,g,x) 
  C_phg_gompertz_makeham_mu(a,b,c,g,x)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
phg_gompertz_makeham_s<-function(a,b,c,g,x) C_phg_gompertz_makeham_s(a,b,c,g,x)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
LL_poisson_gompertz_makeham<-function(par,x,I,E,
                                      lo.par=rep(1e-12,length(par)),
                                      hi.par=rep(0.75,length(par)),
                                      neg=TRUE) 
  C_LL_poisson_gompertz_makeham(par,x,I,E,lo.par,hi.par,neg)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
LL_poisson_gompertz<-function(par,x,I,E,
                              lo.par=rep(1e-12,length(par)),
                              hi.par=rep(0.75,length(par)),
                              neg=TRUE) 
  C_LL_poisson_gompertz(par,x,I,E,lo.par,hi.par,neg)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
LL_poisson_pgh_gompertz<-function(par,x,I,E,
                                  lo.par=rep(1e-12,length(par)),
                                  hi.par=rep(0.75,length(par)),
                                  neg=TRUE) 
  C_LL_poisson_phg_gompertz(par,x,I,E,lo.par,hi.par,neg)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
LL_poisson_phg_gompertz<-function(par,x,I,E,
                                  lo.par=rep(1e-12,length(par)),
                                  hi.par=rep(0.75,length(par)),
                                  neg=TRUE) 
  C_LL_poisson_phg_gompertz(par,x,I,E,lo.par,hi.par,neg)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
LL_poisson_phg_gompertz_makeham<-function(par,x,I,E,
                                          lo.par=rep(1e-12,length(par)),
                                          hi.par=rep(0.75,length(par)),
                                          neg=TRUE) 
  C_LL_poisson_phg_gompertz_makeham(par,x,I,E,lo.par,hi.par,neg)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
LL_poisson_aftg_gompertz<-function(par,x,I,E,
                                   lo.par=rep(1e-12,length(par)),
                                   hi.par=rep(0.75,length(par)),
                                   neg=TRUE,
                                   mzb=1,
                                   lzb=1e-10,
                                   hzb=5,
                                   steps=500) 
  C_LL_poisson_aftg_gompertz(par,x,I,E,lo.par,hi.par,mzb,lzb,hzb,steps,neg)


#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
LL_poisson_aftg_gompertz_makeham<-function(par,x,I,E,
                                           lo.par=rep(1e-12,length(par)),
                                           hi.par=rep(0.75,length(par)),
                                           neg=TRUE,
                                           mzb=1,
                                           lzb=1e-10,
                                           hzb=5,
                                           steps=500) {
  
  C_LL_poisson_aftg_gompertz_makeham(par,x,I,E,lo.par,hi.par,mzb,lzb,hzb,
                                     steps,neg)
}
