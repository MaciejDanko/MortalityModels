#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
aftg_gompertz_makeham_mus <-function(a,b,c,g,x,mzb=1,lzb=1e-10,hzb=5,steps=500)
  C_aftg_gompertz_makeham_mus(a,b,c,g,x,mzb,lzb,hzb,steps)

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
phg_gompertz_makeham_mu<-function(a,b,c,g,x) C_phg_gompertz_makeham_mu(a,b,c,g,x)

#' @useDynLib MortalityModels
#' @importFrom Rcpp evalCpp
#' @export
#' @noRd
phg_gompertz_makeham_s<-function(a,b,c,g,x) C_phg_gompertz_makeham_s(a,b,c,g,x)
