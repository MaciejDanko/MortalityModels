#' @noRd
repmat <- function(X,m,n){
  mx = dim(X)[1]
  nx = dim(X)[2]
  if (is.null(mx)){
    nx=length(X)
    mx=1
  }
  return(base:::matrix(t(base:::matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T))
}