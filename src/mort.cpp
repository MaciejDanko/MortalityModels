#include <stdio.h>
#include <RcppEigen.h>
using namespace Rcpp;

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd gompertz_mu(const double a, 
                            const double b, 
                            const Eigen::VectorXd x){   
  Eigen::VectorXd res = Eigen::VectorXd(x*b);
  return(a*res.array().exp().matrix());
}



// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd gompertz_makeham_mu(const double a, 
                                    const double b, 
                                    const double c, 
                                    const Eigen::VectorXd x){   
  return(gompertz_mu(a,b,x).array() + c);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd gompertz_makeham_mu_par(const Eigen::VectorXd par,
                                    const Eigen::VectorXd x){   
  Eigen::VectorXd res = Eigen::VectorXd(x.array()*par(1));
  return((par(0)*res.array().exp()+par(2)));
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd phg_gompertz_makeham_mu_par(const Eigen::VectorXd par,
                                        const Eigen::VectorXd x){   
  if (par(3)>1e-10){
    Eigen::VectorXd expbx = Eigen::VectorXd((x.array()*par(1)).exp());
    Eigen::VectorXd denom = Eigen::VectorXd(1+par(0)*par(3)*(expbx.array()-1)/par(1));
    return(par(0)*expbx.array()*denom.array().cwiseInverse()+par(2));
  } else {
    Eigen::VectorXd res = Eigen::VectorXd(x.array()*par(1));
    return((par(0)*res.array().exp()+par(2)));
  }
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List aftg_gompertz_makeham_dst_par(const Eigen::VectorXd par,
                                              const Eigen::VectorXd x,
                                              const double mzb,
                                              const double lzb,
                                              const double hzb,
                                              const double steps){   
  
}
// AFT.gam.gomp<-function(a,b,sb,x,mean_zb=1,steps=500,lo.zb=10^-10,hi.zb=5){
//   gomp.mu<-function(A,B,x) A*exp(B*x)
//   gomp.sx<-function(A,B,x) exp(-A*(exp(B*x)-1)/B)
//   if (sb<=0){
//     MU=gomp.mu(a,b*mean_zb,x)
//     SX=gomp.sx(a,b*mean_zb,x)
//     Correction=1
//   } else {
//     dgam.a=(mean_zb*mean_zb)/(sb^2) 
//     dgam.b=(sb^2)/mean_zb
//     zb=seq(from=lo.zb,to=hi.zb,length.out=steps)
//     db=diff(zb)[1]
//     gpdf=dgamma(x=zb,shape=dgam.a,scale=dgam.b)*db
//     Correction=sum(gpdf)  #should be 1 anyway
//     gpdf=gpdf/Correction
//     zb=zb[gpdf>0]
//     gpdf=gpdf[gpdf>0]
//     A=zb*a
//     B=zb*b
//     A=repmat(X=A,m=length(x),n=1)
//     B=repmat(X=B,m=length(x),n=1)
//     gpdf=repmat(X=gpdf,m=length(x),n=1)
//     X=t(repmat(X=x,m=length(zb),n=1))
//     H=gomp.mu(A,B,x=X)
//     S=gomp.sx(A,B,x=X)
//     SX=rowSums(S*gpdf) ///
//     Fx=rowSums(S*H*gpdf) ///
//     MU=Fx/SX
//   }
//   return(list(MuBar=MU,SBar=SX,DistrMissed=1-Correction))
// }
// 
// 
