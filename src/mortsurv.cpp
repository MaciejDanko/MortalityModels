#include <stdio.h>
#include <RcppEigen.h>
using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
SEXP aftg_gompertz_makeham_mus(const double a, 
                               const double b, 
                               const double c,
                               const double s,
                               const Eigen::VectorXd x,
                               const double mzb,
                               const double lzb,
                               const double hzb,
                               const unsigned int steps){
  
  Eigen::VectorXd MU(x.size());
  Eigen::VectorXd SX(x.size());
  Eigen::ArrayXd gpdf(steps);
  Eigen::ArrayXd zb(steps);
  double bmzb = b*mzb;
  if(s<1e-11)  {
    MU = a*(x*bmzb).array().exp();
    SX = (-a*((x.array()*bmzb).exp() - 1)/bmzb).exp();
    gpdf = NA_REAL;
    zb = NA_REAL;
  } else {
    double shape = mzb*mzb/(s*s);
    double scale = s*s/mzb;
    double db = (hzb-lzb) / (steps-1);
    for (int i = 0L; i < steps; ++i) zb(i) = lzb + db*i;
    gpdf = db*zb.pow(shape-1).cwiseProduct((-zb/scale).exp());
    gpdf = gpdf / gpdf.sum();
    Eigen::ArrayXd H(zb.size());
    Eigen::ArrayXd S(zb.size());
    Eigen::ArrayXd bzba = b*zb;
    Eigen::ArrayXd azba = a*zb;
    for (int i = 0L; i < x.size(); ++i) {
      H = azba.cwiseProduct((x(i)*bzba).exp()); //a*zb * exp(x * b * zb)
      S = (-c*x(i)-azba*((x(i)*bzba).exp() - 1)/bzba/mzb).exp(); // 
      SX(i) = (gpdf.cwiseProduct(S)).sum();
      MU(i) = (gpdf.cwiseProduct(H.cwiseProduct(S))).sum()/SX(i)+c;
    }
  }
  return List::create(
    Named("mu") = MU,
    Named("s") = SX,
    Named("pdf") = gpdf.matrix(),
    Named("z") = zb.matrix()
  );
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd gompertz_mu(const double a, 
                            const double b, 
                            const Eigen::VectorXd x){   
  return(a*(x*b).array().exp());
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd gompertz_s(const double a, 
                           const double b, 
                           const Eigen::VectorXd x){   
  return((-a*((x.array()*b).exp() - 1)/b).exp());
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd gompertz_makeham_mu(const double a, 
                                    const double b, 
                                    const double c, 
                                    const Eigen::VectorXd x){   
  return(a*(x*b).array().exp()+c);
}
 
 
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd gompertz_makeham_s(const double a, 
                           const double b, 
                           const double c, 
                           const Eigen::VectorXd x){   
  return((-c*x.array()-a*((x.array()*b).exp() - 1)/b).exp());
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd phg_gompertz_mu(const double a, 
                                const double b, 
                                const double s,
                                const Eigen::VectorXd x){   
  if (s>1e-11){
    Eigen::VectorXd expbx = Eigen::VectorXd((x.array()*b).exp());
    return(a*expbx.array()*(1+a*s*(expbx.array()-1)/b).cwiseInverse());
  } else {
    return(a*(x*b).array().exp());
  }
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd phg_gompertz_s(const double a, 
                                const double b, 
                                const double s,
                                const Eigen::VectorXd x){   
  if (s>1e-11){
    Eigen::VectorXd expbx = Eigen::VectorXd((x.array()*b).exp());
    return((1+a*s*(expbx.array()-1)/b).pow(-1/s));
  } else {
    return((-a*((x.array()*b).exp() - 1)/b).exp());
  }
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd phg_gompertz_makeham_mu(const double a, 
                                        const double b, 
                                        const double c,
                                        const double s,
                                        const Eigen::VectorXd x){   
  if (s>1e-11){
    Eigen::VectorXd expbx = Eigen::VectorXd((x.array()*b).exp());
    return(a*expbx.array()*(1+a*s*(expbx.array()-1)/b).cwiseInverse()+c);
  } else {
    return(a*(x*b).array().exp()+c);
  }
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd phg_gompertz_makeham_s(const double a, 
                               const double b, 
                               const double c,
                               const double s,
                               const Eigen::VectorXd x){   
  if (s>1e-11){
    return((-c*x.array()).exp().cwiseProduct((1+
           a*s*((x.array()*b).exp()-1)/b).pow(-1/s)));
  } else {
    return((-c*x.array()-a*((x.array()*b).exp() - 1)/b).exp());
  }
}


