#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::MatrixXd;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
SEXP C_aftg_gompertz_makeham_mus(const double a, 
                                 const double b, 
                                 const double c,
                                 const double s,
                                 const Eigen::VectorXd x,
                                 const double mzb,
                                 const double lzb,
                                 const double hzb,
                                 const unsigned int steps){
  
  VectorXd MU(x.size());
  VectorXd SX(x.size());
  ArrayXd gpdf(steps);
  ArrayXd zb(steps);
  double bmzb = b*mzb;
  if(s<1e-11)  {
    MU = a*(x*bmzb).array().exp()+c;
    SX = (-c*x.array()-a*((x.array()*bmzb).exp() - 1)/bmzb).exp();
    gpdf = NA_REAL;
    zb = NA_REAL;
  } else {
    double shape = mzb*mzb/(s*s);
    double scale = s*s/mzb;
    double db = (hzb-lzb) / (steps-1);
    zb = VectorXd::LinSpaced(steps,lzb,hzb); 
    gpdf = db*zb.pow(shape-1)*(-zb/scale).exp();
    gpdf = gpdf / gpdf.sum();
    ArrayXd H(zb.size());
    ArrayXd S(zb.size());
    ArrayXd bzba = b*zb;
    ArrayXd azba = a*zb;
    for (int i = 0L; i < x.size(); ++i) {
      H = azba*(x(i)*bzba).exp(); //a*zb * exp(x * b * zb)
      S = (-c*x(i)-azba*((x(i)*bzba).exp() - 1)/bzba/mzb).exp(); // 
      SX(i) = (gpdf*S).sum();
      MU(i) = (gpdf*H*S).sum()/SX(i)+c;
    }
  }
  return List::create(
    Named("mu") = MU,
    Named("s") = SX
  //  Named("pdf") = gpdf.matrix(),
  //  Named("z") = zb.matrix()
  );
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
SEXP C_aftg_gompertz_mus(const double a, 
                         const double b, 
                         const double s,
                         const Eigen::VectorXd x,
                         const double mzb,
                         const double lzb,
                         const double hzb,
                         const unsigned int steps){
  
  VectorXd MU(x.size());
  VectorXd SX(x.size());
  ArrayXd gpdf(steps);
  ArrayXd zb(steps);
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
    zb = VectorXd::LinSpaced(steps,lzb,hzb);
    gpdf = db*zb.pow(shape-1)*(-zb/scale).exp();
    gpdf = gpdf / gpdf.sum();
    ArrayXd H(zb.size());
    ArrayXd S(zb.size());
    ArrayXd bzba = b*zb;
    ArrayXd azba = a*zb;
    for (int i = 0L; i < x.size(); ++i) {
      H = azba*(x(i)*bzba).exp(); //a*zb * exp(x * b * zb)
      S = (-azba*((x(i)*bzba).exp() - 1)/bzba/mzb).exp(); // 
      SX(i) = (gpdf*S).sum();
      MU(i) = (gpdf*H*S).sum()/SX(i);
    }
  }
  return List::create(
    Named("mu") = MU,
    Named("s") = SX
    //Named("pdf") = gpdf.matrix(),
    //Named("z") = zb.matrix()
  );
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd C_gompertz_mu(const double a, 
                              const double b, 
                              const Eigen::VectorXd x){   
  return(a*(x*b).array().exp());
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd C_d_gompertz_mu(const double a, 
                              const double b, 
                              const Eigen::VectorXd x){   
  MatrixXd res(2,x.size());
  ArrayXd dmu_da = (x*b).array().exp();
  ArrayXd dmu_db = a*dmu_da.cwiseProduct(x.array());  
  res.row(0) = dmu_da;
  res.row(1) = dmu_db;
  return(res);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd C_gompertz_s(const double a, 
                             const double b, 
                             const Eigen::VectorXd x){   
  return((-a*((x.array()*b).exp() - 1)/b).exp());
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd C_gompertz_makeham_mu(const double a, 
                                      const double b, 
                                      const double c, 
                                      const Eigen::VectorXd x){   
  return(a*(x*b).array().exp()+c);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd C_d_gompertz_makeham_mu(const double a, 
                                const double b, 
                                const double c, 
                                const Eigen::VectorXd x){   
  MatrixXd res(3,x.size());
  ArrayXd dmu_da = (x*b).array().exp();
  ArrayXd dmu_db = a*dmu_da.cwiseProduct(x.array());  
  res.row(0) = dmu_da;
  res.row(1) = dmu_db;
  res.row(2) = VectorXd::Ones(x.size());
  return(res);
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd C_gompertz_makeham_s(const double a, 
                                     const double b, 
                                     const double c, 
                                     const Eigen::VectorXd x){   
  return((-c*x.array()-a*((x.array()*b).exp() - 1)/b).exp());
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd C_phg_gompertz_mu(const double a, 
                                  const double b, 
                                  const double s,
                                  const Eigen::VectorXd x){   
  if (s>1e-11){
    VectorXd expbx = VectorXd((x.array()*b).exp());
    return(a*expbx.array()*(1+a*s*(expbx.array()-1)/b).cwiseInverse());
  } else {
    return(a*(x*b).array().exp());
  }
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd C_d_phg_gompertz_mu(const double a, 
                                        const double b, 
                                        const double s, 
                                        const Eigen::VectorXd x){   
  MatrixXd res(3,x.size()); 
  ArrayXd tmp1 = (x*b).array().exp(); //2
  ArrayXd tmp2 = tmp1-1; //4
  ArrayXd tmp3 = tmp2 * s; //5
  ArrayXd tmp4 = 1/a + tmp2/b; //7
  ArrayXd tmp5 = tmp3.square(); //12
  ArrayXd tmp6 = tmp1 * x.array(); //14
  res.row(0) = (1/(a*a)) * tmp1 / tmp5; //dmu_da
  res.row(1) = tmp6/tmp4 - tmp1 * (s*tmp6/b - tmp3/(b*b)) / tmp5; //dmu_db
  res.row(2) = -(tmp1*(tmp2/b)/tmp5); //dmu_ds
  return(res);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd C_phg_gompertz_s(const double a, 
                                 const double b, 
                                 const double s,
                                 const Eigen::VectorXd x){   
  if (s>1e-11){
    VectorXd expbx = VectorXd((x.array()*b).exp());
    return((1+a*s*(expbx.array()-1)/b).pow(-1/s));
  } else {
    return((-a*((x.array()*b).exp() - 1)/b).exp());
  }
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd C_phg_gompertz_makeham_mu(const double a, 
                                          const double b, 
                                          const double c,
                                          const double s,
                                          const Eigen::VectorXd x){   
  if (s>1e-11){
    VectorXd expbx = VectorXd((x.array()*b).exp());
    return(a*expbx.array()*(1+a*s*(expbx.array()-1)/b).cwiseInverse()+c);
  } else {
    return(a*(x*b).array().exp()+c);
  }
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd C_d_phg_gompertz_makeham_mu(const double a, 
                                            const double b, 
                                            const double c, 
                                            const double s, 
                                            const Eigen::VectorXd x){   
  MatrixXd res(4,x.size()); 
  ArrayXd tmp1 = (x*b).array().exp(); //2
  ArrayXd tmp2 = tmp1-1; //4
  ArrayXd tmp3 = tmp2 * s; //5
  ArrayXd tmp4 = 1/a + tmp2/b; //7
  ArrayXd tmp5 = tmp3.square(); //12
  ArrayXd tmp6 = tmp1 * x.array(); //14
  res.row(0) = (1/(a*a)) * tmp1 / tmp5; //dmu_da
  res.row(1) = tmp6/tmp4 - tmp1 * (s*tmp6/b - tmp3/(b*b)) / tmp5; //dmu_db
  res.row(2) = VectorXd::Ones(x.size());
  res.row(3) = -(tmp1*(tmp2/b)/tmp5); //dmu_ds
  return(res);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd C_phg_gompertz_makeham_s(const double a, 
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


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double C_LL_poisson_aftg_gompertz_makeham(const Eigen::VectorXd par,
                                          const Eigen::VectorXd x,
                                          const Eigen::VectorXd I,
                                          const Eigen::VectorXd E,
                                          const Eigen::VectorXd lo_par,
                                          const Eigen::VectorXd hi_par,
                                          const double mzb,
                                          const double lzb,
                                          const double hzb,
                                          const unsigned int steps,
                                          const unsigned int neg){
  VectorXd PAR = par;
  for (int i = 0L; i < par.size(); ++i) {
    if (par(i)<lo_par(i)) PAR(i)=lo_par(i);
    if (par(i)>hi_par(i)) PAR(i)=hi_par(i);
  }
  List tmp =
    C_aftg_gompertz_makeham_mus(PAR(0),PAR(1),PAR(2),PAR(3),x,mzb,lzb,hzb,steps);
  ArrayXd mu = as<ArrayXd>(tmp["mu"]);
  double out = (I.array().cwiseProduct(mu.log())-E.array().cwiseProduct(mu)).sum();
  if (traits::is_nan<REALSXP>(out) || 
      NumericVector::is_na(out) || 
      traits::is_infinite<REALSXP>(out)) out = -1e200;
  if (neg) out = -out;
  return(out);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double C_LL_poisson_aftg_gompertz(const Eigen::VectorXd par,
                                  const Eigen::VectorXd x,
                                  const Eigen::VectorXd I,
                                  const Eigen::VectorXd E,
                                  const Eigen::VectorXd lo_par,
                                  const Eigen::VectorXd hi_par,
                                  const double mzb,
                                  const double lzb,
                                  const double hzb,
                                  const unsigned int steps,
                                  const unsigned int neg){
  VectorXd PAR = par;
  for (int i = 0L; i < par.size(); ++i) {
    if (par(i)<lo_par(i)) PAR(i)=lo_par(i);
    if (par(i)>hi_par(i)) PAR(i)=hi_par(i);
  }
  List tmp =
    C_aftg_gompertz_mus(PAR(0),PAR(1),PAR(2),x,mzb,lzb,hzb,steps);
  ArrayXd mu = as<ArrayXd>(tmp["mu"]);
  double out = (I.array().cwiseProduct(mu.log())-E.array().cwiseProduct(mu)).sum();
  if (traits::is_nan<REALSXP>(out) || 
      NumericVector::is_na(out) || 
      traits::is_infinite<REALSXP>(out)) out = -1e200;
  if (neg) out = -out;
  return(out);
}



// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double C_LL_poisson_phg_gompertz_makeham(const Eigen::VectorXd par, 
                                         const Eigen::VectorXd x,
                                         const Eigen::VectorXd I,
                                         const Eigen::VectorXd E,
                                         const Eigen::VectorXd lo_par, 
                                         const Eigen::VectorXd hi_par,
                                         const unsigned int neg){
  VectorXd PAR = par;
  for (int i = 0L; i < par.size(); ++i) {
    if (par(i)<lo_par(i)) PAR(i)=lo_par(i);
    if (par(i)>hi_par(i)) PAR(i)=hi_par(i);
  }
  ArrayXd mu = 
    C_phg_gompertz_makeham_mu(PAR(0),PAR(1),PAR(2),PAR(3),x);
  double out = (I.array().cwiseProduct(mu.log())-E.array().cwiseProduct(mu)).sum(); 
  if (traits::is_nan<REALSXP>(out) || 
      NumericVector::is_na(out) || 
      traits::is_infinite<REALSXP>(out)) out = -1e200;
  if (neg) out = -out;
  return(out);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double C_LL_poisson_phg_gompertz(const Eigen::VectorXd par, 
                                 const Eigen::VectorXd x,
                                 const Eigen::VectorXd I,
                                 const Eigen::VectorXd E,
                                 const Eigen::VectorXd lo_par, 
                                 const Eigen::VectorXd hi_par,
                                 const unsigned int neg){
  VectorXd PAR = par;
  for (int i = 0L; i < par.size(); ++i) {
    if (par(i)<lo_par(i)) PAR(i)=lo_par(i);
    if (par(i)>hi_par(i)) PAR(i)=hi_par(i);
  }
  ArrayXd mu = 
    C_phg_gompertz_mu(PAR(0),PAR(1),PAR(2),x);
  double out = (I.array().cwiseProduct(mu.log())-E.array().cwiseProduct(mu)).sum(); 
  if (traits::is_nan<REALSXP>(out) || 
      NumericVector::is_na(out) || 
      traits::is_infinite<REALSXP>(out)) out = -1e200;
  if (neg) out = -out;
  return(out);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double C_LL_poisson_gompertz(const Eigen::VectorXd par, 
                             const Eigen::VectorXd x,
                             const Eigen::VectorXd I,
                             const Eigen::VectorXd E,
                             const Eigen::VectorXd lo_par, 
                             const Eigen::VectorXd hi_par,
                             const unsigned int neg){
  VectorXd PAR = par;
  for (int i = 0L; i < par.size(); ++i) {
    if (par(i)<lo_par(i)) PAR(i)=lo_par(i);
    if (par(i)>hi_par(i)) PAR(i)=hi_par(i);
  }
  ArrayXd mu = 
    C_gompertz_mu(PAR(0),PAR(1),x);
  double out = (I.array().cwiseProduct(mu.log())-E.array().cwiseProduct(mu)).sum(); 
  if (traits::is_nan<REALSXP>(out) || 
      NumericVector::is_na(out) || 
      traits::is_infinite<REALSXP>(out)) out = -1e200;
  if (neg) out = -out;
  return(out);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double C_LL_poisson_gompertz_makeham(const Eigen::VectorXd par, 
                                     const Eigen::VectorXd x,
                                     const Eigen::VectorXd I,
                                     const Eigen::VectorXd E,
                                     const Eigen::VectorXd lo_par, 
                                     const Eigen::VectorXd hi_par,
                                     const unsigned int neg){
  VectorXd PAR = par;
  for (int i = 0L; i < par.size(); ++i) {
    if (par(i)<lo_par(i)) PAR(i)=lo_par(i);
    if (par(i)>hi_par(i)) PAR(i)=hi_par(i);
  }
  ArrayXd mu = 
    C_gompertz_makeham_mu(PAR(0),PAR(1),PAR(2),x);
  double out = (I.array().cwiseProduct(mu.log())-E.array().cwiseProduct(mu)).sum(); 
  if (traits::is_nan<REALSXP>(out) || 
      NumericVector::is_na(out) || 
      traits::is_infinite<REALSXP>(out)) out = -1e200;
  if (neg) out = -out;
  return(out);
}
