// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// C_aftg_gompertz_makeham_mus
SEXP C_aftg_gompertz_makeham_mus(const double a, const double b, const double c, const double s, const Eigen::VectorXd x, const double mzb, const double lzb, const double hzb, const unsigned int steps);
RcppExport SEXP _MortalityModels_C_aftg_gompertz_makeham_mus(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP sSEXP, SEXP xSEXP, SEXP mzbSEXP, SEXP lzbSEXP, SEXP hzbSEXP, SEXP stepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type mzb(mzbSEXP);
    Rcpp::traits::input_parameter< const double >::type lzb(lzbSEXP);
    Rcpp::traits::input_parameter< const double >::type hzb(hzbSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type steps(stepsSEXP);
    rcpp_result_gen = Rcpp::wrap(C_aftg_gompertz_makeham_mus(a, b, c, s, x, mzb, lzb, hzb, steps));
    return rcpp_result_gen;
END_RCPP
}
// C_gompertz_mu
Eigen::VectorXd C_gompertz_mu(const double a, const double b, const Eigen::VectorXd x);
RcppExport SEXP _MortalityModels_C_gompertz_mu(SEXP aSEXP, SEXP bSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_gompertz_mu(a, b, x));
    return rcpp_result_gen;
END_RCPP
}
// C_gompertz_s
Eigen::VectorXd C_gompertz_s(const double a, const double b, const Eigen::VectorXd x);
RcppExport SEXP _MortalityModels_C_gompertz_s(SEXP aSEXP, SEXP bSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_gompertz_s(a, b, x));
    return rcpp_result_gen;
END_RCPP
}
// C_gompertz_makeham_mu
Eigen::VectorXd C_gompertz_makeham_mu(const double a, const double b, const double c, const Eigen::VectorXd x);
RcppExport SEXP _MortalityModels_C_gompertz_makeham_mu(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_gompertz_makeham_mu(a, b, c, x));
    return rcpp_result_gen;
END_RCPP
}
// C_gompertz_makeham_s
Eigen::VectorXd C_gompertz_makeham_s(const double a, const double b, const double c, const Eigen::VectorXd x);
RcppExport SEXP _MortalityModels_C_gompertz_makeham_s(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_gompertz_makeham_s(a, b, c, x));
    return rcpp_result_gen;
END_RCPP
}
// C_phg_gompertz_mu
Eigen::VectorXd C_phg_gompertz_mu(const double a, const double b, const double s, const Eigen::VectorXd x);
RcppExport SEXP _MortalityModels_C_phg_gompertz_mu(SEXP aSEXP, SEXP bSEXP, SEXP sSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_phg_gompertz_mu(a, b, s, x));
    return rcpp_result_gen;
END_RCPP
}
// C_phg_gompertz_s
Eigen::VectorXd C_phg_gompertz_s(const double a, const double b, const double s, const Eigen::VectorXd x);
RcppExport SEXP _MortalityModels_C_phg_gompertz_s(SEXP aSEXP, SEXP bSEXP, SEXP sSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_phg_gompertz_s(a, b, s, x));
    return rcpp_result_gen;
END_RCPP
}
// C_phg_gompertz_makeham_mu
Eigen::VectorXd C_phg_gompertz_makeham_mu(const double a, const double b, const double c, const double s, const Eigen::VectorXd x);
RcppExport SEXP _MortalityModels_C_phg_gompertz_makeham_mu(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP sSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_phg_gompertz_makeham_mu(a, b, c, s, x));
    return rcpp_result_gen;
END_RCPP
}
// C_phg_gompertz_makeham_s
Eigen::VectorXd C_phg_gompertz_makeham_s(const double a, const double b, const double c, const double s, const Eigen::VectorXd x);
RcppExport SEXP _MortalityModels_C_phg_gompertz_makeham_s(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP, SEXP sSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const double >::type c(cSEXP);
    Rcpp::traits::input_parameter< const double >::type s(sSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_phg_gompertz_makeham_s(a, b, c, s, x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MortalityModels_C_aftg_gompertz_makeham_mus", (DL_FUNC) &_MortalityModels_C_aftg_gompertz_makeham_mus, 9},
    {"_MortalityModels_C_gompertz_mu", (DL_FUNC) &_MortalityModels_C_gompertz_mu, 3},
    {"_MortalityModels_C_gompertz_s", (DL_FUNC) &_MortalityModels_C_gompertz_s, 3},
    {"_MortalityModels_C_gompertz_makeham_mu", (DL_FUNC) &_MortalityModels_C_gompertz_makeham_mu, 4},
    {"_MortalityModels_C_gompertz_makeham_s", (DL_FUNC) &_MortalityModels_C_gompertz_makeham_s, 4},
    {"_MortalityModels_C_phg_gompertz_mu", (DL_FUNC) &_MortalityModels_C_phg_gompertz_mu, 4},
    {"_MortalityModels_C_phg_gompertz_s", (DL_FUNC) &_MortalityModels_C_phg_gompertz_s, 4},
    {"_MortalityModels_C_phg_gompertz_makeham_mu", (DL_FUNC) &_MortalityModels_C_phg_gompertz_makeham_mu, 5},
    {"_MortalityModels_C_phg_gompertz_makeham_s", (DL_FUNC) &_MortalityModels_C_phg_gompertz_makeham_s, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_MortalityModels(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
