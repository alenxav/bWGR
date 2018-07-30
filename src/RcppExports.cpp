// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// KMUP
SEXP KMUP(NumericMatrix X, NumericVector b, NumericVector d, NumericVector xx, NumericVector e, NumericVector L, double Ve, double pi);
RcppExport SEXP _bWGR_KMUP(SEXP XSEXP, SEXP bSEXP, SEXP dSEXP, SEXP xxSEXP, SEXP eSEXP, SEXP LSEXP, SEXP VeSEXP, SEXP piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type e(eSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type Ve(VeSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    rcpp_result_gen = Rcpp::wrap(KMUP(X, b, d, xx, e, L, Ve, pi));
    return rcpp_result_gen;
END_RCPP
}
// KMUP2
SEXP KMUP2(NumericMatrix X, NumericVector Use, NumericVector b, NumericVector d, NumericVector xx, NumericVector E, NumericVector L, double Ve, double pi);
RcppExport SEXP _bWGR_KMUP2(SEXP XSEXP, SEXP UseSEXP, SEXP bSEXP, SEXP dSEXP, SEXP xxSEXP, SEXP ESEXP, SEXP LSEXP, SEXP VeSEXP, SEXP piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Use(UseSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d(dSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type E(ESEXP);
    Rcpp::traits::input_parameter< NumericVector >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type Ve(VeSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    rcpp_result_gen = Rcpp::wrap(KMUP2(X, Use, b, d, xx, E, L, Ve, pi));
    return rcpp_result_gen;
END_RCPP
}
// emBA
SEXP emBA(NumericVector y, NumericMatrix gen, double df, double R2);
RcppExport SEXP _bWGR_emBA(SEXP ySEXP, SEXP genSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(emBA(y, gen, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// emBB
SEXP emBB(NumericVector y, NumericMatrix gen, double df, double R2, double Pi);
RcppExport SEXP _bWGR_emBB(SEXP ySEXP, SEXP genSEXP, SEXP dfSEXP, SEXP R2SEXP, SEXP PiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< double >::type Pi(PiSEXP);
    rcpp_result_gen = Rcpp::wrap(emBB(y, gen, df, R2, Pi));
    return rcpp_result_gen;
END_RCPP
}
// emBC
SEXP emBC(NumericVector y, NumericMatrix gen, double df, double R2, double Pi);
RcppExport SEXP _bWGR_emBC(SEXP ySEXP, SEXP genSEXP, SEXP dfSEXP, SEXP R2SEXP, SEXP PiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< double >::type Pi(PiSEXP);
    rcpp_result_gen = Rcpp::wrap(emBC(y, gen, df, R2, Pi));
    return rcpp_result_gen;
END_RCPP
}
// emRR
SEXP emRR(NumericVector y, NumericMatrix gen, double df, double R2);
RcppExport SEXP _bWGR_emRR(SEXP ySEXP, SEXP genSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(emRR(y, gen, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// emBL
SEXP emBL(NumericVector y, NumericMatrix gen, double R2, double alpha);
RcppExport SEXP _bWGR_emBL(SEXP ySEXP, SEXP genSEXP, SEXP R2SEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(emBL(y, gen, R2, alpha));
    return rcpp_result_gen;
END_RCPP
}
// emDE
SEXP emDE(NumericVector y, NumericMatrix gen, double R2);
RcppExport SEXP _bWGR_emDE(SEXP ySEXP, SEXP genSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(emDE(y, gen, R2));
    return rcpp_result_gen;
END_RCPP
}
// emEN
SEXP emEN(NumericVector y, NumericMatrix gen, double R2, double alpha);
RcppExport SEXP _bWGR_emEN(SEXP ySEXP, SEXP genSEXP, SEXP R2SEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type gen(genSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(emEN(y, gen, R2, alpha));
    return rcpp_result_gen;
END_RCPP
}
// BayesA
SEXP BayesA(NumericVector y, NumericMatrix X, double it, double bi, double df, double R2);
RcppExport SEXP _bWGR_BayesA(SEXP ySEXP, SEXP XSEXP, SEXP itSEXP, SEXP biSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type bi(biSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(BayesA(y, X, it, bi, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// BayesB
SEXP BayesB(NumericVector y, NumericMatrix X, double it, double bi, double pi, double df, double R2);
RcppExport SEXP _bWGR_BayesB(SEXP ySEXP, SEXP XSEXP, SEXP itSEXP, SEXP biSEXP, SEXP piSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type bi(biSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(BayesB(y, X, it, bi, pi, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// BayesC
SEXP BayesC(NumericVector y, NumericMatrix X, double it, double bi, double pi, double df, double R2);
RcppExport SEXP _bWGR_BayesC(SEXP ySEXP, SEXP XSEXP, SEXP itSEXP, SEXP biSEXP, SEXP piSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type bi(biSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(BayesC(y, X, it, bi, pi, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// BayesL
SEXP BayesL(NumericVector y, NumericMatrix X, double it, double bi, double df, double R2);
RcppExport SEXP _bWGR_BayesL(SEXP ySEXP, SEXP XSEXP, SEXP itSEXP, SEXP biSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type bi(biSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(BayesL(y, X, it, bi, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// BayesRR
SEXP BayesRR(NumericVector y, NumericMatrix X, double it, double bi, double df, double R2);
RcppExport SEXP _bWGR_BayesRR(SEXP ySEXP, SEXP XSEXP, SEXP itSEXP, SEXP biSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type bi(biSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(BayesRR(y, X, it, bi, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// BayesA2
SEXP BayesA2(NumericVector y, NumericMatrix X1, NumericMatrix X2, double it, double bi, double df, double R2);
RcppExport SEXP _bWGR_BayesA2(SEXP ySEXP, SEXP X1SEXP, SEXP X2SEXP, SEXP itSEXP, SEXP biSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< double >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type bi(biSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(BayesA2(y, X1, X2, it, bi, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// BayesB2
SEXP BayesB2(NumericVector y, NumericMatrix X1, NumericMatrix X2, double it, double bi, double pi, double df, double R2);
RcppExport SEXP _bWGR_BayesB2(SEXP ySEXP, SEXP X1SEXP, SEXP X2SEXP, SEXP itSEXP, SEXP biSEXP, SEXP piSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< double >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type bi(biSEXP);
    Rcpp::traits::input_parameter< double >::type pi(piSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(BayesB2(y, X1, X2, it, bi, pi, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// BayesRR2
SEXP BayesRR2(NumericVector y, NumericMatrix X1, NumericMatrix X2, double it, double bi, double df, double R2);
RcppExport SEXP _bWGR_BayesRR2(SEXP ySEXP, SEXP X1SEXP, SEXP X2SEXP, SEXP itSEXP, SEXP biSEXP, SEXP dfSEXP, SEXP R2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X1(X1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X2(X2SEXP);
    Rcpp::traits::input_parameter< double >::type it(itSEXP);
    Rcpp::traits::input_parameter< double >::type bi(biSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    Rcpp::traits::input_parameter< double >::type R2(R2SEXP);
    rcpp_result_gen = Rcpp::wrap(BayesRR2(y, X1, X2, it, bi, df, R2));
    return rcpp_result_gen;
END_RCPP
}
// CNT
void CNT(NumericMatrix X);
RcppExport SEXP _bWGR_CNT(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    CNT(X);
    return R_NilValue;
END_RCPP
}
// IMP
void IMP(NumericMatrix X);
RcppExport SEXP _bWGR_IMP(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    IMP(X);
    return R_NilValue;
END_RCPP
}
// GAU
NumericMatrix GAU(NumericMatrix X);
RcppExport SEXP _bWGR_GAU(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(GAU(X));
    return rcpp_result_gen;
END_RCPP
}
// SPC
NumericVector SPC(NumericVector y, NumericVector blk, NumericVector row, NumericVector col, int rN, int cN);
RcppExport SEXP _bWGR_SPC(SEXP ySEXP, SEXP blkSEXP, SEXP rowSEXP, SEXP colSEXP, SEXP rNSEXP, SEXP cNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type blk(blkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type row(rowSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type col(colSEXP);
    Rcpp::traits::input_parameter< int >::type rN(rNSEXP);
    Rcpp::traits::input_parameter< int >::type cN(cNSEXP);
    rcpp_result_gen = Rcpp::wrap(SPC(y, blk, row, col, rN, cN));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bWGR_KMUP", (DL_FUNC) &_bWGR_KMUP, 8},
    {"_bWGR_KMUP2", (DL_FUNC) &_bWGR_KMUP2, 9},
    {"_bWGR_emBA", (DL_FUNC) &_bWGR_emBA, 4},
    {"_bWGR_emBB", (DL_FUNC) &_bWGR_emBB, 5},
    {"_bWGR_emBC", (DL_FUNC) &_bWGR_emBC, 5},
    {"_bWGR_emRR", (DL_FUNC) &_bWGR_emRR, 4},
    {"_bWGR_emBL", (DL_FUNC) &_bWGR_emBL, 4},
    {"_bWGR_emDE", (DL_FUNC) &_bWGR_emDE, 3},
    {"_bWGR_emEN", (DL_FUNC) &_bWGR_emEN, 4},
    {"_bWGR_BayesA", (DL_FUNC) &_bWGR_BayesA, 6},
    {"_bWGR_BayesB", (DL_FUNC) &_bWGR_BayesB, 7},
    {"_bWGR_BayesC", (DL_FUNC) &_bWGR_BayesC, 7},
    {"_bWGR_BayesL", (DL_FUNC) &_bWGR_BayesL, 6},
    {"_bWGR_BayesRR", (DL_FUNC) &_bWGR_BayesRR, 6},
    {"_bWGR_BayesA2", (DL_FUNC) &_bWGR_BayesA2, 7},
    {"_bWGR_BayesB2", (DL_FUNC) &_bWGR_BayesB2, 8},
    {"_bWGR_BayesRR2", (DL_FUNC) &_bWGR_BayesRR2, 7},
    {"_bWGR_CNT", (DL_FUNC) &_bWGR_CNT, 1},
    {"_bWGR_IMP", (DL_FUNC) &_bWGR_IMP, 1},
    {"_bWGR_GAU", (DL_FUNC) &_bWGR_GAU, 1},
    {"_bWGR_SPC", (DL_FUNC) &_bWGR_SPC, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_bWGR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
