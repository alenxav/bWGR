#include <RcppEigen.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
SEXP GS2EIGEN(Eigen::Map<Eigen::VectorXd> e,
              Eigen::MappedSparseMatrix<double> X,
              Eigen::Map<Eigen::VectorXd> b,
              Eigen::MappedSparseMatrix<double> XX,
              double Lmb){
  int P = X.cols();
  int N = X.rows();
  Eigen::VectorXd Y(N);
  Eigen::VectorXd r(P);
  Y = X * b + e;
  r = X.transpose() * Y;
  double b0;
  Eigen::VectorXd Xi;
  for(int i=0; i<P; i++){
    b0 = b(i);
    Xi = XX.col(i);
    b(i) = ( r(i) - Xi.transpose()*b + Xi(i)*b0  ) / (Xi(i)+Lmb);
  }
  e = Y - X * b;
  return List::create(Named("b")=b,Named("e")=e);
}

// [[Rcpp::export]]
NumericMatrix NNSEARCH(NumericVector blk, NumericVector row, NumericVector col, int rN, int cN){
  int n = blk.size(); NumericMatrix X(n,(rN*2+1)*(cN*2+1)); NumericVector Obs(n);
  for(int i=0; i<n; i++){ for(int j=0; j<n; j++){
    if( (i>j) && (blk[i]==blk[j]) && (abs(row[i]-row[j])<=rN) && (abs(col[i]-col[j])<=cN) ){
      X(i,Obs[i]) = j+1;  X(j,Obs[j]) = i+1; Obs[i] = Obs[i]+1; Obs[j] = Obs[j]+1; }}}
  return X;}

// [[Rcpp::export]]
SEXP GSFLM(NumericVector y,
           NumericVector e,
           NumericMatrix gen,
           NumericVector b, 
           NumericVector Lmb,
           NumericVector xx,
           double cxx,
           int maxit = 50){
  double tol = 10e-8;
  double phi = cxx;
  NumericVector e0 = e+0.0;
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Initial parameters
  double vy = var(y);
  double vna = sum(y*e)/(n-1);
  // Beta, mu and epsilon
  double b0,eM;
  double mu = mean(e);
  e = e-mu;
  // Regulation coefficients
  NumericVector Vb(p);
  double b1;
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1.0;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0.0;
    for(int j=0; j<p; j++){
      // Gauss-Seidel
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(Lmb(j)+xx(j)+0.01);
      b[j] = b1;
      // Residuals update
      e = e-gen(_,j)*(b1-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Variance components
    vna = sum(e*e0)/n;
    Vb = b*b+(vna/(xx+Lmb));
    Lmb = sqrt(phi*vna/Vb);
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}
  }
  // Output
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("h2")=1-vna/vy,
                      Named("e")=e,
                      Named("Lmb")=Lmb,
                      Named("vb")=Vb);
}

// [[Rcpp::export]]
SEXP GSRR(NumericVector y,
          NumericVector e,
          NumericMatrix gen,
          NumericVector b,
          NumericVector Lmb,
          NumericVector xx,
          double cxx,
          int maxit = 50){
  double tol = 10e-8;
  double phi = cxx;
  NumericVector e0 = e+0.0;
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Initial parameters
  double vy = var(y);
  double vna = sum(y*e)/(n-1);
  // Beta, mu and epsilon
  double b0,eM,vg,LmbTmp;
  double mu = mean(e);
  e = e-mu;
  // Regulation coefficients
  NumericVector Vb(p);
  double b1;
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0.0;
    for(int j=0; j<p; j++){
      // Gauss-Seidel
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(Lmb(j)+xx(j)+0.01);
      b[j] = b1;
      // Residuals update
      e = e-gen(_,j)*(b1-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Variance components
    vna = sum(e*e0)/n;
    vg = (vy-vna)/phi;
    LmbTmp = vna/vg;
    for(int j=0; j<p; j++){
      Vb[j] = vg+0.0;
      Lmb[j] = LmbTmp+0.0;
    }
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}
  }
  // Output
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("h2")=1-vna/vy,
                      Named("e")=e,
                      Named("Lmb")=Lmb,
                      Named("vb")=Vb);
}
