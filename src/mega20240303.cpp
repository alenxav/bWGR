// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <random>

Eigen::VectorXd solver1x(Eigen::VectorXd Y, Eigen::MatrixXd X,
                         int maxit = 100, double tol = 10e-7, double df0 = 20.0){
  int n = X.rows(), p = X.cols(), numit = 0, J;
  double mu = Y.mean(), mu0;
  Eigen::VectorXd y = Y.array()-mu;
  Eigen::VectorXd tilde = X.transpose() * y;
  for(int i=0; i<p; i++){ X.col(i) = X.col(i).array() - X.col(i).mean(); }
  Eigen::VectorXd XX = X.colwise().squaredNorm().array();
  double TrXSX = XX.sum();
  double MSx = TrXSX/(n-1), vy = y.transpose()*Y; vy = vy/(n-1);
  double ve = vy*0.5, vb=(vy*0.5)/(MSx);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(p), beta0(p);
  Eigen::VectorXd e = y*1.0;
  double b0, b1, lambda=ve/vb, vb0=vb*df0, ve0=ve*df0, cnv = 10.0, logtol = log10(tol);
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  while(numit<maxit){
    beta0 = b*1.0;
    std::shuffle(RGSvec.begin(),RGSvec.end(), g);
    for(int j=0; j<p; j++){
      J = RGSvec[j]; b0 = b[J]*1.0;
      b1 = (e.transpose()*X.col(J)+XX(J)*b0)/(XX[J]+lambda);
      e = e - X.col(J)*(b1-b0); b[J] = b1*1.0;}
    mu0 = e.array().mean(); mu+=mu0; e=e.array()-mu0;
    ve = e.transpose()*y;
    ve += e.transpose()*e; 
    ve = (ve+ve0)/(2*n-1+df0);
    vb = b.transpose()*b;
    vb += tilde.transpose()*b;
    vb = (vb+vb0)/(TrXSX+p+df0);  lambda = ve/vb;
    cnv = log10((beta0.array()-b.array()).square().sum());
    ++numit; if( cnv<logtol || numit == maxit || std::isnan(cnv) ) break;}
  return b;
}

Eigen::VectorXd solver2x(Eigen::VectorXd Y, Eigen::MatrixXd X1, Eigen::MatrixXd X2,
                         int maxit = 100, double tol = 10e-7, double df0 = 20.0){
  int n = X1.rows(), p1 = X1.cols(), p2 = X2.cols(), numit = 0, J;
  double mu = Y.mean(), mu0;
  Eigen::VectorXd y = Y.array()-mu;
  Eigen::VectorXd tilde1 = X1.transpose() * y, tilde2 = X2.transpose() * y;
  for(int i=0; i<p1; i++){ X1.col(i) = X1.col(i).array()-X1.col(i).mean();}
  for(int i=0; i<p2; i++){ X2.col(i) = X2.col(i).array()-X2.col(i).mean();}
  Eigen::VectorXd XX1 = X1.colwise().squaredNorm().array(), XX2 = X2.colwise().squaredNorm().array();
  double TrXSX1 = XX1.sum(), TrXSX2 = XX2.sum();
  double MSx1 = TrXSX1/(n-1), MSx2 = TrXSX2/(n-1),  vy=y.transpose()*Y; vy = vy/(n-1);
  double ve = vy*0.5, vb1=(vy*0.5)/(MSx1), vb2=(vy*0.5)/(MSx2), h2=0.5;
  Eigen::VectorXd b_1 = Eigen::VectorXd::Zero(p1), beta01(p1);
  Eigen::VectorXd b_2 = Eigen::VectorXd::Zero(p2), beta02(p2);
  Eigen::VectorXd e = y*1.0;
  double b0, b1, lambda1=ve/vb1, lambda2=ve/vb2, cnv=10.0, logtol=log10(tol);
  double vb01 = vb1*df0, vb02 = vb2*df0, ve0 = ve*df0;
  std::vector<int> RGSvec1(p1), RGSvec2(p2);
  for(int j=0; j<p1; j++){RGSvec1[j]=j;}
  for(int j=0; j<p2; j++){RGSvec2[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  while(numit<maxit){
    beta01 = b_1*1.0; beta02 = b_2*1.0;
    std::shuffle(RGSvec1.begin(),RGSvec1.end(), g);
    std::shuffle(RGSvec2.begin(),RGSvec2.end(), g);
    for(int j=0; j<p1; j++){
      J = RGSvec1[j]; b0 = b_1[J]*1.0;
      b1 = (e.transpose()*X1.col(J)+XX1(J)*b0)/(XX1[J]+lambda1);
      e = e - X1.col(J)*(b1-b0); b_1[J] = b1*1.0;}
    for(int j=0; j<p2; j++){
      J = RGSvec2[j]; b0 = b_2[J]*1.0;
      b1 = (e.transpose()*X2.col(J)+XX2(J)*b0)/(XX2[J]+lambda2);
      e = e - X2.col(J)*(b1-b0); b_2[J] = b1*1.0;}
    mu0=e.array().mean(); mu+=mu0; e=e.array()-mu0;
    ve = e.transpose()*e;
    ve += e.transpose()*y; 
    ve = (ve+ve0)/(2*n-1+df0);
    vb1 = tilde1.transpose()*b_1; vb1+=b_1.transpose()*b_1; vb1+=vb01;  
    vb2 = tilde2.transpose()*b_2; vb2+=b_2.transpose()*b_2; vb2+=vb02;
    vb1 = vb1/(TrXSX1+p1+df0); vb2 = vb2/(TrXSX2+p2+df0);
    lambda1 = ve/vb1; lambda2 = ve/vb2;
    cnv = log10((beta01.array()-b_1.array()).square().sum()+(beta02.array()-b_2.array()).square().sum());
    ++numit; if( cnv<logtol || numit == maxit || std::isnan(cnv) ) break;  }
  Eigen::VectorXd xxx(1+p1+p2);
  xxx(0) = mu;
  for(int j=0; j<p1 ; j++){xxx(1+j)=b_1(j);}
  for(int j=0; j<p2 ; j++){xxx(1+p1+j)=b_2(j);}
  return xxx;
}

Eigen::MatrixXd submat_f(Eigen::MatrixXd X, Eigen::VectorXi w){
  int n=w.sum(), N=X.rows(), p=X.cols(), n0=0; Eigen::MatrixXd XX(n,p);
  for(int i=0; i<N; i++){ if(w[i]==1){ XX.row(n0) = X.row(i).array(); n0+=1;}}
  return XX;}

Eigen::VectorXd subvec_f(Eigen::VectorXd X, Eigen::VectorXi w){
  int n=w.sum(), N=X.size(), n0=0; Eigen::VectorXd XX(n);
  for(int i=0; i<N; i++){ if(w[i]==1){ XX[n0] = X[i]; n0+=1;}}
  return XX;}

Eigen::MatrixXd UVBETA(Eigen::MatrixXd Y, Eigen::MatrixXd X){
  int n0=Y.rows(), p=X.cols(), k=Y.cols(); Eigen::MatrixXd BETA(p,k); Eigen::MatrixXi W(n0,k);
  for(int i=0;i<n0;i++){for(int j=0;j<k;j++){if(std::isnan(Y(i,j))){W(i,j)=0;}else{W(i,j)=1;}}}
  for(int i=0;i<k;i++){
    if(W.col(i).array().sum()>0){
      BETA.col(i) = solver1x(
        subvec_f( Y.col(i).array(), W.col(i).array()),
        submat_f( X, W.col(i).array())).array();}else{
          BETA.col(i) = Eigen::VectorXd::Zero(p);}}
  return BETA;}

Eigen::MatrixXd GetImputedY(Eigen::MatrixXd Y, Eigen::MatrixXd X, Eigen::MatrixXd BETA){
  int n0=Y.rows(),k=Y.cols();
  Eigen::VectorXd Mu = Eigen::VectorXd::Zero(k), N = Eigen::VectorXd::Zero(k);
  for(int j=0;j<k;j++){for(int i=0;i<n0;i++){
    if(!std::isnan(Y(i,j))){N(j)+=1.0;Mu(j)+=Y(i,j);}}}
  Mu = Mu.array() / N.array();
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(!std::isnan(Y(i,j))){
        Y(i,j) -= Mu(j);}else{
          Y(i,j) = X.row(i)*BETA.col(j);}}}
  return Y;}

Eigen::MatrixXd LatentSpaces(Eigen::MatrixXd Y, Eigen::MatrixXd X, Eigen::MatrixXd BETA){
  int n=Y.rows(),k=Y.cols();
  Eigen::MatrixXd Y2 = GetImputedY(Y,X,BETA);
  Eigen::VectorXd SD = Y2.colwise().squaredNorm().array(); SD = (SD.array()/(n-1)).sqrt();
  for(int i=0; i<k; i++){ Y2.col(i) /= SD(i);};
  Eigen::BDCSVD<Eigen::MatrixXd> svd(Y2, Eigen::ComputeThinU | Eigen::ComputeThinV );
  return svd.matrixU() * svd.singularValues().matrix().asDiagonal();
}

// [[Rcpp::export]]
SEXP MEGA(Eigen::MatrixXd Y, Eigen::MatrixXd X){
  int n0=Y.rows(), p1=X.cols(), k=Y.cols(); Eigen::MatrixXi W(n0,k);
  for(int i=0;i<n0;i++){for(int j=0;j<k;j++){if(std::isnan(Y(i,j))){W(i,j)=0;}else{W(i,j)=1;}}}
  Eigen::MatrixXd BETA = UVBETA(Y,X);
  Eigen::MatrixXd LS = LatentSpaces(Y,X,BETA);
  Eigen::MatrixXd LS_BETA = UVBETA(LS,X);
  int p2 = LS.cols();
  Eigen::VectorXd xxx(1+p1+p2);
  // store outputs
  Eigen::VectorXd mu(k), h2(k);
  Eigen::MatrixXd b1(p2,k), b2(p1,k);
  for(int i=0; i<k; i++){
    xxx = solver2x(
      subvec_f( Y.col(i).array(), W.col(i).array()),
      submat_f( LS, W.col(i).array()),
      submat_f( X, W.col(i).array())).array();
    mu(i) = xxx(0);
    for(int j=0; j<p2 ; j++){b1(j,i) = xxx(1+j);}
    for(int j=0; j<p1 ; j++){b2(j,i) = xxx(1+p2+j);}
  }
  // Fitted values
  Eigen::MatrixXd end_beta = LS_BETA * b1 + b2;
  Eigen::MatrixXd hat = LS*b1+X*b2;
  Eigen::MatrixXd gebv = X*end_beta;
  for(int i=0; i<k; i++){
    hat.col(i) = hat.col(i).array() + mu(i);
    gebv.col(i) = gebv.col(i).array() + mu(i);
  }
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=end_beta,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("LS")=LS,
                            Rcpp::Named("LS_BETA")=LS_BETA,
                            Rcpp::Named("BETA1")=b1,
                            Rcpp::Named("BETA2")=b2,
                            Rcpp::Named("gebv")=gebv);
}

// [[Rcpp::export]]
SEXP GSEM(Eigen::MatrixXd Y, Eigen::MatrixXd X){
  int n0=Y.rows(), p1=X.cols(), k=Y.cols(); Eigen::MatrixXi W(n0,k);
  for(int i=0;i<n0;i++){for(int j=0;j<k;j++){if(std::isnan(Y(i,j))){W(i,j)=0;}else{W(i,j)=1;}}}
  Eigen::MatrixXd BETA = UVBETA(Y,X);
  Eigen::BDCSVD<Eigen::MatrixXd> svd(X*BETA, Eigen::ComputeThinU | Eigen::ComputeThinV );
  Eigen::MatrixXd LS = svd.matrixU() * svd.singularValues().matrix().asDiagonal();
  int p2 = LS.cols();
  Eigen::VectorXd xxx(1+p1+p2);
  // store outputs
  Eigen::VectorXd mu(k), h2(k);
  Eigen::MatrixXd b1(p2,k), b2(p1,k);
  for(int i=0; i<k; i++){
    xxx = solver2x(
      subvec_f( Y.col(i).array(), W.col(i).array()),
      submat_f( LS, W.col(i).array()),
      submat_f( X, W.col(i).array())).array();
    mu(i) = xxx(0);
    for(int j=0; j<p2 ; j++){b1(j,i) = xxx(1+j);}
    for(int j=0; j<p1 ; j++){b2(j,i) = xxx(1+p2+j);}
  }
  // Fitted values
  Eigen::MatrixXd hat = LS*b1+X*b2;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=BETA*svd.matrixV()*b1+b2,
                            Rcpp::Named("hat")=hat);
}
