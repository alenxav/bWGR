// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <random>

// [[Rcpp::export]]
Eigen::MatrixXd EigenARC(Eigen::MatrixXd X, bool centralizeX = true, int cores = 1){
  // cseweb.ucsd.edu/~saul/papers/nips09_kernel.pdf
  if(cores!=1) Eigen::setNbThreads(cores);
  int p = X.cols(), n = X.rows(); 
  double tmp, Npi=3.1416, theta, J1, Kij, Norm;
  if(centralizeX){
    for(int i=0; i<p; i++){
      tmp = (X.col(i).array()).mean();
      X.col(i) = X.col(i).array()-tmp;}}
  Eigen::MatrixXd XXp = X*X.transpose();
  tmp = 1/(XXp.diagonal().mean()); XXp *= tmp;
  Eigen::VectorXd DiagXXp = XXp.diagonal().array();
  for(int i=0; i<n; i++){ for(int j=i; j<n; j++){ 
    Norm = sqrt(DiagXXp(i)*DiagXXp(j)*1.001);
    theta = acos( XXp(i,j)/Norm);
    J1 = sin(theta) + (Npi-theta)*cos(theta);
    Kij = Norm/Npi*J1;
    XXp(i,j) = Kij*1.0; XXp(j,i) = Kij*1.0;}}
  return XXp;}

// [[Rcpp::export]]
Eigen::MatrixXd EigenGAU(Eigen::MatrixXd X, double phi = 1.0, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores);
  int n = X.rows(); double tmp;
  Eigen::MatrixXd XXp = X*X.transpose();
  for(int i=0; i<n; i++){ for(int j=0; j<n; j++){ if(i>j){
    tmp = sqrt(XXp(i,i) + XXp(j,j) - 2*XXp(i,j));
    XXp(i,j) = tmp*1.0; XXp(j,i) = tmp*1.0;}}};
  for(int i=0; i<n; i++){XXp(i,i) = 0.0;}
  tmp = phi * (-n*(n-1)) / (XXp.colwise().sum()).sum();
  XXp *= tmp; return exp(XXp.array());}

// [[Rcpp::export]]
Eigen::MatrixXd EigenGRM(Eigen::MatrixXd X, bool centralizeZ = true, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores); 
  int p = X.cols(); double tmp;
  if(centralizeZ){
    for(int i=0; i<p; i++){
      tmp = (X.col(i).array()).mean();
      X.col(i) = X.col(i).array()-tmp;}}
  Eigen::MatrixXd XXp = X*X.transpose();
  tmp = 1/(XXp.diagonal().mean());
  XXp *= tmp; return XXp;}

// [[Rcpp::export]]
Eigen::MatrixXd EigenCNT(Eigen::MatrixXd X, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores); 
  int p = X.cols();
  Eigen::VectorXd xx = X.colwise().mean();
  for(int i=0; i<p; i++){ X.col(i) = X.col(i).array() - xx(i);}
  return X;}

// [[Rcpp::export]]
SEXP EigenEVD(Eigen::MatrixXd A, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores); 
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  return Rcpp::List::create(Rcpp::Named("U")=es.eigenvectors(),
                            Rcpp::Named("D")=es.eigenvalues());}

// [[Rcpp::export]]
SEXP EigenBDCSVD(Eigen::MatrixXd X, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores);
  Eigen::BDCSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV );
  return Rcpp::List::create(Rcpp::Named("U")=svd.matrixU(),
                            Rcpp::Named("D")=svd.singularValues(),
                            Rcpp::Named("V")=svd.matrixV());}

// [[Rcpp::export]]
SEXP EigenJacobiSVD(Eigen::MatrixXd X, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV );
  return Rcpp::List::create(Rcpp::Named("U")=svd.matrixU(),
                            Rcpp::Named("D")=svd.singularValues(),
                            Rcpp::Named("V")=svd.matrixV());}

// [[Rcpp::export]]
Eigen::VectorXd EigenAcc(Eigen::MatrixXd X1, Eigen::MatrixXd X2, double h2 = 0.5, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores);
  Eigen::MatrixXd X1X1 = X1*X1.transpose(), X1X2 = X1*X2.transpose();
  double Ve = (1.0-h2)/h2, alpha = 1.0/(X1X1.diagonal().array()).mean();
  Eigen::MatrixXd V = X1X1*alpha; V.diagonal() = V.diagonal().array() + Ve;
  return sqrt( alpha * (X1X2.transpose()*(V.llt().solve(X1X2))).diagonal().array()/X2.rowwise().squaredNorm().array());}


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
  return Rcpp::List::create(Rcpp::Named("b")=b,Rcpp::Named("e")=e);
}

// [[Rcpp::export]]
SEXP mrr(Eigen::MatrixXd Y, Eigen::MatrixXd X){
  
  // Basic info
  int k = Y.cols(), n0 = Y.rows(), p = X.cols();
  int maxit = 200;
  
  // Incidence matrix Z
  Eigen::MatrixXd Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXd n = Z.colwise().sum();
  Eigen::VectorXd iN = n.array().inverse();
  
  // Centralize y
  Eigen::VectorXd mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXd y(n0,k);
  for(int i=0; i<k; i++){
    y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  
  // Sum of squares of X
  Eigen::MatrixXd XX(p,k);
  for(int i=0; i<p; i++){
    XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  
  // Compute Tr(XSX);
  Eigen::MatrixXd XSX(p,k);
  for(int i=0; i<p; i++){
    XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
      ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXd MSx = XSX.colwise().sum();
  Eigen::VectorXd TrXSX = n.array()*MSx.array();
  
  // Variances
  iN = (n.array()-1).inverse();
  Eigen::VectorXd vy = y.colwise().squaredNorm(); vy = vy.array() * iN.array();
  Eigen::VectorXd ve = vy * 0.5;
  Eigen::VectorXd iVe = ve.array().inverse();
  Eigen::MatrixXd vb(k,k), TildeHat(k,k);
  vb = (ve.array()/MSx.array()).matrix().asDiagonal();
  Eigen::MatrixXd iG = vb.inverse();
  Eigen::VectorXd h2 = 1 - ve.array()/vy.array();
  
  // Beta tilde;
  Eigen::MatrixXd tilde = X.transpose() * y;
  
  // Initialize coefficient matrices
  Eigen::MatrixXd LHS(k,k);
  Eigen::VectorXd RHS(k);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(p,k);
  Eigen::VectorXd b0(k), b1(k);
  Eigen::MatrixXd e(n0,k); e = y*1.0;
  
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Convergence control
  Eigen::MatrixXd A = vb*1.0; double MinDVb, inflate;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> EVDofA(A);
  Eigen::MatrixXd beta0(p,k);
  double cnv = 10.0, logtol = -8.0;
  int numit = 0;
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    
    // Randomized Gauss-Seidel loop
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    for(int j=0; j<p; j++){
      J = RGSvec[j];
      // Update coefficient
      b0 = b.row(J)*1.0;
      LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVe.array()).matrix();
      RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
      RHS = RHS.array() *iVe.array();
      b1 = LHS.llt().solve(RHS);
      b.row(J) = b1;
      // Update residuals
      e = (e-(X.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = ve.array() * iN.array();
    iVe = ve.array().inverse();
    
    // Genetic variance
    TildeHat = b.transpose()*tilde;
    for(int i=0; i<k; i++){for(int j=0; j<k; j++){
      if(i==j){ vb(i,i) = TildeHat(i,i)/TrXSX(i); }else{
        vb(i,j) = (TildeHat(i,j)+TildeHat(j,i))/(TrXSX(i)+TrXSX(j));}}}
    
    // Bending
    A = vb*1.0;
    EVDofA.compute(A); MinDVb = EVDofA.eigenvalues().minCoeff();
    if( MinDVb < 0.0 ){ inflate = abs(MinDVb*1.1);
      A.diagonal().array()+=inflate; vb=A*1.0;}
    iG = vb.completeOrthogonalDecomposition().pseudoInverse();
    
    // Print status
    cnv = log10((beta0.array()-b.array()).square().sum());  ++numit;
    if( numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; }
    if( cnv<logtol ){break;}
    if(std::isnan(cnv)){ break;}
    
  }
  
  // Fitting the model
  h2 = 1 - ve.array()/vy.array();
  Eigen::MatrixXd hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  
  // Genetic correlations
  Eigen::MatrixXd GC(k,k);
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC,
                            Rcpp::Named("Vb")=vb,
                            Rcpp::Named("Ve")=ve,
                            Rcpp::Named("MSx")=MSx,
                            Rcpp::Named("cnv")=cnv);
  
}

// [[Rcpp::export]]
SEXP mrr_float(Eigen::MatrixXf Y, Eigen::MatrixXf X){
  
  // Basic info
  int k = Y.cols(), n0 = Y.rows(), p = X.cols();
  int maxit = 200;
  
  // Incidence matrix Z
  Eigen::MatrixXf Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXf n = Z.colwise().sum();
  Eigen::VectorXf iN = n.array().inverse();
  
  // Centralize y
  Eigen::VectorXf mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXf y(n0,k);
  for(int i=0; i<k; i++){
    y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  
  // Sum of squares of X
  Eigen::MatrixXf XX(p,k);
  for(int i=0; i<p; i++){
    XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  
  // Compute Tr(XSX);
  Eigen::MatrixXf XSX(p,k);
  for(int i=0; i<p; i++){
    XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
      ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXf MSx = XSX.colwise().sum();
  Eigen::VectorXf TrXSX = n.array()*MSx.array();
  
  // Variances
  iN = (n.array()-1).inverse();
  Eigen::VectorXf vy = y.colwise().squaredNorm(); vy = vy.array() * iN.array();
  Eigen::VectorXf ve = vy * 0.5;
  Eigen::VectorXf iVe = ve.array().inverse();
  Eigen::MatrixXf vb(k,k), TildeHat(k,k);
  vb = (ve.array()/MSx.array()).matrix().asDiagonal();
  Eigen::VectorXf h2 = 1 - ve.array()/vy.array();
  
  // Bending
  Eigen::MatrixXf A = vb*1.0;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> EVDofA(A); float MinDVb, inflate;
  Eigen::MatrixXf iG = vb.completeOrthogonalDecomposition().pseudoInverse();
  
  // Beta tilde;
  Eigen::MatrixXf tilde = X.transpose() * y;
  
  // Initialize coefficient matrices
  Eigen::MatrixXf LHS(k,k);
  Eigen::VectorXf RHS(k);
  Eigen::MatrixXf b = Eigen::MatrixXf::Zero(p,k);
  Eigen::VectorXf b0(k), b1(k);
  Eigen::MatrixXf e(n0,k); e = y*1.0;
  
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Inner RGS
  std::vector<int> InnerRGSvec(k);
  for(int j=0; j<k; j++){InnerRGSvec[j]=j;}
  std::random_device rd2;
  std::mt19937 g2(rd2());
  int ri;
  
  // Convergence control
  Eigen::MatrixXf beta0(p,k);
  float cnv = 10.0, logtol = -8.0;
  int numit = 0;
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    
    // Randomized Gauss-Seidel loop
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    std::shuffle(InnerRGSvec.begin(), InnerRGSvec.end(), g2);
    
    for(int j=0; j<p; j++){
      
      J = RGSvec[j];
      // Update coefficient
      b0 = b.row(J)*1.0;
      LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVe.array()).matrix();
      RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
      RHS = RHS.array() *iVe.array();
      // Inner GS
      b1 = b.row(J)*1.0;
      for(int i=0; i<k; i++){
        ri = InnerRGSvec[i];
        b1(ri) = (RHS(ri)-(LHS.col(ri).array()*b1.array()).sum()+LHS(ri,ri)*b1(ri))/LHS(ri,ri);}
      // Update residuals
      b.row(J) = b1;
      e = (e-(X.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = ve.array() * iN.array();
    iVe = ve.array().inverse();
    
    // Genetic variance
    TildeHat = b.transpose()*tilde;
    for(int i=0; i<k; i++){for(int j=0; j<k; j++){
      if(i==j){ vb(i,i) = TildeHat(i,i)/TrXSX(i); }else{
        vb(i,j) = (TildeHat(i,j)+TildeHat(j,i))/(TrXSX(i)+TrXSX(j));}}}
    
    // Bending and inverse of vb
    A = vb*1.0;
    EVDofA.compute(A); MinDVb = EVDofA.eigenvalues().minCoeff();
    if( MinDVb < 0.0 ){ inflate = abs(MinDVb*1.1);
      A.diagonal().array()+=inflate; vb=A*1.0;}
    iG = vb.completeOrthogonalDecomposition().pseudoInverse();
    
    // Print status
    cnv = log10((beta0.array()-b.array()).square().colwise().sum().maxCoeff());
    ++numit; if( cnv<logtol ){break;}
    if(std::isnan(cnv)){ break;}
    
  }
  
  // Fitting the model
  h2 = 1 - ve.array()/vy.array();
  Eigen::MatrixXf hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  
  // Genetic correlations
  Eigen::MatrixXf GC(k,k);
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC,
                            Rcpp::Named("Vb")=vb,
                            Rcpp::Named("Ve")=ve,
                            Rcpp::Named("MSx")=MSx,
                            Rcpp::Named("cnv")=cnv);
  
}


// [[Rcpp::export]]
SEXP mrr2X(Eigen::MatrixXd Y, Eigen::MatrixXd X1, Eigen::MatrixXd X2){
  
  // Basic info
  int maxit = 1000;
  int k = Y.cols(), n0 = Y.rows(), p1 = X1.cols(), p2 = X2.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXd Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXd n = Z.colwise().sum();
  Eigen::VectorXd iN = n.array().inverse();
  
  // Centralize y
  Eigen::VectorXd mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXd y(n0,k);
  for(int i=0; i<k; i++){
    y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  
  // Sum of squares of X1
  Eigen::MatrixXd XX1(p1,k);
  for(int i=0; i<p1; i++){
    XX1.row(i) = X1.col(i).array().square().matrix().transpose() * Z;}
  
  // Sum of squares of X2
  Eigen::MatrixXd XX2(p2,k);
  for(int i=0; i<p2; i++){
    XX2.row(i) = X2.col(i).array().square().matrix().transpose() * Z;}
  
  // Compute Tr(XSX)1;
  Eigen::MatrixXd XSX1(p1,k);
  for(int i=0; i<p1; i++){
    XSX1.row(i) = XX1.row(i).transpose().array()*iN.array() - 
      ((X1.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXd MSx1 = XSX1.colwise().sum();
  Eigen::VectorXd TrXSX1 = n.array()*MSx1.array();
  
  // Compute Tr(XSX)2;
  Eigen::MatrixXd XSX2(p2,k);
  for(int i=0; i<p2; i++){
    XSX2.row(i) = XX2.row(i).transpose().array()*iN.array() - 
      ((X2.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXd MSx2 = XSX2.colwise().sum();
  Eigen::VectorXd TrXSX2 = n.array()*MSx2.array();
  
  // Variances
  iN = (n.array()-1).inverse();
  // VE
  Eigen::VectorXd vy = y.colwise().squaredNorm(); vy = vy.array()*iN.array();
  Eigen::VectorXd ve = vy * 0.5;
  Eigen::VectorXd iVe = ve.array().inverse();
  // VB1
  Eigen::MatrixXd vb1(k,k), TildeHat1(k,k);
  vb1 = (ve.array()/MSx1.array()).matrix().asDiagonal();
  Eigen::MatrixXd iG1 = vb1.inverse();
  // VB2
  Eigen::MatrixXd vb2(k,k), TildeHat2(k,k);
  vb2 = (ve.array()/MSx2.array()).matrix().asDiagonal();
  Eigen::MatrixXd iG2 = vb2.inverse();
  
  // Beta tilde;
  Eigen::MatrixXd tilde1 = X1.transpose() * y;
  Eigen::MatrixXd tilde2 = X2.transpose() * y;
  
  // Initialize coefficient matrices
  Eigen::VectorXd b0(k), b1(k);
  Eigen::MatrixXd e(n0,k); e = y*1.0;
  Eigen::MatrixXd LHS(k,k);
  Eigen::VectorXd RHS(k);
  Eigen::MatrixXd bA = Eigen::MatrixXd::Zero(p1,k);
  Eigen::MatrixXd bB = Eigen::MatrixXd::Zero(p2,k);
  
  // RGS
  std::vector<int> RGSvec1(p1);
  std::vector<int> RGSvec2(p2);
  for(int j=0; j<p1; j++){RGSvec1[j]=j;}
  for(int j=0; j<p2; j++){RGSvec2[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Convergence control
  Eigen::MatrixXd beta01(p1,k), beta02(p2,k), A(k,k);
  double cnv = 10.0, logtol = -10.0, MinDVb, inflate; int numit = 0;
  A = vb1*1.0; Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> EVDofA(A);
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta01 = bA*1.0;
    beta02 = bB*1.0;
    
    // Randomized Gauss-Seidel loop 1
    std::shuffle(RGSvec1.begin(), RGSvec1.end(), g);
    for(int j=0; j<p1; j++){
      J = RGSvec1[j];
      // Update coefficient
      b0 = bA.row(J)*1.0;
      LHS = iG1;
      LHS.diagonal() += (XX1.row(J).transpose().array() * iVe.array()).matrix();
      RHS = (X1.col(J).transpose()*e).array() + XX1.row(J).array()*b0.transpose().array();
      RHS = RHS.array() *iVe.array();
      b1 = LHS.llt().solve(RHS);
      bA.row(J) = b1;
      // Update residuals
      e = (e-(X1.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Randomized Gauss-Seidel loop 2
    std::shuffle(RGSvec2.begin(), RGSvec2.end(), g);
    for(int j=0; j<p2; j++){
      J = RGSvec2[j];
      // Update coefficient
      b0 = bB.row(J)*1.0;
      LHS = iG2;
      LHS.diagonal() += (XX2.row(2).transpose().array() * iVe.array()).matrix();
      RHS = (X2.col(J).transpose()*e).array() + XX2.row(J).array()*b0.transpose().array();
      RHS = RHS.array() *iVe.array();
      b1 = LHS.llt().solve(RHS);
      bB.row(J) = b1;
      // Update residuals
      e = (e-(X2.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = ve.array() * iN.array();
    iVe = ve.array().inverse();
    
    // Genetic variance 1
    TildeHat1 = bA.transpose()*tilde1;
    for(int i=0; i<k; i++){for(int j=0; j<k; j++){
      if(i==j){ vb1(i,i) = TildeHat1(i,i)/TrXSX1(i); }else{
        vb1(i,j) = (TildeHat1(i,j)+TildeHat1(j,i))/(TrXSX1(i)+TrXSX1(j));}}}
    // Bending 1
    A = vb1*1.0;
    EVDofA.compute(A); MinDVb = EVDofA.eigenvalues().minCoeff();
    if( MinDVb < 0.0 ){ inflate = abs(MinDVb*1.1);
      A.diagonal().array()+=inflate; vb1=A*1.0;}
    iG1 = vb1.completeOrthogonalDecomposition().pseudoInverse();
    
    // Genetic variance 1
    TildeHat2 = bB.transpose()*tilde2;
    for(int i=0; i<k; i++){for(int j=0; j<k; j++){
      if(i==j){ vb2(i,i) = TildeHat2(i,i)/TrXSX2(i); }else{
        vb2(i,j) = (TildeHat2(i,j)+TildeHat2(j,i))/(TrXSX2(i)+TrXSX2(j));}}}
    // Bending 2
    A = vb2*1.0;
    EVDofA.compute(A); MinDVb = EVDofA.eigenvalues().minCoeff();
    if( MinDVb < 0.0 ){ inflate = abs(MinDVb*1.1);
      A.diagonal().array()+=inflate; vb2=A*1.0;}
    iG2 = vb2.completeOrthogonalDecomposition().pseudoInverse();
    
    // Print status
    cnv = log10((beta01.array()-bA.array()).square().sum()) + log10((beta02.array()-bB.array()).square().sum());
    ++numit;
    if( numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; } 
    if( cnv<logtol ){break;}
    if(std::isnan(cnv)){ break;}
    
  }
  
  // Fitting the model
  Eigen::MatrixXd hat = X1 * bA + X2 * bB;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  
  // Genetic correlations and genetic variance 1
  Eigen::MatrixXd GC1(k,k), va1(k,k);
  va1.diagonal() = vb1.diagonal().array() * MSx1.array();
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){ GC1(i,j)=vb1(i,j)/(sqrt(vb1(i,i)*vb1(j,j)));}}
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){ if(i!=j){ va1(i,j)= GC1(i,j)*sqrt(va1(i,i)*va1(j,j));}}}
  
  // Genetic correlations and genetic variance 2
  Eigen::MatrixXd GC2(k,k), va2(k,k);
  va2.diagonal() = vb2.diagonal().array() * MSx2.array();
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){ GC2(i,j)=vb2(i,j)/(sqrt(vb2(i,i)*vb2(j,j)));}}
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){ if(i!=j){ va2(i,j)= GC2(i,j)*sqrt(va2(i,i)*va2(j,j));}}}
  
  // Heritability
  Eigen::VectorXd h2 = 1 - ve.array()/vy.array();
  Eigen::VectorXd h2A = va1.diagonal().array() / ( va1.diagonal().array() + ve.array()).array();
  Eigen::VectorXd h2B = va2.diagonal().array() / ( va2.diagonal().array() + ve.array()).array();
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b1")=bA, Rcpp::Named("b2")=bB,
                            Rcpp::Named("hat")=hat, Rcpp::Named("h2")=h2,
                            Rcpp::Named("h2_1")=h2A, Rcpp::Named("h2_2")=h2B,
                            Rcpp::Named("GC1")=GC1, Rcpp::Named("GC2")=GC2,
                            Rcpp::Named("VE")=ve,
                            Rcpp::Named("VB1")=vb1,Rcpp::Named("VB2")=vb2,
                            Rcpp::Named("VA1")=va1, Rcpp::Named("VA2")=va2,
                            Rcpp::Named("cnv")=cnv);}

// [[Rcpp::export]]
SEXP mrr_svd(Eigen::MatrixXd Y, Eigen::MatrixXd W){
  
  // Start setup
  int maxit = 500;
  double tol = 10e-9;
  
  // Gather basic info
  int k = Y.cols(), n0 = Y.rows(), m = W.cols();
  
  // Center X
  Rcpp::Rcout << "Centering marker score matrix\n";
  Eigen::VectorXd xx = W.colwise().mean();
  for(int i=0; i<m; i++){ W.col(i) = W.col(i).array() - xx(i);}
  
  // Single value decomposition
  Rcpp::Rcout << "SVD of marker scores\n";
  Eigen::BDCSVD<Eigen::MatrixXd> svd(W, Eigen::ComputeThinU | Eigen::ComputeThinV );
  Eigen::MatrixXd V = svd.matrixV();
  Eigen::MatrixXd X = svd.matrixU() * svd.singularValues().array().matrix().asDiagonal();
  int p = X.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXd Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXd n = Z.colwise().sum();
  Eigen::VectorXd iN = n.array().inverse();
  
  // Centralize y
  Rcpp::Rcout << "Centering Y\n";
  Eigen::VectorXd mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXd y(n0,k);
  for(int i=0; i<k; i++){
    y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();
  }
  
  // Sum of squares of X
  Rcpp::Rcout << "Computing diagonal elements of Z'Z\n";
  Eigen::MatrixXd XX(p,k);
  for(int i=0; i<p; i++){ XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  // Compute Tr(XSX);
  Eigen::MatrixXd XSX(p,k);
  for(int i=0; i<p; i++){ XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
      ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXd MSx = XSX.colwise().sum();
  Eigen::VectorXd TrXSX = n.array()*MSx.array();
  
  Rcpp::Rcout << "Set starting values for coefficients and variances\n";
  // Variances
  iN = (n.array()-1).inverse();
  Eigen::VectorXd vy = y.colwise().squaredNorm(); vy = vy.array() * iN.array();
  Eigen::VectorXd ve = vy * 0.5;
  Eigen::VectorXd iVe = ve.array().inverse();
  Eigen::MatrixXd vb(k,k), TildeHat(k,k);
  vb = (ve.array()/MSx.array()).matrix().asDiagonal();
  Eigen::MatrixXd iG = vb.inverse();
  Eigen::VectorXd h2 = 1 - ve.array()/vy.array();
  // Beta tilde;
  Eigen::MatrixXd tilde = X.transpose() * y;
  Eigen::VectorXd TrDinvXSX(k);
  Eigen::MatrixXd Dinv(p,k);
  for(int i=0; i<k; i++){ XSX.col(i) = XSX.col(i).array() * n(i); }
  // Initialize coefficient matrices
  Eigen::MatrixXd LHS(k,k);
  Eigen::VectorXd RHS(k);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(p,k);
  Eigen::VectorXd b0(k), b1(k);
  Eigen::MatrixXd e(n0,k); e = y*1.0;
  // Bending and convergence control
  Eigen::MatrixXd A = vb*1.0, GC(k,k);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> EVDofA(A);
  Eigen::MatrixXd beta0(p,k), vb0(k,k);
  Eigen::VectorXd CNV1(maxit),CNV2(maxit),CNV3(maxit), ve0(k), h20(k);
  double cnv = 10.0, MinDVb, inflate;
  int numit = 0;
  double logtol = log10(tol);
  
  // Loop
  Rcpp::Rcout << "Starting Gauss-Seidel\n";
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    vb0 = vb*1.0;
    ve0 = ve*1.0;
    h20 = h2*1.0;
    
    // Randomized Gauss-Seidel loop
    for(int J=0; J<p; J++){
      // Update coefficient
      b0 = b.row(J)*1.0;
      LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVe.array()).matrix();
      RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
      RHS = RHS.array() * iVe.array();
      b1 = LHS.llt().solve(RHS);
      b.row(J) = b1;
      // Update residuals
      e = (e-(X.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = ve.array() * iN.array();
    iVe = ve.array().inverse();
    h2 = 1 - ve.array()/vy.array();
    
    // Get tilde-hat
    for(int i=0; i<k; i++){
        Dinv.col(i) = (XSX.col(i).array()/ve(i) + iG(i,i)).inverse().array();
        TrDinvXSX(i)  = (XSX.col(i).transpose() * Dinv.col(i));}
    TildeHat = b.transpose()* Dinv.cwiseProduct(tilde);
    
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i==j){ // Variances
          vb(i,i) = TildeHat(i,i)/TrDinvXSX(i);
        }else{ // Covariances
          vb(i,j) = (TildeHat(i,j)+TildeHat(j,i))/(TrDinvXSX(i)+TrDinvXSX(j));
        }}}
    
    // Bending
    A = vb*1.0;
    EVDofA.compute(A); MinDVb = EVDofA.eigenvalues().minCoeff();
    if( MinDVb < 0.0 ){ inflate = abs(MinDVb*1.1);
      A.diagonal().array()+=inflate; vb=A*1.0;}
    iG = vb.completeOrthogonalDecomposition().pseudoInverse();
    
    // Covariances
    cnv = log10((beta0.array()-b.array()).square().sum());  CNV1(numit) = cnv;
    CNV2(numit) = log10((h20.array()-h2.array()).square().sum());
    CNV3(numit) = log10((vb0.array()-vb.array()).square().sum());
    
    // Print
    ++numit;
    if( numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; } 
    if( cnv<logtol ){ Rcpp::Rcout << "Model coverged in "<< numit << " iterations\n"; break; }
    if(std::isnan(cnv)){ break;}
    
  }
  
  Rcpp::Rcout << "Fitting final model\n";
  // Fitting the model
  Eigen::MatrixXd hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  Eigen::MatrixXd beta = V*b;
  
  // Correlations
  Rcpp::Rcout << "Estimating correlations\n";
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  
  // Resize convergence vectors
  Rcpp::Rcout << "Convergence statistics\n";
  Eigen::VectorXd CNV1b(numit),CNV2b(numit),CNV3b(numit);
  for(int i=0; i<numit; i++){ CNV1b(i)=CNV1(i);CNV2b(i)=CNV2(i);CNV3b(i)=CNV3(i);}
  
  // Null model Output
  Rcpp::List NullModelOutput = Rcpp::List::create(Rcpp::Named("Intercepts")=mu,
                                                  Rcpp::Named("MarkerEffects")=beta,
                                                  Rcpp::Named("FittedValues")=hat,
                                                  Rcpp::Named("Heritability")=h2,
                                                  Rcpp::Named("WCorrelations")=GC,
                                                  Rcpp::Named("VarBeta")=vb,
                                                  Rcpp::Named("VarResiduals")=ve,
                                                  Rcpp::Named("ConvergenceBeta")=CNV1b,
                                                  Rcpp::Named("ConvergenceH2")=CNV2b,
                                                  Rcpp::Named("ConvergenceVar")=CNV3b,
                                                  Rcpp::Named("NumOfIterations")=numit);
  NullModelOutput.attr("class") = "WModel";
  return NullModelOutput;
  
}

// [[Rcpp::export]]
SEXP MRR3(Eigen::MatrixXd Y,
          Eigen::MatrixXd X,
          int maxit = 500,
          double tol = 10e-9,
          int cores = 1,
          bool TH = false,
          double NLfactor = 0.0,
          bool InnerGS = false,
          bool NoInv = false,
          bool HCS = false,
          bool XFA = false,
          int NumXFA = 3,
          double R2 = 0.5,
          double gc0 = 0.5, 
          double df0 = 1.0, 
          double weight_prior_h2 = 0.01,
          double weight_prior_gc = 0.01,
          double PenCor = 0.0,
          double MinCor = 1.0,
          double uncorH2below = 0.0,
          double roundGCupFrom = 1.0,
          double roundGCupTo = 1.0,
          double roundGCdownFrom = 1.0,
          double roundGCdownTo = 0.0,
          double bucketGCfrom = 1.0,
          double bucketGCto = 1.0,
          double DeflateMax = 0.9,
          double DeflateBy = 0.0,
          bool OneVarB = false,
          bool OneVarE = false,
          bool verbose = false){
  
  //Set multi-core processing
  if(cores!=1) Eigen::setNbThreads(cores);
  
  // Gather basic info
  int k = Y.cols(), n0 = Y.rows(), p = X.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXd Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXd n = Z.colwise().sum();
  Eigen::VectorXd iN = n.array().inverse();
  
  // Centralize y
  Eigen::VectorXd mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXd y(n0,k);
  for(int i=0; i<k; i++){y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  
  // Center X
  Eigen::VectorXd xx = X.colwise().mean();
  for(int i=0; i<p; i++){ X.col(i) = X.col(i).array() - xx(i);}
  
  // Sum of squares of X
  Eigen::MatrixXd XX(p,k);
  for(int i=0; i<p; i++){
    XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  
  // Compute Tr(XSX);
  Eigen::MatrixXd XSX(p,k);
  for(int i=0; i<p; i++){
    XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
      ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXd MSx = XSX.colwise().sum();
  Eigen::VectorXd TrXSX = n.array()*MSx.array();
  
  // Variances
  iN = (n.array()-1).inverse();
  Eigen::VectorXd vy = y.colwise().squaredNorm(); vy = vy.array() * iN.array();
  
  Eigen::VectorXd ve = vy * (1-R2);
  Eigen::VectorXd iVe = ve.array().inverse();
  Eigen::MatrixXd vb(k,k), TildeHat(k,k);
  Eigen::VectorXd vbInit = ((vy*R2).array()/MSx.array());
  Eigen::VectorXd veInit = ve*1.0;
  vb = vbInit.array().matrix().asDiagonal();
  Eigen::MatrixXd iG = vb.inverse();
  Eigen::VectorXd h2 = 1 - ve.array()/vy.array();
  
  // Starting covariance values
  double tmp;
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      if(i>j){
        tmp = gc0 * sqrt(vb(i,i)*vb(j,j));
        vb(i,j) = tmp;
        vb(j,i) = tmp;
      }
    }
  }
  
  // Beta tilde;
  Eigen::MatrixXd tilde = X.transpose() * y;
  Eigen::VectorXd TrDinvXSX(k);
  Eigen::MatrixXd Dinv(p,k);
  if(TH){
    for(int i=0; i<k; i++){
      XSX.col(i) = XSX.col(i).array() * n(i);
    }
  }
  
  // Prior shape
  Eigen::MatrixXd Sb = vb*df0;
  Eigen::VectorXd Se = ve*df0;
  Eigen::VectorXd iNp = (n.array()+df0-1).inverse();
  
  // Initialize coefficient matrices
  Eigen::MatrixXd LHS(k,k);
  Eigen::VectorXd RHS(k);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(p,k);
  Eigen::VectorXd b0(k), b1(k);
  Eigen::MatrixXd e(n0,k); e = y*1.0;
  
  // Bending and convergence control
  Eigen::MatrixXd A = vb*1.0, GC(k,k);
  double bucketMean = 0.5*(bucketGCfrom+bucketGCto);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> EVDofA(A);
  double MinDVb, inflate = 0.0, Deflate = 1.0;
  Eigen::MatrixXd beta0(p,k), vb0(k,k);
  Eigen::VectorXd CNV1(maxit),CNV2(maxit),CNV3(maxit), ve0(k), h20(k);
  double cnv = 10.0;
  int numit = 0;
  double logtol = log10(tol);
  
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Inner RGS
  std::vector<int> InnerRGSvec(k);
  for(int j=0; j<k; j++){InnerRGSvec[j]=j;}
  std::random_device rd2;
  std::mt19937 g2(rd2());
  int ri;
  
  // Non-Linear weights for marker effects
  bool NonLinear = NLfactor!=0.0;
  Eigen::MatrixXd W(p,k);
  for(int i=0; i<p; i++){ for(int j=0; j<k; j++){  W(i,j) = 1.0; }}
  Eigen::VectorXd iVeWj = iVe*1.0;
  Eigen::VectorXd tmpW(p);
  double maxW, minW;
  
  // Objects for other variance structures
  double gs;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(vb);
  Eigen::MatrixXd UDU(k,k);
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    vb0 = vb*1.0;
    ve0 = ve*1.0;
    h20 = h2*1.0;
    
    // Randomized Gauss-Seidel loop
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    std::shuffle(InnerRGSvec.begin(), InnerRGSvec.end(), g2);
    
    for(int j=0; j<p; j++){
      J = RGSvec[j];
      
      // System of equations - Traditional vs Stranden and Garrick 2009
      if(NoInv){
        LHS = vb * (XX.row(J).transpose().array() * iVeWj.array()).matrix().asDiagonal(); 
        for(int i=0; i<k; i++){ LHS(i,i) += 1.0; }
        RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
        RHS = (vb * (RHS.array() * iVeWj.array()).matrix()).array();
      }else{
        LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVeWj.array()).matrix();
        RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
        RHS = RHS.array() * iVeWj.array();
      }
      
      // Update coefficient
      b0 = b.row(J)*1.0;
      for(int i=0; i<k; i++){ iVeWj(i) = iVe(i)*W(J,i); }
      LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVeWj.array()   ).matrix();
      RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
      RHS = RHS.array() * iVeWj.array();
      
      // Inner GS
      if(InnerGS){
        b1 = b.row(J)*1.0;
        for(int i=0; i<k; i++){
          ri = InnerRGSvec[i];
          b1(ri) = (RHS(ri)-(LHS.col(ri).array()*b1.array()).sum()+LHS(ri,ri)*b1(ri))/LHS(ri,ri);}
      }else{
        b1 = LHS.llt().solve(RHS); 
      }
      
      // Update residuals
      b.row(J) = b1;
      e = (e-(X.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Update marker weights
    if(NonLinear){
      W = b.cwiseAbs();
      for(int j=0; j<k; j++){
        maxW = W.col(j).maxCoeff(); minW = W.col(j).minCoeff();
        tmpW = NLfactor * (W.col(j).array()-minW)/(maxW-minW) + (1.0-NLfactor);
        tmpW = tmpW.array() + (1.0-tmpW.mean());
        W.col(j) = tmpW.array();
      }
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = (ve.array()+Se.array()) * iNp.array();
    h2 = 1 - ve.array()/vy.array();
    // Proportion-based prior
    if(weight_prior_h2>0){for(int i=0; i<k; i++){gs = ve(i)*(1-weight_prior_h2) + weight_prior_h2*veInit(i); ve(i) = gs*1.0;}}
    // Single variance
    if(OneVarE){tmp = ve.array().mean(); for(int i=0; i<k; i++) ve(i) = tmp*1.0;}
    iVe = ve.array().inverse();
    iVeWj = iVe*1.0;
    
    //Genetic variance
    
    // Get tilde-hat
    if(TH){
      for(int i=0; i<k; i++){
        Dinv.col(i) = (XSX.col(i).array()/ve(i) + iG(i,i)).inverse().array();
        TrDinvXSX(i)  = (XSX.col(i).transpose() * Dinv.col(i));
      }
      TildeHat = b.transpose()* Dinv.cwiseProduct(tilde);
    }else{
      TildeHat = b.transpose()*tilde;
    }
    
    // Estimate variances and covariance components
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i==j){ // Variances
          if(TH){
            vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrDinvXSX(i)+df0);
          }else{
            vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrXSX(i)+df0);
          }
        }else{ // Covariances
          if(TH){
            vb(i,j) = (TildeHat(i,j)+TildeHat(j,i)+Sb(i,j))/(TrDinvXSX(i)+TrDinvXSX(j)+df0);
          }else{
            vb(i,j) = (TildeHat(i,j)+TildeHat(j,i)+Sb(i,j))/(TrXSX(i)+TrXSX(j)+df0);
          }
        }}}
    
    
    if(weight_prior_h2>0){ // Proportion-based prior H2
      for(int i=0; i<k; i++){gs = vb(i,i)*(1-weight_prior_h2) + weight_prior_h2*vbInit(i); vb(i,i) = gs*1.0;}}
    if(weight_prior_gc>0){ // Proportion-based prior GC
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){
        if(i!=j){ GC(i,j) = (1.0-weight_prior_gc)*vb(i,j)/(sqrt(vb(i,i)*vb(j,j))) + gc0*weight_prior_gc;}else{GC(i,j) = 1.0;}}}
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){ if(i!=j){ vb(i,j) = GC(i,j)*sqrt(vb(i,i)*vb(j,j));}}}}else{
        // Once calculation of GC without prior
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}}
      
      // Heterogeneous Compound Symmetry
      if(HCS){
        gs = 0.0;
        for(int i=0; i<k; i++){
          for(int j=0; j<k; j++){
            if(i>j){gs += GC(i,j);}}}
        gs = gs/((k*(k-1))/2);
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){ 
          if(i!=j){ GC(i,j) =  gs*1.0;}else{ GC(i,j) = 1.0; }}}
        // Extended Factor Analytics
      }else if(XFA){
        es.compute(GC);
        UDU = es.eigenvalues()[k] * es.eigenvectors().col(k) * es.eigenvectors().col(k).transpose();
        for(int i=1; i<NumXFA; i++) UDU += es.eigenvalues()[k-i] * es.eigenvectors().col(k-i) * es.eigenvectors().col(k-i).transpose();
        GC = UDU * 1.0; for(int i=0; i<k; i++){ GC(i,i)=1.0; };
      }
      
      // Monkeying with the correlations
      for(int i=0; i<k; i++){
        for(int j=0; j<k; j++){
          if(i!=j){
            // Zero'ing  Correlations
            if(MinCor<1.0){ if(GC(i,j)<MinCor){ GC(i,j) = 0.0; }}
            // Penalize Correlations
            if(PenCor>0.0){  GC(i,j) = tanh(PenCor*abs(GC(i,j)))*GC(i,j);} 
            // Round Down
            if(roundGCdownFrom<1.0){ if(GC(i,j)<roundGCdownFrom){ GC(i,j) = roundGCdownTo*1.0; }}
            // Round Up
            if(roundGCupFrom<1.0){ if(GC(i,j)>roundGCupFrom){ GC(i,j) = roundGCupTo*1.0; }}
            // Bucket round
            if(bucketGCfrom<1.0){ if(GC(i,j)>bucketGCfrom && GC(i,j)<bucketGCto  ){ GC(i,j) =  bucketMean*1.0; }}
            // Min H2
            if(uncorH2below>0.0){ if(h2(i)<uncorH2below || h2(j)<uncorH2below  ){ GC(i,j) = 0.0; }}
          }}}
      
      // BEND AND RECONSTRUCT COVARIANCE HERE AND ONLY ONCE
      if(!NoInv||TH){ 
        A = GC*1.0;
        // Deflate
        if(DeflateBy>0){
          A = GC*Deflate; for(int i=0; i<k; i++){ A(i,i)=1.0; }
          if(A.llt().info()==Eigen::NumericalIssue && Deflate>DeflateMax){
            Deflate -= DeflateBy; 
            if(verbose){Rcpp::Rcout << "Deflate GC " <<  Deflate << '\n';}
            A = GC*Deflate; for(int i=0; i<k; i++){ A(i,i)=1.0;}}}
        // Bend
        EVDofA.compute(A); MinDVb = EVDofA.eigenvalues().minCoeff();
        if( MinDVb < 0.0 ){ inflate = abs(MinDVb*1.1);
          if(verbose) Rcpp::Rcout << "Inflate " << inflate << "\n";
          A.diagonal().array()+=inflate; A/=(1.0+inflate); GC=A*1.0;}}
      // Cor to Cov
      if(OneVarB){tmp = TildeHat.diagonal().mean(); vb=GC*tmp; }else{
        for(int i=0; i<k; i++){ for(int j=0; j<k; j++){
          vb(i,j) = GC(i,j)*sqrt(vb(i,i)*vb(j,j));}}}
      if(!NoInv||TH){ 
        iG=vb.completeOrthogonalDecomposition().pseudoInverse();}
      
      // Compute convergence and print status
      
      //cnv = log10((beta0.array()-b.array()).square().sum());
      cnv = log10((beta0.array()-b.array()).square().colwise().sum().maxCoeff());
      CNV1(numit) = cnv; if(std::isnan(cnv)){ if(verbose){Rcpp::Rcout << "Numerical issue! Job aborted (it=" << numit << ")\n";} break;}
      CNV2(numit) = log10((h20.array()-h2.array()).square().sum());
      CNV3(numit) = log10((vb0.array()-vb.array()).square().sum());
      
      // Print
      ++numit;
      if( verbose && numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; } 
      if( cnv<logtol ){ if(verbose){Rcpp::Rcout << "Model coverged in "<< numit << " iterations\n";} break;}
      if( numit == maxit && verbose){ Rcpp::Rcout << "Model did not converge\n"; }    
    
  }
  
  // Fitting the model
  Eigen::MatrixXd hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  
  // Resize convergence vectors
  Eigen::VectorXd CNV1b(numit),CNV2b(numit),CNV3b(numit);
  for(int i=0; i<numit; i++){
    CNV1b(i)=CNV1(i);
    CNV2b(i)=CNV2(i);
    CNV3b(i)=CNV3(i);
  }
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC,
                            Rcpp::Named("vb")=vb,
                            Rcpp::Named("ve")=ve,
                            Rcpp::Named("MSx")=MSx,
                            Rcpp::Named("cnvB")=CNV1b,
                            Rcpp::Named("cnvH2")=CNV2b,
                            Rcpp::Named("cnvV")=CNV3b,
                            Rcpp::Named("b_Weights")=W,
                            Rcpp::Named("Its")=numit);
}

// [[Rcpp::export]]
SEXP MRR3F(Eigen::MatrixXf Y,
          Eigen::MatrixXf X,
          int maxit = 500,
          float tol = 10e-9,
          int cores = 1,
          bool TH = false,
          float NonLinearFactor = 0.0,
          bool InnerGS = false,
          bool NoInv = false,
          bool HCS = false,
          bool XFA = false,
          int NumXFA = 3,
          float R2 = 0.5,
          float gc0 = 0.5, 
          float df0 = 1.0, 
          float weight_prior_h2 = 0.01,
          float weight_prior_gc = 0.01,
          float PenCor = 0.0,
          float MinCor = 1.0,
          float uncorH2below = 0.0,
          float roundGCupFrom = 1.0,
          float roundGCupTo = 1.0,
          float roundGCdownFrom = 1.0,
          float roundGCdownTo = 0.0,
          float bucketGCfrom = 1.0,
          float bucketGCto = 1.0,
          float DeflateMax = 0.9,
          float DeflateBy = 0.0,
          bool OneVarB = false,
          bool OneVarE = false,
          bool verbose = false){
  
  //Set multi-core processing
  if(cores!=1) Eigen::setNbThreads(cores);
  
  // Gather basic info
  int k = Y.cols(), n0 = Y.rows(), p = X.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXf Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXf n = Z.colwise().sum();
  Eigen::VectorXf iN = n.array().inverse();
  
  // Centralize y
  Eigen::VectorXf mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXf y(n0,k);
  for(int i=0; i<k; i++){y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  
  // Center X
  Eigen::VectorXf xx = X.colwise().mean();
  for(int i=0; i<p; i++){ X.col(i) = X.col(i).array() - xx(i);}
  
  // Sum of squares of X
  Eigen::MatrixXf XX(p,k);
  for(int i=0; i<p; i++){
    XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  
  // Compute Tr(XSX);
  Eigen::MatrixXf XSX(p,k);
  for(int i=0; i<p; i++){
    XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
      ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXf MSx = XSX.colwise().sum();
  Eigen::VectorXf TrXSX = n.array()*MSx.array();
  
  // Variances
  iN = (n.array()-1).inverse();
  Eigen::VectorXf vy = y.colwise().squaredNorm(); vy = vy.array() * iN.array();
  
  Eigen::VectorXf ve = vy * (1-R2);
  Eigen::VectorXf iVe = ve.array().inverse();
  Eigen::MatrixXf vb(k,k), TildeHat(k,k);
  Eigen::VectorXf vbInit = ((vy*R2).array()/MSx.array());
  Eigen::VectorXf veInit = ve*1.0;
  vb = vbInit.array().matrix().asDiagonal();
  Eigen::MatrixXf iG = vb.inverse();
  Eigen::VectorXf h2 = 1 - ve.array()/vy.array();
  
  // Starting covariance values
  float tmp;
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      if(i>j){
        tmp = gc0 * sqrt(vb(i,i)*vb(j,j));
        vb(i,j) = tmp;
        vb(j,i) = tmp;
      }
    }
  }
  
  // Beta tilde;
  Eigen::MatrixXf tilde = X.transpose() * y;
  Eigen::VectorXf TrDinvXSX(k);
  Eigen::MatrixXf Dinv(p,k);
  if(TH){
    for(int i=0; i<k; i++){
      XSX.col(i) = XSX.col(i).array() * n(i);
    }
  }
  
  // Prior shape
  Eigen::MatrixXf Sb = vb*df0;
  Eigen::VectorXf Se = ve*df0;
  Eigen::VectorXf iNp = (n.array()+df0-1).inverse();
  
  // Initialize coefficient matrices
  Eigen::MatrixXf LHS(k,k);
  Eigen::VectorXf RHS(k);
  Eigen::MatrixXf b = Eigen::MatrixXf::Zero(p,k);
  Eigen::VectorXf b0(k), b1(k);
  Eigen::MatrixXf e(n0,k); e = y*1.0;
  
  // Bending and convergence control
  Eigen::MatrixXf A = vb*1.0, GC(k,k);
  float bucketMean = 0.5*(bucketGCfrom+bucketGCto);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> EVDofA(A);
  float MinDVb, inflate = 0.0, Deflate = 1.0;
  Eigen::MatrixXf beta0(p,k), vb0(k,k);
  Eigen::VectorXf CNV1(maxit),CNV2(maxit),CNV3(maxit), ve0(k), h20(k);
  float cnv = 10.0;
  int numit = 0;
  float logtol = log10(tol);
  
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Inner RGS
  std::vector<int> InnerRGSvec(k);
  for(int j=0; j<k; j++){InnerRGSvec[j]=j;}
  std::random_device rd2;
  std::mt19937 g2(rd2());
  int ri;
  
  // Non-Linear weights for marker effects
  bool NonLinear = NonLinearFactor!=0.0;
  Eigen::MatrixXf W(p,k);
  for(int i=0; i<p; i++){ for(int j=0; j<k; j++){  W(i,j) = 1.0; }}
  Eigen::VectorXf iVeWj = iVe*1.0;
  Eigen::VectorXf tmpW(p);
  float maxW, minW;
  
  // Objects for other variance structures
  float gs;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(vb);
  Eigen::MatrixXf UDU(k,k);
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    vb0 = vb*1.0;
    ve0 = ve*1.0;
    h20 = h2*1.0;
    
    // Randomized Gauss-Seidel loop
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    std::shuffle(InnerRGSvec.begin(), InnerRGSvec.end(), g2);
    
    for(int j=0; j<p; j++){
      
      J = RGSvec[j];
      b0 = b.row(J)*1.0;
      
      // System of equations - Traditional vs Stranden and Garrick 2009
      if(NoInv){
        LHS = vb * (XX.row(J).transpose().array() * iVeWj.array()).matrix().asDiagonal(); 
        for(int i=0; i<k; i++){ LHS(i,i) += 1.0; }
        RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
        RHS = (vb * (RHS.array() * iVeWj.array()).matrix()).array();
      }else{
        LHS = iG;  LHS.diagonal() += (XX.row(J).transpose().array() * iVeWj.array()).matrix();
        RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
        RHS = RHS.array() * iVeWj.array();
      }
      
      // Inner GS
      if(InnerGS){
        b1 = b.row(J)*1.0;
        for(int i=0; i<k; i++){
          ri = InnerRGSvec[i];
          b1(ri) = (RHS(ri)-(LHS.col(ri).array()*b1.array()).sum()+LHS(ri,ri)*b1(ri))/LHS(ri,ri);}
      }else{
        b1 = LHS.llt().solve(RHS); 
      }
      
      // Update residuals
      b.row(J) = b1;
      e = (e-(X.col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
    }
    
    // Update marker weights
    if(NonLinear){
      W = b.cwiseAbs();
      for(int j=0; j<k; j++){
        maxW = W.col(j).maxCoeff(); minW = W.col(j).minCoeff();
        tmpW = NonLinearFactor * (W.col(j).array()-minW)/(maxW-minW) + (1.0-NonLinearFactor);
        tmpW = tmpW.array() + (1.0-tmpW.mean());
        W.col(j) = tmpW.array();
      }
    }
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = (ve.array()+Se.array()) * iNp.array();
    h2 = 1 - ve.array()/vy.array();
    // Proportion-based prior
    if(weight_prior_h2>0){for(int i=0; i<k; i++){gs = ve(i)*(1-weight_prior_h2) + weight_prior_h2*veInit(i); ve(i) = gs*1.0;}}
    // Single variance
    if(OneVarE){tmp = ve.array().mean(); for(int i=0; i<k; i++) ve(i) = tmp*1.0;}
    iVe = ve.array().inverse();
    iVeWj = iVe*1.0;
    
    //Genetic variance
    
    // Get tilde-hat
    if(TH){
      for(int i=0; i<k; i++){
        Dinv.col(i) = (XSX.col(i).array()/ve(i) + iG(i,i)).inverse().array();
        TrDinvXSX(i)  = (XSX.col(i).transpose() * Dinv.col(i));
      }
      TildeHat = b.transpose()* Dinv.cwiseProduct(tilde);
    }else{
      TildeHat = b.transpose()*tilde;
    }
    
    // Estimate variances and covariance components
    for(int i=0; i<k; i++){
      for(int j=0; j<k; j++){
        if(i==j){ // Variances
          if(TH){
            vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrDinvXSX(i)+df0);
          }else{
            vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrXSX(i)+df0);
          }
        }else{ // Covariances
          if(TH){
            vb(i,j) = (TildeHat(i,j)+TildeHat(j,i)+Sb(i,j))/(TrDinvXSX(i)+TrDinvXSX(j)+df0);
          }else{
            vb(i,j) = (TildeHat(i,j)+TildeHat(j,i)+Sb(i,j))/(TrXSX(i)+TrXSX(j)+df0);
          }
        }}}
    
    
    if(weight_prior_h2>0){ // Proportion-based prior H2
      for(int i=0; i<k; i++){gs = vb(i,i)*(1-weight_prior_h2) + weight_prior_h2*vbInit(i); vb(i,i) = gs*1.0;}}
    if(weight_prior_gc>0){ // Proportion-based prior GC
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){
        if(i!=j){ GC(i,j) = (1.0-weight_prior_gc)*vb(i,j)/(sqrt(vb(i,i)*vb(j,j))) + gc0*weight_prior_gc;}else{GC(i,j) = 1.0;}}}
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){ if(i!=j){ vb(i,j) = GC(i,j)*sqrt(vb(i,i)*vb(j,j));}}}}else{
        // Once calculation of GC without prior
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}}
      
      // Heterogeneous Compound Symmetry
      if(HCS){
        gs = 0.0;
        for(int i=0; i<k; i++){
          for(int j=0; j<k; j++){
            if(i>j){gs += GC(i,j);}}}
        gs = gs/((k*(k-1))/2);
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){ 
          if(i!=j){ GC(i,j) =  gs*1.0;}else{ GC(i,j) = 1.0; }}}
        // Extended Factor Analytics
      }else if(XFA){
        es.compute(GC);
        UDU = es.eigenvalues()[k] * es.eigenvectors().col(k) * es.eigenvectors().col(k).transpose();
        for(int i=1; i<NumXFA; i++) UDU += es.eigenvalues()[k-i] * es.eigenvectors().col(k-i) * es.eigenvectors().col(k-i).transpose();
        GC = UDU * 1.0; for(int i=0; i<k; i++){ GC(i,i)=1.0; };
      }
      
      // Monkeying with the correlations
      for(int i=0; i<k; i++){
        for(int j=0; j<k; j++){
          if(i!=j){
            // Zero'ing  Correlations
            if(MinCor<1.0){ if(GC(i,j)<MinCor){ GC(i,j) = 0.0; }}
            // Penalize Correlations
            if(PenCor>0.0){  GC(i,j) = tanh(PenCor*abs(GC(i,j)))*GC(i,j);} 
            // Round Down
            if(roundGCdownFrom<1.0){ if(GC(i,j)<roundGCdownFrom){ GC(i,j) = roundGCdownTo*1.0; }}
            // Round Up
            if(roundGCupFrom<1.0){ if(GC(i,j)>roundGCupFrom){ GC(i,j) = roundGCupTo*1.0; }}
            // Bucket round
            if(bucketGCfrom<1.0){ if(GC(i,j)>bucketGCfrom && GC(i,j)<bucketGCto  ){ GC(i,j) =  bucketMean*1.0; }}
            // Min H2
            if(uncorH2below>0.0){ if(h2(i)<uncorH2below || h2(j)<uncorH2below  ){ GC(i,j) = 0.0; }}
          }}}
      
      // BEND AND RECONSTRUCT COVARIANCE HERE AND ONLY ONCE
      if(!NoInv||TH){ 
        A = GC*1.0;
        // Deflate
        if(DeflateBy>0){
          A = GC*Deflate; for(int i=0; i<k; i++){ A(i,i)=1.0; }
          if(A.llt().info()==Eigen::NumericalIssue && Deflate>DeflateMax){
            Deflate -= DeflateBy; 
            if(verbose){Rcpp::Rcout << "Deflate GC " <<  Deflate << '\n';}
            A = GC*Deflate; for(int i=0; i<k; i++){ A(i,i)=1.0;}}}
        // Bend
        EVDofA.compute(A); MinDVb = EVDofA.eigenvalues().minCoeff();
        if( MinDVb < 0.0 ){ inflate = abs(MinDVb*1.1);
          if(verbose) Rcpp::Rcout << "Inflate " << inflate << "\n";
          A.diagonal().array()+=inflate; A/=(1.0+inflate); GC=A*1.0;}}
      // Cor to Cov
      if(OneVarB){tmp = TildeHat.diagonal().mean(); vb=GC*tmp; }else{
        for(int i=0; i<k; i++){ for(int j=0; j<k; j++){
          vb(i,j) = GC(i,j)*sqrt(vb(i,i)*vb(j,j));}}}
      if(!NoInv||TH){ 
        iG=vb.completeOrthogonalDecomposition().pseudoInverse();}
      
      // Compute convergence and print status
      
      //cnv = log10((beta0.array()-b.array()).square().sum());
      cnv = log10((beta0.array()-b.array()).square().colwise().sum().maxCoeff());
      CNV1(numit) = cnv; if(std::isnan(cnv)){ if(verbose){Rcpp::Rcout << "Numerical issue! Job aborted (it=" << numit << ")\n";} break;}
      CNV2(numit) = log10((h20.array()-h2.array()).square().sum());
      CNV3(numit) = log10((vb0.array()-vb.array()).square().sum());
      
      // Print
      ++numit;
      if( verbose && numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; } 
      if( cnv<logtol ){ if(verbose){Rcpp::Rcout << "Model coverged in "<< numit << " iterations\n";} break;}
      if( numit == maxit && verbose){ Rcpp::Rcout << "Model did not converge\n"; }
      
  }
  
  // Fitting the model
  Eigen::MatrixXf hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  
  // Resize convergence vectors
  Eigen::VectorXf CNV1b(numit),CNV2b(numit),CNV3b(numit);
  for(int i=0; i<numit; i++){
    CNV1b(i)=CNV1(i);
    CNV2b(i)=CNV2(i);
    CNV3b(i)=CNV3(i);
  }
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC,
                            Rcpp::Named("vb")=vb,
                            Rcpp::Named("ve")=ve,
                            Rcpp::Named("MSx")=MSx,
                            Rcpp::Named("cnvB")=CNV1b,
                            Rcpp::Named("cnvH2")=CNV2b,
                            Rcpp::Named("cnvV")=CNV3b,
                            Rcpp::Named("b_Weights")=W,
                            Rcpp::Named("Its")=numit);
}


// [[Rcpp::export]]
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

// [[Rcpp::export]]
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

// [[Rcpp::export]]
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

Eigen::MatrixXd LatentSpaces(Eigen::MatrixXd Y, Eigen::MatrixXd X, Eigen::MatrixXd BETA, int NPC = 0){
  int n=Y.rows(),k=Y.cols();
  Eigen::MatrixXd Y2 = GetImputedY(Y,X,BETA);
  Eigen::VectorXd SD = Y2.colwise().squaredNorm().array(); SD = (SD.array()/(n-1)).sqrt();
  for(int i=0; i<k; i++){ Y2.col(i) /= SD(i);};
  Eigen::BDCSVD<Eigen::MatrixXd> svd(Y2, Eigen::ComputeThinU | Eigen::ComputeThinV );
  Eigen::MatrixXd LS = svd.matrixU() * svd.singularValues().matrix().asDiagonal();
  if(NPC<0) NPC = round(2*sqrt(svd.matrixU().cols()));
  if(NPC==0) NPC += svd.matrixU().cols();
  return LS.leftCols(NPC);}

// [[Rcpp::export]]
SEXP MEGA(Eigen::MatrixXd Y, Eigen::MatrixXd X, int npc = -1){
  int n0=Y.rows(), p1=X.cols(), k=Y.cols(); Eigen::MatrixXi W(n0,k);
  for(int i=0;i<n0;i++){for(int j=0;j<k;j++){if(std::isnan(Y(i,j))){W(i,j)=0;}else{W(i,j)=1;}}}
  Eigen::MatrixXd BETA = UVBETA(Y,X);
  Eigen::MatrixXd LS = LatentSpaces(Y,X,BETA,npc);
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
SEXP GSEM(Eigen::MatrixXd Y, Eigen::MatrixXd X, int npc = -1){
  int n0=Y.rows(), p1=X.cols(), k=Y.cols(); Eigen::MatrixXi W(n0,k);
  for(int i=0;i<n0;i++){for(int j=0;j<k;j++){if(std::isnan(Y(i,j))){W(i,j)=0;}else{W(i,j)=1;}}}
  Eigen::MatrixXd BETA = UVBETA(Y,X);
  Eigen::BDCSVD<Eigen::MatrixXd> svd(X*BETA, Eigen::ComputeThinU | Eigen::ComputeThinV );
  if(npc<0) npc = round(2*sqrt(svd.matrixU().cols()));
  if(npc==0) npc += svd.matrixU().cols();
  Eigen::MatrixXd LS = (svd.matrixU() * svd.singularValues().matrix().asDiagonal()).leftCols(npc);
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

// [[Rcpp::export]]
Eigen::VectorXf solver1xF(Eigen::VectorXf Y, Eigen::MatrixXf X,
                         int maxit = 100, float tol = 10e-7, float df0 = 20.0){
  int n = X.rows(), p = X.cols(), numit = 0, J;
  float mu = Y.mean(), mu0;
  Eigen::VectorXf y = Y.array()-mu;
  Eigen::VectorXf tilde = X.transpose() * y;
  for(int i=0; i<p; i++){ X.col(i) = X.col(i).array() - X.col(i).mean(); }
  Eigen::VectorXf XX = X.colwise().squaredNorm().array();
  float TrXSX = XX.sum();
  float MSx = TrXSX/(n-1), vy = y.transpose()*Y; vy = vy/(n-1);
  float ve = vy*0.5, vb=(vy*0.5)/(MSx);
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p), beta0(p);
  Eigen::VectorXf e = y*1.0;
  float b0, b1, lambda=ve/vb, vb0=vb*df0, ve0=ve*df0, cnv = 10.0, logtol = log10(tol);
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

// [[Rcpp::export]]
Eigen::VectorXf solver2xF(Eigen::VectorXf Y, Eigen::MatrixXf X1, Eigen::MatrixXf X2,
                         int maxit = 100, float tol = 10e-7, float df0 = 20.0){
  int n = X1.rows(), p1 = X1.cols(), p2 = X2.cols(), numit = 0, J;
  float mu = Y.mean(), mu0;
  Eigen::VectorXf y = Y.array()-mu;
  Eigen::VectorXf tilde1 = X1.transpose() * y, tilde2 = X2.transpose() * y;
  for(int i=0; i<p1; i++){ X1.col(i) = X1.col(i).array()-X1.col(i).mean();}
  for(int i=0; i<p2; i++){ X2.col(i) = X2.col(i).array()-X2.col(i).mean();}
  Eigen::VectorXf XX1 = X1.colwise().squaredNorm().array(), XX2 = X2.colwise().squaredNorm().array();
  float TrXSX1 = XX1.sum(), TrXSX2 = XX2.sum();
  float MSx1 = TrXSX1/(n-1), MSx2 = TrXSX2/(n-1),  vy=y.transpose()*Y; vy = vy/(n-1);
  float ve = vy*0.5, vb1=(vy*0.5)/(MSx1), vb2=(vy*0.5)/(MSx2), h2=0.5;
  Eigen::VectorXf b_1 = Eigen::VectorXf::Zero(p1), beta01(p1);
  Eigen::VectorXf b_2 = Eigen::VectorXf::Zero(p2), beta02(p2);
  Eigen::VectorXf e = y*1.0;
  float b0, b1, lambda1=ve/vb1, lambda2=ve/vb2, cnv=10.0, logtol=log10(tol);
  float vb01 = vb1*df0, vb02 = vb2*df0, ve0 = ve*df0;
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
  Eigen::VectorXf xxx(1+p1+p2);
  xxx(0) = mu;
  for(int j=0; j<p1 ; j++){xxx(1+j)=b_1(j);}
  for(int j=0; j<p2 ; j++){xxx(1+p1+j)=b_2(j);}
  return xxx;
}

Eigen::MatrixXf submat_fF(Eigen::MatrixXf X, Eigen::VectorXi w){
  int n=w.sum(), N=X.rows(), p=X.cols(), n0=0; Eigen::MatrixXf XX(n,p);
  for(int i=0; i<N; i++){ if(w[i]==1){ XX.row(n0) = X.row(i).array(); n0+=1;}}
  return XX;}

Eigen::VectorXf subvec_fF(Eigen::VectorXf X, Eigen::VectorXi w){
  int n=w.sum(), N=X.size(), n0=0; Eigen::VectorXf XX(n);
  for(int i=0; i<N; i++){ if(w[i]==1){ XX[n0] = X[i]; n0+=1;}}
  return XX;}

// [[Rcpp::export]]
Eigen::MatrixXf FUVBETA(Eigen::MatrixXf Y, Eigen::MatrixXf X){
  int n0=Y.rows(), p=X.cols(), k=Y.cols(); Eigen::MatrixXf BETA(p,k); Eigen::MatrixXi W(n0,k);
  for(int i=0;i<n0;i++){for(int j=0;j<k;j++){if(std::isnan(Y(i,j))){W(i,j)=0;}else{W(i,j)=1;}}}
  for(int i=0;i<k;i++){
    if(W.col(i).array().sum()>0){
      BETA.col(i) = solver1xF(
        subvec_fF( Y.col(i).array(), W.col(i).array()),
        submat_fF( X, W.col(i).array())).array();}else{
          BETA.col(i) = Eigen::VectorXf::Zero(p);}}
  return BETA;}

Eigen::MatrixXf GetImputedYF(Eigen::MatrixXf Y, Eigen::MatrixXf X, Eigen::MatrixXf BETA){
  int n0=Y.rows(),k=Y.cols();
  Eigen::VectorXf Mu = Eigen::VectorXf::Zero(k), N = Eigen::VectorXf::Zero(k);
  for(int j=0;j<k;j++){for(int i=0;i<n0;i++){
    if(!std::isnan(Y(i,j))){N(j)+=1.0;Mu(j)+=Y(i,j);}}}
  Mu = Mu.array() / N.array();
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(!std::isnan(Y(i,j))){
        Y(i,j) -= Mu(j);}else{
          Y(i,j) = X.row(i)*BETA.col(j);}}}
  return Y;}

Eigen::MatrixXf LatentSpaces(Eigen::MatrixXf Y, Eigen::MatrixXf X, Eigen::MatrixXf BETA, int NPC = 0){
  int n=Y.rows(),k=Y.cols();
  Eigen::MatrixXf Y2 = GetImputedYF(Y,X,BETA);
  Eigen::VectorXf SD = Y2.colwise().squaredNorm().array(); SD = (SD.array()/(n-1)).sqrt();
  for(int i=0; i<k; i++){ Y2.col(i) /= SD(i);};
  Eigen::BDCSVD<Eigen::MatrixXf> svd(Y2, Eigen::ComputeThinU | Eigen::ComputeThinV );
  Eigen::MatrixXf LS = svd.matrixU() * svd.singularValues().matrix().asDiagonal();
  if(NPC<0) NPC = round(2*sqrt(svd.matrixU().cols()));
  if(NPC==0) NPC += svd.matrixU().cols();
  return LS.leftCols(NPC);}

// [[Rcpp::export]]
SEXP MEGAF(Eigen::MatrixXf Y, Eigen::MatrixXf X, int npc = -1){
  int n0=Y.rows(), p1=X.cols(), k=Y.cols(); Eigen::MatrixXi W(n0,k);
  for(int i=0;i<n0;i++){for(int j=0;j<k;j++){if(std::isnan(Y(i,j))){W(i,j)=0;}else{W(i,j)=1;}}}
  Eigen::MatrixXf BETA = FUVBETA(Y,X);
  Eigen::MatrixXf LS = LatentSpaces(Y,X,BETA,npc);
  Eigen::MatrixXf LS_BETA = FUVBETA(LS,X);
  int p2 = LS.cols();
  Eigen::VectorXf xxx(1+p1+p2);
  // store outputs
  Eigen::VectorXf mu(k), h2(k);
  Eigen::MatrixXf b1(p2,k), b2(p1,k);
  for(int i=0; i<k; i++){
    xxx = solver2xF(
      subvec_fF( Y.col(i).array(), W.col(i).array()),
      submat_fF( LS, W.col(i).array()),
      submat_fF( X, W.col(i).array())).array();
    mu(i) = xxx(0);
    for(int j=0; j<p2 ; j++){b1(j,i) = xxx(1+j);}
    for(int j=0; j<p1 ; j++){b2(j,i) = xxx(1+p2+j);}
  }
  // Fitted values
  Eigen::MatrixXf end_beta = LS_BETA * b1 + b2;
  Eigen::MatrixXf hat = LS*b1+X*b2;
  Eigen::MatrixXf gebv = X*end_beta;
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
SEXP GSEMF(Eigen::MatrixXf Y, Eigen::MatrixXf X, int npc = -1){
  int n0=Y.rows(), p1=X.cols(), k=Y.cols(); Eigen::MatrixXi W(n0,k);
  for(int i=0;i<n0;i++){for(int j=0;j<k;j++){if(std::isnan(Y(i,j))){W(i,j)=0;}else{W(i,j)=1;}}}
  Eigen::MatrixXf BETA = FUVBETA(Y,X);
  Eigen::BDCSVD<Eigen::MatrixXf> svd(X*BETA, Eigen::ComputeThinU | Eigen::ComputeThinV );
  if(npc<0) npc = round(2*sqrt(svd.matrixU().cols()));
  if(npc==0) npc += svd.matrixU().cols();
  Eigen::MatrixXf LS = (svd.matrixU() * svd.singularValues().matrix().asDiagonal()).leftCols(npc);
  int p2 = LS.cols();
  Eigen::VectorXf xxx(1+p1+p2);
  // store outputs
  Eigen::VectorXf mu(k), h2(k);
  Eigen::MatrixXf b1(p2,k), b2(p1,k);
  for(int i=0; i<k; i++){
    xxx = solver2xF(
      subvec_fF( Y.col(i).array(), W.col(i).array()),
      submat_fF( LS, W.col(i).array()),
      submat_fF( X, W.col(i).array())).array();
    mu(i) = xxx(0);
    for(int j=0; j<p2 ; j++){b1(j,i) = xxx(1+j);}
    for(int j=0; j<p1 ; j++){b2(j,i) = xxx(1+p2+j);}
  }
  // Fitted values
  Eigen::MatrixXf hat = LS*b1+X*b2;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=BETA*svd.matrixV().leftCols(npc)*b1+b2,
                            Rcpp::Named("hat")=hat);
}

//// SIMPLIFIED
Eigen::VectorXf xsolver1xF(Eigen::VectorXf Y, Eigen::MatrixXf X){
  int maxit = 100, n = X.rows(), p = X.cols(), numit = 0, J;
  float mu = Y.mean(), mu0, tol = 10e-7;
  Eigen::VectorXf y = Y.array()-mu;
  Eigen::VectorXf tilde = X.transpose() * y;
  for(int i=0; i<p; i++){ X.col(i) = X.col(i).array() - X.col(i).mean(); }
  Eigen::VectorXf XX = X.colwise().squaredNorm().array();
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p), beta0(p);
  Eigen::VectorXf e = y*1.0;
  float b0, b1, lambda= XX.mean(), cnv = 10.0, logtol = log10(tol);
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  while(numit<maxit){  beta0 = b*1.0;
    std::shuffle(RGSvec.begin(),RGSvec.end(), g);
    for(int j=0; j<p; j++){
      J = RGSvec[j]; b0 = b[J]*1.0;
      b1 = (e.transpose()*X.col(J)+XX(J)*b0)/(XX[J]+lambda);
      e = e - X.col(J)*(b1-b0); b[J] = b1*1.0;}
    mu0 = e.array().mean(); mu+=mu0; e=e.array()-mu0;
    cnv = log10((beta0.array()-b.array()).square().sum());
    ++numit; if( cnv<logtol || numit == maxit || std::isnan(cnv) ) break;}
  return b;
}

// [[Rcpp::export]]
Eigen::MatrixXf XFUVBETA(Eigen::MatrixXf Y, Eigen::MatrixXf X){
  int n0=Y.rows(), p=X.cols(), k=Y.cols(); Eigen::MatrixXf BETA(p,k); Eigen::MatrixXi W(n0,k);
  for(int i=0;i<n0;i++){for(int j=0;j<k;j++){if(std::isnan(Y(i,j))){W(i,j)=0;}else{W(i,j)=1;}}}
  for(int i=0;i<k;i++){
      BETA.col(i) = xsolver1xF(
        subvec_fF( Y.col(i).array(), W.col(i).array()),
        submat_fF( X, W.col(i).array())).array();}
  return BETA;}

// [[Rcpp::export]]
SEXP XSEMF(Eigen::MatrixXf Y, Eigen::MatrixXf X, int npc = 0){
  Eigen::MatrixXf BETA = XFUVBETA(Y,X);
  Eigen::MatrixXf G = X*BETA;
  Eigen::BDCSVD<Eigen::MatrixXf> svd(G, Eigen::ComputeThinU | Eigen::ComputeThinV );
  if(npc<0) npc = round(2*sqrt(svd.matrixU().cols()));
  if(npc==0) npc += svd.matrixU().cols();
  Eigen::MatrixXf Z = (svd.matrixU() * svd.singularValues().matrix().asDiagonal()).leftCols(npc);
  Eigen::MatrixXf ALPHA = XFUVBETA(Y,Z);  
  Eigen::MatrixXf b = BETA * svd.matrixV().leftCols(npc) * ALPHA;  G = X*b;
  for(int i=0; i<k; i++){ G.col(i) = G.col(i).array() - G.col(i).mean(); }
  Eigen::VectorXf vg=G.colwise().squaredNorm(); vg/=(Y.rows()); vg=vg.array().sqrt();
  for(int i=0; i<k; i++){ G.col(i) = G.col(i).array()/vg(i); }
  Eigen::MatrixXf GC = (G.transpose()*G)/(Y.rows());
  return Rcpp::List::create(Rcpp::Named("b")=b,Rcpp::Named("GC")=GC,Rcpp::Named("hat")=G);}

Eigen::VectorXf zsolver1xF(Eigen::VectorXf Y, Eigen::MatrixXf X){
  int maxit = 100; float tol = 10e-7; float df0 = 20.0;
  int n = X.rows(), p = X.cols(), numit = 0, J;
  float mu = Y.mean(), mu0;
  Eigen::VectorXf y = Y.array()-mu;
  Eigen::VectorXf tilde = X.transpose() * y;
  for(int i=0; i<p; i++){ X.col(i) = X.col(i).array() - X.col(i).mean(); }
  Eigen::VectorXf XX = X.colwise().squaredNorm().array();
  float TrXSX = XX.sum();
  float MSx = TrXSX/(n-1), vy = y.transpose()*Y; vy = vy/(n-1);
  float ve = vy*0.5, vb=(vy*0.5)/(MSx);
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p), beta0(p);
  Eigen::VectorXf e = y*1.0;
  float b0, b1, lambda=ve/vb, vb0=vb*df0, ve0=ve*df0, cnv = 10.0, logtol = log10(tol);
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
    ve = (ve+ve0)/(n+df0);
    vb = tilde.transpose()*b;
    vb = (vb+vb0)/(TrXSX+df0);
    lambda = ve/vb;
    cnv = log10((beta0.array()-b.array()).square().sum());
    ++numit; if( cnv<logtol || numit == maxit || std::isnan(cnv) ) break;}
  Eigen::VectorXf xxx(p+2); xxx(0)=1-ve/vy; xxx(1) = mu;
  for(int j=0; j<p ; j++){xxx(2+j)=b(j);}
  return xxx;}

Eigen::MatrixXf ZFUVBETA(Eigen::MatrixXf Y, Eigen::MatrixXf X){
  int n0=Y.rows(), p=X.cols(), k=Y.cols(); Eigen::MatrixXf BETA(p+2,k); Eigen::MatrixXi W(n0,k);
  for(int i=0;i<n0;i++){for(int j=0;j<k;j++){if(std::isnan(Y(i,j))){W(i,j)=0;}else{W(i,j)=1;}}}
  for(int i=0;i<k;i++){
    if(W.col(i).array().sum()>0){
      BETA.col(i) = zsolver1xF(
        subvec_fF( Y.col(i).array(), W.col(i).array()),
        submat_fF( X, W.col(i).array())).array();}else{
          BETA.col(i) = Eigen::VectorXf::Zero(p+2);}}
  return BETA;}

// [[Rcpp::export]]
SEXP ZSEMF(Eigen::MatrixXf Y, Eigen::MatrixXf X){
  int k = Y.cols(), N = Y.rows();
  Eigen::MatrixXf BETA = ZFUVBETA(Y,X);
  Eigen::MatrixXf G = X*BETA.bottomRows(X.cols());
  Eigen::VectorXf h2 = BETA.row(0).array();
  Eigen::BDCSVD<Eigen::MatrixXf> svd(G, Eigen::ComputeThinU | Eigen::ComputeThinV );
  Eigen::MatrixXf Z = (svd.matrixU() * svd.singularValues().matrix().asDiagonal());
  Eigen::MatrixXf Coef = ZFUVBETA(Y,Z);
  // Hat and H2
  Eigen::MatrixXf beta_final = BETA.bottomRows(X.cols())*Coef.bottomRows( Z.cols())*svd.matrixV();
  G = X*beta_final; Eigen::MatrixXf hat = G * 1.0;
  Eigen::VectorXf mu = Coef.row(1).array();
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i); }
  h2 = Coef.row(0).array();
  // GC
  for(int i=0; i<k; i++){ G.col(i) = G.col(i).array() - G.col(i).mean(); }
  Eigen::VectorXf vg = G.colwise().squaredNorm(); vg /= N; vg = vg.array().sqrt();
  for(int i=0; i<k; i++){ G.col(i) = G.col(i).array() / vg(i); }
  Eigen::MatrixXf GC = (G.transpose()*G)/N;
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=beta_final,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC);}



// [[Rcpp::export]]
Eigen::MatrixXd EigenArcZ( Eigen::MatrixXd Zfndr, Eigen::MatrixXd Zsamp, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores);  
  int p = Zfndr.cols(), nf = Zfndr.rows(), ns = Zsamp.rows();
  // Centralize matrices to create relationship matrix
  Eigen::VectorXd MeanColumnZfndr = Zfndr.colwise().mean();
  for(int i=0; i<p; i++){
    Zfndr.col(i) = Zfndr.col(i).array()-MeanColumnZfndr(i);
    Zsamp.col(i) = Zsamp.col(i).array()-MeanColumnZfndr(i);}
  Eigen::MatrixXd Kff = Zfndr * Zfndr.transpose();
  Eigen::MatrixXd Kfs = Zfndr * Zsamp.transpose();
  double Kscalar, tmp, NormProd, Npi = 3.14159;
  Eigen::VectorXd DiagKff = Kff.diagonal().array();
  Eigen::VectorXd DiagKss = (Zsamp.cwiseProduct(Zsamp)).rowwise().sum();
  // Relationship(Founder)
  for(int i=0; i<nf; i++){  for(int j=0; j<nf; j++){ if(j>=i){
    NormProd = sqrt(DiagKff(i)*DiagKff(j)*1.001);
    tmp = acos( Kff(i,j)/ NormProd);
    tmp = NormProd*(sin(tmp)+(Npi-tmp)*cos(tmp));
    tmp /= Npi; Kff(i,j) = tmp*1.0; Kff(j,i) = tmp*1.0;}}}
  Kscalar = 1 / Kff.diagonal().mean(); Kff *= Kscalar; 
  // Relationship(Founder,Sample)
  for(int i=0; i<nf; i++){ for(int j=0; j<ns; j++){
    NormProd = sqrt(DiagKff(i)*DiagKss(j)*1.001);
    tmp = acos( Kfs(i,j)/NormProd );
    tmp = NormProd*(sin(tmp)+(Npi-tmp)*cos(tmp));
    tmp /= Npi; Kfs(i,j) = tmp*Kscalar;}}
  // Spectral decomposition
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Kff);
  Eigen::MatrixXd L = es.eigenvectors() * es.eigenvalues().array().rsqrt().matrix().asDiagonal();
  return Kfs.transpose() * L;
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenGauZ( Eigen::MatrixXd Zfndr, Eigen::MatrixXd Zsamp, double phi = 1.0, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores);  
  int p = Zfndr.cols(), nf = Zfndr.rows(), ns = Zsamp.rows();
  // Centralize matrices to create relationship matrix
  Eigen::MatrixXd Kff = Zfndr * Zfndr.transpose();
  Eigen::MatrixXd Kfs = Zfndr * Zsamp.transpose();
  Eigen::VectorXd DiagKff = Kff.diagonal().array();
  Eigen::VectorXd DiagKss = (Zsamp.cwiseProduct(Zsamp)).rowwise().sum();
  double tmp;
  // Relationship(Founder,Sample)
  for(int i=0; i<nf; i++){ for(int j=0; j<ns; j++){
    tmp = sqrt(DiagKff(i) + DiagKss(j) - 2*Kfs(i,j));
    Kfs(i,j) = tmp*1.0;}}
  // Relationship(Founder)
  for(int i=0; i<nf; i++){for(int j=i; j<nf; j++){
    tmp = sqrt(DiagKff(i) + DiagKff(j) - 2*Kff(i,j));
    Kff(i,j) = tmp*1.0; Kff(j,i) = tmp*1.0;}}
  // Scaler
  for(int i=0; i<nf; i++){Kff(i,i) = 0.0;}
  tmp = phi * (-nf*(nf-1)) / (Kff.colwise().sum()).sum();
  //Kff *= tmp; Kfs *= tmp;  Kff = exp(Kff); Kfs = exp(Kfs);
  for(int i=0; i<nf; i++){
    Kff.row(i) = exp(Kff.row(i).array()*tmp);
    Kfs.row(i) = exp(Kfs.row(i).array()*tmp);
  }
  // Spectral decomposition
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Kff);
  Eigen::MatrixXd L = es.eigenvectors() * es.eigenvalues().array().rsqrt().matrix().asDiagonal();
  return Kfs.transpose() * L;
}

// [[Rcpp::export]]
Eigen::MatrixXd K2X(Eigen::MatrixXd K, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(K);
  Eigen::BDCSVD<Eigen::MatrixXd> svd(K, Eigen::ComputeThinU | Eigen::ComputeThinV );
  return svd.matrixU() * svd.singularValues().matrix().asDiagonal();
}

Eigen::MatrixXd GetL(Eigen::MatrixXd A){
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  return es.eigenvectors() * es.eigenvalues().array().sqrt().matrix().asDiagonal();
}

// [[Rcpp::export]]
SEXP MvSimY(
    Eigen::MatrixXd Ufndr,
    Eigen::MatrixXd Zfndr,
    Eigen::MatrixXd Zsamp,
    Eigen::VectorXd GxY,
    Eigen::VectorXd GxL,
    Eigen::VectorXd H2plot,
    int nLoc = 20,
    int Seed = 123){
  
  int p = Zfndr.cols(), nf = Zfndr.rows(), ns = Zsamp.rows(), k = Ufndr.cols();
  
  // Normal sampler
  std::mt19937 gen(Seed);
  std::normal_distribution<> d{0,1};
  
  // Centralize matrices to create relationship matrix
  Eigen::VectorXd MeanColumnZfndr = Zfndr.colwise().mean();
  for(int i=0; i<p; i++){
    Zfndr.col(i) = Zfndr.col(i).array()-MeanColumnZfndr(i);
    Zsamp.col(i) = Zsamp.col(i).array()-MeanColumnZfndr(i);}
  Eigen::MatrixXd Kff = Zfndr * Zfndr.transpose();
  Eigen::MatrixXd Kfs = Zfndr * Zsamp.transpose();
  double Kscalar, tmp, NormProd, Npi = 3.14159;
  Eigen::VectorXd DiagKff = Kff.diagonal().array();
  Eigen::VectorXd DiagKss = (Zsamp.cwiseProduct(Zsamp)).rowwise().sum();
  // Relationship(Founder)
  for(int i=0; i<nf; i++){  for(int j=i; j<nf; j++){ if(j>=i){
    NormProd = sqrt(DiagKff(i)*DiagKff(j));
    tmp = acos( Kff(i,j)/ NormProd);
    tmp = NormProd*(sin(tmp)+(Npi-tmp)*cos(tmp));
    tmp /= Npi; Kff(i,j) = tmp*1.0; Kff(j,i) = tmp*1.0;}}}
  Kscalar = 1 / Kff.diagonal().mean(); Kff *= Kscalar; 
  // Relationship(Founder,Sample)
  for(int i=0; i<nf; i++){ for(int j=0; j<ns; j++){
    NormProd = sqrt(DiagKff(i)*DiagKss(j));
    tmp = acos( Kfs(i,j)/NormProd );
    tmp = NormProd*(sin(tmp)+(Npi-tmp)*cos(tmp));
    tmp /= Npi; Kfs(i,j) = tmp*Kscalar;}}
 
  // Intercept
  Eigen::VectorXd Mu = Ufndr.colwise().mean();
  for(int i=0; i<k; i++){ Ufndr.col(i) = Ufndr.col(i).array()-Mu(i); }
  
  // Covariances
  Eigen::MatrixXd iKu = Kff.llt().solve(Ufndr);
  Eigen::MatrixXd Vg = Ufndr.transpose() * iKu / nf;
  Eigen::MatrixXd GC = Vg*1.0;
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=Vg(i,j)/(sqrt(Vg(i,i)*Vg(j,j)));}}
  
  Eigen::VectorXd TmpVec(k);
  for(int i=0; i<k; i++){ TmpVec(i) =(1-GxY(i))/GxY(i); }
  
  Eigen::VectorXd VecVy = Vg.diagonal().array() * TmpVec.array(); 
  for(int i=0; i<k; i++){ TmpVec(i) =(1-GxL(i))/GxL(i)/nLoc; }
  
  Eigen::VectorXd VecVl = VecVy.array() * TmpVec.array(); 
  Eigen::MatrixXd Vy(k,k), Vl(k,k);
  
  Eigen::VectorXd StdDev = sqrt(VecVy.array());
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){Vy(i,j)= GC(i,j)*StdDev(i)*StdDev(j);}}
  
  StdDev = sqrt(VecVl.array());
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){Vl(i,j)= GC(i,j)*StdDev(i)*StdDev(j);}}
  for(int i=0; i<k; i++){ TmpVec(i) =(1-H2plot(i))/H2plot(i); }
  //
  Eigen::MatrixXd Ve = (Vy.diagonal().array() * TmpVec.array() / nLoc).matrix().asDiagonal(); 
  
  // // Genetic effects
  Eigen::MatrixXd S1(nf,k), S2(nf,k), S3(ns,k);
  for(int i=0; i<nf; i++){for(int j=0; j<k; j++){S1(i,j)= d(gen);S2(i,j)= d(gen);}}
  for(int i=0; i<ns; i++){for(int j=0; j<k; j++){S3(i,j)= d(gen);}}
  Eigen::MatrixXd L = GetL(Kff);
  
  // Sample from founder pop
  Eigen::MatrixXd TmpMat(k,k);
  TmpMat = GetL(Vy);
  Eigen::MatrixXd Ufy = L * S1 * TmpMat.transpose();
  TmpMat = GetL(Vl);
  Eigen::MatrixXd Ufl = L * S2 * TmpMat.transpose();
  TmpMat = GetL(Ve);
  Eigen::MatrixXd E =  S3 * TmpMat.transpose();
  // Conditional expectations
  Eigen::MatrixXd Usg = Kfs * iKu;
  Eigen::MatrixXd Usy = Kfs * Kff.llt().solve(Ufy);
  Eigen::MatrixXd Usl = Kfs * Kff.llt().solve(Ufl);
  //Phenotype
  Eigen::MatrixXd Y = Usg + Usy + Usl + E;
  for(int i=0; i<k; i++){ Y.col(i) = Y.col(i).array() + Mu(i); }
  
  // Output
  return Rcpp::List::create(
    Rcpp::Named("Kfs")=Kfs,
    Rcpp::Named("Kff")=Kff,
    Rcpp::Named("L")=L,
    Rcpp::Named("Mu")=Mu,
    Rcpp::Named("GC")=GC,
    Rcpp::Named("Vg")=Vg,
    Rcpp::Named("Vy")=Vy,
    Rcpp::Named("Vl")=Vl,
    Rcpp::Named("Ve")=Ve,
    Rcpp::Named("Ufy")=Ufy,
    Rcpp::Named("Ufl")=Ufl,
    Rcpp::Named("Usg")=Usg,
    Rcpp::Named("Usy")=Usy,
    Rcpp::Named("Usl")=Usl,
    Rcpp::Named("E")=E,
    Rcpp::Named("Y")=Y
  );
  
}

// [[Rcpp::export]]
SEXP MLM(Eigen::MatrixXd Y, Eigen::MatrixXd X, Eigen::MatrixXd Z,
         int maxit = 500, double logtol = -8, int cores = 1, bool verb = false){
  
  // Basic info
  double df0 = 1.1; 
  if(cores!=1) Eigen::setNbThreads(cores);
  int k = Y.cols(), n0 = Y.rows(), f = X.cols(), p = Z.cols();
  
  // Incidence matrix W
  Eigen::MatrixXd W(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        W(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ W(i,j) = 1.0;}}}
  Eigen::VectorXd n = W.colwise().sum();
  Eigen::VectorXd iN = (n.array()-f).inverse();
  
  // Compute SY, SZ
  Eigen::MatrixXd y(n0,k), WX(n0,f), MU(f,k), iXX(f,f);
  for(int i=0; i<k; i++){
    for(int j=0; j<f; j++){WX.col(j)=X.col(j).array()*W.col(i).array();}
    iXX = (WX.transpose()*WX).inverse();
    MU.col(i) = iXX * WX.transpose()*Y.col(i);
    y.col(i) = (Y.col(i)-WX*MU.col(i) ).array()*W.col(i).array(); }
  iXX = (X.transpose()*X).inverse();
  for(int j=0; j<p; j++){ Z.col(j) = (Z.col(j) - X*(iXX*X.transpose()*Z.col(j))).array(); }
  
  // Sum of squares of Z
  Eigen::MatrixXd ZZ(p,k); 
  for(int i=0; i<p; i++){ ZZ.row(i) = Z.col(i).array().square().matrix().transpose() * W;}
  Eigen::VectorXd TrZSZ = ZZ.colwise().sum().array();
  
  // Initialize coefficient matrices
  Eigen::MatrixXd LHS(k,k);
  Eigen::VectorXd RHS(k);
  Eigen::MatrixXd b = Eigen::MatrixXd::Zero(p,k);
  Eigen::VectorXd b0(k), b1(k);
  Eigen::MatrixXd e(n0,k); e = y*1.0;
  
  // Variances
  Eigen::VectorXd vy = y.colwise().squaredNorm(); vy=vy.array()*iN.array();
  Eigen::VectorXd ve = vy * 0.5;
  Eigen::VectorXd iVe = ve.array().inverse();
  Eigen::MatrixXd vb(k,k), TildeHat(k,k);
  vb = (ve.array()/ (TrZSZ.array()*iN.array())  ).matrix().asDiagonal();
  Eigen::MatrixXd iG = vb.inverse();
  Eigen::VectorXd h2 = 1 - ve.array()/vy.array();
  Eigen::MatrixXd tilde = Z.transpose() * y;
  
  // Bending
  Eigen::MatrixXd A = vb*1.0;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> EVDofA(A);
  double MinDVb, inflate = 0.0;
  
  // Prior for stability
  Eigen::MatrixXd Sb = vb*df0;
  Eigen::VectorXd Se = ve*df0;
  Eigen::VectorXd iNp = (n.array()+df0-f).inverse();
  
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  
  // Convergence control
  Eigen::MatrixXd beta0(p,k);
  double cnv = 10.0;
  int numit = 0;
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    
    // Randomized Gauss-Seidel loop
    std::shuffle(RGSvec.begin(), RGSvec.end(), g);
    for(int j=0; j<p; j++){
      J = RGSvec[j];
      
      // Update coefficient
      b0 = b.row(J)*1.0;
      LHS = iG;  LHS.diagonal() += (ZZ.row(J).transpose().array() * iVe.array()).matrix();
      RHS = (Z.col(J).transpose()*e).array() + ZZ.row(J).array()*b0.transpose().array();
      RHS = RHS.array() *iVe.array();
      b1 = LHS.llt().solve(RHS);
      b.row(J) = b1;
      
      // Update residuals
      e = (e-(Z.col(J)*(b1-b0).transpose()).cwiseProduct(W)).matrix();}
    
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = (ve.array()+Se.array()) * iNp.array();
    iVe = ve.array().inverse();
    
    // Genetic variance
    TildeHat = b.transpose()*tilde;
    for(int i=0; i<k; i++){for(int j=0; j<k; j++){
      if(i==j){
        vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrZSZ(i)+df0);            }else{
          vb(i,j) = (TildeHat(i,j)+TildeHat(j,i))/(TrZSZ(i)+TrZSZ(j)); }}}
    
    // Bending
    A = vb*1.0;
    EVDofA.compute(A); MinDVb = EVDofA.eigenvalues().minCoeff();
    if( MinDVb < 0.0 ){ inflate = abs(MinDVb*1.1);
      A.diagonal().array()+=inflate; vb=A*1.0;}
    iG = vb.completeOrthogonalDecomposition().pseudoInverse();
    
    // Print status
    ++numit;
    cnv = log10((beta0.array()-b.array()).square().sum());
    
    if(verb){
      if( std::isnan(cnv) ){Rcpp::Rcout << "Numerical issue! Job aborted (it=" << numit << ")\n"; break;}
      if( numit % 100 == 0 ){ Rcpp::Rcout << "Iter: "<< numit << " || log10 Conv: "<< cnv << "\n"; } 
      if(  cnv<logtol ){ Rcpp::Rcout << "Model coverged in "<< numit << " iterations\n"; break; }
      if( numit == maxit ){ Rcpp::Rcout << "Model did not converge\n";}}else if( std::isnan(cnv) ){ break;}
    
  }
  
  // Fitting the model
  h2 = 1 - ve.array()/vy.array();
  Eigen::MatrixXd hat = Z*b;
  for(int i=0; i<k; i++){ hat.col(i) = ( X * MU.col(i) + hat.col(i)).array(); }
  
  // Genetic correlations
  Eigen::MatrixXd GC(k,k);
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  
  // Name and create outputs
  Rcpp::List OutputList = Rcpp::List::create(Rcpp::Named("b")=MU,
                                             Rcpp::Named("u")=b,
                                             Rcpp::Named("hat")=hat,
                                             Rcpp::Named("h2")=h2,
                                             Rcpp::Named("GC")=GC,
                                             Rcpp::Named("vb")=vb,
                                             Rcpp::Named("ve")=ve,
                                             Rcpp::Named("cnv")=cnv,
                                             Rcpp::Named("its")=numit);
  
  // Output
  OutputList.attr("class") = "PEGSmodel";
  return OutputList;
  
}
