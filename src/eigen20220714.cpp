// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// [[Rcpp::export]]
Eigen::MatrixXd EigenARC(Eigen::MatrixXd X, bool centralizeZ = true, int cores = 2){
   // cseweb.ucsd.edu/~saul/papers/nips09_kernel.pdf
  Eigen::setNbThreads(cores); int p = X.cols(), n = X.rows(); 
  double tmp, Npi=3.1416, theta, J1, Kij, Norm;
  if(centralizeX){
    for(int i=0; i<p; i++){
      tmp = (X.col(i).array()).mean();
      X.col(i) = X.col(i).array()-tmp;}}
  Eigen::MatrixXd XXp = X*X.transpose();
  tmp = 1/(XXp.diagonal().mean()); XXp *= tmp;
  Eigen::VectorXd DiagXXp = XXp.diagonal().array();
  for(int i=0; i<n; i++){ for(int j=0; j<n; j++){ if(i>=j){ 
    Norm = sqrt(DiagXXp(i)*DiagXXp(j));
    theta = acos( XXp(i,j)/Norm);
    J1 = sin(theta) + (Npi-theta)*cos(theta);
    Kij = Norm/Npi*J1;
    XXp(i,j) = Kij*1.0; XXp(j,i) = Kij*1.0;}}}
  return XXp;}

// [[Rcpp::export]]
Eigen::MatrixXd EigenGAU(Eigen::MatrixXd X, double phi = 1.0, int cores = 2){
  Eigen::setNbThreads(cores);
  int n = X.rows(); double tmp;
  Eigen::MatrixXd XXp = X*X.transpose();
  for(int i=0; i<n; i++){ for(int j=0; j<n; j++){ if(i>j){
    tmp = sqrt(XXp(i,i) + XXp(j,j) - 2*XXp(i,j));
    XXp(i,j) = tmp*1.0; XXp(j,i) = tmp*1.0;}}};
  for(int i=0; i<n; i++){XXp(i,i) = 0.0;}
  tmp = phi * (-n*(n-1)) / (XXp.colwise().sum()).sum();
  XXp *= tmp; return exp(XXp.array());}

// [[Rcpp::export]]
Eigen::MatrixXd EigenGRM(Eigen::MatrixXd X, bool centralizeZ = true, int cores = 2){
  Eigen::setNbThreads(cores); int p = X.cols(); double tmp;
  if(centralizeZ){
    for(int i=0; i<p; i++){
      tmp = (X.col(i).array()).mean();
      X.col(i) = X.col(i).array()-tmp;}}
  Eigen::MatrixXd XXp = X*X.transpose();
  tmp = 1/(XXp.diagonal().mean());
  XXp *= tmp; return XXp;}

// [[Rcpp::export]]
Eigen::MatrixXd EigenCNT(Eigen::MatrixXd X, int cores = 2){
  Eigen::setNbThreads(cores); int p = X.cols();
  Eigen::VectorXd xx = X.colwise().mean();
  for(int i=0; i<p; i++){ X.col(i) = X.col(i).array() - xx(i);}
  return X;}

// [[Rcpp::export]]
SEXP EigenEVD(Eigen::MatrixXd A, int cores = 2){
  Eigen::setNbThreads(cores); 
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
  return Rcpp::List::create(Rcpp::Named("U")=es.eigenvectors(),
                            Rcpp::Named("D")=es.eigenvalues());}

// [[Rcpp::export]]
SEXP EigenBDCSVD(Eigen::MatrixXd X, int cores = 2){
  Eigen::setNbThreads(cores);
  Eigen::BDCSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV );
  return Rcpp::List::create(Rcpp::Named("U")=svd.matrixU(),
                            Rcpp::Named("D")=svd.singularValues(),
                            Rcpp::Named("V")=svd.matrixV());}

// [[Rcpp::export]]
SEXP EigenJacobiSVD(Eigen::MatrixXd X, int cores = 2){
  Eigen::setNbThreads(cores);
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV );
  return Rcpp::List::create(Rcpp::Named("U")=svd.matrixU(),
                            Rcpp::Named("D")=svd.singularValues(),
                            Rcpp::Named("V")=svd.matrixV());}

// [[Rcpp::export]]
Eigen::VectorXd EigenAcc(Eigen::MatrixXd X1, Eigen::MatrixXd X2, double h2 = 0.5, int cores = 2){
  Eigen::setNbThreads(cores);
  Eigen::MatrixXd X1X1 = X1*X1.transpose(), X1X2 = X1*X2.transpose();
  double Ve = (1.0-h2)/h2, alpha = 1.0/(X1X1.diagonal().array()).mean();
  Eigen::MatrixXd V = X1X1*alpha; V.diagonal() = V.diagonal().array() + Ve;
  return sqrt( alpha * (X1X2.transpose()*(V.llt().solve(X1X2))).diagonal().array()/X2.rowwise().squaredNorm().array());}
