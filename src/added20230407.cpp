// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <random>

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
  double deflate = 1.0, deflateMax = 0.9;
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
    //ve = ve.array() * iN.array();
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
    for(int i=0; i<k; i++){ for(int j=0; j<k; j++){if(i!=j){A(i,j) *= deflate;}}}
    if(A.llt().info()==Eigen::NumericalIssue && deflate>deflateMax){
      deflate -= 0.005; Rcpp::Rcout <<  deflate;
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){A(i,j) = vb(i,j)*deflate;}}}}
    EVDofA.compute(A);
    MinDVb = EVDofA.eigenvalues().minCoeff();
    if( MinDVb < 0.0 ){ 
      Rcpp::Rcout << ".";
      inflate = inflate - MinDVb*1.001;
      A.diagonal().array() += inflate;}
    iG = A.inverse();
    
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
  Rcpp::List OutputList = Rcpp::List::create(Rcpp::Named("b")=MU,Rcpp::Named("u")=b,
                                             Rcpp::Named("hat")=hat,Rcpp::Named("h2")=h2,
                                             Rcpp::Named("GC")=GC,Rcpp::Named("bend")=deflate,
                                             Rcpp::Named("vb")=vb,Rcpp::Named("ve")=ve,
                                             Rcpp::Named("cnv")=cnv,Rcpp::Named("its")=numit);
  
  // Output
  OutputList.attr("class") = "PEGSmodel";
  return OutputList;
  
}