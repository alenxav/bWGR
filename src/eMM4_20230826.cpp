// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include <random>
using Eigen::MatrixXd;  
using Eigen::VectorXd;
using Rcpp::List;
using Rcpp::Named;
using Rcpp::ListOf;
using Rcpp::Rcout;

List PrepInput( MatrixXd Z, MatrixXd y, MatrixXd W, MatrixXd X, MatrixXd iXX,
                VectorXd iN, VectorXd ve, int k, double df0 ){
  int p = Z.cols();
  // SZ
  for(int j=0; j<p; j++){ Z.col(j) = (Z.col(j) - X*(iXX*X.transpose()*Z.col(j))).array(); }
  // beta tilde
  MatrixXd tilde = Z.transpose() * y;
  // Tr(Z'SZ)
  MatrixXd ZZ(p,k); 
  for(int i=0; i<p; i++){ ZZ.row(i) = Z.col(i).array().square().matrix().transpose() * W;}
  VectorXd TrZSZ = ZZ.colwise().sum().array();
  // Starting values for covariances and coefficients
  MatrixXd b = MatrixXd::Zero(p,k);
  MatrixXd vb(k,k), GC=MatrixXd::Zero(k,k), beta0 = MatrixXd::Zero(p,k);
  vb = (ve.array()/ (TrZSZ.array()*iN.array())).matrix().asDiagonal();
  MatrixXd iG = vb.inverse();
  double deflate = 1.0, inflate=0.0;
  MatrixXd Zb = Z*b;
  // Name and create outputs
  return List::create(Named("Zb") = Zb,
                      Named("p") = p,
                      Named("Z") = Z,
                      Named("ZZ") = ZZ,
                      Named("TrZSZ") = TrZSZ,
                      Named("b") = b,
                      Named("tilde") = tilde,
                      Named("deflate") = deflate,
                      Named("inflate") = inflate,
                      Named("vb") = vb,
                      Named("iG") = iG,
                      Named("Sb") = vb*df0,
                      Named("GC") = GC);}

void UpdateRE( MatrixXd e, MatrixXd W, List LL, int which, int k, VectorXd iVe,
               VectorXd ve, VectorXd iN, double df0, double deflateMax, VectorXd CnvB){
  // Take stuff out of the list
  List L = LL[which];
  double deflate=L["deflate"],inflate=L["inflate"],p=L["p"];//, CnvB=L["CnvB"];
  MatrixXd tilde=L["tilde"],vb=L["vb"],b=L["b"],iG=L["iG"],ZZ=L["ZZ"],Z=L["Z"],Zb=L["Zb"],GC=L["GC"],Sb=L["Sb"];
  VectorXd TrZSZ=L["TrZSZ"];
  Sb.diagonal() = (df0 * ve.array()/ (TrZSZ.array()*iN.array())).array();
  // Initialize coefficient matrices
  MatrixXd LHS(k,k);
  VectorXd RHS(k), b0(k), b1(k);
  // RGS
  std::vector<int> RGSvec(p);
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  std::random_device rd;
  std::mt19937 g(rd());
  int J;
  // Store coefficients pre-iteration
  MatrixXd beta0 = b*1.0;
  // Randomized Gauss-Seidel loop
  std::shuffle(RGSvec.begin(), RGSvec.end(), g);
  for(int j=0; j<p; j++){
    J = RGSvec[j];
    // Update coefficient
    b0 = b.row(J)*1.0;
    LHS = iG;  LHS.diagonal() += (ZZ.row(J).transpose().array() * iVe.array()).matrix();
    RHS = (Z.col(J).transpose()*e).array() + ZZ.row(J).array()*b0.transpose().array();
    RHS = RHS.array() * iVe.array();
    b1 = LHS.llt().solve(RHS);
    b.row(J) = b1;
    // Update residuals
    e = (e-(Z.col(J)*(b1-b0).transpose()).cwiseProduct(W)).matrix();}
  // Genetic variance
  Zb = Z*b;
  MatrixXd TildeHat = b.transpose()*tilde;
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){ if(i==j){
        vb(i,i) = (TildeHat(i,i)+Sb(i,i))/(TrZSZ(i)+df0); }else{
        vb(i,j) = (TildeHat(i,j)+TildeHat(j,i))/(TrZSZ(i)+TrZSZ(j));}}}
  // Bending procedure 1
  MatrixXd A = vb*1.0;
  for(int i=0; i<k; i++){ for(int j=0; j<k; j++){if(i!=j){A(i,j) *= deflate;}}}
  if(A.llt().info()==Eigen::NumericalIssue && deflate>deflateMax){
    deflate -= 0.005; Rcout <<  deflate;
    for(int i=0; i<k; i++){for(int j=0; j<k; j++){if(i!=j){A(i,j) = vb(i,j)*deflate;}}}}
  // Bending procedure 2
  Eigen::SelfAdjointEigenSolver<MatrixXd> EVDofA(A);
  EVDofA.compute(A);
  double MinDVb = EVDofA.eigenvalues().minCoeff();
  if( MinDVb < 0.0 ){ //  Rcout << ".";
    inflate = 0.0 - MinDVb*1.001;
    A.diagonal().array() += inflate;}
  // Inverse G, GC and convergence
  iG = A.inverse();
  for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  CnvB[which] = log10((beta0.array()-b.array()).square().sum());
  // Updates
  L["b"] = b;
  L["vb"] = vb;
  L["iG"] = iG;
  L["GC"] = GC;
  L["CnvB"] = CnvB;
  L["deflate"] = deflate;
  L["inflate"] = inflate;
  L["Sb"] = Sb;
  L["Zb"] = Zb;
  LL[which] = L;
}

MatrixXd GetZb( List L, int WhichTerm){
  List SubList = L(WhichTerm);
  MatrixXd Zb = SubList[0];//1];
  return Zb;
}

// [[Rcpp::export]]
SEXP MLMX(MatrixXd Y, MatrixXd X, List Zlist,
         int maxit = 500, double logtol = -8, int cores = 1,
         double df0 = 0.1, double deflateMax = 0.9, bool CHECK1=false){
  // Basic info
  if(cores!=1) Eigen::setNbThreads(cores);
  int k = Y.cols(), n0 = Y.rows(), f = X.cols(), re = Zlist.size();
  VectorXd CnvB(re);
  // Incidence matrix W
  MatrixXd W(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        W(i,j) = 0.0; Y(i,j) = 0.0;
      }else{ W(i,j) = 1.0;}}}
  VectorXd n = W.colwise().sum();
  VectorXd iN = (n.array()-f).inverse();
  // Compute SY, e
  MatrixXd y(n0,k), WX(n0,f), MU(f,k), iXX(f,f);
  for(int i=0; i<k; i++){
    for(int j=0; j<f; j++){WX.col(j)=X.col(j).array()*W.col(i).array();}
    iXX = (WX.transpose()*WX).inverse();
    MU.col(i) = iXX * WX.transpose()*Y.col(i);
    y.col(i) = (Y.col(i)-WX*MU.col(i) ).array()*W.col(i).array(); }
  iXX = (X.transpose()*X).inverse();
  MatrixXd e = y*1.0;
  // Variances
  VectorXd vy = y.colwise().squaredNorm(); vy=vy.array()*iN.array();
  VectorXd ve = vy * 0.5;
  VectorXd iVe = ve.array().inverse();
  VectorXd iNp = (n.array()+df0-f).inverse();
  VectorXd Se = ve*df0;
  // Random effects
  ListOf<List> RndEff(re);
  for(int i=0; i<re; i++){RndEff[i] = PrepInput(  Zlist[i], y, W, X, iXX, iN, ve, k, df0 );}
  if(CHECK1){return(RndEff);} 
  // Convergence control
  double cnv = 10.0;
  int numit = 0;
  List tmpList;
  // Loop
  //Rcout << "Start loop \n";
  while(numit<maxit){
    // Update random effects
    for(int j=0; j<re; j++){UpdateRE( e, W, RndEff, j, k, iVe, ve, iN, df0, deflateMax, CnvB );}
    Rcout << "Update residual term\n";
    // Residual variance
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = (ve.array()+Se.array()) * iNp.array();
    iVe = ve.array().inverse();
    // Print status
    ++numit;
    cnv = CnvB.maxCoeff();
    if( std::isnan(cnv) ){ Rcout << "Numerical issue! Job aborted (it=" << numit << ")\n"; break;}
    if( numit % 100 == 0 ){ Rcout << "Iter: "<< numit << " || log10 Conv: "<< cnv << "\n"; } 
    if( cnv<logtol ){ Rcout << "Model coverged in "<< numit << " iterations\n"; break; }
  }
  // Fitting the model
  Rcout << "Fit model\n";
  MatrixXd hat = MatrixXd::Zero(n0,k), Zb(n0,k);
  for(int q=0; q<re; q++){
    Zb = GetZb(RndEff,q);
    hat = hat + Zb;}
  for(int i=0; i<k; i++){ hat.col(i) = ( X * MU.col(i) + hat.col(i)).array(); }
  // Name and create outputs
  List OutputList = List::create(
    Named("InputPrepRF") = List::create(Named("y")=y,Named("W")= W,Named("X")= X, Named("iXX")=iXX,Named("iN")= iN, Named("ve")=ve, Named("k")=k, Named("df0")=df0),
    Named("Fits") = List::create(Named("ve")=ve, Named("hat")=hat,Named("cnv")=cnv, Named("its")=numit),
    Named("RndEff")=RndEff);
  OutputList.attr("class") = "MLMX";
  return OutputList;}
