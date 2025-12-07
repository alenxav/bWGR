// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <random>
#include <algorithm>  // for std::shuffle
#include <cmath>      // for std::isnan

// [[Rcpp::export]]
SEXP PEGS(Eigen::MatrixXf Y, // matrix response variables
          Eigen::MatrixXf X, // design matrix of random effects
          int maxit = 100, // maximum number of iterations
          float logtol = -4.0, // convergence tolerance
          float covbend = 1.1, // covariance bending factor
          float covMinEv = 10e-6, // minimum eigenvalue to bend covariance
          int XFA = -1, // number of principal components to fit
          bool NNC = true){ // non-negative correlations
  
  // Get input dimensions
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
  for(int i=0; i<k; i++){
    y.col(i) = (Y.col(i).array()-mu(i)).array()*Z.col(i).array();}
  
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
  Eigen::MatrixXf iG = vb.inverse();
  Eigen::VectorXf h2 = 1 - ve.array()/vy.array();
  
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
  
  // Convergence control
  Eigen::MatrixXf beta0(p,k);
  float cnv = 10.0, MinDVb = 0.0, inflate = 0.0;
  int numit = 0;
  
  // Bending objects
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> EVDofA(vb);
  Eigen::VectorXf std_dev = vb.array().sqrt();
  Eigen::VectorXf inv_std_dev = std_dev.array().inverse();
  Eigen::MatrixXf GC = inv_std_dev.asDiagonal() * vb * inv_std_dev.asDiagonal();
  
  // XFA
  if(XFA<0) XFA = k;
  Eigen::VectorXf sd = vb.diagonal().array().sqrt();
  Eigen::VectorXf inv_sd = sd.array().inverse();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigen_solver(GC);
  Eigen::MatrixXf V_reduced = eigen_solver.eigenvectors().rightCols(XFA);
  Eigen::VectorXf D_reduced_diag = eigen_solver.eigenvalues().tail(XFA);
  
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
    
    // XFA
    if(XFA == 0){
      // If XFA is set to zero, make traits independent
      sd = vb.diagonal().array();
      vb.setZero();
      vb.diagonal() = sd.array();
    }else if(XFA>0){
      // Compute GC
      sd = vb.diagonal().array().sqrt();
      for (int t = 0; t < k; ++t) sd(t) = std::max(sd(t), 1e-12f);
      inv_sd = sd.array().inverse();
      GC = inv_sd.asDiagonal() * vb * inv_sd.asDiagonal();
      // Decompose and reconstruct GC
      eigen_solver.compute(GC);
      V_reduced = eigen_solver.eigenvectors().rightCols(XFA);
      D_reduced_diag = eigen_solver.eigenvalues().tail(XFA);
      GC = V_reduced * D_reduced_diag.asDiagonal() * V_reduced.transpose();
      GC.diagonal().setOnes();
      // Scale correlations back to covariances
      vb = sd.asDiagonal() * GC * sd.asDiagonal();
    }
    
    // Bending
    if(NNC) vb = vb.array().cwiseMax(0.0).matrix();
    EVDofA.compute(vb); MinDVb = EVDofA.eigenvalues().minCoeff();
    if( MinDVb < covMinEv ){if(abs(MinDVb*covbend)>inflate) inflate = abs(MinDVb*covbend);}
    if( k>=5 || MinDVb < covMinEv ){ vb.diagonal().array() += inflate; }
    iG = vb.completeOrthogonalDecomposition().pseudoInverse();
    
    // Update intercept
    b0 = e.colwise().sum();
    b0 = b0.array() * iN.array();
    for(int i=0; i<k; i++){ mu(i) += b0(i);
      e.col(i) = (e.col(i).array()-b0(i)).array() * Z.col(i).array();}
    
    // Print status
    cnv = log10((beta0.array()-b.array()).square().sum());  ++numit;
    if( numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; } 
    if( cnv<logtol ){break;}
  }
  
  // Fitting the model
  h2 = 1 - ve.array()/vy.array();
  Eigen::MatrixXf hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}

  // GC
  sd = vb.diagonal().array().sqrt();
  for (int t = 0; t < k; ++t) sd(t) = std::max(sd(t), 1e-12f);
  inv_sd = sd.array().inverse();
  GC = inv_sd.asDiagonal() * vb * inv_sd.asDiagonal();
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC,
                            Rcpp::Named("bend")=inflate,
                            Rcpp::Named("numit")=numit,
                            Rcpp::Named("cnv")=cnv);
}

// [[Rcpp::export]]
SEXP PEGSX(Eigen::MatrixXf Y,
           Eigen::MatrixXf X,
           Rcpp::List Z_list, 
           int maxit = 500,
           float logtol = -8.0,
           float covbend = 1.1,
           float covMinEv = 10e-6,
           int cores = 1,
           bool verbose = false,
           float df0 = 1.1,
           bool NNC = false,
           bool InnerGS = false,
           bool NoInv = false,
           bool XFA = false,
           int NumXFA = 3){
          
  if (cores != 1) Eigen::setNbThreads(cores);
  
  // Dimensions
  const int n = Y.rows();
  const int k = Y.cols();
  const int f = X.cols();
  const int R = Z_list.size();
  if (verbose) Rcpp::Rcout << "Rout: Start PEGSX\n n=" << n << " k=" << k << " f=" << f << " R=" << R
                           << " \n InnerGS=" << InnerGS
                           << " \n NoInv=" << NoInv
                           << " \n XFA=" << XFA
                           << " \n NumXFA=" << NumXFA
                           << " \n NNC=" << NNC << "\n";
  
  // Build mask W and zero-out NA in Y
  Eigen::MatrixXf W(n, k);
  for (int i = 0; i < n; ++i) {
    for (int t = 0; t < k; ++t) {
      if (std::isnan(Y(i, t))) {
        W(i, t) = 0.0f;
        Y(i, t) = 0.0f;
      } else {
        W(i, t) = 1.0f;
      }
    }
  }
  Eigen::VectorXf n_each = W.colwise().sum();
  if (verbose) Rcpp::Rcout << "Rout: Built mask W and counts per trait\n";
  
  // Masked X per trait (WX) and iXX = (WX'WX)^+
  std::vector<Eigen::MatrixXf> WX_list(k);
  std::vector<Eigen::MatrixXf> iXX_list(k);
  Eigen::MatrixXf b = Eigen::MatrixXf::Zero(f, k); // fixed-effect coefficients
  for (int t = 0; t < k; ++t) {
    WX_list[t].resize(n, f);
    for (int j = 0; j < f; ++j) WX_list[t].col(j) = X.col(j).array() * W.col(t).array();
    Eigen::MatrixXf XX = WX_list[t].transpose() * WX_list[t];
    iXX_list[t] = XX.completeOrthogonalDecomposition().pseudoInverse();
    Eigen::VectorXf RHS = WX_list[t].transpose() * Y.col(t);
    b.col(t).noalias() = iXX_list[t] * RHS;
  }
  if (verbose) Rcpp::Rcout << "Rout: Estimated fixed effects b\n";
  
  // Residuals: y = (Y - WX*b) masked by W
  Eigen::MatrixXf y(n, k);
  for (int t = 0; t < k; ++t){
    y.col(t) = (Y.col(t) - WX_list[t] * b.col(t)).array() * W.col(t).array();
  }
  if (verbose) Rcpp::Rcout << "Rout: Computed masked residuals y\n";
  
  // Precompute p_r and RGS indices
  std::vector<int> pR(R, 0);
  for (int r = 0; r < R; ++r) {
    Eigen::MatrixXf Zr_once = Rcpp::as<Eigen::MatrixXf>(Z_list[r]);
    pR[r] = Zr_once.cols();
    if (verbose) Rcpp::Rcout << "Rout: Z[" << r << "] cols=" << pR[r] << "\n";
  }
  std::vector<std::vector<int>> RGS_index(R);
  for (int r = 0; r < R; ++r) {
    RGS_index[r].resize(pR[r]);
    for (int j = 0; j < pR[r]; ++j) RGS_index[r][j] = j;
  }
  
  // === Precomputations (formerly helper functions) ===
  Rcpp::List ZZ_list_R(R);
  for (int r = 0; r < R; ++r) {
    Eigen::MatrixXf Zs = Rcpp::as<Eigen::MatrixXf>(Z_list[r]);
    if (Zs.rows() != n) Rcpp::stop("Inlined Compute_ZpZ: each Z must have n rows.");
    const int P_r = Zs.cols();
    Eigen::MatrixXf ZpZ(P_r, k);
    for (int c = 0; c < k; ++c) {
      for (int p = 0; p < P_r; ++p) {
        Eigen::VectorXf wz = Zs.col(p).array() * W.col(c).array();
        ZpZ(p, c) = wz.squaredNorm();
      }
    }
    ZZ_list_R[r] = ZpZ;
  }
  
  Rcpp::List TrZSZ_list_R(R);
  for (int r = 0; r < R; ++r) {
    Eigen::MatrixXf Zs = Rcpp::as<Eigen::MatrixXf>(Z_list[r]);
    if (Zs.rows() != n) Rcpp::stop("Inlined Compute_TrZSZ: each Z must have n rows equal to nrow(X).");
    const int q = Zs.cols();
    Eigen::VectorXf trait_traces(k);
    for (int c = 0; c < k; ++c) {
      Eigen::VectorXf w_c = W.col(c);
      if (w_c.size() != n) Rcpp::stop("Inlined Compute_TrZSZ: length(w) != nrow(X).");
      Eigen::MatrixXf M = X.array().colwise() * w_c.array();
      Eigen::MatrixXf InvMpM;
      Eigen::MatrixXf Sym = M.transpose() * M;
      Eigen::LLT<Eigen::MatrixXf> llt(Sym);
      if (llt.info() == Eigen::Success) {
        Eigen::MatrixXf I = Eigen::MatrixXf::Identity(Sym.rows(), Sym.cols());
        InvMpM = llt.solve(I);
      } else {
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXf> cod(Sym);
        InvMpM = cod.pseudoInverse();
      }
      float current_trace = 0.0f;
      Eigen::MatrixXf Mt = M.transpose();
      for (Eigen::Index j = 0; j < q; ++j) {
        Eigen::VectorXf z_trace = Zs.col(j).array() * w_c.array();
        Eigen::VectorXf beta_trace = InvMpM * (Mt * z_trace);
        Eigen::VectorXf sz = (z_trace - M * beta_trace).array() * w_c.array();
        current_trace += z_trace.dot(sz);
      }
      trait_traces(c) = current_trace;
    }
    TrZSZ_list_R[r] = trait_traces;
  }
  
  Rcpp::List tilde_list_R(R);
  for (int r = 0; r < R; ++r) {
    Eigen::MatrixXf Zs = Rcpp::as<Eigen::MatrixXf>(Z_list[r]);
    if (Zs.rows() != n) Rcpp::stop("Inlined Compute_beta_tilde: each Z must have n rows.");
    tilde_list_R[r] = Zs.transpose() * y;
  }
  if (verbose) Rcpp::Rcout << "Rout: Precomputed ZZ/TrZSZ/tilde\n";
  
  std::vector<Eigen::MatrixXf> u_list(R);
  for (int r = 0; r < R; ++r) u_list[r] = Eigen::MatrixXf::Zero(pR[r], k);
  
  Eigen::VectorXf ssy = y.colwise().squaredNorm();
  Eigen::VectorXf denom = (n_each.array() - f).matrix();
  for (int t = 0; t < k; ++t) denom(t) = std::max(denom(t), 1.0f);
  Eigen::VectorXf iN_mlm = denom.array().inverse().matrix();
  Eigen::VectorXf ve = (ssy.array() * iN_mlm.array()).matrix() * 0.5f;
  ve = ve.array().max(1e-8f);
  Eigen::VectorXf iVe = ve.array().inverse().matrix();
  if (verbose) Rcpp::Rcout << "Rout: Initialized ve\n";
  
  std::vector<Eigen::MatrixXf> vb_list(R);
  std::vector<Eigen::MatrixXf> iG_list(R);
  std::vector<float> bend_inflate(R, 0.0f);
  for (int r = 0; r < R; ++r) {
    vb_list[r].resize(k, k);
    vb_list[r].setZero();
    Eigen::VectorXf Tr_r = Rcpp::as<Eigen::VectorXf>(TrZSZ_list_R[r]);
    for (int t = 0; t < k; ++t) {
      float denom_r = std::max(Tr_r(t) * iN_mlm(t), 1e-8f);
      vb_list[r](t, t) = ve(t) / denom_r;
    }
    iG_list[r] = vb_list[r].completeOrthogonalDecomposition().pseudoInverse();
  }
  if (verbose) Rcpp::Rcout << "Rout: Initialized vb and iG for all effects\n";
  
  std::vector<Eigen::MatrixXf> Sb_list(R);
  for (int r = 0; r < R; ++r) Sb_list[r] = vb_list[r] * df0;
  Eigen::VectorXf Se = ve * df0;
  Eigen::VectorXf iNp = (n_each.array() + df0 - f).matrix();
  for (int t = 0; t < k; ++t) iNp(t) = 1.0f / std::max(iNp(t), 1.0f);
  if (verbose) Rcpp::Rcout << "Rout: Set priors\n";
  
  Eigen::MatrixXf e = y;
  
  std::vector<Eigen::MatrixXf> u0_list(R);
  float cnv = 10.0f;
  int numit = 0;
  std::random_device rd_dev;
  std::mt19937 gen(rd_dev());
  if (verbose) Rcpp::Rcout << "Rout: Starting GS iterations\n";
  
  // Iteration loop
  while (numit < maxit) {
    for (int r = 0; r < R; ++r) u0_list[r] = u_list[r];
    
    // Randomized Gauss–Seidel across effects
    for (int r = 0; r < R; ++r) {
      std::shuffle(RGS_index[r].begin(), RGS_index[r].end(), gen);
      
      Eigen::MatrixXf Zr        = Rcpp::as<Eigen::MatrixXf>(Z_list[r]);
      Eigen::MatrixXf ZZr       = Rcpp::as<Eigen::MatrixXf>(ZZ_list_R[r]);
      Eigen::MatrixXf tilde_r   = Rcpp::as<Eigen::MatrixXf>(tilde_list_R[r]);
      Eigen::VectorXf Tr_r      = Rcpp::as<Eigen::VectorXf>(TrZSZ_list_R[r]);
      
      Eigen::MatrixXf &ur  = u_list[r];
      Eigen::MatrixXf &iGr = iG_list[r];
      Eigen::MatrixXf &vbr = vb_list[r];
      
      const int p_r = Zr.cols();
      for (int jj = 0; jj < p_r; ++jj) {
        int J = RGS_index[r][jj];
        Eigen::VectorXf u0 = ur.row(J).transpose();
        Eigen::VectorXf diagVec = (ZZr.row(J).transpose().array() * iVe.array()).matrix();
        
        Eigen::MatrixXf LHS(k, k);
        Eigen::VectorXf RHS(k);
        if (NoInv) {
          LHS.noalias() = vbr * diagVec.asDiagonal();
          for (int i = 0; i < k; ++i) LHS(i, i) += 1.0f;
          Eigen::VectorXf base = (Zr.col(J).transpose() * e).transpose();
          base.array() += ZZr.row(J).transpose().array() * u0.array();
          base.array() *= iVe.array();
          RHS.noalias() = vbr * base;
        } else {
          LHS = iGr;
          LHS.diagonal().array() += diagVec.array();
          RHS = (Zr.col(J).transpose() * e).transpose();
          RHS.array() += ZZr.row(J).transpose().array() * u0.array();
          RHS.array() *= iVe.array();
        }
        
        Eigen::VectorXf u1(k);
        if (InnerGS) {
          u1 = ur.row(J).transpose();
          for (int i = 0; i < k; ++i) {
            float diag_ii = LHS(i, i);
            if (std::abs(diag_ii) < 1e-12f) continue;
            float accum = 0.0f;
            for (int l = 0; l < k; ++l) {
              if (l == i) continue;
              accum += LHS(i, l) * u1(l);
            }
            u1(i) = (RHS(i) - accum) / diag_ii;
          }
        } else {
          u1 = LHS.llt().solve(RHS);
        }
        ur.row(J) = u1.transpose();
        e = (e - (Zr.col(J) * (u1 - u0).transpose()).cwiseProduct(W)).eval();
      }
      
      // ---- vb update for effect r ----
      Eigen::MatrixXf TildeHat(k, k);
      TildeHat.noalias() = ur.transpose() * tilde_r;
      
      for (int i = 0; i < k; ++i) {
        float denom_i = std::max(Tr_r(i), 1e-8f);
        vbr(i,i) = (TildeHat(i,i) + Sb_list[r](i,i)) / (denom_i + df0);
      }
      for (int i = 0; i < k; ++i) {
        for (int j = i + 1; j < k; ++j) {
          float denom_ij = std::max(Tr_r(i) + Tr_r(j), 1e-8f);
          float cov_ij   = (TildeHat(i,j) + TildeHat(j,i)) / denom_ij;
          vbr(i,j) = cov_ij;
          vbr(j,i) = cov_ij;
        }
      }
      
      if (NNC) {
        vbr = vbr.array().cwiseMax(0.0f).matrix();
        vbr = 0.5f * (vbr + vbr.transpose());
        if (verbose) Rcpp::Rcout << "Rout: Applied NNC to vb for effect r=" << r << "\n";
      }
      
      if (XFA && NumXFA > 0) {
        Eigen::VectorXf sd = vbr.diagonal().array().sqrt();
        for (int t = 0; t < k; ++t) sd(t) = std::max(sd(t), 1e-12f);
        Eigen::VectorXf inv_sd = sd.array().inverse();
        Eigen::MatrixXf GC = inv_sd.asDiagonal() * vbr * inv_sd.asDiagonal();
        
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> esGC(GC);
        Eigen::VectorXf evals = esGC.eigenvalues();
        Eigen::MatrixXf evecs = esGC.eigenvectors();
        
        int useF = std::min(NumXFA, k);
        Eigen::MatrixXf GC_red = Eigen::MatrixXf::Zero(k, k);
        for (int ii = 0; ii < useF; ++ii) {
          int idx = k - 1 - ii;
          float lam = evals(idx);
          Eigen::VectorXf v = evecs.col(idx);
          GC_red.noalias() += lam * (v * v.transpose());
        }
        for (int i = 0; i < k; ++i) GC_red(i, i) = 1.0f;
        
        for (int i = 0; i < k; ++i) {
          for (int j = 0; j < k; ++j) {
            vbr(i, j) = GC_red(i, j) * std::sqrt(std::max(vbr(i, i), 1e-12f) * std::max(vbr(j, j), 1e-12f));
          }
        }
        if (verbose) Rcpp::Rcout << "Rout: Applied XFA with NumXFA=" << useF << " for effect r=" << r << "\n";
      }
      
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(vbr);
      float min_ev = es.eigenvalues().minCoeff();
      if (min_ev < covMinEv ) {
        float need = std::abs(min_ev * covbend);
        bend_inflate[r] = std::max(bend_inflate[r], need);
      }

      if( k>=5 || MinDVb < covMinEv ){ vbr.diagonal().array() += bend_inflate[r]; }
      iG_list[r] = vbr.completeOrthogonalDecomposition().pseudoInverse();

    } // end effects loop
    
    // Update fixed effects
    for (int t = 0; t < k; ++t) {
      Eigen::VectorXf RHS_tmp = WX_list[t].transpose() * e.col(t);
      Eigen::VectorXf b_tmp   = iXX_list[t] * RHS_tmp;
      b.col(t).noalias() += b_tmp;
      Eigen::VectorXf delta = WX_list[t] * b_tmp;
      e.col(t) = (e.col(t) - delta).array() * W.col(t).array();
    }
    
    // Update residual variance
    Eigen::VectorXf new_ve = (e.cwiseProduct(y)).colwise().sum();
    new_ve = (new_ve.array() + Se.array()) * iNp.array();
    ve  = new_ve.array().max(1e-8f);
    iVe = ve.array().inverse();
    if (verbose) Rcpp::Rcout << "Rout: Updated residual variances ve\n";
    
    // Convergence
    double ss_max = 0.0;
    for (int r = 0; r < R; ++r) {
      double diff = (u0_list[r].array() - u_list[r].array()).square().sum();
      ss_max = std::max(ss_max, diff);
    }
    cnv = std::log10(static_cast<float>(std::max(ss_max, 1e-20)));
    ++numit;
    if (verbose) Rcpp::Rcout << "Rout: Iter " << numit << " \n log10(cnv)=" << cnv << " \n ss_max=" << ss_max << "\n";
    if (numit >= 5 && (cnv < logtol || std::isnan(cnv))) {
      if (verbose) Rcpp::Rcout << "Rout: Convergence reached or NaN detected, breaking loop\n";
      break;
    }
  } // end while
  
  if (verbose) Rcpp::Rcout << "Rout: Iterations finished. numit=" << numit << " cnv=" << cnv << "\n";
  
  // Fitted values
  Eigen::MatrixXf hat = X * b;
  for (int r = 0; r < R; ++r) {
    Eigen::MatrixXf Zr = Rcpp::as<Eigen::MatrixXf>(Z_list[r]);
    hat.noalias() += Zr * u_list[r];
  }
  if (verbose) Rcpp::Rcout << "Rout: Computed fitted values hat\n";
  
  // Prepare outputs
  Rcpp::List GC_out(R);
  Rcpp::List vb_out(R);
  Rcpp::List u_out(R);
  Rcpp::NumericVector bend_out(R);
  for (int r = 0; r < R; ++r) {
    Eigen::VectorXf sd = vb_list[r].diagonal().array().sqrt();
    for (int t = 0; t < k; ++t) sd(t) = std::max(sd(t), 1e-12f);
    Eigen::VectorXf inv_sd = sd.array().inverse();
    Eigen::MatrixXf GC = inv_sd.asDiagonal() * vb_list[r] * inv_sd.asDiagonal();
    vb_out[r]   = vb_list[r];
    GC_out[r]   = GC;
    u_out[r]    = u_list[r];
    bend_out[r] = bend_inflate[r];
  }
  if (verbose) Rcpp::Rcout << "Rout: Prepared outputs vb_list, GC_list, u, bend\n";
  
  // Return list
  return Rcpp::List::create(
    Rcpp::Named("TrZSZ")   = TrZSZ_list_R,
    Rcpp::Named("b")       = b,
    Rcpp::Named("u")       = u_out,
    Rcpp::Named("vb_list") = vb_out,
    Rcpp::Named("GC_list") = GC_out,
    Rcpp::Named("ve")      = ve,
    Rcpp::Named("hat")     = hat,
    Rcpp::Named("its")     = numit,
    Rcpp::Named("cnv")     = cnv,
    Rcpp::Named("bend")    = bend_out
  );
}

// [[Rcpp::export]]
SEXP PEGSZ(Eigen::MatrixXf Y, // matrix response variables
           Rcpp::List X_list, // LIST of design matrices of random effects
           int maxit = 100, // maximum number of iterations
           float logtol = -4.0, // convergence tolerance
           float covbend = 1.1, // covariance bending factor
           float covMinEv = 10e-6,  // minimum eigenvalue to bend
           int XFA = -1, // number of principal components to fit
           bool NNC = true){ // non-negative correlations
  
  // Get input dimensions
  int k = Y.cols(), n0 = Y.rows();
  int n_effects = X_list.size();
  
  std::vector<Eigen::MatrixXf> X_mats;
  Eigen::VectorXi p_vec(n_effects);
  for(int i=0; i<n_effects; ++i){
    X_mats.push_back(Rcpp::as<Eigen::MatrixXf>(X_list[i]));
    p_vec(i) = X_mats[i].cols();
  }
  
  // Incidence matrix Z
  Eigen::MatrixXf Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait and get inverses
  Eigen::VectorXf n = Z.colwise().sum();
  Eigen::VectorXf iN_mu = n.array().inverse(); // for mean calculation
  
  // Centralize y
  Eigen::VectorXf mu = Y.colwise().sum();
  mu = mu.array() * iN_mu.array();
  Eigen::MatrixXf y(n0,k);
  for(int i=0; i<k; i++){
    y.col(i) = (Y.col(i).array()-mu(i)).array()*Z.col(i).array();}
  
  // Pre-compute for each effect: Sum of squares of X and Tr(XSX)
  std::vector<Eigen::MatrixXf> XX_list;
  Eigen::MatrixXf TrXSX(n_effects, k);
  Eigen::MatrixXf MSx_mat(n_effects, k);
  
  for(int eff=0; eff<n_effects; ++eff){
    int p = p_vec(eff);
    Eigen::MatrixXf XX(p,k);
    for(int i=0; i<p; i++){
      XX.row(i) = X_mats[eff].col(i).array().square().matrix().transpose() * Z;
    }
    XX_list.push_back(XX);
    
    Eigen::MatrixXf XSX(p,k);
    for(int i=0; i<p; i++){
      XSX.row(i) = XX.row(i).transpose().array()*iN_mu.array() - 
        ((X_mats[eff].col(i).transpose()*Z).transpose().array()*iN_mu.array()).square();
    }
    MSx_mat.row(eff) = XSX.colwise().sum();
    TrXSX.row(eff) = n.transpose().array() * MSx_mat.row(eff).array();
  }
  
  // Variances
  Eigen::VectorXf iN_var = (n.array()-1).inverse(); // for variance calculation
  Eigen::VectorXf vy = y.colwise().squaredNorm(); vy = vy.array() * iN_var.array();
  Eigen::VectorXf ve = vy * 0.5;
  Eigen::VectorXf iVe = ve.array().inverse();
  
  std::vector<Eigen::MatrixXf> vb_list(n_effects, Eigen::MatrixXf(k,k));
  std::vector<Eigen::MatrixXf> iG_list(n_effects, Eigen::MatrixXf(k,k));
  for(int eff=0; eff<n_effects; ++eff){
    vb_list[eff] = (ve.array() / MSx_mat.row(eff).transpose().array()).matrix().asDiagonal();
    iG_list[eff] = vb_list[eff].completeOrthogonalDecomposition().pseudoInverse();
  }
  
  // Beta tilde for each effect
  std::vector<Eigen::MatrixXf> tilde_list;
  for(int eff=0; eff<n_effects; ++eff){
    tilde_list.push_back(X_mats[eff].transpose() * y);
  }
  
  // Initialize coefficient matrices and residuals
  Eigen::MatrixXf LHS(k,k);
  Eigen::VectorXf RHS(k);
  std::vector<Eigen::MatrixXf> b_list(n_effects);
  for(int eff=0; eff<n_effects; ++eff){
    b_list[eff] = Eigen::MatrixXf::Zero(p_vec(eff), k);
  }
  Eigen::VectorXf b0(k), b1(k);
  Eigen::MatrixXf e = y;
  
  // RGS
  std::random_device rd;
  std::mt19937 g(rd());
  
  // Convergence control
  std::vector<Eigen::MatrixXf> beta0_list(n_effects);
  float cnv = 10.0;
  Eigen::VectorXf inflate(n_effects); inflate.setZero();
  int numit = 0;
  
  // Bending & XFA objects
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> EVDofA(k);
  if(XFA<0) XFA = k;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigen_solver(k);
  
  // Main iterative loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0_list = b_list;
    
    // Randomized Gauss-Seidel loop for each effect
    for(int eff=0; eff<n_effects; ++eff){
      int p = p_vec(eff);
      std::vector<int> RGSvec(p);
      for(int j=0; j<p; j++){RGSvec[j]=j;}
      std::shuffle(RGSvec.begin(), RGSvec.end(), g);
      
      for(int j=0; j<p; j++){
        int J = RGSvec[j];
        // Update coefficient
        b0 = b_list[eff].row(J);
        LHS = iG_list[eff];
        LHS.diagonal() += (XX_list[eff].row(J).transpose().array() * iVe.array()).matrix();
        RHS = (X_mats[eff].col(J).transpose()*e).array() + XX_list[eff].row(J).array()*b0.transpose().array();
        RHS = RHS.array() *iVe.array();
        b1 = LHS.llt().solve(RHS);
        b_list[eff].row(J) = b1;
        // Update total residuals sequentially
        e = (e - (X_mats[eff].col(J)*(b1-b0).transpose()).cwiseProduct(Z)).matrix();
      }
    }
    
    // Update residual variance (using total residual)
    ve = (e.cwiseProduct(y)).colwise().sum();
    ve = ve.array() * iN_var.array();
    iVe = ve.array().inverse();
    
    // Update genetic variance, XFA, and Bending for each effect
    for(int eff=0; eff<n_effects; ++eff){
      Eigen::MatrixXf TildeHat = b_list[eff].transpose() * tilde_list[eff];
      Eigen::MatrixXf vb(k,k);
      for(int r=0; r<k; r++){
        for(int c=0; c<k; c++){
          if(r==c){ 
            if(TrXSX(eff,r) != 0) vb(r,c) = TildeHat(r,c)/TrXSX(eff,r); else vb(r,c) = 0;
          }else{
            if((TrXSX(eff,r)+TrXSX(eff,c)) != 0) vb(r,c) = (TildeHat(r,c)+TildeHat(c,r))/(TrXSX(eff,r)+TrXSX(eff,c)); else vb(r,c) = 0;
          }
        }
      }
      
      // XFA
      if(XFA == 0){
        Eigen::VectorXf sd_diag = vb.diagonal();
        vb.setZero();
        vb.diagonal() = sd_diag;
      }else if(XFA>0 && XFA < k){
        Eigen::VectorXf sd = vb.diagonal().array().sqrt();
        for (int t = 0; t < k; ++t) sd(t) = std::max(sd(t), 1e-12f);
        Eigen::VectorXf inv_sd = sd.array().inverse();
        Eigen::MatrixXf GC = inv_sd.asDiagonal() * vb * inv_sd.asDiagonal();
        
        eigen_solver.compute(GC);
        Eigen::MatrixXf V_reduced = eigen_solver.eigenvectors().rightCols(XFA);
        Eigen::VectorXf D_reduced_diag = eigen_solver.eigenvalues().tail(XFA);
        GC = V_reduced * D_reduced_diag.asDiagonal() * V_reduced.transpose();
        GC.diagonal().setOnes();
        
        vb = sd.asDiagonal() * GC * sd.asDiagonal();
      }
      
      // Bending
      if(NNC) vb = vb.array().cwiseMax(0.0).matrix();
      EVDofA.compute(vb); 
      float MinDVb = EVDofA.eigenvalues().minCoeff();
      if( MinDVb < covMinEv ){if(abs(MinDVb*covbend) > inflate(eff)) inflate(eff) = abs(MinDVb*covbend);}
      if( k>=5 || MinDVb < covMinEv ){ vb.diagonal().array() += inflate(eff); }      
      vb_list[eff] = vb;
      iG_list[eff] = vb.completeOrthogonalDecomposition().pseudoInverse();
    }
    // Update intercept
    b0 = e.colwise().sum();
    b0 = b0.array() * iN_var.array(); // This uses n-1, replicating original code's behavior
    for(int i=0; i<k; i++){ 
      mu(i) += b0(i);
      e.col(i) = (e.col(i).array()-b0(i)).array() * Z.col(i).array();
    }
    // Check for convergence
    cnv = 0.0;
    for(int eff=0; eff<n_effects; ++eff){
      cnv += (beta0_list[eff].array() - b_list[eff].array()).square().sum();
    }
    cnv = log10(cnv);
    ++numit;
    if( numit % 100 == 0){ Rcpp::Rcout << "Iter: "<< numit << " || Conv: "<< cnv << "\n"; } 
    if( cnv<logtol ){break;}
  }
  
  // Fit final predicted values
  Eigen::MatrixXf hat = Eigen::MatrixXf::Zero(n0,k);
  for(int eff=0; eff<n_effects; ++eff){
    hat += X_mats[eff] * b_list[eff];
  }
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  
  // Heritability (total, as per original formula)
  Eigen::VectorXf h2 = 1.0 - ve.array()/vy.array();
  
  // Prepare output lists for b and GC
  Rcpp::List b_out(n_effects);
  Rcpp::List GC_out(n_effects);
  for(int eff=0; eff<n_effects; ++eff){
    b_out[eff] = b_list[eff];
    
    Eigen::VectorXf sd = vb_list[eff].diagonal().array().sqrt();
    for (int t = 0; t < k; ++t) sd(t) = std::max(sd(t), 1e-12f);
    Eigen::VectorXf inv_sd = sd.array().inverse();
    GC_out[eff] = inv_sd.asDiagonal() * vb_list[eff] * inv_sd.asDiagonal();
  }
  
  // Output
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b_out,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC_out,
                            Rcpp::Named("bend")=inflate,
                            Rcpp::Named("numit")=numit,
                            Rcpp::Named("cnv")=cnv);

}


