
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <random>
#include <algorithm>  // for std::shuffle
#include <cmath>      // for std::isnan

using namespace Rcpp;
using namespace Eigen;

// Helper: M = X .* w (column-wise scaling by w)
static inline MatrixXf build_M(const VectorXf &w, const MatrixXf &X) {
  if (X.rows() != w.size()) {
    stop("build_M: length(w) must equal nrow(X).");
  }
  MatrixXf M(X.rows(), X.cols());
  M = X.array().colwise() * w.array();
  return M;
}

// Compute inverse (or pseudo-inverse) of t(M) %*% M once, outside the Z loop.
// Set `force_explicit_inverse = true` to mirror R's `solve(t(M) %*% M)` semantics.
static inline MatrixXf inv_tMtM(const MatrixXf &M, bool force_explicit_inverse = true) {
  MatrixXf Sym = M.transpose() * M; // p x p
  if (force_explicit_inverse) {
    LLT<MatrixXf> llt(Sym);
    if (llt.info() == Success) {
      MatrixXf I = MatrixXf::Identity(Sym.rows(), Sym.cols());
      return llt.solve(I);
    }
    CompleteOrthogonalDecomposition<MatrixXf> cod(Sym);
    return cod.pseudoInverse();
  } else {
    MatrixXf I = MatrixXf::Identity(Sym.rows(), Sym.cols());
    LLT<MatrixXf> llt(Sym);
    if (llt.info() == Success) {
      return llt.solve(I);
    }
    CompleteOrthogonalDecomposition<MatrixXf> cod(Sym);
    return cod.pseudoInverse();
  }
}

// --------------------------- TrZSZ (float) ---------------------------
float Compute_TrZSZ(const VectorXf &w,
                    const MatrixXf &X,
                    const MatrixXf &Z,
                    const bool force_explicit_inverse = true) {
  const Index n = X.rows();
  const Index q = Z.cols();
  if (w.size() != n) stop("Compute_TrZSZ: length(w) != nrow(X).");
  if (Z.rows() != n) stop("Compute_TrZSZ: nrow(Z) must equal nrow(X).");
  
  MatrixXf M  = build_M(w, X);         // n x p
  MatrixXf Mt = M.transpose();         // p x n
  MatrixXf InvMpM = inv_tMtM(M, force_explicit_inverse);  // p x p
  
  float TrZSZ = 0.0f;
  for (Index j = 0; j < q; ++j) {
    VectorXf z = Z.col(j).array() * w.array();            // n
    VectorXf beta = InvMpM * (Mt * z);                    // p
    VectorXf sz = (z - M * beta).array() * w.array();     // n
    TrZSZ += z.dot(sz);
  }
  return TrZSZ;
}

VectorXf Compute_All_TrZSZ(const MatrixXf &W,
                           const MatrixXf &X,
                           const MatrixXf &Z,
                           const bool force_explicit_inverse = true) {
  const Index n = X.rows();
  if (W.rows() != n) stop("Compute_All_TrZSZ: nrow(W) must equal nrow(X).");
  if (Z.rows() != n) stop("Compute_All_TrZSZ: nrow(Z) must equal nrow(X).");
  const Index K = W.cols();
  VectorXf out(K);
  for (Index k = 0; k < K; ++k) {
    out[k] = Compute_TrZSZ(W.col(k), X, Z, force_explicit_inverse);
  }
  return out;
}

List Compute_TrZSZ_for_All_Zs(const MatrixXf &W,
                              const MatrixXf &X,
                              const List &list_of_Z,
                              const bool force_explicit_inverse = true) {
  const Index n = X.rows();
  if (W.rows() != n) stop("Compute_TrZSZ_for_All_Zs: nrow(W) must equal nrow(X).");
  const int R = list_of_Z.size();
  List out(R);
  for (int r = 0; r < R; ++r) {
    MatrixXf Zs = as<MatrixXf>(list_of_Z[r]);
    if (Zs.rows() != n) {
      stop("Compute_TrZSZ_for_All_Zs: each Z must have n rows equal to nrow(X).");
    }
    out[r] = Compute_All_TrZSZ(W, X, Zs, force_explicit_inverse);
  }
  return out;
}

// --------------------------- ZpZ (float) ---------------------------
MatrixXf Compute_ZpZ(const MatrixXf &W,
                     const MatrixXf &Z) {
  const Index n = Z.rows();
  const Index K = W.cols();
  const Index P = Z.cols();
  if (W.rows() != n) stop("Compute_ZpZ: nrow(W) must equal nrow(Z).");
  
  MatrixXf ZpZ(P, K);  // rows=P, cols=K (matches R)
  for (Index k = 0; k < K; ++k) {
    for (Index p = 0; p < P; ++p) {
      VectorXf wz = Z.col(p).array() * W.col(k).array();
      ZpZ(p, k) = wz.squaredNorm();
    }
  }
  return ZpZ;
}

List Compute_ZpZ_for_All_Zs(const MatrixXf &W,
                            const List &list_of_Z) {
  const Index n = W.rows();
  const int R = list_of_Z.size();
  List out(R);
  for (int r = 0; r < R; ++r) {
    MatrixXf Zs = as<MatrixXf>(list_of_Z[r]);
    if (Zs.rows() != n) stop("Compute_ZpZ_for_All_Zs: each Z must have n rows.");
    out[r] = Compute_ZpZ(W, Zs);
  }
  return out;
}

// --------------------------- beta_tilde (float) ---------------------------
List Compute_beta_tilde_for_All_Zs(const MatrixXf &SY,
                                   const List &list_of_Z) {
  const Index n = SY.rows();
  const int R = list_of_Z.size();
  List out(R);
  for (int r = 0; r < R; ++r) {
    MatrixXf Zs = as<MatrixXf>(list_of_Z[r]);  // n x P_r
    if (Zs.rows() != n) stop("Compute_beta_tilde_for_All_Zs: each Z must have n rows.");
    MatrixXf beta_tilde = Zs.transpose() * SY; // (P_r x K)
    out[r] = beta_tilde;
  }
  return out;
}

// ======================= MLMX (float) =======================
// [[Rcpp::export]]
SEXP PEGSX(MatrixXf Y,     // n x k responses
          MatrixXf X,     // n x f fixed-effects design
          List Z_list,    // list of random-effects designs (n x p_r)
          int maxit = 500,
          float logtol = -8.0,
          int cores = 1,
          bool verbose = false,
          float df0 = 1.1,
          bool NonNegativeCorr = false,
          bool InnerGS = false,
          bool NoInv = false,
          bool XFA = false,
          int NumXFA = 3)
{
  if (cores != 1) Eigen::setNbThreads(cores);
  
  // Dimensions
  const int n = Y.rows();
  const int k = Y.cols();
  const int f = X.cols();
  const int R = Z_list.size();
  if (verbose) Rcout << "Rout: Start PEGSX\n n=" << n << " k=" << k << " f=" << f << " R=" << R
                     << " \n InnerGS=" << InnerGS
                     << " \n NoInv=" << NoInv
                     << " \n XFA=" << XFA
                     << " \n NumXFA=" << NumXFA
                     << " \n NonNegativeCorr=" << NonNegativeCorr << "\n";
  
  // Build mask W and zero-out NA in Y
  MatrixXf W(n, k);
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
  VectorXf n_each = W.colwise().sum();
  if (verbose) Rcout << "Rout: Built mask W and counts per trait\n";
  
  // Masked X per trait (WX) and iXX = (WX'WX)^+
  std::vector<MatrixXf> WX_list(k);
  std::vector<MatrixXf> iXX_list(k);
  MatrixXf b = MatrixXf::Zero(f, k); // fixed-effect coefficients
  for (int t = 0; t < k; ++t) {
    WX_list[t].resize(n, f);
    for (int j = 0; j < f; ++j) WX_list[t].col(j) = X.col(j).array() * W.col(t).array();
    MatrixXf XX = WX_list[t].transpose() * WX_list[t];
    iXX_list[t] = XX.completeOrthogonalDecomposition().pseudoInverse();
    VectorXf RHS = WX_list[t].transpose() * Y.col(t);
    b.col(t).noalias() = iXX_list[t] * RHS;
  }
  if (verbose) Rcout << "Rout: Estimated fixed effects b\n";
  
  // Residuals: y = (Y - WX*b) masked by W
  MatrixXf y(n, k);
  for (int t = 0; t < k; ++t){
    y.col(t) = (Y.col(t) - WX_list[t] * b.col(t)).array() * W.col(t).array();
  }
  if (verbose) Rcout << "Rout: Computed masked residuals y\n";
  
  // Precompute p_r and RGS indices
  std::vector<int> pR(R, 0);
  for (int r = 0; r < R; ++r) {
    MatrixXf Zr_once = as<MatrixXf>(Z_list[r]);
    pR[r] = Zr_once.cols();
    if (verbose) Rcout << "Rout: Z[" << r << "] cols=" << pR[r] << "\n";
  }
  std::vector<std::vector<int>> RGS_index(R);
  for (int r = 0; r < R; ++r) {
    RGS_index[r].resize(pR[r]);
    for (int j = 0; j < pR[r]; ++j) RGS_index[r][j] = j;
  }
  
  // === Precomputations via helpers (no duplication) ===
  List ZZ_list_R    = Compute_ZpZ_for_All_Zs(W, Z_list);        // p_r x k
  List TrZSZ_list_R = Compute_TrZSZ_for_All_Zs(W, X, Z_list);   // k-vector per effect
  List tilde_list_R = Compute_beta_tilde_for_All_Zs(y, Z_list); // p_r x k
  if (verbose) Rcout << "Rout: Precomputed ZZ/TrZSZ/tilde\n";
  
  // Initialize random coefficients u_r
  std::vector<MatrixXf> u_list(R);
  for (int r = 0; r < R; ++r) u_list[r] = MatrixXf::Zero(pR[r], k);
  
  // Residual variances ve initialization
  VectorXf ssy = y.colwise().squaredNorm();
  VectorXf denom = (n_each.array() - f).matrix();
  for (int t = 0; t < k; ++t) denom(t) = std::max(denom(t), 1.0f);
  VectorXf iN_mlm = denom.array().inverse().matrix();
  VectorXf ve = (ssy.array() * iN_mlm.array()).matrix() * 0.5f;
  ve = ve.array().max(1e-8f);
  VectorXf iVe = ve.array().inverse().matrix();
  if (verbose) Rcout << "Rout: Initialized ve\n";
  
  // vb and iG per effect (MLM-like init)
  std::vector<MatrixXf> vb_list(R);
  std::vector<MatrixXf> iG_list(R);
  std::vector<float> bend_inflate(R, 0.0f);
  for (int r = 0; r < R; ++r) {
    vb_list[r].resize(k, k);
    vb_list[r].setZero();
    
    VectorXf Tr_r = as<VectorXf>(TrZSZ_list_R[r]);
    
    for (int t = 0; t < k; ++t) {
      float denom_r = std::max(Tr_r(t) * iN_mlm(t), 1e-8f);
      vb_list[r](t, t) = ve(t) / denom_r;
    }
    iG_list[r] = vb_list[r].completeOrthogonalDecomposition().pseudoInverse();
  }
  if (verbose) Rcout << "Rout: Initialized vb and iG for all effects\n";
  
  // Priors
  std::vector<MatrixXf> Sb_list(R);
  for (int r = 0; r < R; ++r) Sb_list[r] = vb_list[r] * df0;
  VectorXf Se = ve * df0;
  VectorXf iNp = (n_each.array() + df0 - f).matrix();
  for (int t = 0; t < k; ++t) iNp(t) = 1.0f / std::max(iNp(t), 1.0f);
  if (verbose) Rcout << "Rout: Set priors\n";
  
  // Working residuals (start with y)
  MatrixXf e = y;
  
  // Convergence
  std::vector<MatrixXf> u0_list(R);
  float cnv = 10.0f;
  int numit = 0;
  std::random_device rd_dev;
  std::mt19937 gen(rd_dev());
  if (verbose) Rcout << "Rout: Starting GS iterations\n";
  
  // Iteration loop
  while (numit < maxit) {
    // Store previous u
    for (int r = 0; r < R; ++r) u0_list[r] = u_list[r];
    
    // -------------------------
    // Randomized Gauss–Seidel across effects
    // -------------------------
    for (int r = 0; r < R; ++r) {
      std::shuffle(RGS_index[r].begin(), RGS_index[r].end(), gen);
      
      MatrixXf Zr        = as<MatrixXf>(Z_list[r]);         // n x p_r
      MatrixXf ZZr       = as<MatrixXf>(ZZ_list_R[r]);      // p_r x k
      MatrixXf tilde_r   = as<MatrixXf>(tilde_list_R[r]);   // p_r x k
      VectorXf Tr_r      = as<VectorXf>(TrZSZ_list_R[r]);   // k
      
      MatrixXf &ur  = u_list[r];
      MatrixXf &iGr = iG_list[r];
      MatrixXf &vbr = vb_list[r];
      
      const int p_r = Zr.cols();
      for (int jj = 0; jj < p_r; ++jj) {
        int J = RGS_index[r][jj];
        VectorXf u0 = ur.row(J).transpose();
        
        // diag terms from ZZ(J,:) .* iVe
        VectorXf diagVec = (ZZr.row(J).transpose().array() * iVe.array()).matrix();
        
        MatrixXf LHS(k, k);
        VectorXf RHS(k);
        if (NoInv) {
          LHS.noalias() = vbr * diagVec.asDiagonal();
          for (int i = 0; i < k; ++i) LHS(i, i) += 1.0f;
          
          VectorXf base = (Zr.col(J).transpose() * e).transpose();
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
        
        VectorXf u1(k);
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
        
        // Update residuals with mask
        e = (e - (Zr.col(J) * (u1 - u0).transpose()).cwiseProduct(W)).eval();
      }
      
      // ---- vb update for effect r ----
      MatrixXf TildeHat(k, k);
      TildeHat.noalias() = ur.transpose() * tilde_r;
      
      for (int i = 0; i < k; ++i) {
        float denom_i = std::max(Tr_r(i), 1e-8f);
        vbr(i,i) = (TildeHat(i,i) + Sb_list[r](i,i)) / (denom_i + df0);
      }
      for (int i = 0; i < k; ++i)
        for (int j = i + 1; j < k; ++j) {
          float denom_ij = std::max(Tr_r(i) + Tr_r(j), 1e-8f);
          float cov_ij   = (TildeHat(i,j) + TildeHat(j,i)) / denom_ij;
          vbr(i,j) = cov_ij;
          vbr(j,i) = cov_ij;
        }
        
        if (NonNegativeCorr) {
          vbr = vbr.array().cwiseMax(0.0f).matrix();
          vbr = 0.5f * (vbr + vbr.transpose());
          if (verbose) Rcout << "Rout: Applied NonNegativeCorr to vb for effect r=" << r << "\n";
        }
        
        if (XFA && NumXFA > 0) {
          VectorXf sd = vbr.diagonal().array().sqrt();
          for (int t = 0; t < k; ++t) sd(t) = std::max(sd(t), 1e-12f);
          VectorXf inv_sd = sd.array().inverse();
          MatrixXf GC = inv_sd.asDiagonal() * vbr * inv_sd.asDiagonal();
          
          SelfAdjointEigenSolver<MatrixXf> esGC(GC);
          VectorXf evals = esGC.eigenvalues();
          MatrixXf evecs = esGC.eigenvectors();
          
          int useF = std::min(NumXFA, k);
          MatrixXf GC_red = MatrixXf::Zero(k, k);
          for (int ii = 0; ii < useF; ++ii) {
            int idx = k - 1 - ii;
            float lam = evals(idx);
            VectorXf v = evecs.col(idx);
            GC_red.noalias() += lam * (v * v.transpose());
          }
          for (int i = 0; i < k; ++i) GC_red(i, i) = 1.0f;
          
          for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k; ++j) {
              vbr(i, j) = GC_red(i, j) * std::sqrt(std::max(vbr(i, i), 1e-12f) * std::max(vbr(j, j), 1e-12f));
            }
          }
          if (verbose) Rcout << "Rout: Applied XFA with NumXFA=" << useF << " for effect r=" << r << "\n";
        }
        
        SelfAdjointEigenSolver<MatrixXf> es(vbr);
        float min_ev = es.eigenvalues().minCoeff();
        if (min_ev < 0.001f) {
          float need = std::abs(min_ev * 1.1f);
          bend_inflate[r] = std::max(bend_inflate[r], need);
        }
        vbr.diagonal().array() += bend_inflate[r];
        
        iG_list[r] = vbr.completeOrthogonalDecomposition().pseudoInverse();
    } // end effects loop
    
    // ------------------------------------------------------
    // Update fixed effects (MLM-style incremental correction)
    // ------------------------------------------------------
    for (int t = 0; t < k; ++t) {
      // RHS_tmp = WX' * e_t
      VectorXf RHS_tmp = WX_list[t].transpose() * e.col(t);           // f
      // b_tmp = (WX'WX)^+ * RHS_tmp
      VectorXf b_tmp   = iXX_list[t] * RHS_tmp;                       // f
      // accumulate into b(:,t)
      b.col(t).noalias() += b_tmp;                                    // f
      // update residuals e_t := e_t - WX * b_tmp, then re-mask
      VectorXf delta = WX_list[t] * b_tmp;                            // n
      e.col(t) = (e.col(t) - delta).array() * W.col(t).array();       // n
    }
    
    // Update residual variance ve with prior: ve = (sum(e*y) + Se) * iNp
    VectorXf new_ve = (e.cwiseProduct(y)).colwise().sum();
    new_ve = (new_ve.array() + Se.array()) * iNp.array();
    ve  = new_ve.array().max(1e-8f);
    iVe = ve.array().inverse();
    if (verbose) Rcout << "Rout: Updated residual variances ve\n";
    
    // Convergence on u changes (max across effects)
    double ss_max = 0.0;
    for (int r = 0; r < R; ++r) {
      double diff = (u0_list[r].array() - u_list[r].array()).square().sum();
      ss_max = std::max(ss_max, diff);
    }
    cnv = std::log10(static_cast<float>(std::max(ss_max, 1e-20)));
    ++numit;
    if (verbose) Rcout << "Rout: Iter " << numit
                       << " \n log10(cnv)=" << cnv
                       << " \n ss_max=" << ss_max << "\n";
    if (numit >= 5 && (cnv < logtol || std::isnan(cnv))) {
      if (verbose) Rcout << "Rout: Convergence reached or NaN detected, breaking loop\n";
      break;
    }
  } // end while
  
  if (verbose) Rcout << "Rout: Iterations finished. numit=" << numit << " cnv=" << cnv << "\n";
  
  // Fitted values: hat = X*b + sum_r Z_r * u_r (no duplication of Z)
  MatrixXf hat = X * b;
  for (int r = 0; r < R; ++r) {
    MatrixXf Zr = as<MatrixXf>(Z_list[r]);
    hat.noalias() += Zr * u_list[r];
  }
  if (verbose) Rcout << "Rout: Computed fitted values hat\n";
  
  // Genetic correlations per effect (standardize vb)
  List GC_out(R);
  List vb_out(R);
  List u_out(R);
  NumericVector bend_out(R);
  for (int r = 0; r < R; ++r) {
    VectorXf sd = vb_list[r].diagonal().array().sqrt();
    for (int t = 0; t < k; ++t) sd(t) = std::max(sd(t), 1e-12f);
    VectorXf inv_sd = sd.array().inverse();
    MatrixXf GC = inv_sd.asDiagonal() * vb_list[r] * inv_sd.asDiagonal();
    vb_out[r]   = vb_list[r];
    GC_out[r]   = GC;
    u_out[r]    = u_list[r];
    bend_out[r] = bend_inflate[r];
  }
  if (verbose) Rcout << "Rout: Prepared outputs vb_list, GC_list, u, bend\n";
  
  // Return TrZSZ as the original List, and other outputs
  return List::create(
    Named("TrZSZ")   = TrZSZ_list_R,   // list of k-vectors (float), returned directly
    Named("b")       = b,              // f x k fixed-effect coefficients
    Named("u")       = u_out,          // list of p_r x k random-effect coefficients
    Named("vb_list") = vb_out,         // list of k x k covariances per random effect
    Named("GC_list") = GC_out,         // list of k x k genetic correlations per random effect
    Named("ve")      = ve,             // k residual variances
    Named("hat")     = hat,            // n x k fitted values
    Named("its")     = numit,          // iterations
    Named("cnv")     = cnv,            // log10 convergence score (max across effects)
    Named("bend")    = bend_out        // per-effect bending magnitude
  );
}
