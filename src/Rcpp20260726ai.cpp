// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <random>
using namespace Rcpp;

static float fvar(const Eigen::VectorXf& x){
  return (x.array()-x.mean()).square().sum()/(float)(x.size()-1);
}

// [[Rcpp::export]]
SEXP KMUP(Eigen::MatrixXf X, Eigen::VectorXf b, Eigen::VectorXf d, Eigen::VectorXf xx, Eigen::VectorXf e, Eigen::VectorXf L, float Ve, float pi){
  int p = X.cols();
  Eigen::VectorXf e1 = e;
  Eigen::VectorXf e2 = e;
  float b0,b1,b2,cj,dj,pj;
  float C = -0.5f/std::sqrt(Ve);
  for(int j=0; j<p; j++){
    b0 = b[j];
    b1 = R::rnorm( (X.col(j).dot(e)+xx[j]*b0)/(xx[j]+L[j]), std::sqrt(Ve/(xx[j]+L[j])) );
    b2 = R::rnorm( 0, std::sqrt(Ve/(xx[j]+L[j])) );
    e1 = e - X.col(j)*(b1-b0);
    if(pi>0){
      e2 = e - X.col(j)*(b2-b0);
      cj = (1-pi)*std::exp(C*e1.squaredNorm());
      dj = (pi)*std::exp(C*e2.squaredNorm());
      pj = cj/(cj+dj);
      if(R::rbinom(1,pj)==1){
        b[j] = b1; d[j] = 1; e = e1;
      }else{
        b[j] = b2; d[j] = 0; e = e2;
      }
    }else{
      d[j] = 1; b[j] = b1; e = e1;
    }
  }
  return List::create(Named("b") = b, Named("d") = d, Named("e") = e);
}

// [[Rcpp::export]]
SEXP KMUP2(Eigen::MatrixXf X, Eigen::VectorXf Use, Eigen::VectorXf b, Eigen::VectorXf d, Eigen::VectorXf xx, Eigen::VectorXf E, Eigen::VectorXf L, float Ve, float pi){
  int p = X.cols();
  int n0 = X.rows();
  int n = Use.size();
  float b0,b1,b2,cj,dj,pj;
  float C = -0.5f/std::sqrt(Ve);
  float bg = (float)n0/(float)n;
  Eigen::VectorXf e0(n), H(n);
  for(int k=0; k<n; k++){
    e0[k] = E[(int)Use[k]];
  }
  Eigen::VectorXf e1 = e0;
  Eigen::VectorXf e2 = e0;
  for(int j=0; j<p; j++){
    for(int x=0; x<n; x++){
      H[x] = X((int)Use[x],j);
    }
    b0 = b[j];
    b1 = R::rnorm((H.dot(e0)+b0)/(xx(j)*bg+L(j)), std::sqrt(Ve/(xx(j)*bg+L(j))));
    b2 = R::rnorm(0, std::sqrt(Ve/(xx(j)*bg+L(j))));
    e1 = e0 - H*(b1-b0);
    if(pi>0){
      e2 = e0 - H*(b2-b0);
      cj = (1-pi)*std::exp(C*e1.squaredNorm());
      dj = (pi)*std::exp(C*e2.squaredNorm());
      pj = cj/(cj+dj);
      if(R::rbinom(1,pj)==1){
        b[j] = b1; d[j] = 1; e0 = e1;
      }else{
        b[j] = b2; d[j] = 0; e0 = e2;
      }
    }else{
      d[j] = 1; b[j] = b1; e0 = e1;
    }
  }
  return List::create(Named("b") = b, Named("d") = d, Named("e") = e0);
}

// [[Rcpp::export]]
SEXP emBA(Eigen::VectorXf y, Eigen::MatrixXf gen, float df = 10, float R2 = 0.5){
  int it = 200;
  int p = gen.cols();
  int n = gen.rows();
  float ve = 1;
  Eigen::VectorXf d = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf vb = Eigen::VectorXf::Ones(p);
  Eigen::VectorXf Lmb = ve * vb.cwiseInverse();
  float vy = fvar(y);
  Eigen::VectorXf xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = gen.col(i).squaredNorm();
    vx[i] = fvar(gen.col(i));
  }
  float MSx = vx.sum();
  float Sb = R2*(df+2)*vy/MSx;
  float Se = (1-R2)*(df+2)*vy;
  float mu = y.mean();
  Eigen::VectorXf e = y.array()-mu;
  float b0,b1,eM,h2;
  std::vector<int> order(p); for(int j=0;j<p;j++) order[j]=j;
  for(int i=0; i<it; i++){
    std::shuffle(order.begin(), order.end(), std::mt19937(i));
    for(int jj=0; jj<p; jj++){
      int j = order[jj];
      b0 = b[j];
      b1 = (gen.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb[j]);
      e -= gen.col(j)*(b1-b0);
      b[j] = b1;
      vb[j] = (Sb+b[j]*b[j])/(df+1);
      e -= gen.col(j)*(b1-b0);
    }
    ve = (e.squaredNorm()+Se)/(n+df);
    Lmb = ve * vb.cwiseInverse();
    eM = e.mean();
    mu += eM;
    e = e.array()-eM;
  }
  h2 = 1-ve/vy;
  Eigen::VectorXf fit = gen*b;
  fit = fit.array()+mu;
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("Vb") = vb,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emBB(Eigen::VectorXf y, Eigen::MatrixXf gen, float df = 10, float R2 = 0.5, float Pi = 0.75){
  int it = 200;
  int p = gen.cols();
  int n = gen.rows();
  float ve = 1;
  Eigen::VectorXf d = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf vb = Eigen::VectorXf::Ones(p);
  Eigen::VectorXf Lmb = ve * vb.cwiseInverse();
  float vy = fvar(y);
  if(Pi>0.5){Pi = 1-Pi;}
  Eigen::VectorXf xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = gen.col(i).squaredNorm();
    vx[i] = fvar(gen.col(i));
  }
  float MSx = vx.sum()*Pi;
  float Sb = R2*(df+2)*vy/MSx;
  float Se = (1-R2)*(df+2)*vy;
  float mu = y.mean();
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf e1(n), e2(n);
  float b0,b1,LR,eM,h2,C;
  float Pi0 = (1-Pi)/Pi;
  std::vector<int> order(p); for(int j=0;j<p;j++) order[j]=j;
  for(int i=0; i<it; i++){
    C = -0.5f/std::sqrt(ve);
    std::shuffle(order.begin(), order.end(), std::mt19937(i));
    for(int jj=0; jj<p; jj++){
      int j = order[jj];
      b0 = b[j];
      b1 = (gen.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb[j]);
      e1 = e - gen.col(j)*(b1-b0);
      e2 = e - gen.col(j)*(0-b0);
      LR = Pi0*std::exp(C*(e2.squaredNorm()-e1.squaredNorm()));
      d[j] = (1/(1+LR));
      b[j] = b1*d[j];
      vb[j] = (Sb+b[j]*b[j])/(df+1);
      e -= gen.col(j)*(b[j]-b0);
    }
    ve = (e.squaredNorm()+Se)/(n+df);
    Lmb = ve * vb.cwiseInverse();
    eM = e.mean();
    mu += eM;
    e = e.array()-eM;
  }
  h2 = 1-ve/vy;
  Eigen::VectorXf fit = gen*b;
  fit = fit.array()+mu;
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("d") = d,
                      Named("hat") = fit,
                      Named("Vb") = vb,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emBC(Eigen::VectorXf y, Eigen::MatrixXf gen, float df = 10, float R2 = 0.5, float Pi = 0.75){
  int it = 200;
  int p = gen.cols();
  int n = gen.rows();
  Eigen::VectorXf d = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  float vy = fvar(y);
  if(Pi>0.5){Pi = 1-Pi;}
  Eigen::VectorXf xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = gen.col(i).squaredNorm();
    vx[i] = fvar(gen.col(i));
  }
  float MSx = vx.sum()*Pi*(1-Pi);
  float Sa = R2*(df+2)*vy/MSx;
  float Se = (1-R2)*(df+2)*vy;
  float mu = y.mean();
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf e1(n), e2(n);
  float ve = Sa;
  float va = Se;
  float Lmb = ve/va;
  float b0,b1,LR,eM,h2,C;
  float Pi0 = (1-Pi)/Pi;
  std::vector<int> order(p); for(int j=0;j<p;j++) order[j]=j;
  for(int i=0; i<it; i++){
    C = -0.5f/std::sqrt(ve);
    std::shuffle(order.begin(), order.end(), std::mt19937(i));
    for(int jj=0; jj<p; jj++){
      int j = order[jj];
      b0 = b[j];
      b1 = (gen.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb);
      e1 = e - gen.col(j)*(b1-b0);
      e2 = e - gen.col(j)*(0-b0);
      LR = Pi0*std::exp(C*(e2.squaredNorm()-e1.squaredNorm()));
      d[j] = (1/(1+LR));
      b[j] = b1*d[j];
      e -= gen.col(j)*(b[j]-b0);
    }
    ve = (e.squaredNorm()+Se)/(n+df);
    va = (b.squaredNorm()+Sa)/(p+df)/(d.mean()-Pi);
    Lmb = ve/va;
    eM = e.mean();
    mu += eM;
    e = e.array()-eM;
  }
  h2 = 1-ve/vy;
  Eigen::VectorXf fit = gen*b;
  fit = fit.array()+mu;
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("d") = d,
                      Named("hat") = fit,
                      Named("Vg") = va*MSx,
                      Named("Va") = va,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emDE(Eigen::VectorXf y, Eigen::MatrixXf gen, float R2 = 0.5){
  int maxit = 300;
  float tol = 10e-6f;
  int p = gen.cols();
  int n = gen.rows();
  float b0,eM;
  float mu = y.mean();
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf xx(p);
  for(int k=0; k<p; k++){
    xx[k] = gen.col(k).squaredNorm();
    if(xx[k]==0) xx[k]=0.1f;
  }
  Eigen::VectorXf vx(p);
  for(int k=0; k<p; k++){vx[k] = fvar(gen.col(k));}
  float cxx = vx.sum()*(1-R2)/R2;
  float Ve;
  Eigen::VectorXf Vb(p);
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf Lmb = Eigen::VectorXf::Constant(p, (float)p+cxx);
  float b1;
  Eigen::VectorXf bc(p);
  int numit = 0;
  float cnv = 1;
  std::vector<int> order(p); for(int j=0;j<p;j++) order[j]=j;
  while(numit<maxit){
    bc = b;
    std::shuffle(order.begin(), order.end(), std::mt19937(numit));
    for(int jj=0; jj<p; jj++){
      int j = order[jj];
      b0 = b[j];
      b1 = (gen.col(j).dot(e)+xx[j]*b0)/(Lmb(j)+xx(j));
      b[j] = b1;
      e -= gen.col(j)*(b1-b0);
    }
    eM = e.mean();
    mu += eM;
    e = e.array()-eM;
    Ve = e.dot(y)/(n-1);
    Vb = (b.array()*b.array() + Ve/(xx.array()+Lmb.array()+0.0001f)).matrix();
    for(int j=0; j<p; j++){
      Lmb[j] = std::sqrt(cxx*Ve/Vb[j]);
    }
    ++numit;
    cnv = (bc-b).cwiseAbs().sum();
    if(cnv<tol){break;}
  }
  Eigen::VectorXf fit = gen*b;
  fit = fit.array()+mu;
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("Vb")=Vb,
                      Named("Ve")=Ve,
                      Named("h2")=Vb.sum()/(Vb.sum()+Ve));
}

// [[Rcpp::export]]
SEXP emRR(Eigen::VectorXf y, Eigen::MatrixXf gen, float df = 10, float R2 = 0.5){
  int it = 200;
  int p = gen.cols();
  int n = gen.rows();
  Eigen::VectorXf xx(p), vx(p);
  for(int k=0; k<p; k++){
    xx[k] = gen.col(k).squaredNorm();
    vx[k] = fvar(gen.col(k));
  }
  float MSx = vx.sum();
  float Lmb = MSx;
  float Rho = MSx*(1-R2)/R2;
  float vy = fvar(y);
  float ve = 0.5f*vy;
  float vb = ve/MSx;
  float Se = (1-R2)*(df+2)*vy;
  float Sb = R2*(df+2)*vy/MSx;
  float mu = y.mean();
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf e = y.array()-mu;
  float b0,eM,h2;
  std::vector<int> order(p); for(int j=0;j<p;j++) order[j]=j;
  for(int i=0; i<it; i++){
    std::shuffle(order.begin(), order.end(), std::mt19937(i));
    for(int jj=0; jj<p; jj++){
      int j = order[jj];
      b0 = b[j];
      b[j] = (gen.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb);
      e -= gen.col(j)*(b[j]-b0);
    }
    vb = (b.squaredNorm()+Sb)/(p+df);
    ve = (e.squaredNorm()+Se)/(n+df);
    Lmb = std::sqrt(Rho*ve/vb);
    eM = e.mean();
    mu += eM;
    e = e.array()-eM;
  }
  h2 = 1-ve/vy;
  Eigen::VectorXf fit = gen*b;
  fit = fit.array()+mu;
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("Va") = vb,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emBL(Eigen::VectorXf y, Eigen::MatrixXf gen, float R2 = 0.5, float alpha = 0.02){
  int it = 200;
  float h2 = R2;
  int p = gen.cols();
  int n = gen.rows();
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  float b0,eM,Half_L2;
  float mu = y.mean();
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf xx(p);
  for(int k=0; k<p; k++){xx[k] = gen.col(k).squaredNorm();}
  float cxx = xx.mean();
  float Lmb1 = cxx*((1-h2)/h2)*alpha*0.5f;
  float Lmb2 = cxx*((1-h2)/h2)*(1-alpha);
  float OLS, G;
  std::vector<int> order(p); for(int j=0;j<p;j++) order[j]=j;
  for(int i=0; i<it; i++){
    std::shuffle(order.begin(), order.end(), std::mt19937(i));
    for(int jj=0; jj<p; jj++){
      int j = order[jj];
      b0 = b[j];
      OLS = (gen.col(j).dot(e)+xx[j]*b0);
      Half_L2 = 0.5f*OLS/(xx[j]+cxx);
      if(OLS>0){
        G = 0.5f*(OLS-Lmb1)/(Lmb2+xx(j));
        if(G>0){b[j] = G+Half_L2;}else{b[j] = Half_L2;}
      }else{
        G = 0.5f*(OLS+Lmb1)/(Lmb2+xx(j));
        if(G<0){b[j] = G+Half_L2;}else{b[j] = Half_L2;}
      }
      e -= gen.col(j)*(b[j]-b0);
    }
    eM = e.mean();
    mu += eM;
    e = e.array()-eM;
  }
  Eigen::VectorXf fit = gen*b;
  fit = fit.array()+mu;
  h2 = 1 - fvar(e)/fvar(y);
  return List::create(Named("mu") = mu, Named("b") = b, Named("hat") = fit, Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emEN(Eigen::VectorXf y, Eigen::MatrixXf gen, float R2 = 0.5, float alpha = 0.02){
  int maxit = 300;
  float tol = 10e-11f;
  int p = gen.cols();
  int n = gen.rows();
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  float b0,eM;
  float mu = y.mean();
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf xx(p), vx(p);
  for(int k=0; k<p; k++){xx[k] = gen.col(k).squaredNorm();}
  for(int k=0; k<p; k++){vx[k] = fvar(gen.col(k));}
  float cxx = vx.sum()*(1-R2)/R2;
  float Ve, Va;
  float Sy = std::sqrt(fvar(y));
  float Lmb = cxx;
  float Lmb1 = 0.5f*Lmb*alpha*Sy;
  float Lmb2 = Lmb*(1-alpha);
  float trAC22 = 0;
  for(int k=0; k<p; k++) trAC22 += 1.0f/(xx[k]+Lmb);
  float OLS, b1;
  Eigen::VectorXf bc(p);
  int numit = 0;
  float cnv = 1;
  std::vector<int> order(p); for(int j=0;j<p;j++) order[j]=j;
  while(numit<maxit){
    bc = b;
    std::shuffle(order.begin(), order.end(), std::mt19937(numit));
    for(int jj=0; jj<p; jj++){
      int j = order[jj];
      b0 = b[j];
      OLS = (gen.col(j).dot(e)+xx[j]*b0);
      if(OLS>0){
        b1 = (OLS-Lmb1)/(Lmb2+xx(j)); if(b1<0){b1=0;}
      }else{
        b1 = (OLS+Lmb1)/(Lmb2+xx(j)); if(b1>0){b1=0;}
      }
      b[j] = b1;
      e -= gen.col(j)*(b1-b0);
    }
    eM = e.mean();
    mu += eM;
    e = e.array()-eM;
    Ve = e.dot(y)/(n-1);
    Va = (b.squaredNorm()+trAC22*Ve)/p;
    Lmb = Ve/Va;
    Lmb1 = 0.5f*Lmb*alpha*Sy;
    Lmb2 = Lmb*(1-alpha);
    ++numit;
    cnv = (bc-b).cwiseAbs().sum();
    if(cnv<tol){break;}
  }
  Eigen::VectorXf fit = gen*b;
  fit = fit.array()+mu;
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("Va")=Va*cxx,
                      Named("Ve")=Ve,
                      Named("h2")=Va*cxx/(Va*cxx+Ve));
}

// [[Rcpp::export]]
SEXP emML(Eigen::VectorXf y, Eigen::MatrixXf gen,
          Rcpp::Nullable<Rcpp::NumericVector> D = R_NilValue){
  int maxit = 300;
  float tol = 10e-8f;
  int p = gen.cols();
  int n = gen.rows();
  bool P_WEIGHTS = false;
  Eigen::VectorXf d = Eigen::VectorXf::Ones(p);
  if(D.isNotNull()){
    P_WEIGHTS = true;
    Rcpp::NumericVector Dtmp(D);
    for(int i=0;i<p;i++) d[i]=(float)Dtmp[i];
  }
  float b0, eM, ve, vb, h2, mu = y.mean();
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = gen.col(i).squaredNorm();
    vx[i] = fvar(gen.col(i));
  }
  float MSx = vx.sum(), Lmb=MSx;
  Eigen::VectorXf bc(p);
  int numit = 0;
  float cnv = 1;
  std::vector<int> order(p); for(int j=0;j<p;j++) order[j]=j;
  while(numit<maxit){
    bc = b;
    std::shuffle(order.begin(), order.end(), std::mt19937(numit));
    for(int jj=0; jj<p; jj++){
      int j = order[jj];
      b0 = b[j];
      if(P_WEIGHTS){
        b[j] = (gen.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb/d[j]);
      }else{
        b[j] = (gen.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb);
      }
      e -= gen.col(j)*(b[j]-b0);
    }
    eM = e.mean();
    mu += eM;
    e = e.array()-eM;
    ve = (y.array()-mu).matrix().dot(e)/(float)n;
    vb = (y.array()-mu).matrix().dot((y.array()-mu).matrix()-e)/(float)(n*MSx);
    Lmb = ve/vb;
    ++numit;
    cnv = (bc-b).cwiseAbs().sum();
    if(cnv<tol){break;}
  }
  Eigen::VectorXf fit = y - e;
  h2 = vb*MSx/(vb*MSx+ve);
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("h2")=h2,
                      Named("Vb")=vb,
                      Named("Va")=vb*MSx,
                      Named("Ve")=ve);
}

// [[Rcpp::export]]
SEXP emGWA(Eigen::VectorXf y, Eigen::MatrixXf gen){
  int maxit = 500;
  float tol = 10e-8f;
  int p = gen.cols();
  int n = gen.rows();
  float b0, eM, ve, vb, h2, mu = y.mean(), vy = fvar(y);
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = gen.col(i).squaredNorm();
    vx[i] = fvar(gen.col(i));
  }
  float MSx = vx.sum(), Lmb=MSx;
  Eigen::VectorXf bc(p);
  int numit = 0;
  float cnv = 1;
  while(numit<maxit){
    bc = b;
    for(int j=0; j<p; j++){
      b0 = b[j];
      b[j] = (gen.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb);
      e -= gen.col(j)*(b[j]-b0);
    }
    ve = e.dot(y)/(n-1);
    vb = (vy-ve)/MSx;
    Lmb = ve/vb;
    eM = e.mean();
    mu += eM;
    e = e.array()-eM;
    ++numit;
    cnv = (bc-b).cwiseAbs().sum();
    if(cnv<tol){break;}
  }
  Eigen::VectorXf fit = gen*b;
  fit = fit.array()+mu;
  h2 = vb*MSx/(vb*MSx+ve);
  // Genome-wide screening
  Eigen::VectorXf LRT(p), b_ols(p);
  Eigen::VectorXf y0(n), e0(n), e1(n);
  float ve0, ve1, L0, L1;
  for(int j=0; j<p; j++){
    y0 = e + gen.col(j)*b[j];
    b_ols[j] = gen.col(j).dot(y0)/xx[j];
    e0 = y0.array()-y0.mean();
    e1 = y0 - gen.col(j)*b_ols[j];
    e1 = e1.array()-e1.mean();
    ve0 = y0.dot(e0)/(n-1);
    ve1 = y0.dot(e1)/(n-1);
    L0 = -e0.squaredNorm()/(2*ve0)-0.5f*n*std::log(6.28f*ve0);
    L1 = -e1.squaredNorm()/(2*ve1)-0.5f*n*std::log(6.28f*ve1);
    LRT[j] = 2*(L1-L0);
  }
  Rcpp::NumericVector LRT_rcpp(p);
  for(int j=0;j<p;j++) LRT_rcpp[j]=(double)LRT[j];
  Rcpp::NumericVector PVAL_rcpp = -log10(1-Rcpp::pchisq(LRT_rcpp,1,true,false));
  Eigen::VectorXf PVAL(p);
  for(int j=0;j<p;j++) PVAL[j]=(float)PVAL_rcpp[j];
  return List::create(Named("mu")=mu, Named("b")=b, Named("b_LS")=b_ols,
                      Named("h2")=h2, Named("hat")=fit,
                      Named("Vb")=vb, Named("Ve")=ve,
                      Named("LRT")=LRT, Named("PVAL")=PVAL);
}

// [[Rcpp::export]]
SEXP BayesA(Eigen::VectorXf y, Eigen::MatrixXf X,
            float it = 1500, float bi = 500,
            float df = 5, float R2 = 0.5){
  int p = X.cols(), n = X.rows();
  Eigen::VectorXf xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = X.col(i).squaredNorm();
    vx[i] = fvar(X.col(i));
  }
  float MSx = vx.sum();
  float vy = fvar(y);
  float Sb = (R2)*df*vy/MSx;
  float Se = (1-R2)*df*vy;
  float mu = y.mean();
  float b0,b1,eM,h2,MU=0,VE=0,vg,ve=vy;
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf B = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf VB = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf fit(n);
  Eigen::VectorXf vb = Eigen::VectorXf::Constant(p, Sb);
  Eigen::VectorXf Lmb = ve * vb.cwiseInverse();
  Eigen::VectorXf e = y.array()-mu;
  int iit = (int)it, ibi = (int)bi;
  for(int i=0; i<iit; i++){
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = R::rnorm((X.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb[j]), std::sqrt(ve/(xx[j]+Lmb[j])));
      b[j] = b1;
      vb[j] = (Sb+b1*b1)/R::rchisq(df+1);
      e -= X.col(j)*(b1-b0);
    }
    eM = R::rnorm(e.mean(), std::sqrt(ve/n));
    mu += eM; e = e.array()-eM;
    ve = (e.squaredNorm()+Se)/R::rchisq(n+df);
    Lmb = ve * vb.cwiseInverse();
    if(i>ibi){MU+=mu; B+=b; VB+=vb; VE+=ve;}
  }
  float MCMC = it-bi;
  MU/=MCMC; B/=MCMC; VB/=MCMC; VE/=MCMC;
  vg = VB.sum(); h2 = vg/(vg+VE);
  fit = X*B;
  fit = fit.array()+MU;
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("hat") = fit,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("h2") = h2, Named("MSx") = MSx);
}

// [[Rcpp::export]]
SEXP BayesB(Eigen::VectorXf y, Eigen::MatrixXf X,
            float it = 1500, float bi = 500,
            float pi = 0.95, float df = 5, float R2 = 0.5){
  int p = X.cols(), n = X.rows();
  int iit=(int)it, ibi=(int)bi;
  float MCMC = it-bi;
  Eigen::VectorXf xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = X.col(i).squaredNorm();
    vx[i] = fvar(X.col(i));
  }
  float MSx = vx.sum();
  float vy = fvar(y);
  float Sb = (R2)*df*vy/MSx;
  float Se = (1-R2)*df*vy;
  Eigen::VectorXf d = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf D = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf B = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf VB = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf vb = Eigen::VectorXf::Constant(p, Sb);
  float ve = vy;
  Eigen::VectorXf Lmb = ve * vb.cwiseInverse();
  float mu = y.mean();
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf e1(n), e2(n);
  float b0,b1,eM,h2,C,MU=0,VE=0,cj,dj,pj;
  for(int i=0; i<iit; i++){
    C = -0.5f/std::sqrt(ve);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = R::rnorm((X.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb[j]), std::sqrt(ve/(xx[j]+Lmb[j])));
      e1 = e - X.col(j)*(b1-b0);
      e2 = e - X.col(j)*(0-b0);
      cj = (1-pi)*std::exp(C*e1.squaredNorm());
      dj = (pi)*std::exp(C*e2.squaredNorm());
      pj = cj/(cj+dj);
      if(R::rbinom(1,pj)==1){
        b[j] = b1; d[j] = 1;
      }else{
        b[j] = R::rnorm(0, std::sqrt(ve/(xx[j]+Lmb[j]))); d[j] = 0;
      }
      vb[j] = (Sb+b[j]*b[j])/R::rchisq(df+1);
      e -= X.col(j)*(b[j]-b0);
    }
    eM = R::rnorm(e.mean(), std::sqrt(ve/n));
    mu += eM; e = e.array()-eM;
    ve = (e.squaredNorm()+Se)/R::rchisq(n+df);
    Lmb = ve * vb.cwiseInverse();
    if(i>ibi){MU+=mu; B+=b; VB+=vb; VE+=ve; D+=d;}
  }
  MU/=MCMC; B/=MCMC; VB/=MCMC; VE/=MCMC; D/=MCMC;
  float vg = VB.sum();
  h2 = vg/(vg+VE);
  Eigen::VectorXf fit = X*B;
  fit = fit.array()+MU;
  return List::create(Named("mu") = MU,
                      Named("b") = B, Named("d") = D,
                      Named("hat") = fit,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("h2") = h2, Named("MSx") = MSx);
}

// [[Rcpp::export]]
SEXP BayesC(Eigen::VectorXf y, Eigen::MatrixXf X,
            float it = 1500, float bi = 500,
            float pi = 0.95, float df = 5, float R2 = 0.5){
  int p = X.cols(), n = X.rows();
  int iit=(int)it, ibi=(int)bi;
  Eigen::VectorXf xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = X.col(i).squaredNorm();
    vx[i] = fvar(X.col(i));
  }
  float MSx = vx.sum();
  float vy = fvar(y);
  float Sb = df*(R2)*vy/MSx/(1-pi);
  float Se = df*(1-R2)*vy;
  float mu = y.mean();
  float b0,b1,eM,h2,C,MU=0,VB=0,VE=0,cj,dj,pj,vg,ve=vy,vb=Sb;
  Eigen::VectorXf d = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf D = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf B = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf fit(n);
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf e1(n), e2(n);
  float Lmb=ve/vb;
  for(int i=0; i<iit; i++){
    C = -0.5f/std::sqrt(ve);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = R::rnorm((X.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb), std::sqrt(ve/(xx[j]+Lmb)));
      e1 = e - X.col(j)*(b1-b0);
      e2 = e - X.col(j)*(0-b0);
      cj = (1-pi)*std::exp(C*e1.squaredNorm());
      dj = (pi)*std::exp(C*e2.squaredNorm());
      pj = cj/(cj+dj);
      if(R::rbinom(1,pj)==1){
        b[j] = b1; d[j] = 1;
      }else{
        b[j] = R::rnorm(0, std::sqrt(ve/(xx[j]+Lmb))); d[j] = 0;
      }
      e -= X.col(j)*(b[j]-b0);
    }
    eM = R::rnorm(e.mean(), std::sqrt(ve/n));
    mu += eM; e = e.array()-eM;
    vb = (b.squaredNorm()+Sb)/R::rchisq(df+p);
    ve = (e.squaredNorm()+Se)/R::rchisq(n+df);
    Lmb = ve/vb;
    if(i>ibi){MU+=mu; B+=b; D+=d; VB+=vb; VE+=ve;}
  }
  float MCMC = it-bi;
  MU/=MCMC; B/=MCMC; D/=MCMC; VB/=MCMC; VE/=MCMC;
  vg = VB*MSx; h2 = vg/(vg+VE);
  fit = X*B;
  fit = fit.array()+MU;
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("d") = D, Named("hat") = fit,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("h2") = h2, Named("MSx") = MSx);
}

// [[Rcpp::export]]
SEXP BayesL(Eigen::VectorXf y, Eigen::MatrixXf X,
            float it = 1500, float bi = 500,
            float df = 5, float R2 = 0.5){
  int p = X.cols(), n = X.rows();
  int iit=(int)it, ibi=(int)bi;
  Eigen::VectorXf xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = X.col(i).squaredNorm();
    vx[i] = fvar(X.col(i));
  }
  float MSx = vx.sum();
  float Phi = MSx*(1-R2)/R2;
  float vy = fvar(y);
  float Sb = (R2)*df*vy/MSx;
  float Se = (1-R2)*df*vy;
  float mu = y.mean();
  float b0,b1,eM,h2,MU=0,VE=0,vg,ve=vy;
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf B = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf VB = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf fit(n);
  Eigen::VectorXf vb = Eigen::VectorXf::Constant(p, Sb);
  Eigen::VectorXf Lmb = ve * vb.cwiseInverse();
  Eigen::VectorXf e = y.array()-mu;
  for(int i=0; i<iit; i++){
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = R::rnorm((X.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb[j]), std::sqrt(ve/(xx[j]+Lmb[j])));
      b[j] = b1;
      vb[j] = (Sb+b1*b1)/R::rchisq(df+1);
      e -= X.col(j)*(b1-b0);
    }
    eM = R::rnorm(e.mean(), std::sqrt(ve/n));
    mu += eM; e = e.array()-eM;
    ve = (e.squaredNorm()+Se)/R::rchisq(n+df);
    for(int j=0;j<p;j++) Lmb[j] = std::sqrt(Phi*ve/vb[j]);
    if(i>ibi){MU+=mu; B+=b; VB+=vb; VE+=ve;}
  }
  float MCMC = it-bi;
  MU/=MCMC; B/=MCMC; VB/=MCMC; VE/=MCMC;
  vg = VB.sum(); h2 = vg/(vg+VE);
  fit = X*B;
  fit = fit.array()+MU;
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("hat") = fit,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("h2") = h2, Named("MSx") = MSx);
}

// [[Rcpp::export]]
SEXP BayesRR(Eigen::VectorXf y, Eigen::MatrixXf X,
             float it = 1500, float bi = 500,
             float df = 5, float R2 = 0.5){
  int p = X.cols(), n = X.rows();
  int iit=(int)it, ibi=(int)bi;
  Eigen::VectorXf xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = X.col(i).squaredNorm();
    vx[i] = fvar(X.col(i));
  }
  float MSx = vx.sum();
  float vy = fvar(y);
  float Sb = (R2)*df*vy/MSx;
  float Se = (1-R2)*df*vy;
  float mu = y.mean();
  float b0,b1,eM,h2,MU=0,VE=0,VB=0,vg,ve=vy,vb=Sb,Lmb=ve/vb;
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf B = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf fit(n);
  Eigen::VectorXf e = y.array()-mu;
  for(int i=0; i<iit; i++){
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = R::rnorm((X.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb), std::sqrt(ve/(xx[j]+Lmb)));
      e -= X.col(j)*(b1-b0);
      b[j] = b1;
    }
    eM = R::rnorm(e.mean(), std::sqrt(ve/n));
    mu += eM; e = e.array()-eM;
    ve = (e.squaredNorm()+Se)/R::rchisq(n+df);
    vb = (b.squaredNorm()+Sb)/R::rchisq(p+df);
    Lmb = ve/vb;
    if(i>ibi){MU+=mu; B+=b; VB+=vb; VE+=ve;}
  }
  float MCMC = it-bi;
  MU/=MCMC; B/=MCMC; VB/=MCMC; VE/=MCMC;
  vg = VB*MSx; h2 = vg/(vg+VE);
  fit = X*B;
  fit = fit.array()+MU;
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("hat") = fit,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("h2") = h2, Named("MSx") = MSx);
}

// [[Rcpp::export]]
SEXP BayesCpi(Eigen::VectorXf y, Eigen::MatrixXf X,
          float it = 1500, float bi = 500,
          float df = 5, float R2 = 0.5){
  int p = X.cols(), n = X.rows();
  int iit=(int)it, ibi=(int)bi;
  Eigen::VectorXf xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = X.col(i).squaredNorm();
    vx[i] = fvar(X.col(i));
  }
  float MSx = vx.sum();
  float pi = 0.5f;
  float vy = fvar(y);
  float Sb = df*(R2)*vy/MSx/(1-pi);
  float Se = df*(1-R2)*vy;
  float mu = y.mean();
  float b0,b1,b2,eM,h2,C,MU=0,VB=0,VE=0,Pi=0,pj,vg,ve=vy,vb=Sb;
  Eigen::VectorXf d = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf D = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf B = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf fit(n);
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf e1(n), e2(n);
  float Lmb=ve/vb;
  for(int i=0; i<iit; i++){
    C = -0.5f/std::sqrt(ve);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = R::rnorm((X.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb), std::sqrt(ve/(xx[j]+Lmb)));
      b2 = R::rnorm(0, std::sqrt(ve/(xx[j]+Lmb)));
      e1 = e - X.col(j)*(b1-b0);
      e2 = e - X.col(j)*(0-b0);
      pj = (1-pi)*std::exp(C*(e1.squaredNorm()-e2.squaredNorm()));
      if(pj>1) pj = 1;
      if(R::rbinom(1,pj)==1){
        b[j] = b1; d[j] = 1;
      }else{
        b[j] = b2; d[j] = 0;
      }
      e -= X.col(j)*(b[j]-b0);
    }
    eM = R::rnorm(e.mean(), std::sqrt(ve/n));
    mu += eM; e = e.array()-eM;
    vb = (b.squaredNorm()+Sb)/R::rchisq(p+df);
    ve = (e.squaredNorm()+Se)/R::rchisq(n+df);
    Lmb = ve/vb;
    pi = d.mean();
    Sb = df*(R2)*vy/MSx/(1-pi);
    if(i>ibi){MU+=mu; B+=b; D+=d; VB+=vb; VE+=ve; Pi+=pi;}
  }
  float MCMC = it-bi;
  MU/=MCMC; B/=MCMC; D/=MCMC; VB/=MCMC; VE/=MCMC; Pi=1-Pi/MCMC;
  Eigen::VectorXf PVAL = (-1.0f*(1.0f-D.array()).log()).matrix();
  vg = VB*MSx/Pi; h2 = vg/(vg+VE);
  fit = X*B;
  fit = fit.array()+MU;
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("d") = D, Named("pi") = Pi,
                      Named("hat") = fit, Named("h2") = h2,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("PVAL") = PVAL);
}

// [[Rcpp::export]]
SEXP BayesDpi(Eigen::VectorXf y, Eigen::MatrixXf X,
          float it = 1500, float bi = 500,
          float df = 5, float R2 = 0.5){
  int p = X.cols(), n = X.rows();
  int iit=(int)it, ibi=(int)bi;
  Eigen::VectorXf xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = X.col(i).squaredNorm();
    vx[i] = fvar(X.col(i));
  }
  float MSx = vx.sum();
  float pi = 0.5f;
  float vy = fvar(y);
  float Sb = (R2)*df*vy/MSx;
  float Se = (1-R2)*df*vy;
  float mu = y.mean();
  float b0,b1,b2,eM,h2,C,MU=0,VE=0,Pi=0,pj,vg,ve=vy;
  Eigen::VectorXf d = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf D = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf B = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf VB = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf fit(n);
  Eigen::VectorXf vb = Eigen::VectorXf::Constant(p, Sb);
  Eigen::VectorXf Lmb = ve * vb.cwiseInverse();
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf e1(n), e2(n);
  for(int i=0; i<iit; i++){
    C = -0.5f/std::sqrt(ve);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = R::rnorm((X.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb[j]), std::sqrt(ve/(xx[j]+Lmb[j])));
      b2 = R::rnorm(0, std::sqrt(ve/(xx[j]+Lmb[j])));
      e1 = e - X.col(j)*(b1-b0);
      e2 = e - X.col(j)*(b2-b0);
      pj = (1-pi)*std::exp(C*(e1.squaredNorm()-e2.squaredNorm()));
      if(pj>1) pj = 1;
      if(R::rbinom(1,pj)==1){
        b[j] = b1; d[j] = 1;
      }else{
        b[j] = b2; d[j] = 0;
      }
      vb[j] = (Sb+b[j]*b[j])/R::rchisq(df+1);
      e -= X.col(j)*(b[j]-b0);
    }
    eM = R::rnorm(e.mean(), std::sqrt(ve/n));
    mu += eM; e = e.array()-eM;
    ve = (e.squaredNorm()+Se)/R::rchisq(n+df);
    Lmb = ve * vb.cwiseInverse();
    pi = d.mean();
    if(i>ibi){MU+=mu; B+=b; D+=d; VB+=vb; VE+=ve; Pi+=pi;}
  }
  float MCMC = it-bi;
  MU/=MCMC; B/=MCMC; D/=MCMC; VB/=MCMC; VE/=MCMC; Pi=1-Pi/MCMC;
  Eigen::VectorXf PVAL = (-1.0f*(1.0f-D.array()).log()).matrix();
  vg = VB.sum(); h2 = vg/(vg+VE);
  fit = X*B;
  fit = fit.array()+MU;
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("d") = D, Named("pi") = Pi,
                      Named("hat") = fit, Named("h2") = h2,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("PVAL") = PVAL);
}

// [[Rcpp::export]]
SEXP BayesA2(Eigen::VectorXf y, Eigen::MatrixXf X1, Eigen::MatrixXf X2,
             float it = 1500, float bi = 500,
             float df = 5, float R2 = 0.5){
  int n = X1.rows();
  int p1 = X1.cols();
  int p2 = X2.cols();
  int iit=(int)it, ibi=(int)bi;
  Eigen::VectorXf xx1(p1), vx1(p1);
  for(int i=0; i<p1; i++){
    xx1[i] = X1.col(i).squaredNorm();
    vx1[i] = fvar(X1.col(i));
  }
  float MSx1 = vx1.sum();
  Eigen::VectorXf xx2(p2), vx2(p2);
  for(int i=0; i<p2; i++){
    xx2[i] = X2.col(i).squaredNorm();
    vx2[i] = fvar(X2.col(i));
  }
  float MSx2 = vx2.sum();
  float vy = fvar(y);
  float Sb1 = (R2)*df*vy/MSx1;
  float Sb2 = (R2)*df*vy/MSx2;
  float Se = (1-R2)*df*vy;
  float mu = y.mean();
  float b_t0,b_t1,eM,h2,MU=0,VE=0,vg,ve=vy;
  Eigen::VectorXf b1 = Eigen::VectorXf::Zero(p1);
  Eigen::VectorXf B1 = Eigen::VectorXf::Zero(p1);
  Eigen::VectorXf VB1 = Eigen::VectorXf::Zero(p1);
  Eigen::VectorXf b2 = Eigen::VectorXf::Zero(p2);
  Eigen::VectorXf B2 = Eigen::VectorXf::Zero(p2);
  Eigen::VectorXf VB2 = Eigen::VectorXf::Zero(p2);
  Eigen::VectorXf vb1 = Eigen::VectorXf::Constant(p1, Sb1);
  Eigen::VectorXf vb2 = Eigen::VectorXf::Constant(p2, Sb2);
  Eigen::VectorXf Lmb1 = ve * vb1.cwiseInverse();
  Eigen::VectorXf Lmb2 = ve * vb2.cwiseInverse();
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf fit(n);
  for(int i=0; i<iit; i++){
    for(int j=0; j<p1; j++){
      b_t0 = b1[j];
      b_t1 = R::rnorm((X1.col(j).dot(e)+xx1[j]*b_t0)/(xx1[j]+Lmb1[j]), std::sqrt(ve/(xx1[j]+Lmb1[j])));
      b1[j] = b_t1;
      vb1[j] = (Sb1+b1[j]*b1[j])/R::rchisq(df+1);
      e -= X1.col(j)*(b_t1-b_t0);
    }
    for(int j=0; j<p2; j++){
      b_t0 = b2[j];
      b_t1 = R::rnorm((X2.col(j).dot(e)+xx2[j]*b_t0)/(xx2[j]+Lmb2[j]), std::sqrt(ve/(xx2[j]+Lmb2[j])));
      b2[j] = b_t1;
      vb2[j] = (Sb2+b2[j]*b2[j])/R::rchisq(df+1);
      e -= X2.col(j)*(b_t1-b_t0);
    }
    eM = R::rnorm(e.mean(), std::sqrt(ve/n));
    mu += eM; e = e.array()-eM;
    ve = (e.squaredNorm()+Se)/R::rchisq(n+df);
    Lmb1 = ve * vb1.cwiseInverse();
    Lmb2 = ve * vb2.cwiseInverse();
    if(i>ibi){MU+=mu; VE+=ve; B1+=b1; VB1+=vb1; B2+=b2; VB2+=vb2;}
  }
  float MCMC = it-bi;
  MU/=MCMC; VE/=MCMC; B1/=MCMC; VB1/=MCMC; B2/=MCMC; VB2/=MCMC;
  vg = VB1.sum()+VB2.sum(); h2 = vg/(vg+VE);
  fit = X1*B1 + X2*B2;
  fit = fit.array()+MU;
  return List::create(Named("hat") = fit, Named("mu") = MU,
                      Named("b1") = B1, Named("b2") = B2,
                      Named("vb1") = VB1, Named("vb2") = VB2,
                      Named("ve") = VE, Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP BayesB2(Eigen::VectorXf y, Eigen::MatrixXf X1, Eigen::MatrixXf X2,
             float it = 1500, float bi = 500,
             float pi = 0.95, float df = 5, float R2 = 0.5){
  int n = X1.rows();
  int p1 = X1.cols();
  int p2 = X2.cols();
  int iit=(int)it, ibi=(int)bi;
  Eigen::VectorXf xx1(p1), vx1(p1);
  for(int i=0; i<p1; i++){
    xx1[i] = X1.col(i).squaredNorm();
    vx1[i] = fvar(X1.col(i));
  }
  float MSx1 = vx1.sum();
  Eigen::VectorXf xx2(p2), vx2(p2);
  for(int i=0; i<p2; i++){
    xx2[i] = X2.col(i).squaredNorm();
    vx2[i] = fvar(X2.col(i));
  }
  float MSx2 = vx2.sum();
  float vy = fvar(y);
  float Sb1 = (R2)*df*vy/MSx1;
  float Sb2 = (R2)*df*vy/MSx2;
  float Se = (1-R2)*df*vy;
  float mu = y.mean();
  float b_t0,b_t1,b_t2,eM,h2,C,MU=0,VE=0,cj,dj,pj,vg,ve=vy;
  Eigen::VectorXf d1 = Eigen::VectorXf::Zero(p1);
  Eigen::VectorXf b1 = Eigen::VectorXf::Zero(p1);
  Eigen::VectorXf D1 = Eigen::VectorXf::Zero(p1);
  Eigen::VectorXf B1 = Eigen::VectorXf::Zero(p1);
  Eigen::VectorXf VB1 = Eigen::VectorXf::Zero(p1);
  Eigen::VectorXf d2 = Eigen::VectorXf::Zero(p2);
  Eigen::VectorXf b2 = Eigen::VectorXf::Zero(p2);
  Eigen::VectorXf D2 = Eigen::VectorXf::Zero(p2);
  Eigen::VectorXf B2 = Eigen::VectorXf::Zero(p2);
  Eigen::VectorXf VB2 = Eigen::VectorXf::Zero(p2);
  Eigen::VectorXf vb1 = Eigen::VectorXf::Constant(p1, Sb1);
  Eigen::VectorXf vb2 = Eigen::VectorXf::Constant(p2, Sb2);
  Eigen::VectorXf Lmb1 = ve * vb1.cwiseInverse();
  Eigen::VectorXf Lmb2 = ve * vb2.cwiseInverse();
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf e1(n), e2(n), fit(n);
  for(int i=0; i<iit; i++){
    C = -0.5f/std::sqrt(ve);
    for(int j=0; j<p1; j++){
      b_t0 = b1[j];
      b_t1 = R::rnorm((X1.col(j).dot(e)+xx1[j]*b_t0)/(xx1[j]+Lmb1[j]), std::sqrt(ve/(xx1[j]+Lmb1[j])));
      b_t2 = R::rnorm(0, std::sqrt(ve/(xx1[j]+Lmb1[j])));
      e1 = e - X1.col(j)*(b_t1-b_t0);
      e2 = e - X1.col(j)*(b_t2-b_t0);
      cj = (1-pi)*std::exp(C*e1.squaredNorm());
      dj = (pi)*std::exp(C*e2.squaredNorm());
      pj = cj/(cj+dj);
      if(R::rbinom(1,pj)==1){
        b1[j] = b_t1; d1[j] = 1;
      }else{
        b1[j] = b_t2; d1[j] = 0;
      }
      vb1[j] = (Sb1+b1[j]*b1[j])/R::rchisq(df+1);
      e -= X1.col(j)*(b1[j]-b_t0);
    }
    for(int j=0; j<p2; j++){
      b_t0 = b2[j];
      b_t1 = R::rnorm((X2.col(j).dot(e)+xx2[j]*b_t0)/(xx2[j]+Lmb2[j]), std::sqrt(ve/(xx2[j]+Lmb2[j])));
      b_t2 = R::rnorm(0, std::sqrt(ve/(xx2[j]+Lmb2[j])));
      e1 = e - X2.col(j)*(b_t1-b_t0);
      e2 = e - X2.col(j)*(b_t2-b_t0);
      cj = (1-pi)*std::exp(C*e1.squaredNorm());
      dj = (pi)*std::exp(C*e2.squaredNorm());
      pj = cj/(cj+dj);
      if(R::rbinom(1,pj)==1){
        b2[j] = b_t1; d2[j] = 1;
      }else{
        b2[j] = b_t2; d2[j] = 0;
      }
      vb2[j] = (Sb2+b2[j]*b2[j])/R::rchisq(df+1);
      e -= X2.col(j)*(b2[j]-b_t0);
    }
    eM = R::rnorm(e.mean(), std::sqrt(ve/n));
    mu += eM; e = e.array()-eM;
    ve = (e.squaredNorm()+Se)/R::rchisq(n+df);
    Lmb1 = ve * vb1.cwiseInverse();
    Lmb2 = ve * vb2.cwiseInverse();
    if(i>ibi){MU+=mu; VE+=ve; B1+=b1; D1+=d1; VB1+=vb1; B2+=b2; D2+=d2; VB2+=vb2;}
  }
  float MCMC = it-bi;
  MU/=MCMC; VE/=MCMC; B1/=MCMC; D1/=MCMC; VB1/=MCMC; B2/=MCMC; D2/=MCMC; VB2/=MCMC;
  vg = VB1.sum()+VB2.sum(); h2 = vg/(vg+VE);
  fit = X1*B1 + X2*B2;
  fit = fit.array()+MU;
  return List::create(Named("mu") = MU,
                      Named("b1") = B1, Named("d1") = D1, Named("vb1") = VB1,
                      Named("b2") = B2, Named("d2") = D2, Named("vb2") = VB2,
                      Named("ve") = VE, Named("hat") = fit, Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP BayesRR2(Eigen::VectorXf y, Eigen::MatrixXf X1, Eigen::MatrixXf X2,
              float it = 1500, float bi = 500,
              float df = 5, float R2 = 0.5){
  int n = X1.rows();
  int p1 = X1.cols();
  int p2 = X2.cols();
  int iit=(int)it, ibi=(int)bi;
  Eigen::VectorXf xx1(p1), vx1(p1);
  for(int i=0; i<p1; i++){
    xx1[i] = X1.col(i).squaredNorm();
    vx1[i] = fvar(X1.col(i));
  }
  float MSx1 = vx1.sum();
  Eigen::VectorXf xx2(p2), vx2(p2);
  for(int i=0; i<p2; i++){
    xx2[i] = X2.col(i).squaredNorm();
    vx2[i] = fvar(X2.col(i));
  }
  float MSx2 = vx2.sum();
  float vy = fvar(y);
  float Sb1 = (R2)*df*vy/MSx1;
  float Sb2 = (R2)*df*vy/MSx2;
  float Se = (1-R2)*df*vy;
  float mu = y.mean();
  float b_t0,b_t1,eM,h2,MU=0,VE=0,vg,vb1,vb2,VB1=0,VB2=0,Lmb1=MSx1,Lmb2=MSx2,ve=vy;
  Eigen::VectorXf b1 = Eigen::VectorXf::Zero(p1);
  Eigen::VectorXf B1 = Eigen::VectorXf::Zero(p1);
  Eigen::VectorXf b2 = Eigen::VectorXf::Zero(p2);
  Eigen::VectorXf B2 = Eigen::VectorXf::Zero(p2);
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf fit(n);
  for(int i=0; i<iit; i++){
    for(int j=0; j<p1; j++){
      b_t0 = b1[j];
      b_t1 = R::rnorm((X1.col(j).dot(e)+xx1[j]*b_t0)/(xx1[j]+Lmb1), std::sqrt(ve/(xx1[j]+Lmb1)));
      b1[j] = b_t1;
      e -= X1.col(j)*(b_t1-b_t0);
    }
    for(int j=0; j<p2; j++){
      b_t0 = b2[j];
      b_t1 = R::rnorm((X2.col(j).dot(e)+xx2[j]*b_t0)/(xx2[j]+Lmb2), std::sqrt(ve/(xx2[j]+Lmb2)));
      b2[j] = b_t1;
      e -= X2.col(j)*(b_t1-b_t0);
    }
    eM = R::rnorm(e.mean(), std::sqrt(ve/n));
    mu += eM; e = e.array()-eM;
    ve = (e.squaredNorm()+Se)/R::rchisq(n+df);
    vb1 = (Sb1+b1.squaredNorm())/R::rchisq(df+p1);
    vb2 = (Sb2+b2.squaredNorm())/R::rchisq(df+p2);
    Lmb1 = ve/vb1; Lmb2 = ve/vb2;
    if(i>ibi){MU+=mu; VE+=ve; B1+=b1; VB1+=vb1; B2+=b2; VB2+=vb2;}
  }
  float MCMC = it-bi;
  MU/=MCMC; VE/=MCMC; B1/=MCMC; VB1/=MCMC; B2/=MCMC; VB2/=MCMC;
  vg = (VB1*MSx1+VB2*MSx2); h2 = vg/(vg+VE);
  fit = X1*B1 + X2*B2;
  fit = fit.array()+MU;
  return List::create(Named("hat") = fit, Named("mu") = MU,
                      Named("b1") = B1, Named("b2") = B2,
                      Named("vb1") = VB1, Named("vb2") = VB2,
                      Named("ve") = VE, Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emML2(Eigen::VectorXf y, Eigen::MatrixXf X1, Eigen::MatrixXf X2,
           Rcpp::Nullable<Rcpp::NumericVector> D1 = R_NilValue,
           Rcpp::Nullable<Rcpp::NumericVector> D2 = R_NilValue){
  int maxit = 350; float tol = 10e-8f;
  int p1 = X1.cols();
  int p2 = X2.cols();
  int n = X1.rows();
  bool P1_WEIGHTS = false;
  bool P2_WEIGHTS = false;
  Eigen::VectorXf d1 = Eigen::VectorXf::Ones(p1);
  Eigen::VectorXf d2 = Eigen::VectorXf::Ones(p2);
  if(D1.isNotNull()){
    P1_WEIGHTS = true;
    Rcpp::NumericVector Dtmp(D1);
    for(int i=0;i<p1;i++) d1[i]=(float)Dtmp[i];
  }
  if(D2.isNotNull()){
    P2_WEIGHTS = true;
    Rcpp::NumericVector Dtmp(D2);
    for(int i=0;i<p2;i++) d2[i]=(float)Dtmp[i];
  }
  float b0, eM, ve, vb1, vb2, h2, mu = y.mean();
  Eigen::VectorXf b1 = Eigen::VectorXf::Zero(p1);
  Eigen::VectorXf b2 = Eigen::VectorXf::Zero(p2);
  Eigen::VectorXf u1(n), u2(n), cY(n);
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf x1x1(p1), vx1(p1);
  for(int i=0; i<p1; i++){
    x1x1[i] = X1.col(i).squaredNorm();
    vx1[i] = fvar(X1.col(i));
  }
  float MSx1 = vx1.sum(), Lmb1=MSx1;
  Eigen::VectorXf x2x2(p2), vx2(p2);
  for(int i=0; i<p2; i++){
    x2x2[i] = X2.col(i).squaredNorm();
    vx2[i] = fvar(X2.col(i));
  }
  float MSx2 = vx2.sum(), Lmb2=MSx2;
  Eigen::VectorXf bc1(p1), bc2(p2);
  int numit = 0;
  float cnv = 1;
  while(numit<maxit){
    bc1 = b1; bc2 = b2;
    for(int j=0; j<p1; j++){
      b0 = b1[j];
      if(P1_WEIGHTS){
        b1[j] = (X1.col(j).dot(e)+x1x1[j]*b0)/(x1x1[j]+Lmb1/d1[j]);
      }else{
        b1[j] = (X1.col(j).dot(e)+x1x1[j]*b0)/(x1x1[j]+Lmb1);
      }
      e -= X1.col(j)*(b1[j]-b0);
    }
    for(int j=0; j<p2; j++){
      b0 = b2[j];
      if(P2_WEIGHTS){
        b2[j] = (X2.col(j).dot(e)+x2x2[j]*b0)/(x2x2[j]+Lmb2/d2[j]);
      }else{
        b2[j] = (X2.col(j).dot(e)+x2x2[j]*b0)/(x2x2[j]+Lmb2);
      }
      e -= X2.col(j)*(b2[j]-b0);
    }
    u1 = X1*b1;
    u2 = X2*b2;
    eM = e.mean();
    mu += eM;
    e = e.array()-eM;
    cY = u1+u2+e;
    ve = e.dot(cY)/(float)n;
    vb1 = (u1.dot(cY)/(float)n)/MSx1;
    vb2 = (u2.dot(cY)/(float)n)/MSx2;
    Lmb1 = ve/vb1;
    Lmb2 = ve/vb2;
    ++numit;
    cnv = (bc1-b1).cwiseAbs().sum()+(bc2-b2).cwiseAbs().sum();
    if(cnv<tol){break;}
  }
  Eigen::VectorXf fit = Eigen::VectorXf::Constant(n, mu) + u1 + u2;
  h2 = 1-ve/fvar(y);
  return List::create(Named("mu")=mu,
                      Named("b1")=b1, Named("b2")=b2,
                      Named("Vb1")=vb1, Named("Vb2")=vb2, Named("Ve")=ve,
                      Named("u1")=u1, Named("u2")=u2,
                      Named("MSx1")=MSx1, Named("MSx2")=MSx2,
                      Named("h2")=h2, Named("hat")=fit);
}

// [[Rcpp::export]]
Eigen::MatrixXf CNT(Eigen::MatrixXf X){
  for(int j=0; j<X.cols(); j++){
    X.col(j) = X.col(j).array() - X.col(j).mean();
  }
  return X;
}

// [[Rcpp::export]]
Eigen::MatrixXf IMP(Eigen::MatrixXf X){
  int p = X.cols(); int n = X.rows();
  for(int j=0; j<p; j++){
    bool hasna = false;
    for(int i=0; i<n; i++){
      if(std::isnan(X(i,j))){ hasna=true; break; }
    }
    if(hasna){
      float sum=0; int cnt=0;
      for(int i=0; i<n; i++){
        if(!std::isnan(X(i,j))){ sum+=X(i,j); cnt++; }
      }
      float EXP = (cnt>0) ? sum/(float)cnt : 0.0f;
      for(int i=0; i<n; i++){
        if(std::isnan(X(i,j))) X(i,j) = EXP;
      }
    }
  }
  return X;
}

// [[Rcpp::export]]
Eigen::MatrixXf GAU(Eigen::MatrixXf X){
  int n = X.rows();
  Eigen::MatrixXf K(n,n);
  float d2, md=0;
  for(int i=0; i<n; i++){
    K(i,i) = 0;
    for(int j=i+1; j<n; j++){
      Eigen::VectorXf D = (X.row(i) - X.row(j)).transpose();
      d2 = D.squaredNorm();
      K(i,j) = d2; K(j,i) = d2;
    }
  }
  // mean of off-diagonal
  float tot=0; int cnt=0;
  for(int i=0; i<n; i++) for(int j=0; j<n; j++) if(i!=j){ tot+=K(i,j); cnt++; }
  md = tot/(float)cnt;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      K(i,j) = std::exp(-K(i,j)/md);
    }
  }
  return K;
}

// [[Rcpp::export]]
Eigen::MatrixXf GRM(Eigen::MatrixXf X, bool Code012 = false){
  int n = X.rows(), p = X.cols();
  Eigen::MatrixXf K(n,n);
  Eigen::VectorXf xx(p);
  float zz, Sum2pq=0.0f;
  for(int i=0; i<p; i++){ xx[i] = X.col(i).mean(); }
  if(Code012){
    for(int i=0; i<p; i++){ Sum2pq = Sum2pq + xx[i]*xx[i]/2.0f; }
  }else{
    for(int i=0; i<p; i++){ Sum2pq = Sum2pq + fvar(X.col(i)); }
  }
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(i<=j){
        zz = (X.row(i)-xx.transpose()).dot(X.row(j)-xx.transpose());
        K(i,j) = zz; K(j,i) = zz;
      }
    }
  }
  return K/Sum2pq;
}

// [[Rcpp::export]]
Eigen::VectorXf SPC(Eigen::VectorXf y, Eigen::VectorXf blk, Eigen::VectorXf row, Eigen::VectorXf col, float rN=3, float cN=1){
  int n = y.size(), t1=0, t2=0;
  Eigen::VectorXf Cov = Eigen::VectorXf::Zero(n);
  Eigen::VectorXf Phe = Eigen::VectorXf::Zero(n);
  Eigen::VectorXf Obs = Eigen::VectorXf::Zero(n);
  for(int i=0; i<n; i++){ for(int j=0; j<n; j++){
    t1 = (int)(row[i]-row[j]); if(t1<0){t1=-t1;}
    t2 = (int)(col[i]-col[j]); if(t2<0){t2=-t2;}
    if( (i>j) && (blk[i]==blk[j]) && (t1<=rN) && (t2<=cN) ){
      Phe[i] += y[j]; Obs[i] += 1; Phe[j] += y[i]; Obs[j] += 1; }}}
  for(int i=0; i<n; i++){ if(Obs[i]>0) Cov[i] = Phe[i]/Obs[i]; }
  return Cov;
}

// [[Rcpp::export]]
Eigen::MatrixXf SPM(Eigen::VectorXf blk, Eigen::VectorXf row, Eigen::VectorXf col, float rN=3, float cN=1){
  int n=blk.size(), t1=0, t2=0;
  Eigen::MatrixXf X = Eigen::MatrixXf::Zero(n,n);
  for(int i=0; i<n; i++){ for(int j=0; j<n; j++){
    t1 = (int)(row[i]-row[j]); if(t1<0){t1=-t1;}
    t2 = (int)(col[i]-col[j]); if(t2<0){t2=-t2;}
    if( (blk[i]==blk[j]) && (i>j) && (t1<=rN) && (t2<=cN) ){ X(i,j)=1.0f; X(j,i)=1.0f; }
  }}
  return X;
}

// [[Rcpp::export]]
SEXP mtgsru(Eigen::MatrixXf Y, Eigen::MatrixXf X,
            Eigen::MatrixXf b, Eigen::MatrixXf vb, Eigen::VectorXf ve,
            Eigen::MatrixXf iG, int maxit = 50){
  float tol = 10e-8f;
  int k = Y.cols(), p = X.cols(), n0 = X.rows();
  Eigen::MatrixXf fit(n0,k), o(n0,k), y(n0,k), e(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){ o(i,j)=0.0f; y(i,j)=0.0f; }
      else{ o(i,j)=1.0f; y(i,j)=Y(i,j); }}}
  Eigen::VectorXf n = o.colwise().sum().transpose();
  Eigen::MatrixXf xx(p,k), vx(p,k); float tmp;
  for(int i=0; i<p; i++){ for(int j=0; j<k; j++){
    xx(i,j) = (X.col(i).array()*X.col(i).array()*o.col(j).array()).sum();
    tmp = (X.col(i).array()*o.col(j).array()).sum()/n(j);
    vx(i,j) = xx(i,j)/n(j) - tmp*tmp; }}
  Eigen::VectorXf MSx = vx.colwise().sum().transpose();
  Eigen::MatrixXf LHS(k,k);
  Eigen::VectorXf b0(k), b1(k), vy(k), RHS(k);
  for(int i=0; i<k; i++){
    e.col(i) = y.col(i);
    vy(i) = e.col(i).squaredNorm()/n(i); }
  iG = vb.inverse();
  Eigen::MatrixXf bc(p,k);
  int numit=0; float cnv=1;
  while(numit<maxit){
    bc = b;
    for(int j=0; j<p; j++){
      b0 = b.row(j).transpose();
      LHS = iG;
      for(int i=0; i<k; i++){ LHS(i,i) += xx(j,i)/ve(i); }
      for(int i=0; i<k; i++){ RHS(i) = (e.col(i).dot(X.col(j)) + xx(j,i)*b0(i))/ve(i); }
      b1 = LHS.llt().solve(RHS);
      b.row(j) = b1.transpose();
      for(int i=0; i<k; i++){ e.col(i) = (e.col(i) - X.col(j)*(b1(i)-b0(i))).array()*o.col(i).array(); }}
    for(int i=0; i<k; i++){ ve(i) = e.col(i).dot(y.col(i))/(n(i)-1); }
    for(int i=0; i<n0; i++){ for(int j=0; j<k; j++){ fit(i,j) = X.row(i).dot(b.col(j)); }}
    for(int i=0; i<k; i++){ for(int j=0; j<k; j++){
      vb(i,j) = (fit.col(i).dot(y.col(j))+fit.col(j).dot(y.col(i))) /
                ((n(i)*MSx(i))+(n(j)*MSx(j))); }}
    for(int i=0; i<k; i++){ vb(i,i) *= 1.01f; }
    iG = vb.inverse();
    ++numit; cnv = (bc-b).cwiseAbs().sum(); if(cnv<tol){break;}}
  for(int i=0; i<n0; i++){ for(int j=0; j<k; j++){ fit(i,j) = X.row(i).dot(b.col(j)); }}
  Eigen::VectorXf h2(k); h2 = (1.0f - ve.array()/vy.array()).matrix();
  return List::create(Named("b")=b, Named("hat")=fit, Named("e")=fit, Named("MSx")=MSx,
                      Named("vb")=vb, Named("ve")=ve, Named("iG")=iG, Named("h2")=h2);
}

// [[Rcpp::export]]
SEXP lasso(Eigen::VectorXf y, Eigen::MatrixXf gen){
  int maxit = 300;
  float tol = 10e-8f;
  int p = gen.cols(), n = gen.rows();
  float eM, mu = y.mean();
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf e = y.array() - mu;
  Eigen::VectorXf xx(p);
  for(int i=0; i<p; i++){ xx[i] = gen.col(i).squaredNorm(); }
  float Lmb = xx.mean()/p;
  Eigen::VectorXf bc(p), yx(p);
  int numit=0; float cnv=1;
  while(numit<maxit){
    bc = b;
    for(int j=0; j<p; j++){
      e += gen.col(j)*b[j];
      yx[j] = e.dot(gen.col(j));
      if(yx[j]>0){
        b[j] = (yx[j]-Lmb)/xx[j];
        if(b[j]<0) b[j]=0;
      }else{
        b[j] = (yx[j]+Lmb)/xx[j];
        if(b[j]>0) b[j]=0;
      }
      e -= gen.col(j)*b[j]; }
    float tmp = 0.0f;
    for(int j=0; j<p; j++) tmp += std::fabs(yx[j]) - std::fabs(b[j]*xx[j]);
    Lmb = 2.0f * tmp / p;
    Lmb = 2.0f * std::sqrt(std::fabs(Lmb));
    eM = e.mean(); mu += eM; e = e.array()-eM;
    ++numit; cnv = (bc-b).cwiseAbs().sum(); if(cnv<tol){break;}}
  Eigen::VectorXf fit = y - e;
  float h2 = 1.0f - (e.dot(y)/(n-1))/fvar(y);
  return List::create(Named("mu")=mu, Named("b")=b,
                      Named("h2")=h2, Named("hat")=fit,
                      Named("Lmb")=Lmb);
}

// [[Rcpp::export]]
SEXP emBCpi(Eigen::VectorXf y, Eigen::MatrixXf gen, float df = 10, float R2 = 0.5, float Pi = 0.75){
  int it = 200;
  int p = gen.cols(), n = gen.rows();
  Eigen::VectorXf d = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  float vy = fvar(y);
  if(Pi>0.5f) Pi = 1.0f-Pi;
  Eigen::VectorXf xx(p), vx(p);
  for(int i=0; i<p; i++){ xx[i]=gen.col(i).squaredNorm(); vx[i]=fvar(gen.col(i)); }
  float PriorPi = Pi;
  float MSx = vx.sum()*Pi*(1.0f-Pi);
  float Sa = R2*(df+2)*vy/MSx;
  float Se = (1.0f-R2)*(df+2)*vy;
  float mu = y.mean();
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf e1(n), e2(n);
  float ve=Sa, va=Se, Lmb=ve/va;
  float b0,b1,LR,eM,h2,C;
  float Pi0 = (1.0f-Pi)/Pi;
  for(int i=0; i<it; i++){
    C = -0.5f/std::sqrt(ve);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = (gen.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb);
      e1 = e - gen.col(j)*(b1-b0);
      e2 = e - gen.col(j)*(0.0f-b0);
      LR = Pi0*std::exp(C*(e2.squaredNorm()-e1.squaredNorm()));
      d[j] = 1.0f/(1.0f+LR);
      b[j] = b1*d[j];
      e -= gen.col(j)*(b[j]-b0); }
    Pi = ((1.0f-d.mean())*p + PriorPi*df)/(p+df);
    Pi0 = (1.0f-Pi)/Pi;
    MSx = vx.sum()*Pi*(1.0f-Pi);
    Sa = R2*(df+2)*vy/MSx;
    eM = e.mean();
    ve = (e.squaredNorm()+Se)/(n+df);
    va = (b.squaredNorm()+Sa)/(p+df)/(d.mean()-Pi);
    Lmb = ve/va;
    eM = e.mean(); mu += eM; e = e.array()-eM; }
  h2 = 1.0f-ve/vy;
  Eigen::VectorXf fit = gen*b; fit = fit.array()+mu;
  return List::create(Named("mu")=mu, Named("b")=b, Named("d")=d,
                      Named("pi")=Pi, Named("hat")=fit,
                      Named("Vg")=va*MSx, Named("Va")=va,
                      Named("Ve")=ve, Named("h2")=h2);
}

// [[Rcpp::export]]
Eigen::MatrixXf NNSEARCH(Eigen::VectorXf blk, Eigen::VectorXf row, Eigen::VectorXf col, int rN, int cN){
  int n = blk.size();
  Eigen::MatrixXf X = Eigen::MatrixXf::Zero(n,(rN*2+1)*(cN*2+1));
  Eigen::VectorXf Obs = Eigen::VectorXf::Zero(n);
  for(int i=0; i<n; i++){ for(int j=0; j<n; j++){
    if( (i>j) && (blk[i]==blk[j]) &&
        (std::fabs(row[i]-row[j])<=rN) && (std::fabs(col[i]-col[j])<=cN) ){
      X(i,(int)Obs[i]) = (float)(j+1);
      X(j,(int)Obs[j]) = (float)(i+1);
      Obs[i] += 1; Obs[j] += 1; }}}
  return X;
}

// [[Rcpp::export]]
SEXP GSFLM(Eigen::VectorXf y, Eigen::VectorXf e, Eigen::MatrixXf gen,
           Eigen::VectorXf b, Eigen::VectorXf Lmb, Eigen::VectorXf xx,
           float cxx, int maxit = 50){
  float tol = 10e-8f;
  float phi = cxx;
  Eigen::VectorXf e0 = e;
  int p = gen.cols(), n = gen.rows();
  float vy = fvar(y);
  float vna = y.dot(e)/(n-1);
  float b0, eM;
  float mu = e.mean(); e = e.array()-mu;
  Eigen::VectorXf Vb(p);
  float b1;
  Eigen::VectorXf bc(p);
  int numit=0; float cnv=1.0f;
  while(numit<maxit){
    bc = b;
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = (gen.col(j).dot(e)+xx[j]*b0)/(Lmb(j)+xx(j)+0.01f);
      b[j] = b1;
      e -= gen.col(j)*(b1-b0); }
    eM = e.mean(); mu += eM; e = e.array()-eM;
    vna = e.dot(e0)/n;
    Vb = (b.array()*b.array() + vna/(xx.array()+Lmb.array())).matrix();
    for(int j=0; j<p; j++) Lmb[j] = std::sqrt(phi*vna/Vb[j]);
    ++numit; cnv = (bc-b).cwiseAbs().sum(); if(cnv<tol){break;}}
  return List::create(Named("mu")=mu, Named("b")=b,
                      Named("h2")=1.0f-vna/vy, Named("e")=e,
                      Named("Lmb")=Lmb, Named("vb")=Vb);
}

// [[Rcpp::export]]
SEXP GSRR(Eigen::VectorXf y, Eigen::VectorXf e, Eigen::MatrixXf gen,
          Eigen::VectorXf b, Eigen::VectorXf Lmb, Eigen::VectorXf xx,
          float cxx, int maxit = 50){
  float tol = 10e-8f;
  float phi = cxx;
  Eigen::VectorXf e0 = e;
  int p = gen.cols(), n = gen.rows();
  float vy = fvar(y);
  float vna = y.dot(e)/(n-1);
  float b0, eM, vg, LmbTmp;
  float mu = e.mean(); e = e.array()-mu;
  Eigen::VectorXf Vb(p);
  float b1;
  Eigen::VectorXf bc(p);
  int numit=0; float cnv=1.0f;
  while(numit<maxit){
    bc = b;
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = (gen.col(j).dot(e)+xx[j]*b0)/(Lmb(j)+xx(j)+0.01f);
      b[j] = b1;
      e -= gen.col(j)*(b1-b0); }
    eM = e.mean(); mu += eM; e = e.array()-eM;
    vna = e.dot(e0)/n;
    vg = (vy-vna)/phi;
    LmbTmp = vna/vg;
    for(int j=0; j<p; j++){ Vb[j]=vg; Lmb[j]=LmbTmp; }
    ++numit; cnv = (bc-b).cwiseAbs().sum(); if(cnv<tol){break;}}
  return List::create(Named("mu")=mu, Named("b")=b,
                      Named("h2")=1.0f-vna/vy, Named("e")=e,
                      Named("Lmb")=Lmb, Named("vb")=Vb);
}
