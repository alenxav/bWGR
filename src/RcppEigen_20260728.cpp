// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>


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
  return Rcpp::List::create(Rcpp::Named("b") = b, Rcpp::Named("d") = d, Rcpp::Named("e") = e);
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
  return Rcpp::List::create(Rcpp::Named("b") = b, Rcpp::Named("d") = d, Rcpp::Named("e") = e0);
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
  return Rcpp::List::create(Rcpp::Named("mu") = mu,
                            Rcpp::Named("b") = b,
                            Rcpp::Named("hat") = fit,
                            Rcpp::Named("Vb") = vb,
                            Rcpp::Named("Ve") = ve,
                            Rcpp::Named("h2") = h2);
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
  return Rcpp::List::create(Rcpp::Named("mu") = mu,
                            Rcpp::Named("b") = b,
                            Rcpp::Named("d") = d,
                            Rcpp::Named("hat") = fit,
                            Rcpp::Named("Vb") = vb,
                            Rcpp::Named("Ve") = ve,
                            Rcpp::Named("h2") = h2);
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
  return Rcpp::List::create(Rcpp::Named("mu") = mu,
                            Rcpp::Named("b") = b,
                            Rcpp::Named("d") = d,
                            Rcpp::Named("hat") = fit,
                            Rcpp::Named("Vg") = va*MSx,
                            Rcpp::Named("Va") = va,
                            Rcpp::Named("Ve") = ve,
                            Rcpp::Named("h2") = h2);
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
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=fit,
                            Rcpp::Named("Vb")=Vb,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("h2")=Vb.sum()/(Vb.sum()+Ve));
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
  return Rcpp::List::create(Rcpp::Named("mu") = mu,
                            Rcpp::Named("b") = b,
                            Rcpp::Named("hat") = fit,
                            Rcpp::Named("Va") = vb,
                            Rcpp::Named("Ve") = ve,
                            Rcpp::Named("h2") = h2);
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
  return Rcpp::List::create(Rcpp::Named("mu") = mu, Rcpp::Named("b") = b, Rcpp::Named("hat") = fit, Rcpp::Named("h2") = h2);
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
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=fit,
                            Rcpp::Named("Va")=Va*cxx,
                            Rcpp::Named("Ve")=Ve,
                            Rcpp::Named("h2")=Va*cxx/(Va*cxx+Ve));
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
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("hat")=fit,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("Vb")=vb,
                            Rcpp::Named("Va")=vb*MSx,
                            Rcpp::Named("Ve")=ve);
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
  return Rcpp::List::create(Rcpp::Named("mu")=mu, Rcpp::Named("b")=b, Rcpp::Named("b_LS")=b_ols,
                                        Rcpp::Named("h2")=h2, Rcpp::Named("hat")=fit,
                                        Rcpp::Named("Vb")=vb, Rcpp::Named("Ve")=ve,
                                        Rcpp::Named("LRT")=LRT, Rcpp::Named("PVAL")=PVAL);
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
  return Rcpp::List::create(Rcpp::Named("mu") = MU, Rcpp::Named("b") = B,
                            Rcpp::Named("hat") = fit,
                            Rcpp::Named("vb") = VB, Rcpp::Named("ve") = VE,
                            Rcpp::Named("h2") = h2, Rcpp::Named("MSx") = MSx);
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
  float b0,b1,eM,h2,C,MU=0,VE=0,LR,pj;
  float Pi0 = pi/(1.0f-pi);
  for(int i=0; i<iit; i++){
    C = -0.5f/std::sqrt(ve);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = R::rnorm((X.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb[j]), std::sqrt(ve/(xx[j]+Lmb[j])));
      e1 = e - X.col(j)*(b1-b0);
      e2 = e - X.col(j)*(0-b0);
      LR = Pi0*std::exp(C*(e2.squaredNorm()-e1.squaredNorm()));
      pj = 1.0f/(1.0f+LR);
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
  return Rcpp::List::create(Rcpp::Named("mu") = MU,
                            Rcpp::Named("b") = B, Rcpp::Named("d") = D,
                            Rcpp::Named("hat") = fit,
                            Rcpp::Named("vb") = VB, Rcpp::Named("ve") = VE,
                            Rcpp::Named("h2") = h2, Rcpp::Named("MSx") = MSx);
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
  float b0,b1,eM,h2,C,MU=0,VB=0,VE=0,LR,pj,vg,ve=vy,vb=Sb;
  float Pi0 = pi/(1.0f-pi);
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
      LR = Pi0*std::exp(C*(e2.squaredNorm()-e1.squaredNorm()));
      pj = 1.0f/(1.0f+LR);
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
  return Rcpp::List::create(Rcpp::Named("mu") = MU, Rcpp::Named("b") = B,
                            Rcpp::Named("d") = D, Rcpp::Named("hat") = fit,
                            Rcpp::Named("vb") = VB, Rcpp::Named("ve") = VE,
                            Rcpp::Named("h2") = h2, Rcpp::Named("MSx") = MSx);
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
  return Rcpp::List::create(Rcpp::Named("mu") = MU, Rcpp::Named("b") = B,
                            Rcpp::Named("hat") = fit,
                            Rcpp::Named("vb") = VB, Rcpp::Named("ve") = VE,
                            Rcpp::Named("h2") = h2, Rcpp::Named("MSx") = MSx);
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
  return Rcpp::List::create(Rcpp::Named("mu") = MU, Rcpp::Named("b") = B,
                            Rcpp::Named("hat") = fit,
                            Rcpp::Named("vb") = VB, Rcpp::Named("ve") = VE,
                            Rcpp::Named("h2") = h2, Rcpp::Named("MSx") = MSx);
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
  float b0,b1,b2,eM,h2,C,MU=0,VB=0,VE=0,Pi=0,LR,pj,vg,ve=vy,vb=Sb;
  Eigen::VectorXf d = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf b = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf D = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf B = Eigen::VectorXf::Zero(p);
  Eigen::VectorXf fit(n);
  Eigen::VectorXf e = y.array()-mu;
  Eigen::VectorXf e1(n), e2(n);
  float Lmb=ve/vb;
  float Pi0 = pi/(1.0f-pi);
  for(int i=0; i<iit; i++){
    C = -0.5f/std::sqrt(ve);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = R::rnorm((X.col(j).dot(e)+xx[j]*b0)/(xx[j]+Lmb), std::sqrt(ve/(xx[j]+Lmb)));
      b2 = R::rnorm(0, std::sqrt(ve/(xx[j]+Lmb)));
      e1 = e - X.col(j)*(b1-b0);
      e2 = e - X.col(j)*(0-b0);
      LR = Pi0*std::exp(C*(e2.squaredNorm()-e1.squaredNorm()));
      pj = 1.0f/(1.0f+LR);
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
  return Rcpp::List::create(Rcpp::Named("mu") = MU, Rcpp::Named("b") = B,
                            Rcpp::Named("d") = D, Rcpp::Named("pi") = Pi,
                            Rcpp::Named("hat") = fit, Rcpp::Named("h2") = h2,
                            Rcpp::Named("vb") = VB, Rcpp::Named("ve") = VE,
                            Rcpp::Named("PVAL") = PVAL);
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
  return Rcpp::List::create(Rcpp::Named("mu") = MU, Rcpp::Named("b") = B,
                            Rcpp::Named("d") = D, Rcpp::Named("pi") = Pi,
                            Rcpp::Named("hat") = fit, Rcpp::Named("h2") = h2,
                            Rcpp::Named("vb") = VB, Rcpp::Named("ve") = VE,
                            Rcpp::Named("PVAL") = PVAL);
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
  return Rcpp::List::create(Rcpp::Named("hat") = fit, Rcpp::Named("mu") = MU,
                            Rcpp::Named("b1") = B1, Rcpp::Named("b2") = B2,
                            Rcpp::Named("vb1") = VB1, Rcpp::Named("vb2") = VB2,
                            Rcpp::Named("ve") = VE, Rcpp::Named("h2") = h2);
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
  return Rcpp::List::create(Rcpp::Named("mu") = MU,
                            Rcpp::Named("b1") = B1, Rcpp::Named("d1") = D1, Rcpp::Named("vb1") = VB1,
                                        Rcpp::Named("b2") = B2, Rcpp::Named("d2") = D2, Rcpp::Named("vb2") = VB2,
                                                    Rcpp::Named("ve") = VE, Rcpp::Named("hat") = fit, Rcpp::Named("h2") = h2);
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
  return Rcpp::List::create(Rcpp::Named("hat") = fit, Rcpp::Named("mu") = MU,
                            Rcpp::Named("b1") = B1, Rcpp::Named("b2") = B2,
                            Rcpp::Named("vb1") = VB1, Rcpp::Named("vb2") = VB2,
                            Rcpp::Named("ve") = VE, Rcpp::Named("h2") = h2);
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
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b1")=b1, Rcpp::Named("b2")=b2,
                            Rcpp::Named("Vb1")=vb1, Rcpp::Named("Vb2")=vb2, Rcpp::Named("Ve")=ve,
                                        Rcpp::Named("u1")=u1, Rcpp::Named("u2")=u2,
                                        Rcpp::Named("MSx1")=MSx1, Rcpp::Named("MSx2")=MSx2,
                                        Rcpp::Named("h2")=h2, Rcpp::Named("hat")=fit);
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
  return Rcpp::List::create(Rcpp::Named("b")=b, Rcpp::Named("hat")=fit, Rcpp::Named("e")=fit, Rcpp::Named("MSx")=MSx,
                            Rcpp::Named("vb")=vb, Rcpp::Named("ve")=ve, Rcpp::Named("iG")=iG, Rcpp::Named("h2")=h2);
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
  return Rcpp::List::create(Rcpp::Named("mu")=mu, Rcpp::Named("b")=b,
                            Rcpp::Named("h2")=h2, Rcpp::Named("hat")=fit,
                            Rcpp::Named("Lmb")=Lmb);
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
  return Rcpp::List::create(Rcpp::Named("mu")=mu, Rcpp::Named("b")=b, Rcpp::Named("d")=d,
                                        Rcpp::Named("pi")=Pi, Rcpp::Named("hat")=fit,
                                        Rcpp::Named("Vg")=va*MSx, Rcpp::Named("Va")=va,
                                        Rcpp::Named("Ve")=ve, Rcpp::Named("h2")=h2);
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
  return Rcpp::List::create(Rcpp::Named("mu")=mu, Rcpp::Named("b")=b,
                            Rcpp::Named("h2")=1.0f-vna/vy, Rcpp::Named("e")=e,
                            Rcpp::Named("Lmb")=Lmb, Rcpp::Named("vb")=Vb);
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
  return Rcpp::List::create(Rcpp::Named("mu")=mu, Rcpp::Named("b")=b,
                            Rcpp::Named("h2")=1.0f-vna/vy, Rcpp::Named("e")=e,
                            Rcpp::Named("Lmb")=Lmb, Rcpp::Named("vb")=Vb);
}


// [[Rcpp::export]]
SEXP PEGS(Eigen::MatrixXf Y, // matrix response variables
          Eigen::MatrixXf X, // design matrix of random effects
          int maxit = 100, // maximum number of iterations
          float logtol = -4.0, // convergence tolerance
          float covbend = 1.1, // covariance bending factor
          float covMinEv = 10e-4, // minimum eigenvalue to bend covariance
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
           float covMinEv = 10e-4,
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
      
      if( k>=5 || min_ev < covMinEv ){ 
        vbr.diagonal().array() += bend_inflate[r];
      }
      
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
           float covMinEv = 10e-4,  // minimum eigenvalue to bend
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

// [[Rcpp::export]]
SEXP PEGS_sparse(Eigen::MatrixXf Y, // matrix response variables
                 const Eigen::MappedSparseMatrix<double>& S, // SPARSE design matrix of random effects
                 int maxit = 100, // maximum number of iterations
                 float logtol = -4.0, // convergence tolerance
                 float covbend = 1.1, // covariance bending factor
                 float covMinEv = 10e-4, // minimum eigenvalue to bend covariance
                 int XFA = -1, // number of principal components to fit
                 bool NNC = true){ // non-negative correlations
  
  // Cast input sparse matrix to float for consistency
  Eigen::SparseMatrix<float> X = S.cast<float>();
  
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
  
  // Sum of squares of X (sparse matrix implementation)
  Eigen::MatrixXf XX(p,k);
  for(int i=0; i<p; i++){
    // Square the elements of the sparse column vector
    Eigen::SparseVector<float> x_col_sq = X.col(i).cwiseProduct(X.col(i));
    // Efficiently compute row i of XX
    XX.row(i) = x_col_sq.transpose() * Z;
  }
  
  // Compute Tr(XSX);
  Eigen::MatrixXf XSX(p,k);
  for(int i=0; i<p; i++){
    // X.col(i).transpose() * Z is an efficient sparse-vector dense-matrix multiply
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
  
  // Beta tilde; (efficient sparse-matrix dense-matrix multiply)
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
      // X.col(J).transpose()*e is an efficient sparse-vector dense-matrix multiply
      RHS = (X.col(J).transpose()*e).array() + XX.row(J).array()*b0.transpose().array();
      RHS = RHS.array() *iVe.array();
      b1 = LHS.llt().solve(RHS);
      b.row(J) = b1;
      
      // Update residuals (sparse implementation)
      // This avoids forming a large temporary n0 x k matrix
      Eigen::VectorXf delta_b = b1 - b0;
      Eigen::SparseVector<float> x_col_J = X.col(J);
      for (Eigen::SparseVector<float>::InnerIterator it(x_col_J); it; ++it) {
        // e.row(i) -= x_ij * delta_b' .* z_i'
        e.row(it.index()) -= (it.value() * delta_b.transpose()).cwiseProduct(Z.row(it.index()));
      }
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
      sd = vb.diagonal().array();
      vb.setZero();
      vb.diagonal() = sd.array();
    }else if(XFA>0){
      sd = vb.diagonal().array().sqrt();
      for (int t = 0; t < k; ++t) sd(t) = std::max(sd(t), 1e-12f);
      inv_sd = sd.array().inverse();
      GC = inv_sd.asDiagonal() * vb * inv_sd.asDiagonal();
      eigen_solver.compute(GC);
      V_reduced = eigen_solver.eigenvectors().rightCols(XFA);
      D_reduced_diag = eigen_solver.eigenvalues().tail(XFA);
      GC = V_reduced * D_reduced_diag.asDiagonal() * V_reduced.transpose();
      GC.diagonal().setOnes();
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
  
  // Fitting the model (efficient sparse-matrix dense-matrix multiply)
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
SEXP PEGSZ_sparse(Eigen::MatrixXf Y, // matrix response variables
                  Rcpp::List X_list, // LIST of SPARSE design matrices of random effects
                  int maxit = 100, // maximum number of iterations
                  float logtol = -4.0, // convergence tolerance
                  float covbend = 1.1, // covariance bending factor
                  float covMinEv = 10e-4, // minimum eigenvalue to bend
                  int XFA = -1, // number of principal components to fit
                  bool NNC = true){ // non-negative correlations
  
  // Get input dimensions
  int k = Y.cols(), n0 = Y.rows();
  int n_effects = X_list.size();
  
  // Store input sparse matrices (cast from double to float)
  std::vector<Eigen::SparseMatrix<float>> X_mats;
  Eigen::VectorXi p_vec(n_effects);
  for(int i=0; i<n_effects; ++i){
    const Eigen::MappedSparseMatrix<double> X_in(Rcpp::as<Eigen::MappedSparseMatrix<double>>(X_list[i]));
    X_mats.push_back(X_in.cast<float>());
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
  
  // Pre-compute for each effect: Sum of squares of X and Tr(XSX) (sparse implementation)
  std::vector<Eigen::MatrixXf> XX_list;
  Eigen::MatrixXf TrXSX(n_effects, k);
  Eigen::MatrixXf MSx_mat(n_effects, k);
  
  for(int eff=0; eff<n_effects; ++eff){
    int p = p_vec(eff);
    Eigen::MatrixXf XX(p,k);
    for(int i=0; i<p; i++){
      Eigen::SparseVector<float> x_col_sq = X_mats[eff].col(i).cwiseProduct(X_mats[eff].col(i));
      XX.row(i) = x_col_sq.transpose() * Z;
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
  Eigen::VectorXf iN_var = (n.array().max(1.0f)-1).array().inverse(); // for variance calculation, avoid n=1
  Eigen::VectorXf vy = y.colwise().squaredNorm(); vy = vy.array() * iN_var.array();
  Eigen::VectorXf ve = vy * 0.5;
  Eigen::VectorXf iVe = ve.array().inverse();
  
  std::vector<Eigen::MatrixXf> vb_list(n_effects, Eigen::MatrixXf(k,k));
  std::vector<Eigen::MatrixXf> iG_list(n_effects, Eigen::MatrixXf(k,k));
  for(int eff=0; eff<n_effects; ++eff){
    vb_list[eff] = (ve.array() / MSx_mat.row(eff).transpose().array()).matrix().asDiagonal();
    iG_list[eff] = vb_list[eff].completeOrthogonalDecomposition().pseudoInverse();
  }
  
  // Beta tilde for each effect (sparse-matrix dense-matrix multiply)
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
        
        // Update total residuals sequentially (sparse implementation)
        Eigen::VectorXf delta_b = b1 - b0;
        Eigen::SparseVector<float> x_col_J = X_mats[eff].col(J);
        for (Eigen::SparseVector<float>::InnerIterator it(x_col_J); it; ++it) {
          e.row(it.index()) -= (it.value() * delta_b.transpose()).cwiseProduct(Z.row(it.index()));
        }
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
    // Using iN_mu for intercept update to be consistent with mean calculation
    b0 = b0.array() * iN_mu.array(); 
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

// End


// [[Rcpp::export]]
Eigen::MatrixXf EigenARC(Eigen::MatrixXf X, bool centralizeX = true, int cores = 1){
  // cseweb.ucsd.edu/~saul/papers/nips09_kernel.pdf
  if(cores!=1) Eigen::setNbThreads(cores);
  int p = X.cols(), n = X.rows(); 
  float tmp, Npi=3.1416, theta, J1, Kij, Norm;
  if(centralizeX){
    for(int i=0; i<p; i++){
      tmp = (X.col(i).array()).mean();
      X.col(i) = X.col(i).array()-tmp;}}
  Eigen::MatrixXf XXp = X*X.transpose();
  tmp = 1/(XXp.diagonal().mean()); XXp *= tmp;
  Eigen::VectorXf DiagXXp = XXp.diagonal().array();
  for(int i=0; i<n; i++){ for(int j=i; j<n; j++){ 
    Norm = sqrt(DiagXXp(i)*DiagXXp(j)*1.001);
    theta = acos( XXp(i,j)/Norm);
    J1 = sin(theta) + (Npi-theta)*cos(theta);
    Kij = Norm/Npi*J1;
    XXp(i,j) = Kij*1.0; XXp(j,i) = Kij*1.0;}}
  return XXp;}

// [[Rcpp::export]]
Eigen::MatrixXf EigenGAU(Eigen::MatrixXf X, float phi = 1.0, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores);
  int n = X.rows(); float tmp;
  Eigen::MatrixXf XXp = X*X.transpose();
  for(int i=0; i<n; i++){ for(int j=0; j<n; j++){ if(i>j){
    tmp = sqrt(XXp(i,i) + XXp(j,j) - 2*XXp(i,j));
    XXp(i,j) = tmp*1.0; XXp(j,i) = tmp*1.0;}}};
  for(int i=0; i<n; i++){XXp(i,i) = 0.0;}
  tmp = phi * (-n*(n-1)) / (XXp.colwise().sum()).sum();
  XXp *= tmp; return exp(XXp.array());}

// [[Rcpp::export]]
Eigen::MatrixXf EigenGRM(Eigen::MatrixXf X, bool centralizeZ = true, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores); 
  int p = X.cols(); float tmp;
  if(centralizeZ){
    for(int i=0; i<p; i++){
      tmp = (X.col(i).array()).mean();
      X.col(i) = X.col(i).array()-tmp;}}
  Eigen::MatrixXf XXp = X*X.transpose();
  XXp.diagonal() = XXp.diagonal().array() + 1.0;
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
SEXP EigenEVD(Eigen::MatrixXf A, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores); 
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(A);
  return Rcpp::List::create(Rcpp::Named("U")=es.eigenvectors(),
                            Rcpp::Named("D")=es.eigenvalues());}

// [[Rcpp::export]]
SEXP EigenBDCSVD(Eigen::MatrixXf X, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores);
  Eigen::BDCSVD<Eigen::MatrixXf> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV );
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
    std::shuffle(RGSvec1.begin(), RGSvec1.end(), std::mt19937(numit));
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
    std::shuffle(RGSvec2.begin(), RGSvec2.end(), std::mt19937(numit));
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
          bool ACS = false,
          int NumXFA = 3,
          double R2 = 0.5,
          double gc0 = 0.5, 
          double df0 = 1.0, 
          bool updateMu = false,
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
  std::vector<int> RGSvec(p); int J;
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  
  // Inner RGS
  std::vector<int> InnerRGSvec(k); int ri;
  for(int j=0; j<k; j++){InnerRGSvec[j]=j;}
  
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
    std::shuffle(RGSvec.begin(), RGSvec.end(), std::mt19937(numit));
    std::shuffle(InnerRGSvec.begin(), InnerRGSvec.end(), std::mt19937(numit));
    
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
      
      // Update coefficient      
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
    
    // Proportion-base prior
    if(weight_prior_h2>0){ // Proportion-based prior H2
      for(int i=0; i<k; i++){gs = vb(i,i)*(1-weight_prior_h2) + weight_prior_h2*vbInit(i); vb(i,i) = gs*1.0;}}
    if(weight_prior_gc>0){ // Proportion-based prior GC
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){
        if(i!=j){ GC(i,j) = (1.0-weight_prior_gc)*vb(i,j)/(sqrt(vb(i,i)*vb(j,j))) + gc0*weight_prior_gc;}else{GC(i,j) = 1.0;}}}
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){ if(i!=j){ vb(i,j) = GC(i,j)*sqrt(vb(i,i)*vb(j,j));}}}}else{
        // Once calculation of GC without prior
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}}
    
    // Average covariance structure (HCS+XFA)
    if(ACS){
      gs = (GC.array().sum()-k)/((k*(k-1)))/2.0;
      es.compute(GC);
      UDU = es.eigenvalues()[k-1] * es.eigenvectors().col(k-1) * es.eigenvectors().col(k-1).transpose();
      for(int i=1; i<NumXFA; i++) UDU += es.eigenvalues()[k-i-1] * es.eigenvectors().col(k-i-1) * es.eigenvectors().col(k-i-1).transpose();
      UDU.array() += gs; GC = UDU * 0.5;
      for(int i=0; i<k; i++){ GC(i,i)=1.0; };
    }else if(HCS){ // Heterogeneous Compound Symmetry
      gs = 0.0;
      for(int i=0; i<k; i++){
        for(int j=0; j<k; j++){
          if(i>j){gs += GC(i,j);}}}
      gs = gs/((k*(k-1))/2);
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){ 
        if(i!=j){ GC(i,j) =  gs*1.0;}else{ GC(i,j) = 1.0; }}}
    }else if(XFA){ // Extended Factor Analytics
      es.compute(GC);
      UDU = es.eigenvalues()[k-1] * es.eigenvectors().col(k-1) * es.eigenvectors().col(k-1).transpose();
      for(int i=1; i<NumXFA; i++) UDU += es.eigenvalues()[k-i-1] * es.eigenvectors().col(k-i-1) * es.eigenvectors().col(k-i-1).transpose();
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
    
    // Update intercept
    if(updateMu){
      b0 = e.colwise().sum();
      b0 = b0.array() * iN.array();
      for(int i=0; i<k; i++){ mu(i) += b0(i);
        e.col(i) = (e.col(i).array()-b0(i)).array() * Z.col(i).array();} 
    }
    
    
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
           bool ACS = false,
           int NumXFA = 3,
           float R2 = 0.5,
           float gc0 = 0.5, 
           float df0 = 1.0, 
           bool updateMu = false,
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
  std::vector<int> RGSvec(p); int J;
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  
  // Inner RGS
  std::vector<int> InnerRGSvec(k); int ri;
  for(int j=0; j<k; j++){InnerRGSvec[j]=j;}
  
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
    std::shuffle(RGSvec.begin(), RGSvec.end(), std::mt19937(numit));
    std::shuffle(InnerRGSvec.begin(), InnerRGSvec.end(), std::mt19937(numit));
    
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
    
    // Proportion-base prior
    if(weight_prior_h2>0){ // Proportion-based prior H2
      for(int i=0; i<k; i++){gs = vb(i,i)*(1-weight_prior_h2) + weight_prior_h2*vbInit(i); vb(i,i) = gs*1.0;}}
    if(weight_prior_gc>0){ // Proportion-based prior GC
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){
        if(i!=j){ GC(i,j) = (1.0-weight_prior_gc)*vb(i,j)/(sqrt(vb(i,i)*vb(j,j))) + gc0*weight_prior_gc;}else{GC(i,j) = 1.0;}}}
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){ if(i!=j){ vb(i,j) = GC(i,j)*sqrt(vb(i,i)*vb(j,j));}}}}else{
        // Once calculation of GC without prior
        for(int i=0; i<k; i++){for(int j=0; j<k; j++){GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}}
    
    // Average covariance structure (HCS+XFA)
    if(ACS){
      gs = (GC.array().sum()-k)/((k*(k-1)))/2.0;
      es.compute(GC);
      UDU = es.eigenvalues()[k-1] * es.eigenvectors().col(k-1) * es.eigenvectors().col(k-1).transpose();
      for(int i=1; i<NumXFA; i++) UDU += es.eigenvalues()[k-i-1] * es.eigenvectors().col(k-i-1) * es.eigenvectors().col(k-i-1).transpose();
      UDU.array() += gs; GC = UDU * 0.5;
      for(int i=0; i<k; i++){ GC(i,i)=1.0; };
    }else if(HCS){ // Heterogeneous Compound Symmetry
      gs = 0.0;
      for(int i=0; i<k; i++){
        for(int j=0; j<k; j++){
          if(i>j){gs += GC(i,j);}}}
      gs = gs/((k*(k-1))/2);
      for(int i=0; i<k; i++){for(int j=0; j<k; j++){ 
        if(i!=j){ GC(i,j) =  gs*1.0;}else{ GC(i,j) = 1.0; }}}
    }else if(XFA){ // Extended Factor Analytics
      es.compute(GC);
      UDU = es.eigenvalues()[k-1] * es.eigenvectors().col(k-1) * es.eigenvectors().col(k-1).transpose();
      for(int i=1; i<NumXFA; i++) UDU += es.eigenvalues()[k-i-1] * es.eigenvectors().col(k-i-1) * es.eigenvectors().col(k-i-1).transpose();
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
    
    // Update intercept
    if(updateMu){
      b0 = e.colwise().sum();
      b0 = b0.array() * iN.array();
      for(int i=0; i<k; i++){ mu(i) += b0(i);
        e.col(i) = (e.col(i).array()-b0(i)).array() * Z.col(i).array();} 
    }    
    
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
SEXP mrr_svd(Eigen::MatrixXf Y, Eigen::MatrixXf W){
  
  // Start setup
  int maxit = 500;
  float tol = 10e-9;
  
  // Gather basic info
  int k = Y.cols(), n0 = Y.rows(), m = W.cols();
  
  // Center X
  Rcpp::Rcout << "Centering marker score matrix\n";
  Eigen::VectorXf xx = W.colwise().mean();
  for(int i=0; i<m; i++){ W.col(i) = W.col(i).array() - xx(i);}
  
  // Single value decomposition
  Rcpp::Rcout << "SVD of marker scores\n";
  Eigen::BDCSVD<Eigen::MatrixXf> svd(W, Eigen::ComputeThinU | Eigen::ComputeThinV );
  Eigen::MatrixXf V = svd.matrixV();
  Eigen::MatrixXf X = svd.matrixU() * svd.singularValues().array().matrix().asDiagonal();
  int p = X.cols();
  
  // Incidence matrix Z
  Eigen::MatrixXf Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0;Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  
  // Count observations per trait
  Eigen::VectorXf n = Z.colwise().sum();
  Eigen::VectorXf iN = n.array().inverse();
  
  // Centralize y
  Rcpp::Rcout << "Centering Y\n";
  Eigen::VectorXf mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXf y(n0,k);
  for(int i=0; i<k; i++){
    y.col(i) = (Y.col(i).array()-mu(i)).array() * Z.col(i).array();
  }
  
  // Sum of squares of X
  Rcpp::Rcout << "Computing diagonal elements of Z'Z\n";
  Eigen::MatrixXf XX(p,k);
  for(int i=0; i<p; i++){ XX.row(i) = X.col(i).array().square().matrix().transpose() * Z;}
  // Compute Tr(XSX);
  Eigen::MatrixXf XSX(p,k);
  for(int i=0; i<p; i++){ XSX.row(i) = XX.row(i).transpose().array()*iN.array() - 
    ((X.col(i).transpose()*Z).transpose().array()*iN.array()).square();}
  Eigen::VectorXf MSx = XSX.colwise().sum();
  Eigen::VectorXf TrXSX = n.array()*MSx.array();
  
  Rcpp::Rcout << "Set starting values for coefficients and variances\n";
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
  Eigen::VectorXf TrDinvXSX(k);
  Eigen::MatrixXf Dinv(p,k);
  for(int i=0; i<k; i++){ XSX.col(i) = XSX.col(i).array() * n(i); }
  // Initialize coefficient matrices
  Eigen::MatrixXf LHS(k,k);
  Eigen::VectorXf RHS(k);
  Eigen::MatrixXf b = Eigen::MatrixXf::Zero(p,k);
  Eigen::VectorXf b0(k), b1(k);
  Eigen::MatrixXf e(n0,k); e = y*1.0;
  // Bending and convergence control
  Eigen::MatrixXf A = vb*1.0, GC(k,k);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> EVDofA(A);
  Eigen::MatrixXf beta0(p,k), vb0(k,k);
  Eigen::VectorXf CNV1(maxit),CNV2(maxit),CNV3(maxit), ve0(k), h20(k);
  float cnv = 10.0, MinDVb, inflate;
  int numit = 0;
  float logtol = log10(tol);
  
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
  Eigen::MatrixXf hat = X * b;
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i);}
  Eigen::MatrixXf beta = V*b;
  
  // Correlations
  Rcpp::Rcout << "Estimating correlations\n";
  for(int i=0; i<k; i++){
    for(int j=0; j<k; j++){
      GC(i,j)=vb(i,j)/(sqrt(vb(i,i)*vb(j,j)));}}
  
  // Resize convergence vectors
  Rcpp::Rcout << "Convergence statistics\n";
  Eigen::VectorXf CNV1b(numit),CNV2b(numit),CNV3b(numit);
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
  return NullModelOutput;}

// [[Rcpp::export]]
SEXP MLM(Eigen::MatrixXf Y, Eigen::MatrixXf X, Eigen::MatrixXf Z,
         int maxit = 500, float logtol = -8, int cores = 1, 
         bool verb = false, float df0 = 1.1){
  
  // Basic info
  if(cores!=1) Eigen::setNbThreads(cores);
  int k = Y.cols(), n0 = Y.rows(), f = X.cols(), p = Z.cols();
  
  // Incidence matrix W
  Eigen::MatrixXf W(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        W(i,j) = 0.0;
        Y(i,j) = 0.0;
      }else{ W(i,j) = 1.0;}}}
  Eigen::VectorXf n = W.colwise().sum();
  Eigen::VectorXf iN = (n.array()-f).inverse();
  
  // Compute SY, SZ
  Eigen::MatrixXf y(n0,k), WX(n0,f), MU(f,k), iXX(f,f);
  for(int i=0; i<k; i++){
    for(int j=0; j<f; j++){WX.col(j)=X.col(j).array()*W.col(i).array();}
    iXX = (WX.transpose()*WX).inverse();
    MU.col(i) = iXX * WX.transpose()*Y.col(i);
    y.col(i) = (Y.col(i)-WX*MU.col(i) ).array()*W.col(i).array(); }
  iXX = (X.transpose()*X).inverse();
  for(int j=0; j<p; j++){ Z.col(j) = (Z.col(j) - X*(iXX*X.transpose()*Z.col(j))).array(); }
  
  // Sum of squares of Z
  Eigen::MatrixXf ZZ(p,k); 
  for(int i=0; i<p; i++){ ZZ.row(i) = Z.col(i).array().square().matrix().transpose() * W;}
  Eigen::VectorXf TrZSZ = ZZ.colwise().sum().array();
  
  // Initialize coefficient matrices
  Eigen::MatrixXf LHS(k,k);
  Eigen::VectorXf RHS(k);
  Eigen::MatrixXf b = Eigen::MatrixXf::Zero(p,k);
  Eigen::VectorXf b0(k), b1(k);
  Eigen::MatrixXf e(n0,k); e = y*1.0;
  
  // Variances
  Eigen::VectorXf vy = y.colwise().squaredNorm(); vy=vy.array()*iN.array();
  Eigen::VectorXf ve = vy * 0.5;
  Eigen::VectorXf iVe = ve.array().inverse();
  Eigen::MatrixXf vb(k,k), TildeHat(k,k);
  vb = (ve.array()/ (TrZSZ.array()*iN.array())  ).matrix().asDiagonal();
  Eigen::MatrixXf iG = vb.inverse();
  Eigen::VectorXf h2 = 1 - ve.array()/vy.array();
  Eigen::MatrixXf tilde = Z.transpose() * y;
  
  // Bending
  Eigen::MatrixXf A = vb*1.0;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> EVDofA(A);
  float MinDVb, inflate = 0.0;
  
  // Prior for stability
  Eigen::MatrixXf Sb = vb*df0;
  Eigen::VectorXf Se = ve*df0;
  Eigen::VectorXf iNp = (n.array()+df0-f).inverse();
  
  // RGS
  std::vector<int> RGSvec(p); int J;
  for(int j=0; j<p; j++){RGSvec[j]=j;}
  
  // Convergence control
  Eigen::MatrixXf beta0(p,k);
  float cnv = 10.0;
  int numit = 0;
  
  // Loop
  while(numit<maxit){
    
    // Store coefficients pre-iteration
    beta0 = b*1.0;
    
    // Randomized Gauss-Seidel loop
    std::shuffle(RGSvec.begin(), RGSvec.end(), std::mt19937(numit));
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
    if( MinDVb < 0.001 ){if(abs(MinDVb*1.1)>inflate) inflate = abs(MinDVb*1.1);}
    vb.diagonal().array()+=inflate; 
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
  Eigen::MatrixXf hat = Z*b;
  for(int i=0; i<k; i++){ hat.col(i) = ( X * MU.col(i) + hat.col(i)).array(); }
  
  // Genetic correlations
  Eigen::MatrixXf GC(k,k);
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
  return OutputList;
  
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
  while(numit<maxit){
    beta0 = b*1.0;
    std::shuffle(RGSvec.begin(), RGSvec.end(), std::mt19937(numit));
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
  while(numit<maxit){
    beta01 = b_1*1.0; beta02 = b_2*1.0;
    std::shuffle(RGSvec1.begin(), RGSvec1.end(), std::mt19937(numit));
    std::shuffle(RGSvec2.begin(), RGSvec2.end(), std::mt19937(numit));
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
                            Rcpp::Named("hat")=hat);}

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
  while(numit<maxit){
    beta0 = b*1.0;
    std::shuffle(RGSvec.begin(), RGSvec.end(), std::mt19937(numit));
    for(int j=0; j<p; j++){
      J = RGSvec[j]; if(XX[J]>0.00001){ b0 = b[J]*1.0;
        b1 = (e.transpose()*X.col(J)+XX(J)*b0)/(XX[J]+lambda);
        e = e - X.col(J)*(b1-b0); b[J] = b1*1.0;}else{ b[J] = 0.0;}}
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
  while(numit<maxit){
    beta01 = b_1*1.0; beta02 = b_2*1.0;
    std::shuffle(RGSvec1.begin(), RGSvec2.end(), std::mt19937(numit));
    std::shuffle(RGSvec1.begin(), RGSvec2.end(), std::mt19937(numit));
    for(int j=0; j<p1; j++){
      J = RGSvec1[j]; if(XX1[J]>0.00001){ b0 = b_1[J]*1.0;
        b1 = (e.transpose()*X1.col(J)+XX1(J)*b0)/(XX1[J]+lambda1);
        e = e - X1.col(J)*(b1-b0); b_1[J] = b1*1.0; }else{ b_1[J] = 0.0;}}
    for(int j=0; j<p2; j++){
      J = RGSvec2[j]; if(XX2[J]>0.00001){ b0 = b_2[J]*1.0;
        b1 = (e.transpose()*X2.col(J)+XX2(J)*b0)/(XX2[J]+lambda2);
        e = e - X2.col(J)*(b1-b0); b_2[J] = b1*1.0;}else{ b_2[J] = 0.0;}}
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
  while(numit<maxit){  beta0 = b*1.0;
    std::shuffle(RGSvec.begin(), RGSvec.end(), std::mt19937(numit));
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
  for(int i=0; i<(Y.cols()); i++){ G.col(i) = G.col(i).array() - G.col(i).mean(); }
  Eigen::VectorXf vg=G.colwise().squaredNorm(); vg/=(Y.rows()); vg=vg.array().sqrt();
  for(int i=0; i<(Y.cols()); i++){ G.col(i) = G.col(i).array()/vg(i); }
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
  while(numit<maxit){
    beta0 = b*1.0;
    std::shuffle(RGSvec.begin(), RGSvec.end(), std::mt19937(numit));
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

// [[Rcpp::export]]
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
SEXP ZSEMF(Eigen::MatrixXf Y, Eigen::MatrixXf X, int npc = 0){
  int k = Y.cols(), N = Y.rows();
  Eigen::MatrixXf BETA = ZFUVBETA(Y,X);
  Eigen::MatrixXf G = X*BETA.bottomRows(X.cols());
  Eigen::VectorXf h2 = BETA.row(0).array();
  Eigen::BDCSVD<Eigen::MatrixXf> svd(G, Eigen::ComputeThinU | Eigen::ComputeThinV );
  if(npc<0) npc = round(2*sqrt(svd.matrixU().cols()));
  if(npc==0) npc += svd.matrixU().cols();
  Eigen::MatrixXf Z = (svd.matrixU() * svd.singularValues().matrix().asDiagonal()).leftCols(npc);
  Eigen::MatrixXf Coef = ZFUVBETA(Y,Z);
  // Hat and H2
  Eigen::MatrixXf beta_final = BETA.bottomRows(X.cols()) * svd.matrixV().leftCols(npc) * Coef.bottomRows( Z.cols()); 
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
SEXP YSEMF(Eigen::MatrixXf Y, Eigen::MatrixXf X, int npc = -1){
  int k = Y.cols(), N = Y.rows();
  Eigen::MatrixXf BETA = ZFUVBETA(Y,X);
  Eigen::MatrixXf G = X*BETA.bottomRows(X.cols());
  Eigen::BDCSVD<Eigen::MatrixXf> svd(G, Eigen::ComputeThinU | Eigen::ComputeThinV );
  if(npc<0) npc = round(2*sqrt(svd.matrixU().cols()));
  if(npc==0) npc += svd.matrixU().cols();
  Eigen::MatrixXf Z = (svd.matrixU() * svd.singularValues().matrix().asDiagonal()).leftCols(npc);
  Eigen::MatrixXf ALPHA = ZFUVBETA(Y,Z);
  Eigen::MatrixXf beta_Fa = BETA.bottomRows(X.cols()) * svd.matrixV().leftCols(npc) * ALPHA.bottomRows( Z.cols()); 
  G = X*beta_Fa;
  Eigen::MatrixXf beta_Xd = ZFUVBETA(Y-G,X);
  Eigen::MatrixXf beta_FaXd = beta_Fa+beta_Xd.bottomRows( X.cols());
  G = X*beta_FaXd;
  Eigen::MatrixXf hat = G * 1.0;
  Eigen::VectorXf mu = beta_Xd.row(1).array();
  for(int i=0; i<k; i++){ hat.col(i) = hat.col(i).array() + mu(i); }
  Eigen::VectorXf h2 = ALPHA.row(0).array() + beta_Xd.row(0).array();
  for(int i=0; i<k; i++){ G.col(i) = G.col(i).array() - G.col(i).mean(); }
  Eigen::VectorXf vg = G.colwise().squaredNorm(); vg /= N; vg = vg.array().sqrt();
  for(int i=0; i<k; i++){ G.col(i) = G.col(i).array() / vg(i); }
  Eigen::MatrixXf GC = (G.transpose()*G)/N;
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=beta_FaXd,
                            Rcpp::Named("hat")=hat,
                            Rcpp::Named("h2")=h2,
                            Rcpp::Named("GC")=GC);}

// [[Rcpp::export]]
Eigen::MatrixXf EigenArcZ( Eigen::MatrixXf Zfndr, Eigen::MatrixXf Zsamp, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores);  
  int p = Zfndr.cols(), nf = Zfndr.rows(), ns = Zsamp.rows();
  // Centralize matrices to create relationship matrix
  Eigen::VectorXf MeanColumnZfndr = Zfndr.colwise().mean();
  for(int i=0; i<p; i++){
    Zfndr.col(i) = Zfndr.col(i).array()-MeanColumnZfndr(i);
    Zsamp.col(i) = Zsamp.col(i).array()-MeanColumnZfndr(i);}
  Eigen::MatrixXf Kff = Zfndr * Zfndr.transpose();
  Eigen::MatrixXf Kfs = Zfndr * Zsamp.transpose();
  float Kscalar, tmp, NormProd, Npi = 3.14159;
  Eigen::VectorXf DiagKff = Kff.diagonal().array();
  Eigen::VectorXf DiagKss = (Zsamp.cwiseProduct(Zsamp)).rowwise().sum();
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
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(Kff);
  Eigen::MatrixXf L = es.eigenvectors() * es.eigenvalues().array().rsqrt().matrix().asDiagonal();
  return Kfs.transpose() * L;
}

// [[Rcpp::export]]
Eigen::MatrixXf EigenGauZ( Eigen::MatrixXf Zfndr, Eigen::MatrixXf Zsamp, float phi = 1.0, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores);  
  int p = Zfndr.cols(), nf = Zfndr.rows(), ns = Zsamp.rows();
  // Centralize matrices to create relationship matrix
  Eigen::MatrixXf Kff = Zfndr * Zfndr.transpose();
  Eigen::MatrixXf Kfs = Zfndr * Zsamp.transpose();
  Eigen::VectorXf DiagKff = Kff.diagonal().array();
  Eigen::VectorXf DiagKss = (Zsamp.cwiseProduct(Zsamp)).rowwise().sum();
  float tmp;
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
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(Kff);
  Eigen::MatrixXf L = es.eigenvectors() * es.eigenvalues().array().rsqrt().matrix().asDiagonal();
  return Kfs.transpose() * L;
}

// [[Rcpp::export]]
Eigen::MatrixXf K2X(Eigen::MatrixXf K, float MinEV = 1e-8, int cores = 1){
  if(cores!=1) Eigen::setNbThreads(cores);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(K);
  Eigen::BDCSVD<Eigen::MatrixXf> svd(K, Eigen::ComputeThinU | Eigen::ComputeThinV );
  int NPC = 0; Eigen::VectorXf D = svd.singularValues().array();
  for(int i=0; i< D.size(); i++){ if(D[i]>MinEV) NPC +=1; };
  return svd.matrixU().leftCols(NPC) * svd.singularValues().head(NPC).matrix().asDiagonal();}

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

// SUPPORTING FUNCTIONS 10/29/2025

Eigen::VectorXf standardize_vector(Eigen::VectorXf vector) {
  float mean = vector.mean();
  float sum_of_squares = (vector.array() - mean).square().sum();
  float std_dev = std::sqrt(sum_of_squares / (vector.size() - 1));
  Eigen::VectorXf standardized_vector = (vector.array() - mean) / std_dev;
  return standardized_vector;
}

Eigen::MatrixXf standardize_matrix(Eigen::MatrixXf matrix) {
  int k = matrix.cols(), n=matrix.rows();
  Eigen::MatrixXf N_mat(n,k);
  for(int i=0; i<k; i++){ N_mat.col(i) = standardize_vector(matrix.col(i).array()).array();};
  return N_mat;
}

Eigen::MatrixXf MaxRow(Eigen::MatrixXf matrix) {
  int n = matrix.rows(), p = matrix.cols(); 
  float tmp;
  for(int i=0; i<n; i++){
    tmp = matrix.row(i).array().maxCoeff(); 
    for(int j=0; j<p; j++){
      if(matrix(i,j)==tmp){
        matrix(i,j)=1.0;
      }else{
        matrix(i,j)=0.0;}
    }
  }
  return matrix;
}

Eigen::VectorXf uniqueValues(const Eigen::VectorXf& inputVector) {
  std::set<int> uniqueSet(inputVector.data(), inputVector.data() + inputVector.size());
  Eigen::VectorXf uniqueVector(uniqueSet.size());
  int i = 0;
  for (int value : uniqueSet) { uniqueVector(i++) = value;}
  return uniqueVector;
}

Eigen::MatrixXf DropIncZero(Eigen::MatrixXf matrix) {
  int k = matrix.cols(), n=matrix.rows();
  Eigen::VectorXf Counts = matrix.colwise().sum();
  int non_zero_cols = 0;
  for(int i=0; i<k; i++){ if(Counts(i)>0) non_zero_cols ++;  }
  Eigen::VectorXi ind(non_zero_cols);
  int counts_included_columns = 0;
  for(int i=0; i<k; i++){ 
    if(Counts(i)>0){
      ind(counts_included_columns) = i;
      counts_included_columns ++;}}
  return matrix(Eigen::all,ind);
}

// [[Rcpp::export]]
Eigen::MatrixXf ClusterBlup(Eigen::MatrixXf Y, Eigen::MatrixXf X, float lambda = 2.0){
  // Basic info
  int k = Y.cols(), n0 = Y.rows(), p = X.cols();
  // Incidence matrix Z
  Eigen::MatrixXf Z(n0,k);
  for(int i=0; i<n0; i++){
    for(int j=0; j<k; j++){
      if(std::isnan(Y(i,j))){
        Z(i,j) = 0.0; Y(i,j) = 0.0;
      }else{ Z(i,j) = 1.0;}}}
  // Count observations per trait
  Eigen::VectorXf n = Z.colwise().sum();
  Eigen::VectorXf iN = n.array().inverse();
  // Centralize y
  Eigen::VectorXf mu = Y.colwise().sum();
  mu = mu.array() * iN.array();
  Eigen::MatrixXf y(n0,k);
  for(int i=0; i<k; i++){y.col(i)=(Y.col(i).array()-mu(i)).array() * Z.col(i).array();}
  // Sum of squares of X
  Eigen::MatrixXf XY = y * X;
  Eigen::MatrixXf ZX = Z * X;
  ZX = ZX.array() + lambda;
  Eigen::MatrixXf DEV = XY.cwiseProduct(ZX.cwiseInverse());
  Eigen::VectorXf INTERCEPT = ((X.transpose()*X).llt().solve(X.transpose()*mu)).array();
  //Eigen::VectorXf INTERCEPT = (X.transpose()*mu).array() / Z.colwise().sum().array();
  for(int i=0; i<p; i++){ 
    DEV.col(i) = DEV.col(i).array() + INTERCEPT(i); 
    for(int j=0; j<n0; j++){ if(XY(j,i)==0) DEV(j,i) = std::nanf(""); }}
  return(DEV);
}

//  [[Rcpp::export]]
Eigen::MatrixXf IncMatrix(Eigen::VectorXf x){
  int n = x.size();
  Eigen::VectorXf set = uniqueValues(x);
  int p = set.size();
  Eigen::MatrixXf X(n,p);
  for(int i=0; i<n; i++){for(int j=0; j<p; j++){if(x(i)==set(j)){X(i,j)=1.0;}else{X(i,j)=0.0;}}}
  return(X);
}

// [[Rcpp::export]]
Eigen::MatrixXf EM_recluster(Eigen::MatrixXf Y, Eigen::MatrixXf C, int rounds = 3){
  int k = Y.cols(), n = Y.rows(), p = C.cols(); float tmp;
  // Normalize Y
  Y = standardize_matrix(Y);
  // Normalized Centroids
  Eigen::MatrixXf Centroids = standardize_matrix(Y*C);
  // Corr
  Eigen::MatrixXf covar = MaxRow( Y.transpose() * Centroids); 
  // Reclassify
  for(int its=0; its<rounds; its++){
    Centroids = Y*covar;
    covar = MaxRow( Y.transpose() * Centroids);
  }
  return DropIncZero(covar);
}

// [[Rcpp::export]]
Eigen::MatrixXf Get_Cluster_Corr(Eigen::MatrixXf Y, Eigen::MatrixXf C){
  int k = Y.cols(), n = Y.rows(), p = C.cols(); float tmp;
  Y = standardize_matrix(Y);
  Eigen::MatrixXf Centroids = standardize_matrix(Y*C);
  Eigen::MatrixXf corr = Centroids.transpose() * Centroids; 
  corr = corr.array() / n; 
  for(int i=0; i<p; i++){ corr(i,i)=1.0;}
  return corr;
}

