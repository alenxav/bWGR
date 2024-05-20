#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP KMUP(NumericMatrix X, NumericVector b, NumericVector d, NumericVector xx, NumericVector e, NumericVector L, double Ve, double pi){  
  int p = X.ncol();
  NumericVector e1 = e+0;
  NumericVector e2 = e+0;
  double b0,b1,b2,cj,dj,pj;
  double C = -0.5/sqrt(Ve);
  //
  for(int j=0; j<p; j++){
    b0 = b[j];
    b1 = R::rnorm( (sum(X(_,j)*e)+xx[j]*b0)/(xx[j]+L[j]), sqrt(Ve/(xx[j]+L[j]))  );
    b2 = R::rnorm( 0, sqrt(Ve/(xx[j]+L[j]))  );
    e1 = e-X(_,j)*(b1-b0); // with marker
    if(pi>0){
      e2 = e-X(_,j)*(b2-b0); // without marker
      //
      cj = (1-pi)*exp(C*sum(e1*e1));
      dj = (pi)*exp(C*sum(e2*e2));
      pj = cj/(cj+dj);
      //
      if(R::rbinom(1,pj)==1){
        b[j] = b1;
        d[j] = 1;
        e = e1;
      }else{
        b[j] = b2;
        d[j] = 0;
        e = e2;
      }
    }else{
      // WO Variable Selection
      d[j] = 1;
      b[j] = b1;
      e = e1;
    }
  }
 return List::create(Named("b") = b,Named("d") = d,Named("e") = e);
}

// [[Rcpp::export]]
SEXP KMUP2(NumericMatrix X, NumericVector Use, NumericVector b,  NumericVector d, NumericVector xx, NumericVector E, NumericVector L, double Ve, double pi){  
  int p = X.ncol();
  int n0 = X.nrow();
  int n = Use.size();
  NumericVector e1 = E+0;
  NumericVector e2 = E+0;
  double b0,b1,b2,cj,dj,pj;
  double C = -0.5/sqrt(Ve);
  //
  double bg = n0/n;
  NumericVector e0(n);
  NumericVector H(n);
  //
  for(int k=0; k<n; k++){
    e0[k] = E[Use[k]];
  }
  //
  for(int j=0; j<p; j++){
    //
    for(int x=0; x<n; x++){
      H[x] = X(Use[x],j);
    }
    //
    b0 = b[j];
    b1 = R::rnorm((sum(H*e0)+b0)/(xx(j)*bg+L(j)),sqrt(Ve/(xx(j)*bg+L(j))));
    b2 = R::rnorm( 0 ,sqrt(Ve/(xx(j)*bg+L(j))));
    e1 = e0 - H*(b1-b0);
    //
    if(pi>0){
      e2 = e0 - H*(b2-b0);
      //
      cj = (1-pi)*exp(C*sum(e1*e1));
      dj = (pi)*exp(C*sum(e2*e2));
      pj = cj/(cj+dj);
      //
      if(R::rbinom(1,pj)==1){
        b[j] = b1;
        d[j] = 1;
        e0 = e1;
      }else{
        b[j] = b2;
        d[j] = 0;
        e0 = e2;
      }
      //
    }else{
      d[j] = 1;
      b[j] = b1;
      e0 = e1;
    }
    //
  }  
  return List::create(Named("b") = b,
                      Named("d") = d,
                      Named("e") = e0);
}

// [[Rcpp::export]]
SEXP emBA(NumericVector y, NumericMatrix gen, double df = 10, double R2 = 0.5){
  int it = 200;
  int p = gen.ncol();
  int n = gen.nrow();
  double ve = 1;
  NumericVector d(p);
  NumericVector b(p);
  NumericVector vb = b+1;
  NumericVector Lmb = ve/vb;
  double vy = var(y);
  NumericVector xx(p);
  NumericVector vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(gen(_,i)*gen(_,i));
    vx[i] = var(gen(_,i));}
  double MSx = sum(vx);
  double Sb = R2*(df+2)*vy/MSx;
  double Se = (1-R2)*(df+2)*vy;
  double mu = mean(y);
  NumericVector e = y-mu;
  double b0,b1,eM,h2;
  for(int i=0; i<it; i++){
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]);
      e = e-gen(_,j)*(b1-b0);
      b[j] = b1;
      vb[j] = (Sb+b[j]*b[j])/(df+1);
      e = e - gen(_,j)*(b1-b0);}
    ve = (sum(e*e)+Se)/(n+df);
    Lmb = ve/vb;
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;}
  h2 = 1-ve/vy;
  NumericVector fit(n);
  for(int k=0; k<n; k++){fit[k] = sum(gen(k,_)*b)+mu;}
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("Vb") = vb,
                      Named("Ve") = ve,
                      Named("h2") = h2);}

// [[Rcpp::export]]
SEXP emBB(NumericVector y, NumericMatrix gen, double df = 10, double R2 = 0.5, double Pi = 0.75){
  int it = 200;
  int p = gen.ncol();
  int n = gen.nrow();
  double ve = 1;
  NumericVector d(p);
  NumericVector b(p);
  NumericVector vb = b+1;
  NumericVector Lmb = ve/vb;
  double vy = var(y);
  if(Pi>0.5){Pi = 1-Pi;} 
  NumericVector xx(p);
  NumericVector vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(gen(_,i)*gen(_,i));
    vx[i] = var(gen(_,i));}
  double MSx = sum(vx)*Pi;
  double Sb = R2*(df+2)*vy/MSx;
  double Se = (1-R2)*(df+2)*vy;
  double mu = mean(y);
  NumericVector e = y-mu;
  NumericVector e1(n);
  NumericVector e2(n);
  double b0,b1,LR,eM,h2,C;
  double Pi0 = (1-Pi)/Pi;
  for(int i=0; i<it; i++){
    C = -0.5/sqrt(ve);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]);
      e1 = e-gen(_,j)*(b1-b0);
      e2 = e-gen(_,j)*(0-b0);
      LR = Pi0*exp(C*(sum(e2*e2)-sum(e1*e1)));
      d[j] = (1/(1+LR));
      b[j] = b1*d[j];
      vb[j] = (Sb+b[j]*b[j])/(df+1);
      e = e - gen(_,j)*(b[j]-b0);}
    ve = (sum(e*e)+Se)/(n+df);
    Lmb = ve/vb;
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;}
  h2 = 1-ve/vy;
  NumericVector fit(n);
  for(int k=0; k<n; k++){fit[k] = sum(gen(k,_)*b)+mu;}
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("d") = d,
                      Named("hat") = fit,
                      Named("Vb") = vb,
                      Named("Ve") = ve,
                      Named("h2") = h2);}

// [[Rcpp::export]]
SEXP emBC(NumericVector y, NumericMatrix gen, double df = 10, double R2 = 0.5, double Pi = 0.75){
  int it = 200;
  int p = gen.ncol();
  int n = gen.nrow();
  NumericVector d(p);
  NumericVector b(p);
  double vy = var(y);
  if(Pi>0.5){Pi = 1-Pi;} 
  NumericVector xx(p);
  NumericVector vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(gen(_,i)*gen(_,i));
    vx[i] = var(gen(_,i));
  }
  double MSx = sum(vx)*Pi*(1-Pi);
  double Sa = R2*(df+2)*vy/MSx;
  double Se = (1-R2)*(df+2)*vy;
  double mu = mean(y);
  NumericVector e = y-mu;
  NumericVector e1(n);
  NumericVector e2(n);
  double ve = Sa;
  double va = Se;
  double Lmb = ve/va;
  double b0,b1,LR,eM,h2,C;
  double Pi0 = (1-Pi)/Pi;
  for(int i=0; i<it; i++){
    C = -0.5/sqrt(ve);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);
      e1 = e-gen(_,j)*(b1-b0);
      e2 = e-gen(_,j)*(0-b0);
      LR = Pi0*exp(C*(sum(e2*e2)-sum(e1*e1)));
      d[j] = (1/(1+LR));
      b[j] = b1*d[j];
      e = e - gen(_,j)*(b[j]-b0);
    }
    ve = (sum(e*e)+Se)/(n+df);
    va = (sum(b*b)+Sa)/(p+df)/(mean(d)-Pi);
    Lmb = ve/va;
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
  }
  h2 = 1-ve/vy;
  NumericVector fit(n);
  for(int k=0; k<n; k++){fit[k] = sum(gen(k,_)*b)+mu;}
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
SEXP emDE(NumericVector y, NumericMatrix gen, double R2 = 0.5){
  // Convergence criteria
  int maxit = 300;
  double tol = 10e-6;
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  double b0,eM;
  double mu = mean(y);
  NumericVector e = y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){
    xx[k] = sum(gen(_,k)*gen(_,k));
    if(xx[k]==0) xx[k]=0.1;}
  NumericVector vx(p);
  for(int k=0; k<p; k++){vx[k] = var(gen(_,k));}
  double cxx = sum(vx)*(1-R2)/R2;
  // Regulation coefficients
  double Ve;
  NumericVector Vb(p);
  NumericVector b(p);
  NumericVector Lmb = p+cxx;
  double b1;
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      // Ordinary Least Square
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(Lmb(j)+xx(j));
      b[j] = b1;
      // Residuals update
      e = e-gen(_,j)*(b1-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Variance components
    Ve = sum(e*y)/(n-1);
    Vb = b*b+(Ve/(xx+Lmb+0.0001));
    Lmb = sqrt(cxx*Ve/Vb);
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}
  }
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  // Output
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("Vb")=Vb,
                      Named("Ve")=Ve,
                      Named("h2")=sum(Vb)/(sum(Vb)+Ve));
}

// [[Rcpp::export]]
SEXP emRR(NumericVector y, NumericMatrix gen, double df = 10, double R2 = 0.5){
  int it = 200;
  int p = gen.ncol();
  int n = gen.nrow();
  NumericVector xx(p);
  NumericVector vx(p);
  for(int k=0; k<p; k++){
    xx[k] = sum(gen(_,k)*gen(_,k));
    vx[k] = var(gen(_,k));}
  double MSx = sum(vx);
  double Lmb = MSx;
  double Rho = MSx*(1-R2)/R2;
  double vy = var(y);
  double ve = 0.5*vy;
  double vb = ve/MSx;
  double Se = (1-R2)*(df+2)*vy;
  double Sb = R2*(df+2)*vy/MSx;
  double mu = mean(y);
  NumericVector b(p);
  NumericVector e = y-mu;
  double b0,eM,h2;
  for(int i=0; i<it; i++){
    for(int j=0; j<p; j++){
      b0 = b[j];
      b[j] = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);
      e = e-gen(_,j)*(b[j]-b0);
    }
    vb = (sum(b*b)+Sb)/(p+df);
    ve = (sum(e*e)+Se)/(n+df);
    Lmb = sqrt(Rho*ve/vb);
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
  }
  h2 = 1-ve/vy;
  NumericVector fit(n);
  for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu;}
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("hat") = fit,
                      Named("Va") = vb,
                      Named("Ve") = ve,
                      Named("h2") = h2);
}

// [[Rcpp::export]]
SEXP emBL(NumericVector y, NumericMatrix gen, double R2 = 0.5, double alpha = 0.02){
  int it = 200;
  double h2 = R2;  
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  NumericVector b(p);
  double b0,eM,Half_L2;
  double mu = mean(y);
  NumericVector e = y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
  // Regulation coefficients
  double cxx = mean(xx);
  double Lmb1 = cxx*((1-h2)/h2)*alpha*0.5;
  double Lmb2 = cxx*((1-h2)/h2)*(1-alpha);
  double OLS, G;
  // Loop
  for(int i=0; i<it; i++){
    // Regression coefficients loop
    for(int j=0; j<p; j++){
      // Ordinary Least Square
      b0 = b[j];
      OLS = (sum(gen(_,j)*e)+xx[j]*b0);
      Half_L2 = 0.5*OLS/(xx[j]+cxx);
      // Regularization of positive OLS
      if(OLS>0){
        G = 0.5*(OLS-Lmb1)/(Lmb2+xx(j));
        if(G>0){b[j] = G+Half_L2;}else{b[j] = Half_L2;}
      }else{
        // Regularization of negative OLS
        G = 0.5*(OLS+Lmb1)/(Lmb2+xx(j));
        if(G<0){b[j] = G+Half_L2;}else{b[j] = Half_L2;}
      }
      // Residuals update
      e = e-gen(_,j)*(b[j]-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;}
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  h2 = 1 - var(e)/var(y);
  // Output
  return List::create(Named("mu") = mu, Named("b") = b, Named("hat") = fit, Named("h2") = h2); }

// [[Rcpp::export]]
SEXP emEN(NumericVector y, NumericMatrix gen, double R2 = 0.5, double alpha = 0.02){
  // Convergence criteria
  int maxit = 300;
  double tol = 10e-11;
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  NumericVector b(p);
  double b0,eM;
  double mu = mean(y);
  NumericVector e = y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
  NumericVector vx(p);
  for(int k=0; k<p; k++){vx[k] = var(gen(_,k));}
  double cxx = sum(vx)*(1-R2)/R2;
  // Regulation coefficients
  double Ve, Va;
  double Sy = sd(y);
  double Lmb = cxx;
  double Lmb1 = 0.5*Lmb*alpha*Sy;
  double Lmb2 = Lmb*(1-alpha);
  double trAC22 = sum(1/(xx+Lmb));
  double OLS, b1;
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      // Ordinary Least Square
      b0 = b[j];
      OLS = (sum(gen(_,j)*e)+xx[j]*b0);
      // Regularization of positive OLS
      if(OLS>0){
        // Regularization of positive OLS
        b1 = (OLS-Lmb1)/(Lmb2+xx(j));if(b1<0){b1=0;};
      }else{
        // Regularization of negative OLS
        b1 = (OLS+Lmb1)/(Lmb2+xx(j));if(b1>0){b1=0;};
      }
      b[j] = b1;
      // Residuals update
      e = e-gen(_,j)*(b1-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Variance components
    Ve = sum(e*y)/(n-1);
    Va = (sum(b*b)+trAC22*Ve)/p;
    Lmb = Ve/Va;
    Lmb1 = 0.5*Lmb*alpha*Sy;
    Lmb2 = Lmb*(1-alpha);
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}
  }
  // Fitting the model
  NumericVector fit(n); for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  // Output
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("Va")=Va*cxx,
                      Named("Ve")=Ve,
                      Named("h2")=Va*cxx/(Va*cxx+Ve));
}

// [[Rcpp::export]]
SEXP emML(NumericVector y, NumericMatrix gen,
          Rcpp::Nullable<Rcpp::NumericVector> D = R_NilValue){
  // Convergence parameters
  int maxit = 300;
  double tol = 10e-8;
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Weights
  bool P_WEIGHTS = FALSE;
  NumericVector d(p);
  if (D.isNotNull()){P_WEIGHTS = TRUE; d=D;}
  // Beta, mu and epsilon
  double b0, eM, ve, vb, h2, mu = mean(y);
  NumericVector b(p), e = y-mu;
  // Marker variance
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(gen(_,i)*gen(_,i));
    vx[i] = var(gen(_,i));}
  double MSx = sum(vx), Lmb=MSx;
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      b0 = b[j];
      if(P_WEIGHTS){
        b[j] = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb/d[j]);
      }else{
        b[j] = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);}
      e = e-gen(_,j)*(b[j]-b0);}
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Variance components update
    ve = sum((y-mu)*e)/n;
    vb = sum((y-mu)*(y-mu-e))/(n*MSx);
    Lmb = ve/vb;
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}}
  // Fitting the model
  NumericVector fit = y-e;
  h2 = vb*MSx/(vb*MSx+ve);
  // Output
  return List::create(Named("mu")=mu,
                      Named("b")=b,
                      Named("hat")=fit,
                      Named("h2")=h2,
                      Named("Vb")=vb,
                      Named("Va")=vb*MSx,
                      Named("Ve")=ve);}

// [[Rcpp::export]]
SEXP emGWA(NumericVector y, NumericMatrix gen){
  int maxit = 500;
  double tol = 10e-8;
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  double b0, eM, ve, vb, h2, mu = mean(y), vy = var(y);
  NumericVector b(p), e = y-mu;
  // Marker variance
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(gen(_,i)*gen(_,i));
    vx[i] = var(gen(_,i));}
  double MSx = sum(vx), Lmb=MSx;
  // Convergence control
  NumericVector bc(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      b0 = b[j];
      b[j] = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);
      e = e-gen(_,j)*(b[j]-b0);}
    // Variance components update
    ve = sum(e*y)/(n-1);
    vb = (vy-ve)/MSx;
    Lmb = ve/vb;
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}}
  // Fitting the model
  NumericVector fit(n);
  for(int k=0; k<n; k++){ fit[k] = sum(gen(k,_)*b)+mu; }
  h2 = vb*MSx/(vb*MSx+ve);
  // Genome-wide screening
  NumericVector LRT(p),PVAL(p),y0(n),e0(n),e1(n),b_ols(p);
  double ve0, ve1, L0, L1;
  for(int j=0; j<p; j++){
    // Full conditional phenotype
    y0 = e+gen(_,j)*b[j];
    // Fixed effect marker
    b_ols[j] = sum(gen(_,j)*y0)/xx[j];
    // Null model
    e0 = y0-mean(y0);
    // Alternative model
    e1 = y0-gen(_,j)*b_ols[j];
    e1 = e1-mean(e1);
    // ReML variance
    ve0 = sum(y0*e0)/(n-1);
    ve1 = sum(y0*e1)/(n-1);
    // Likelihood ratio
    L0 = -sum(e0*e0)/(2*ve0)-0.5*n*log(6.28*ve0);
    L1 = -sum(e1*e1)/(2*ve1)-0.5*n*log(6.28*ve1);
    LRT[j] = 2*(L1-L0);}
  PVAL = -log10(1-pchisq(LRT,1,true));
  // Output
  return List::create(Named("mu")=mu, Named("b")=b, Named("b_LS")=b_ols,
                      Named("h2")=h2, Named("hat")=fit,
                      Named("Vb")=vb, Named("Ve")=ve,
                      Named("LRT")=LRT, Named("PVAL")=PVAL);}

// [[Rcpp::export]]
SEXP BayesA(NumericVector y, NumericMatrix X,
            double it = 1500, double bi = 500,
            double df = 5, double R2 = 0.5){
  // Get dimensions of X
  int p = X.ncol(), n = X.nrow();
  // Estimate crossproducts and MSx
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(X(_,i)*X(_,i));
    vx[i] = var(X(_,i));}
  double MSx = sum(vx);
  // Get priors
  double vy = var(y);
  double Sb = (R2)*df*vy/MSx;
  double Se = (1-R2)*df*vy;
  double mu = mean(y);
  // Create empty objects
  double b0,b1,eM,h2,MU,VE,vg,ve=vy;
  NumericVector b(p),B(p),VB(p),fit(n);
  NumericVector vb=b+Sb,Lmb=ve/vb,e=y-mu;
  // MCMC loop
  for(int i=0; i<it; i++){
    // Update marker effects
    for(int j=0; j<p; j++){
      b0 = b[j];
      // Sample marker effect
      b1 = R::rnorm((sum(X(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]),sqrt(ve/(xx[j]+Lmb[j])));
      b[j] = b1;
      // Update marker variance and residuals
      vb[j] = (Sb+b1*b1)/R::rchisq(df+1);
      e = e - X(_,j)*(b1-b0);
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update residual variance and lambda
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    Lmb = ve/vb;
    // Store posterior sums
    if(i>bi){MU=MU+mu; B=B+b; VB=VB+vb; VE=VE+ve;}
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC; B = B/MCMC;
  VB = VB/MCMC; VE = VE/MCMC;
  // Get fitted values and h2
  vg = sum(VB); h2 = vg/(vg+VE);
  for(int k=0; k<n; k++){fit[k] = sum(X(k,_)*B)+MU;}
  // Return output
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("hat") = fit,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("h2") = h2, Named("MSx") = MSx);}

// [[Rcpp::export]]
SEXP BayesB(NumericVector y, NumericMatrix X,
            double it = 1500, double bi = 500,
            double pi = 0.95, double df = 5, double R2 = 0.5){
  //
  int p = X.ncol();
  int n = X.nrow();
  int MCMC = it-bi;
  //
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(X(_,i)*X(_,i));
    vx[i] = var(X(_,i));}
  double MSx = sum(vx);
  //
  double vy = var(y);
  double Sb = (R2)*df*vy/MSx;
  double Se = (1-R2)*df*vy;
  //
  NumericVector d(p),b(p),D(p),B(p),VB(p);
  NumericVector vb = b+Sb;
  double ve = vy;
  NumericVector Lmb = ve/vb;
  double mu = mean(y);
  NumericVector e = y-mu;
  NumericVector e1(n);
  NumericVector e2(n);
  double b0,b1,eM,h2,C,MU,VE,cj,dj,pj;
  //
  for(int i=0; i<it; i++){
    //
    C = -0.5/sqrt(ve);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = R::rnorm( (sum(X(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]),
                     sqrt(ve/(xx[j]+Lmb[j]))  );
      e1 = e-X(_,j)*(b1-b0); // with marker
      e2 = e-X(_,j)*(0-b0); // without marker
      //
      cj = (1-pi)*exp(C*sum(e1*e1));
      dj = (pi)*exp(C*sum(e2*e2));
      pj = cj/(cj+dj);
      //
      if(R::rbinom(1,pj)==1){
        b[j] = b1;
        d[j] = 1;
      }else{
        b[j] = R::rnorm(0,sqrt(ve/(xx[j]+Lmb[j])));
        d[j] = 0;
      }
      vb[j] = (Sb+b[j]*b[j])/R::rchisq(df+1);
      e = e - X(_,j)*(b[j]-b0);
    }
    //
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM;
    e = e-eM;
    //
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    Lmb = ve/vb;
    //
    if(i>bi){
      MU=MU+mu;
      B=B+b;
      VB=VB+vb;
      VE=VE+ve;
      D = D+d;
    }
    //
  }
  //
  MU = MU/MCMC; B=B/MCMC;
  VB=VB/MCMC; VE=VE/MCMC;
  D = D/MCMC;
  //
  double vg = sum(VB);
  h2 = vg/(vg+VE);
  NumericVector fit(n);
  for(int k=0; k<n; k++){fit[k] = sum(X(k,_)*B)+MU;}
  //
  return List::create(Named("mu") = MU,
                      Named("b") = B,
                      Named("d") = D,
                      Named("hat") = fit,
                      Named("vb") = VB,
                      Named("ve") = VE,
                      Named("h2") = h2,
                      Named("MSx") = MSx);
}

// [[Rcpp::export]]
SEXP BayesC(NumericVector y, NumericMatrix X,
            double it = 1500, double bi = 500,
            double pi = 0.95, double df = 5, double R2 = 0.5){
  // Get dimensions of X
  int p = X.ncol(), n = X.nrow();
  // Estimate crossproducts and MSx
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(X(_,i)*X(_,i));
    vx[i] = var(X(_,i));}
  double MSx = sum(vx);
  // Get priors
  double vy = var(y);
  double Sb = df*(R2)*vy/MSx/(1-pi);
  double Se = df*(1-R2)*vy;
  double mu = mean(y);
  // Create empty objects
  double b0,b1,eM,h2,C,MU,VB,VE,cj,dj,pj,vg,ve=vy,vb=Sb;
  NumericVector d(p),b(p),D(p),B(p),fit(n);
  NumericVector e=y-mu,e1(n),e2(n);
  double Lmb=ve/vb;
  // MCMC loop
  for(int i=0; i<it; i++){
    C = -0.5/sqrt(ve);
    // Update marker effects
    for(int j=0; j<p; j++){
      b0 = b[j];
      // Sample marker effect
      b1 = R::rnorm((sum(X(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb),
                    sqrt(ve/(xx[j]+Lmb)));
      e1 = e-X(_,j)*(b1-b0); // Pr(with marker)
      e2 = e-X(_,j)*(0-b0); // Pr(without marker)
      // Pr(marker included)
      cj = (1-pi)*exp(C*sum(e1*e1));
      dj = (pi)*exp(C*sum(e2*e2));
      pj = cj/(cj+dj);
      // Smple from Bernoulli
      if(R::rbinom(1,pj)==1){
        b[j] = b1; d[j] = 1;
      }else{
        b[j] = R::rnorm(0,sqrt(ve/(xx[j]+Lmb))); d[j] = 0;
      }
      // Update residuals
      e = e - X(_,j)*(b[j]-b0);
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update variance components and lambda
    vb = (sum(b*b)+Sb)/R::rchisq(df+p);
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    Lmb = ve/vb;
    // Store posterior sums
    if(i>bi){MU=MU+mu; B=B+b; D=D+d; VB=VB+vb; VE=VE+ve;}
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC; B = B/MCMC; D = D/MCMC;
  VB = VB/MCMC; VE = VE/MCMC;
  // Get fitted values and h2
  vg = VB*MSx; h2 = vg/(vg+VE);
  for(int k=0; k<n; k++){fit[k] = sum(X(k,_)*B)+MU;}
  // Return output
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("d") = D,   Named("hat") = fit,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("h2") = h2, Named("MSx") = MSx);}

// [[Rcpp::export]]
SEXP BayesL(NumericVector y, NumericMatrix X,
            double it = 1500, double bi = 500,
            double df = 5, double R2 = 0.5){
  // Get dimensions of X
  int p = X.ncol(), n = X.nrow();
  // Estimate crossproducts and MSx
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(X(_,i)*X(_,i));
    vx[i] = var(X(_,i));}
  double MSx = sum(vx);
  double Phi = MSx*(1-R2)/R2;
  // Get priors
  double vy = var(y);
  double Sb = (R2)*df*vy/MSx;
  double Se = (1-R2)*df*vy;
  double mu = mean(y);
  // Create empty objects
  double b0,b1,eM,h2,MU,VE,vg,ve=vy;
  NumericVector b(p),B(p),VB(p),fit(n);
  NumericVector vb=b+Sb,Lmb=ve/vb,e=y-mu;
  // MCMC loop
  for(int i=0; i<it; i++){
    // Update marker effects
    for(int j=0; j<p; j++){
      b0 = b[j];
      // Sample marker effect
      b1 = R::rnorm((sum(X(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]),sqrt(ve/(xx[j]+Lmb[j])));
      b[j] = b1;
      // Update marker variance and residuals
      vb[j] = (Sb+b1*b1)/R::rchisq(df+1);
      e = e - X(_,j)*(b1-b0);
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update residual variance and lambda
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    Lmb = sqrt(Phi*ve/vb);
    // Store posterior sums
    if(i>bi){MU=MU+mu; B=B+b; VB=VB+vb; VE=VE+ve;}
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC; B = B/MCMC;
  VB = VB/MCMC; VE = VE/MCMC;
  // Get fitted values and h2
  vg = sum(VB); h2 = vg/(vg+VE);
  for(int k=0; k<n; k++){fit[k] = sum(X(k,_)*B)+MU;}
  // Return output
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("hat") = fit,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("h2") = h2, Named("MSx") = MSx);}

// [[Rcpp::export]]
SEXP BayesRR(NumericVector y, NumericMatrix X,
             double it = 1500, double bi = 500,
             double df = 5, double R2 = 0.5){
  // Get dimensions of X
  int p = X.ncol(), n = X.nrow();
  // Estimate crossproducts and MSx
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(X(_,i)*X(_,i));
    vx[i] = var(X(_,i));}
  double MSx = sum(vx);
  // Get priors
  double vy = var(y);
  double Sb = (R2)*df*vy/MSx;
  double Se = (1-R2)*df*vy;
  double mu = mean(y);
  // Create empty objects
  double b0,b1,eM,h2,MU,VE,VB,vg,ve=vy,vb=Sb,Lmb=ve/vb;
  NumericVector b(p),B(p),fit(n),e=y-mu;
  // MCMC loop
  for(int i=0; i<it; i++){
    // Update marker effects
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = R::rnorm((sum(X(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb),sqrt(ve/(xx[j]+Lmb)));
      e = e-X(_,j)*(b1-b0);
      b[j] = b1;
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update variance components and lambda
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    vb = (sum(b*b)+Sb)/R::rchisq(p+df);
    Lmb = ve/vb;
    // Store posterior sums
    if(i>bi){MU=MU+mu; B=B+b; VB=VB+vb; VE=VE+ve;}
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC; B = B/MCMC;
  VB = VB/MCMC; VE = VE/MCMC;
  // Get fitted values and h2
  vg = VB*MSx; h2 = vg/(vg+VE);
  for(int k=0; k<n; k++){fit[k] = sum(X(k,_)*B)+MU;}
  // Return output
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("hat") = fit,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("h2") = h2, Named("MSx") = MSx);}

// [[Rcpp::export]]
SEXP BayesCpi(NumericVector y, NumericMatrix X,
          double it = 1500, double bi = 500,
          double df = 5, double R2 = 0.5){
  // Get dimensions of X
  int p = X.ncol(), n = X.nrow();
  // Estimate crossproducts and MSx
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(X(_,i)*X(_,i));
    vx[i] = var(X(_,i));}
  double MSx = sum(vx);
  // Get priors
  double pi = 0.5;
  double vy = var(y);
  double Sb = df*(R2)*vy/MSx/(1-pi);
  double Se = df*(1-R2)*vy;
  double mu = mean(y);
  // Create empty objects
  double b0,b1,b2,eM,h2,C,MU,VB,VE,Pi,pj,vg,ve=vy,vb=Sb;
  NumericVector d(p),b(p),D(p),B(p),fit(n);
  NumericVector e=y-mu,e1(n),e2(n);
  double Lmb=ve/vb;
  // MCMC loop
  for(int i=0; i<it; i++){
    C = -0.5/sqrt(ve);
    // Update marker effects
    for(int j=0; j<p; j++){
      b0 = b[j];
      // Sample marker effect
      b1 = R::rnorm((sum(X(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb),sqrt(ve/(xx[j]+Lmb)));
      b2 = R::rnorm(0,sqrt(ve/(xx[j]+Lmb)));
      e1 = e-X(_,j)*(b1-b0); // Pr(with marker)
      e2 = e-X(_,j)*(0-b0); // Pr(without marker)
      // Pr(marker included)
      pj = (1-pi)*exp(C*(sum(e1*e1)-sum(e2*e2)));
      if(pj>1) pj = 1;
      // Smple from Bernoulli
      if(R::rbinom(1,pj)==1){
        b[j] = b1; d[j] = 1;
      }else{
        b[j] = b2; d[j] = 0;
      }
      // Update residuals
      e = e - X(_,j)*(b[j]-b0);
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update variance components and lambda
    vb = (sum(b*b)+Sb)/R::rchisq(p+df);
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    Lmb = ve/vb;
    pi = mean(d);
    Sb = df*(R2)*vy/MSx/(1-pi);
    // Store posterior sums
    if(i>bi){MU=MU+mu; B=B+b; D=D+d; VB=VB+vb; VE=VE+ve; Pi = Pi+pi;}
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC; B = B/MCMC; D = D/MCMC;
  VB = VB/MCMC; VE = VE/MCMC; Pi = 1-Pi/MCMC;
  // Getting GWAS results
  NumericVector PVAL = -log(1-D);
  // Get fitted values and h2
  vg = VB*MSx/Pi; h2 = vg/(vg+VE);
  for(int k=0; k<n; k++){fit[k] = sum(X(k,_)*B)+MU;}
  // Return output
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("d") = D, Named("pi") = Pi,
                      Named("hat") = fit, Named("h2") = h2,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("PVAL") = PVAL);}

// [[Rcpp::export]]
SEXP BayesDpi(NumericVector y, NumericMatrix X,
          double it = 1500, double bi = 500,
          double df = 5, double R2 = 0.5){
  // Get dimensions of X
  int p = X.ncol(), n = X.nrow();
  // Estimate crossproducts and MSx
  NumericVector xx(p), vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(X(_,i)*X(_,i));
    vx[i] = var(X(_,i));}
  double MSx = sum(vx);
  // Get priors
  double pi = 0.5;
  double vy = var(y);
  double Sb = (R2)*df*vy/MSx;
  double Se = (1-R2)*df*vy;
  double mu = mean(y);
  // Create empty objects
  double b0,b1,b2,eM,h2,C,MU,VE,Pi,pj,vg,ve=vy;
  NumericVector d(p),b(p),D(p),B(p),VB(p),fit(n);
  NumericVector vb=b+Sb,Lmb=ve/vb,e=y-mu,e1(n),e2(n);
  // MCMC loop
  for(int i=0; i<it; i++){
    C = -0.5/sqrt(ve);
    // Update marker effects
    for(int j=0; j<p; j++){
      b0 = b[j];
      // Sample marker effect
      b1 = R::rnorm((sum(X(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb[j]),sqrt(ve/(xx[j]+Lmb[j])));
      b2 = R::rnorm(0,sqrt(ve/(xx[j]+Lmb[j])));        
      e1 = e-X(_,j)*(b1-b0); // Pr(with marker)
      e2 = e-X(_,j)*(b2-b0); // Pr(without marker)
      // Pr(marker included)
      pj = (1-pi)*exp(C*(sum(e1*e1)-sum(e2*e2)));
      if(pj>1) pj = 1;
      // Smple from Bernoulli
      if(R::rbinom(1,pj)==1){
        b[j] = b1; d[j] = 1;
      }else{
        b[j] = b2; d[j] = 0;
      }
      // Update marker variance and residuals
      vb[j] = (Sb+b[j]*b[j])/R::rchisq(df+1);
      e = e - X(_,j)*(b[j]-b0);
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update residual variance and lambda
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    Lmb = ve/vb;    
    pi = mean(d);
    // Store posterior sums
    if(i>bi){
      MU=MU+mu; B=B+b; D=D+d;
      VB=VB+vb; VE=VE+ve; Pi = Pi+pi;}
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC; B = B/MCMC; D = D/MCMC;
  VB = VB/MCMC; VE = VE/MCMC; Pi = 1-Pi/MCMC;
  // Getting GWAS results
  NumericVector PVAL = -log(1-D);
  // Get fitted values and h2
  vg = sum(VB); h2 = vg/(vg+VE);
  for(int k=0; k<n; k++){fit[k] = sum(X(k,_)*B)+MU;}
  // Return output
  return List::create(Named("mu") = MU, Named("b") = B,
                      Named("d") = D, Named("pi") = Pi, 
                      Named("hat") = fit, Named("h2") = h2,
                      Named("vb") = VB, Named("ve") = VE,
                      Named("PVAL") = PVAL);}

// [[Rcpp::export]]
SEXP BayesA2(NumericVector y, NumericMatrix X1, NumericMatrix X2,
             double it = 1500, double bi = 500,
             double df = 5, double R2 = 0.5){
  // Get dimensions of X
  int n = X1.nrow();
  int p1 = X1.ncol();
  int p2 = X2.ncol();
  // Estimate crossproducts and MSx
  NumericVector xx1(p1), vx1(p1);
  for(int i=0; i<p1; i++){
    xx1[i] = sum(X1(_,i)*X1(_,i));
    vx1[i] = var(X1(_,i));}
  double MSx1 = sum(vx1);
  NumericVector xx2(p2), vx2(p2);
  for(int i=0; i<p2; i++){
    xx2[i] = sum(X2(_,i)*X2(_,i));
    vx2[i] = var(X2(_,i));}
  double MSx2 = sum(vx2);
  // Get priors
  double vy = var(y);
  double Sb1 = (R2)*df*vy/MSx1;
  double Sb2 = (R2)*df*vy/MSx2;
  double Se = (1-R2)*df*vy;
  double mu = mean(y);
  // Create empty objects
  double b_t0,b_t1,eM,h2,MU,VE,vg,ve=vy;
  NumericVector b1(p1),B1(p1),VB1(p1);
  NumericVector b2(p2),B2(p2),VB2(p2);
  NumericVector vb1=b1+Sb1,vb2=b2+Sb2,Lmb1=ve/vb1,Lmb2=ve/vb2,e=y-mu,fit(n);
  // MCMC loop
  for(int i=0; i<it; i++){
    // Update marker effects 1
    for(int j=0; j<p1; j++){
      b_t0 = b1[j];
      // Sample marker effect
      b_t1 = R::rnorm((sum(X1(_,j)*e)+xx1[j]*b_t0)/(xx1[j]+Lmb1[j]),sqrt(ve/(xx1[j]+Lmb1[j])));
      b1[j] = b_t1;
      // Update marker variance and residuals
      vb1[j] = (Sb1+b1[j]*b1[j])/R::rchisq(df+1);
      e = e - X1(_,j)*(b_t1-b_t0);
    }
    // Update marker effects 1
    for(int j=0; j<p2; j++){
      b_t0 = b2[j];
      // Sample marker effect
      b_t1 = R::rnorm((sum(X2(_,j)*e)+xx2[j]*b_t0)/(xx2[j]+Lmb2[j]),sqrt(ve/(xx2[j]+Lmb2[j])));
      b2[j] = b_t1;
      // Update marker variance and residuals
      vb2[j] = (Sb2+b2[j]*b2[j])/R::rchisq(df+1);
      e = e - X2(_,j)*(b_t1-b_t0);
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update residual variance and lambda
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    Lmb1 = ve/vb1;
    Lmb2 = ve/vb2;
    // Store posterior sums
    if(i>bi){
      MU=MU+mu; VE=VE+ve;
      B1=B1+b1; VB1=VB1+vb1; 
      B2=B2+b2; VB2=VB2+vb2; 
    }
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC; VE = VE/MCMC;
  B1 = B1/MCMC; VB1 = VB1/MCMC; 
  B2 = B2/MCMC; VB2 = VB2/MCMC; 
  // Get fitted values and h2
  vg = sum(VB1)+sum(VB2); h2 = vg/(vg+VE);
  for(int k=0; k<n; k++){fit[k] = sum(X1(k,_)*B1)+sum(X2(k,_)*B2)+MU;}
  // Return output
  return List::create(Named("hat") = fit, Named("mu") = MU,
                      Named("b1") = B1, Named("b2") = B2, 
                      Named("vb1") = VB1, Named("vb2") = VB2,
                      Named("ve") = VE, Named("h2") = h2);}

// [[Rcpp::export]]
SEXP BayesB2(NumericVector y, NumericMatrix X1, NumericMatrix X2,
             double it = 1500, double bi = 500,
             double pi = 0.95, double df = 5, double R2 = 0.5){
  // Get dimensions of X
  int n = X1.nrow();
  int p1 = X1.ncol();
  int p2 = X2.ncol();
  // Estimate crossproducts and MSx
  NumericVector xx1(p1), vx1(p1);
  for(int i=0; i<p1; i++){
    xx1[i] = sum(X1(_,i)*X1(_,i));
    vx1[i] = var(X1(_,i));}
  double MSx1 = sum(vx1);
  NumericVector xx2(p2), vx2(p2);
  for(int i=0; i<p2; i++){
    xx2[i] = sum(X2(_,i)*X2(_,i));
    vx2[i] = var(X2(_,i));}
  double MSx2 = sum(vx2);
  // Get priors
  double vy = var(y);
  double Sb1 = (R2)*df*vy/MSx1;
  double Sb2 = (R2)*df*vy/MSx2;
  double Se = (1-R2)*df*vy;
  double mu = mean(y);
  // Create empty objects
  double b_t0,b_t1,b_t2,eM,h2,C,MU,VE,cj,dj,pj,vg,ve=vy;
  NumericVector d1(p1),b1(p1),D1(p1),B1(p1),VB1(p1),fit(n);
  NumericVector d2(p2),b2(p2),D2(p2),B2(p2),VB2(p2);
  NumericVector vb1=b1+Sb1,vb2=b2+Sb2,Lmb1=ve/vb1,Lmb2=ve/vb2,e=y-mu,e1(n),e2(n);
  // MCMC loop
  for(int i=0; i<it; i++){
    C = -0.5/sqrt(ve);
    // Update marker effects 1
    for(int j=0; j<p1; j++){
      b_t0 = b1[j];
      // Sample marker effect
      b_t1 = R::rnorm((sum(X1(_,j)*e)+xx1[j]*b_t0)/(xx1[j]+Lmb1[j]),sqrt(ve/(xx1[j]+Lmb1[j])));
      b_t2 = R::rnorm(0,sqrt(ve/(xx1[j]+Lmb1[j])));
      e1 = e-X1(_,j)*(b_t1-b_t0); // Pr(with marker)
      e2 = e-X1(_,j)*(b_t2-b_t0); // Pr(without marker)
      // Pr(marker included)
      cj = (1-pi)*exp(C*sum(e1*e1));
      dj = (pi)*exp(C*sum(e2*e2));
      pj = cj/(cj+dj);
      // Smple from Bernoulli
      if(R::rbinom(1,pj)==1){
        b1[j] = b_t1; d1[j] = 1;
      }else{
        b1[j] = b_t2; d1[j] = 0;
      }
      // Update marker variance and residuals
      vb1[j] = (Sb1+b1[j]*b1[j])/R::rchisq(df+1);
      e = e - X1(_,j)*(b1[j]-b_t0);
    }
    // Update marker effects 1
    for(int j=0; j<p2; j++){
      b_t0 = b2[j];
      // Sample marker effect
      b_t1 = R::rnorm((sum(X2(_,j)*e)+xx2[j]*b_t0)/(xx2[j]+Lmb2[j]),sqrt(ve/(xx2[j]+Lmb2[j])));
      b_t2 = R::rnorm(0,sqrt(ve/(xx2[j]+Lmb2[j])));
      e1 = e-X2(_,j)*(b_t1-b_t0); // Pr(with marker)
      e2 = e-X2(_,j)*(b_t2-b_t0); // Pr(without marker)
      // Pr(marker included)
      cj = (1-pi)*exp(C*sum(e1*e1));
      dj = (pi)*exp(C*sum(e2*e2));
      pj = cj/(cj+dj);
      // Smple from Bernoulli
      if(R::rbinom(1,pj)==1){
        b2[j] = b_t1; d2[j] = 1;
      }else{
        b2[j] = b_t2; d2[j] = 0;
      }
      // Update marker variance and residuals
      vb2[j] = (Sb2+b2[j]*b2[j])/R::rchisq(df+1);
      e = e - X2(_,j)*(b2[j]-b_t0);
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update residual variance and lambda
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    Lmb1 = ve/vb1;
    Lmb2 = ve/vb2;
    // Store posterior sums
    if(i>bi){
      MU=MU+mu; VE=VE+ve;
      B1=B1+b1; D1=D1+d1; VB1=VB1+vb1; 
      B2=B2+b2; D2=D2+d2; VB2=VB2+vb2; 
    }
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC; VE = VE/MCMC;
  B1 = B1/MCMC; D1 = D1/MCMC; VB1 = VB1/MCMC; 
  B2 = B2/MCMC; D2 = D2/MCMC; VB2 = VB2/MCMC; 
  // Get fitted values and h2
  vg = sum(VB1)+sum(VB2); h2 = vg/(vg+VE);
  for(int k=0; k<n; k++){fit[k] = sum(X1(k,_)*B1)+sum(X2(k,_)*B2)+MU;}
  // Return output
  return List::create(Named("mu") = MU,
                      Named("b1") = B1, Named("d1") = D1, Named("vb1") = VB1,
                            Named("b2") = B2, Named("d2") = D2, Named("vb2") = VB2,
                                  Named("ve") = VE, Named("hat") = fit, Named("h2") = h2);}

// [[Rcpp::export]]
SEXP BayesRR2(NumericVector y, NumericMatrix X1, NumericMatrix X2,
              double it = 1500, double bi = 500,
              double df = 5, double R2 = 0.5){
  // Get dimensions of X
  int n = X1.nrow();
  int p1 = X1.ncol();
  int p2 = X2.ncol();
  // Estimate crossproducts and MSx
  NumericVector xx1(p1), vx1(p1);
  for(int i=0; i<p1; i++){
    xx1[i] = sum(X1(_,i)*X1(_,i));
    vx1[i] = var(X1(_,i));}
  double MSx1 = sum(vx1);
  NumericVector xx2(p2), vx2(p2);
  for(int i=0; i<p2; i++){
    xx2[i] = sum(X2(_,i)*X2(_,i));
    vx2[i] = var(X2(_,i));}
  double MSx2 = sum(vx2);
  // Get priors
  double vy = var(y);
  double Sb1 = (R2)*df*vy/MSx1;
  double Sb2 = (R2)*df*vy/MSx2;
  double Se = (1-R2)*df*vy;
  double mu = mean(y);
  // Create empty objects
  double b_t0,b_t1,eM,h2,MU,VE,vg,vb1,vb2,VB1,VB2,Lmb1=MSx1,Lmb2=MSx2,ve=vy;
  NumericVector b1(p1),B1(p1),b2(p2),B2(p2),e=y-mu,fit(n);
  // MCMC loop
  for(int i=0; i<it; i++){
    // Update marker effects 1
    for(int j=0; j<p1; j++){
      b_t0 = b1[j];
      b_t1 = R::rnorm((sum(X1(_,j)*e)+xx1[j]*b_t0)/(xx1[j]+Lmb1),sqrt(ve/(xx1[j]+Lmb1)));
      b1[j] = b_t1;
      e = e - X1(_,j)*(b_t1-b_t0);
    }
    // Update marker effects 1
    for(int j=0; j<p2; j++){
      b_t0 = b2[j];
      b_t1 = R::rnorm((sum(X2(_,j)*e)+xx2[j]*b_t0)/(xx2[j]+Lmb2),sqrt(ve/(xx2[j]+Lmb2)));
      b2[j] = b_t1;
      e = e - X2(_,j)*(b_t1-b_t0);
    }
    // Update intercept
    eM = R::rnorm(mean(e),sqrt(ve/n));
    mu = mu+eM; e = e-eM;
    // Update residual variance and lambda
    ve = (sum(e*e)+Se)/R::rchisq(n+df);
    vb1 = (Sb1+sum(b1*b1))/R::rchisq(df+p1);
    vb2 = (Sb2+sum(b2*b2))/R::rchisq(df+p2);
    Lmb1 = ve/vb1;
    Lmb2 = ve/vb2;
    // Store posterior sums
    if(i>bi){
      MU=MU+mu; VE=VE+ve;
      B1=B1+b1; VB1=VB1+vb1; 
      B2=B2+b2; VB2=VB2+vb2; 
    }
  }
  // Get posterior means
  double MCMC = it-bi;
  MU = MU/MCMC; VE = VE/MCMC;
  B1 = B1/MCMC; VB1 = VB1/MCMC; 
  B2 = B2/MCMC; VB2 = VB2/MCMC; 
  // Get fitted values and h2
  vg = (VB1*MSx1+VB2*MSx2); h2 = vg/(vg+VE);
  for(int k=0; k<n; k++){fit[k] = sum(X1(k,_)*B1)+sum(X2(k,_)*B2)+MU;}
  // Return output
  return List::create(Named("hat") = fit, Named("mu") = MU,
                      Named("b1") = B1, Named("b2") = B2, 
                      Named("vb1") = VB1, Named("vb2") = VB2,
                      Named("ve") = VE, Named("h2") = h2);}

// [[Rcpp::export]]
SEXP emML2(NumericVector y, NumericMatrix X1, NumericMatrix X2,
           Rcpp::Nullable<Rcpp::NumericVector> D1 = R_NilValue,
           Rcpp::Nullable<Rcpp::NumericVector> D2 = R_NilValue){
  int maxit = 350; double tol = 10e-8;
  // Functions starts here
  int p1 = X1.ncol();
  int p2 = X2.ncol();
  int n = X1.nrow();
  // Weights
  bool P1_WEIGHTS = FALSE;
  bool P2_WEIGHTS = FALSE;
  NumericVector d1(p1), d2(p2);
  if (D1.isNotNull()){P1_WEIGHTS = TRUE; d1=D1;}
  if (D2.isNotNull()){P2_WEIGHTS = TRUE; d2=D2;}
  // Beta, mu and epsilon
  double b0, eM, ve, vb1, vb2, h2, mu = mean(y);
  NumericVector b1(p1), b2(p2), u1(n), u2(n), cY(n), e = y-mu;
  // Marker variance
  NumericVector x1x1(p1), vx1(p1);
  for(int i=0; i<p1; i++){
    x1x1[i] = sum(X1(_,i)*X1(_,i));
    vx1[i] = var(X1(_,i));}
  double MSx1 = sum(vx1), Lmb1=MSx1;
  NumericVector x2x2(p2), vx2(p2);
  for(int i=0; i<p2; i++){
    x2x2[i] = sum(X2(_,i)*X2(_,i));
    vx2[i] = var(X2(_,i));}
  double MSx2 = sum(vx2), Lmb2=MSx2;
  // Convergence control
  NumericVector bc1(p1), bc2(p2);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Save b(t0)
    bc1 = b1+0;
    bc2 = b2+0;
    // Regression coefficients loop 1
    for(int j=0; j<p1; j++){
      b0 = b1[j];
      if(P1_WEIGHTS){
        b1[j] = (sum(X1(_,j)*e)+x1x1[j]*b0)/(x1x1[j]+Lmb1/d1[j]);
      }else{
        b1[j] = (sum(X1(_,j)*e)+x1x1[j]*b0)/(x1x1[j]+Lmb1);}
      e = e-X1(_,j)*(b1[j]-b0);}
    // Regression coefficients loop 2
    for(int j=0; j<p2; j++){
      b0 = b2[j];
      if(P2_WEIGHTS){
        b2[j] = (sum(X2(_,j)*e)+x2x2[j]*b0)/(x2x2[j]+Lmb2/d2[j]);
      }else{
        b2[j] = (sum(X2(_,j)*e)+x2x2[j]*b0)/(x2x2[j]+Lmb2);}
      e = e-X2(_,j)*(b2[j]-b0);}
    // Fitting the model
    for(int k=0; k<n; k++){
      u1[k] = sum(X1(k,_)*b1);
      u2[k] = sum(X2(k,_)*b2);
    }
    // Intercept update
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Variance components update
    cY = u1+u2+e;
    ve = sum(e*cY)/n;
    vb1 = (sum(u1*cY)/n)/MSx1;
    vb2 = (sum(u2*cY)/n)/MSx2;
    Lmb1 = ve/vb1;
    Lmb2 = ve/vb2;
    // Convergence
    ++numit;
    cnv = sum(abs(bc1-b1))+sum(abs(bc2-b2));
    if( cnv<tol ){break;}}
  // Fitting the model
  NumericVector fit = mu+u1+u2;
  h2 = 1-ve/var(y);
  // Output
  return List::create(Named("mu")=mu,
                      Named("b1")=b1, Named("b2")=b2, 
                      Named("Vb1")=vb1, Named("Vb2")=vb2, Named("Ve")=ve,
                            Named("u1")=u1, Named("u2")=u2,
                            Named("MSx1")=MSx1, Named("MSx2")=MSx2,
                            Named("h2")=h2, Named("hat")=fit);}

// [[Rcpp::export]]
NumericMatrix CNT(NumericMatrix X){for(int j=0;j<X.ncol();j++){X(_,j)=X(_,j)-mean(X(_,j));}; return(X);}

// [[Rcpp::export]]
NumericMatrix IMP(NumericMatrix X){
  int p = X.ncol(); int n = X.nrow();
  LogicalVector MIS(n); NumericVector x(n);
  NumericVector z; double EXP;
  for(int j=0; j<p; j++){
    if(is_true(any(is_na(X(_,j))))){
      x = X(_,j); MIS = is_na(x);
      z = x[!MIS]; EXP = mean(z);
      X(_,j) = ifelse(MIS,EXP,x);}
  };return(X);};

// [[Rcpp::export]]
NumericMatrix GAU(NumericMatrix X){
  int n = X.nrow(); NumericVector D; NumericMatrix K(n,n); double d2, md;
  for(int i=0; i<n; i++){; for(int j=0; j<n; j++){
    if(i==j){ K(i,j)=0; }else if(j>i){; D = X(i,_)-X(j,_);
    d2 = sum(D*D); K(i,j)=d2; K(j,i)=d2; }}}; md = mean(K);
    for(int i=0; i<n; i++){K(i,_) = exp(-K(i,_)/md);} return K;}

// [[Rcpp::export]]
NumericMatrix GRM(NumericMatrix X, bool Code012 = false){
  int n = X.nrow(), p = X.ncol();
  NumericMatrix K(n,n); NumericVector xx(p); double zz, Sum2pq=0.0;
  for(int i=0; i<p; i++){ xx[i] = mean(X(_,i)); }
  if(Code012){
    for(int i=0; i<p; i++){ Sum2pq = Sum2pq + xx[i]*xx[i]/2;}
  }else{
    for(int i=0; i<p; i++){Sum2pq = Sum2pq + var(X(_,i));}
  }
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(i<=j ){
        zz = sum( (X(i,_)-xx)*(X(j,_)-xx) );
        K(i,j)=zz; K(j,i)=zz;
      }
    }
  }
  return K/Sum2pq;
}

// [[Rcpp::export]]
NumericVector SPC(NumericVector y, NumericVector blk, NumericVector row, NumericVector col, double rN=3, double cN=1){
  int n = y.size(), t1=0, t2=0; NumericVector Cov(n), Phe(n), Obs(n);
  for(int i=0; i<n; i++){; for(int j=0; j<n; j++){
    t1 = row[i]-row[j]; if(t1<0){t1=-t1;}; t2 = col[i]-col[j]; if(t2<0){t2=-t2;};
    if( (i>j) & (blk[i]==blk[j]) & (t1<=rN) & (t2<=cN) ){
      Phe[i] = Phe[i]+y[j]; Obs[i] = Obs[i]+1; Phe[j] = Phe[j]+y[i]; Obs[j] = Obs[j]+1; }}}
  Cov = Phe/Obs; return Cov;}

// [[Rcpp::export]]
NumericMatrix SPM(NumericVector blk, NumericVector row, NumericVector col, double rN=3, double cN=1){
  int n=blk.size(),t1=0,t2=0; NumericMatrix X(n,n); for(int i=0; i<n; i++){; for(int j=0; j<n; j++){
    t1 = row[i]-row[j]; if(t1<0){t1=-t1;}; t2 = col[i]-col[j]; if(t2<0){t2=-t2;};
    if( (blk[i]==blk[j]) & (i>j)  & (t1<=rN) & (t2<=cN) ){ X(i,j) = 1; X(j,i) = 1;
       }else{ X(i,j) = 0; X(j,i) = 0; }}; X(i,i) = 0;}; return X;}

// [[Rcpp::export]]
SEXP mtgsru(NumericMatrix Y, NumericMatrix X,
          NumericMatrix b, NumericMatrix vb, NumericVector ve,
          NumericMatrix iG, int maxit = 50){
  
  double tol = 10e-8;
  Rcpp::Environment base("package:base");
  Rcpp::Function solve = base["solve"];
  
  int k = Y.ncol(), p = X.ncol(), n0 = X.nrow();
  NumericMatrix fit(n0,k),o(n0,k),y(n0,k),e(n0,k);
  
  for(int i=0; i<k; i++){
    o(_,i) = ifelse(is_na(Y(_,i)),0,1);
    y(_,i) = ifelse(is_na(Y(_,i)),0,Y(_,i));}
  
  NumericVector n = colSums(o);
  NumericMatrix xx(p,k), vx(p,k); double tmp;
  for(int i=0; i<p; i++){ for(int j=0; j<k; j++){
    xx(i,j) = sum(X(_,i)*X(_,i)*o(_,j));
    tmp = sum(X(_,i)*o(_,j))/n(j);
    vx(i,j) = xx(i,j)/n(j)-tmp*tmp;}}
  
  NumericVector MSx = colSums(vx);
  NumericMatrix rho(k,k),LHS(k,k);
  NumericVector b0(k),b1(k),vy(k),RHS(k);
  
  for(int i=0; i<k; i++){ e(_,i) = y(_,i)+0;
    vy(i) = sum(e(_,i)*e(_,i))/n(i); }
  iG = solve(vb); NumericMatrix bc(p,k);
  
  int numit = 0; double cnv = 1; while(numit<maxit){
    bc = b+0; for(int j=0; j<p; j++){ b0 = b(j,_); LHS = iG+0;
    for(int i=0; i<k; i++){ LHS(i,i) = iG(i,i)+(xx(j,i)/ve(i));
      RHS(i) = (sum(e(_,i)*X(_,j))+xx(j,i)*b0(i))/ve(i);}
    b1 = solve(LHS, RHS);  b(j,_) = b1;
    
    for(int i=0; i<k; i++){ e(_,i) = (e(_,i)-X(_,j)*(b1(i)-b0(i)))*o(_,i);} }
    for(int i=0; i<k; i++){ ve(i) = (sum(e(_,i)*y(_,i)))/(n(i)-1);}
    for(int i=0; i<n0; i++){ for(int j=0; j<k; j++){fit(i,j) = sum(X(i,_)*b(_,j));}}
    for(int i=0; i<k; i++){ for(int j=0; j<k; j++){
      vb(i,j) = (sum(fit(_,i)*y(_,j))+sum(fit(_,j)*y(_,i))) / ((n(i)*MSx(i))+(n(j)*MSx(j))); }}
    for(int i=0; i<k; i++){vb(i,i) = vb(i,i)*1.01;}
    iG = solve(vb); ++numit; cnv = sum(abs(bc-b)); if( cnv<tol ){break;}}
  for(int i=0; i<n0; i++){for(int j=0; j<k; j++){fit(i,j) = sum(X(i,_)*b(_,j));}}
  
  NumericVector h2(k); h2 = 1-ve/vy;
  return List::create(Named("b")=b, Named("hat")=fit, Named("e")=fit, Named("MSx")=MSx,
                      Named("vb")=vb, Named("ve")=ve, Named("iG")=iG, Named("h2")=h2);}

// [[Rcpp::export]]
SEXP lasso(NumericVector y, NumericMatrix gen){
  int maxit = 300;
  double tol = 10e-8;
  // Functions starts here
  int p = gen.ncol();
  int n = gen.nrow();
  // Beta, mu and epsilon
  double eM, mu = mean(y);
  NumericVector b(p), e = y-mu;
  // Marker variance
  NumericVector xx(p);
  for(int i=0; i<p; i++){xx[i] = sum(gen(_,i)*gen(_,i));}
  double Lmb=mean(xx)/p;
  // Convergence control
  NumericVector bc(p), yx(p);
  int numit = 0;
  double cnv = 1;
  // Loop
  while(numit<maxit){
    // Regression coefficients loop
    bc = b+0;
    for(int j=0; j<p; j++){
      e = e+gen(_,j)*b[j];
      yx[j] = sum(e*gen(_,j));
      if(yx[j]>0){
        b[j] = (yx[j]-Lmb)/xx[j];
        if(b[j]<0){b[j]=0;}
      }else{
        b[j] = (yx[j]+Lmb)/xx[j];
        if(b[j]>0){b[j]=0;}
      }
      e = e-gen(_,j)*b[j];}
    // Update regularization and intercept
    Lmb = 2*mean(abs(yx)-abs(b*xx));
    Lmb = 2*sqrt(Lmb); // this line was added by ad hoc, without derivations
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
    // Convergence
    ++numit;
    cnv = sum(abs(bc-b));
    if( cnv<tol ){break;}}
  // Fitting the model
  NumericVector fit = y-e;
  double h2 = 1-(sum(e*y)/(n-1))/var(y);
  // Output
  return List::create(Named("mu")=mu, Named("b")=b,
                      Named("h2")=h2, Named("hat")=fit,
                      Named("Lmb")=Lmb);}

// [[Rcpp::export]]
SEXP emBCpi(NumericVector y, NumericMatrix gen, double df = 10, double R2 = 0.5, double Pi = 0.75){
  int it = 200;
  int p = gen.ncol();
  int n = gen.nrow();
  NumericVector d(p);
  NumericVector b(p);
  double vy = var(y);
  if(Pi>0.5){Pi = 1-Pi;} 
  NumericVector xx(p);
  NumericVector vx(p);
  for(int i=0; i<p; i++){
    xx[i] = sum(gen(_,i)*gen(_,i));
    vx[i] = var(gen(_,i));
  }
  double PriorPi = Pi*1.0;
  double MSx = sum(vx)*Pi*(1-Pi);
  double Sa = R2*(df+2)*vy/MSx;
  double Se = (1-R2)*(df+2)*vy;
  double mu = mean(y);
  NumericVector e = y-mu;
  NumericVector e1(n);
  NumericVector e2(n);
  double ve = Sa;
  double va = Se;
  double Lmb = ve/va;
  double b0,b1,LR,eM,h2,C;
  double Pi0 = (1-Pi)/Pi;
  for(int i=0; i<it; i++){
    C = -0.5/sqrt(ve);
    for(int j=0; j<p; j++){
      b0 = b[j];
      b1 = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);
      e1 = e-gen(_,j)*(b1-b0);
      e2 = e-gen(_,j)*(0-b0);
      LR = Pi0*exp(C*(sum(e2*e2)-sum(e1*e1)));
      d[j] = (1/(1+LR));
      b[j] = b1*d[j];
      e = e - gen(_,j)*(b[j]-b0);
    }
    // update pi
    Pi = ((1-mean(d))*p + PriorPi*df)/(p+df);
    Pi0 = (1-Pi)/Pi;
    MSx = sum(vx)*Pi*(1-Pi);
    Sa = R2*(df+2)*vy/MSx;
    // update var
    eM = mean(e);
    ve = (sum(e*e)+Se)/(n+df);
    va = (sum(b*b)+Sa)/(p+df)/(mean(d)-Pi);
    Lmb = ve/va;
    // update mu
    eM = mean(e);
    mu = mu+eM;
    e = e-eM;
  }
  h2 = 1-ve/vy;
  NumericVector fit(n);
  for(int k=0; k<n; k++){fit[k] = sum(gen(k,_)*b)+mu;}
  return List::create(Named("mu") = mu,
                      Named("b") = b,
                      Named("d") = d,
                      Named("pi") = Pi,
                      Named("hat") = fit,
                      Named("Vg") = va*MSx,
                      Named("Va") = va,
                      Named("Ve") = ve,
                      Named("h2") = h2);}

// [[Rcpp::export]]
NumericMatrix NNSEARCH(NumericVector blk, NumericVector row, NumericVector col, int rN, int cN){
  int n = blk.size(); NumericMatrix X(n,(rN*2+1)*(cN*2+1)); NumericVector Obs(n);
  for(int i=0; i<n; i++){ for(int j=0; j<n; j++){
    if( (i>j) && (blk[i]==blk[j]) && (fabs(row[i]-row[j])<=rN) && (fabs(col[i]-col[j])<=cN) ){
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
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("h2")=1-vna/vy,
                            Rcpp::Named("e")=e,
                            Rcpp::Named("Lmb")=Lmb,
                            Rcpp::Named("vb")=Vb);
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
  return Rcpp::List::create(Rcpp::Named("mu")=mu,
                            Rcpp::Named("b")=b,
                            Rcpp::Named("h2")=1-vna/vy,
                            Rcpp::Named("e")=e,
                            Rcpp::Named("Lmb")=Lmb,
                            Rcpp::Named("vb")=Vb);
}