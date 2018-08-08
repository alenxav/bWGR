#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP KMUP(NumericMatrix X, NumericVector b, NumericVector d, NumericVector xx, NumericVector e, NumericVector L, double Ve, double pi){  
  int p = X.ncol();
  NumericVector e1 = e+0;
  NumericVector e2 = e+0;
  double b0,b1,b2,cj,dj,pj;
  double C = -0.5/Ve;
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
  double C = -0.5/Ve;
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
    C = -0.5/ve;
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
    C = -0.5/ve;
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
  for(int k=0; k<p; k++){xx[k] = sum(gen(_,k)*gen(_,k));}
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
    C = -0.5/ve;
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
    C = -0.5/ve;
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
    C = -0.5/ve;
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
SEXP emML(NumericVector y, NumericMatrix gen,
          Rcpp::Nullable<Rcpp::NumericVector> D = R_NilValue){
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
      if(P_WEIGHTS){
        b[j] = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb/d[j]);
      }else{
        b[j] = (sum(gen(_,j)*e)+xx[j]*b0)/(xx[j]+Lmb);}
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
  // Output
  return List::create(Named("mu")=mu, Named("b")=b,
                      Named("h2")=h2, Named("hat")=fit,
                      Named("Vb")=vb, Named("Ve")=ve);}

// [[Rcpp::export]]
void CNT(NumericMatrix X){for(int j=0;j<X.ncol();j++){X(_,j)=X(_,j)-mean(X(_,j));}}

// [[Rcpp::export]]
void IMP(NumericMatrix X){;int p = X.ncol(); int n = X.nrow();
LogicalVector MIS(n); NumericVector x(n); NumericVector z; double EXP;
for(int j=0; j<p; j++){;if(is_true(any(is_na(X(_,j))))){
  x = X(_,j); MIS = is_na(x);z = x[!MIS]; EXP = mean(z);
  X(_,j) = ifelse(MIS,EXP,x);};};};

// [[Rcpp::export]]
NumericMatrix GAU(NumericMatrix X){
  int n = X.nrow(); NumericVector D; NumericMatrix K(n,n); double d2, md;
  for(int i=0; i<n; i++){; for(int j=0; j<n; j++){
    if(i==j){ K(i,j)=0; }else if(j>i){ D = X(i,_)-X(j,_);
    d2 = sum(D*D); d2 = d2*d2; K(i,j)=d2; K(j,i)=d2; }}}; md = mean(K);
    for(int i=0; i<n; i++){K(i,_) = exp(-K(i,_)/md);} return K;}

// [[Rcpp::export]]
NumericVector SPC(NumericVector y, NumericVector blk, NumericVector row, NumericVector col, int rN=3, int cN=1){
  int n = y.size(); NumericVector Cov(n), Phe(n), Obs(n);
  for(int i=0; i<n; i++){; for(int j=0; j<n; j++){
      if( (i>j) & (blk[i]==blk[j]) & (abs(row[i]-row[j])<=rN) & (abs(col[i]-col[j])<=cN) ){
        Phe[i] = Phe[i]+y[j]; Obs[i] = Obs[i]+1; Phe[j] = Phe[j]+y[i]; Obs[j] = Obs[j]+1; }}}
  Cov = Phe/Obs; return Cov;}