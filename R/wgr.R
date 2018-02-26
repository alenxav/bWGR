
wgr = function(y,X,
               it=1500,bi=500,th=1,
               bag=1,rp=FALSE,
               iv=FALSE,de=FALSE,pi=0,
               df=5,R2=0.5,
               eigK=NULL,VarK=0.95,
               verb=FALSE){
  
  if(de) iv=TRUE
  anyNA = function(x) any(is.na(x))
  # Imputation with expectation
  if(anyNA(X)){
    cat('Imputing missing genotypes\n')
    imp = function(x){
      x[is.na(x)] = mean(x,na.rm=TRUE)
      x[is.nan(x)] = 0
      return(x)}
    X = apply(X,2,imp)}
  if(bag!=1) df = df/(bag^2)
  gen0 = X
  # Polygenic term
  if(!is.null(eigK)){
    V = eigK$values
    pk = which.max((cumsum(V)/length(V))>VarK)
    U0 = U = eigK$vectors[,1:pk]
    V = V[1:pk]
    H = h = rep(0,pk)
    dh = rep(0,pk)
    Vk = rep(1,pk)
    xxK = rep(bag,pk)
  }
  # Remove missing y's
  if(anyNA(y)){
    mis = which(is.na(y))
    y = y[-mis]
    X = X[-mis,]
    if(!is.null(eigK)) U = U[-mis,]
  }
  # MCMC settings
  post = seq(bi,it,th)
  mc = length(post)
  # Preparing inputs
  n = nrow(X)
  p = ncol(X)
  xx = apply(X,2,crossprod)*bag
  b = rep(0,p)
  d = rep(1,p)
  mu = mean(y - X%*%b)
  e = y - mu
  MSx = sum(apply(X,2,var,na.rm=T))
  Va = MSx
  Vb = rep(Va,p)
  Ve = 1
  L = Vb/Ve
  vy = var(y,na.rm=T)
  # Priors 10/15/17
  Sb = (R2)*df*vy/MSx;
  Se = (1-R2)*df*vy;
  if(!is.null(eigK)) Sk = R2*var(y,na.rm=T)*(df+2)
  # Storing Posterior
  B0 = VA = VE = VP = S = 0
  VB = D = B = rep(0,p)
  #RUN
  if(verb) pb = txtProgressBar(style = 3)
  for(i in 1:it){
    # Resampling
    if(bag!=1) Use = sort(sample(n,n*bag,rp))-1
    # Update polygenic term and regression coefficients
    if(!is.null(eigK)){
      Lk = Ve/(V*Vk)
      if(bag!=1){
        update = KMUP2(U,Use,h,dh,xxK,e,Lk,Ve,0)
      }else{
        update = KMUP(U,h,dh,xxK,e,Lk,Ve,0)
      }
      h = update[[1]]
      e = update[[3]]
      if(bag!=1){update = KMUP2(X,Use,b,d,xx,e,L,Ve,pi)}else{update = KMUP(X,b,d,xx,e,L,Ve,pi)}
      if(pi>0) d = update[[2]]
      b = update[[1]]
      e = update[[3]]
    }else{
      # Update regression coefficients without polygene
      if(bag!=1){update = KMUP2(X,Use,b,d,xx,e,L,Ve,pi)}else{update = KMUP(X,b,d,xx,e,L,Ve,pi)}
      if(pi>0) d = update[[2]]
      b = update[[1]]
      e = update[[3]]
    }
    # Update genetic variance
    if(iv){
      # Variable selection?
      if(pi>0){
        q = d; q[q==0]=pi
        if(de){
          # Laplace?
          Vb = sqrt( b^2*Ve/MSx )
        }else{
          # T?
          Vb = (Sb+(b)^2)/rchisq(p,df+1)
        }
        # All-in?
      }else{
        # Laplace?
        if(de){
          Vb = sqrt( b^2*Ve/MSx )
        }else{
          # T?
          Vb = (Sb+b^2)/rchisq(p,df+1)
        }
      }
    }else{
      Va = (crossprod(b)+Sb)/rchisq(1,df+p)
      Vb = rep(Va,p)
    }
    if(!is.null(eigK)){
      Vp = (sum(h^2/V)+Sk)/rchisq(1,df+pk)
      Vk = rep(Vp,pk)
    }
    # Residual variance
    Ve = (crossprod(e)+Se)[1,1]/rchisq(1,n*bag+df)
    L = Ve/Vb;
    # Intercept
    if(!is.null(eigK)){e = y-mu-X%*%b-U%*%h}else{e = y-mu-X%*%b}
    mu0 = rnorm(1,mean(e), Ve/n )
    mu = mu+mu0
    e = e-mu0
    # Save posterior
    if(i%in%post){
      B0 = B0+mu;
      B = B+b
      D = D+d
      VE = VE+Ve
      if(iv){VB = VB+Vb}else{VA = VA+Va}
      if(!is.null(eigK)){H = H+h; VP = VP+Vp}
    }
    if(verb) setTxtProgressBar(pb, i/it)
  }
  if(verb) close(pb)
  # Posterior mean
  B0 = B0/mc
  D = D/mc
  B = B/mc/mean(D)
  VE = VE/mc
  if(iv){VB = VB/mc}else{VA = VA/mc}
  if(!is.null(eigK)){H = H/mc; VP = VP/mc}
  # Fitted values
  if(!is.null(eigK)){
    poly = U0 %*% H
    HAT = B0 + gen0 %*% B + poly
  }else{
    HAT = B0 + gen0 %*% B
  }
  # Output
  if(!is.null(eigK)){
    if(iv){
      final = list('mu'=B0,'b'=B,'Vb'=VB,'d'=D,'Ve'=VE,'hat'=HAT,'u'=poly,'Vk'=VP,'cxx'=mean(xx))
    }else{
      final = list('mu'=B0,'b'=B,'Vb'=VA,'d'=D,'Ve'=VE,'hat'=HAT,'u'=poly,'Vk'=VP,'cxx'=mean(xx))
    }
  }else{
    if(iv){
      final = list('mu'=B0,'b'=B,'Vb'=VB,'d'=D,'Ve'=VE,'hat'=HAT,'cxx'=mean(xx))
    }else{
      final = list('mu'=B0,'b'=B,'Vb'=VA,'d'=D,'Ve'=VE,'hat'=HAT,'cxx'=mean(xx))
    }
  }
  return(final)
}
