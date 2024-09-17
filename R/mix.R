
mixed = function(y,random=NULL,fixed=NULL,data=NULL,X=list(),alg=emML,maxit=10,Deregress=FALSE,...){
  
  # Get y vector
  if(!is.null(data)) y = data[[deparse(substitute(y))]]
  Vy = var(y,na.rm=TRUE)
  
  # Remove missing
  if(anyNA(y)){
    wna = which(is.na(y)); y = y[-wna]
    data = droplevels.data.frame(data[-wna,])}
  
  # Random effects
  if(!is.null(random)){
    rnd = attr(terms(random),"term.labels")
    LMB0 = LMB = sapply(data[rnd], function(x) mean(table(x)) )
    df0 = sapply(data[rnd],function(x) ifelse(is.factor(x),length(unique(x)),1))
    RND = TRUE }else{ RND = FALSE }
  
  # Fixed effects or intercept
  if(!is.null(fixed)){
    fxd = attr(terms(fixed),"term.labels")
    df = sum(sapply(data[fxd],function(x) ifelse(is.factor(x),length(unique(x)),1)))
    FIX = TRUE }else{ FIX = FALSE; df=1 }
  B = H = A = list(); n = length(y)
  
  # Fixed Categorical & Continuous Term
  FCT = function(y,x){
    if(is.factor(x)){
      b = tapply(y,x,mean,na.rm=T);
      b[is.na(b)] = 0
      g = b[x]
      e = y-g
    }else{
      b = cov(y,x,use='pair')/var(x,na.rm=T);
      b = c(b);
      g = x*b
      e = y-g }
    return(list(g=g,b=b,e=e))}
  
  # Random Categorical & Continuous Term
  RCT = function(y0,x,lmb0){
    y00 = y0
    x00 = x
    if( anyNA(y0) ){
      x = x[!is.na(y0)]
      y0 = y0[!is.na(y0)]
      }
    Mean = function(x) sum(x)/(length(x)+lmb0)
    if(is.factor(x)){
      b = tapply(y0,x,Mean);
      if(anyNA(b)) b[is.na(b)] = 0
      g = b[x00]
    }else{
      b = cov(y0,x,use='pair')/(var(x,na.rm=T)+lmb0/length(y0));
      b = as.vector(b)
      g = b*x00
    }
    # Deregress
    if(Deregress){
      slope = c(crossprod(y0,g)/(crossprod(g)))
      xb = g*slope
      b0 = y0-xb
      g = b0+xb
    }
    return(list(h=g,b=b,e=y00-g))
  }
  
  # Structured random
  gws = function(e,i,...){
    x = data[[i]]
    if(iter==1){ e00 = e }else{ e00 = e + H[[i]]}
    e00 = e00 - mean(e00,na.rm = TRUE)
    e0 = tapply(e00,x,mean,na.rm=TRUE)[rownames(X[[i]])]
    # Fit WGR
    comn = intersect(names(e0),rownames(X[[i]]))
    h = alg(e0[comn],X[[i]][comn,],...)
    fit = c(h$hat)
    names(fit) = comn
    # Deregress
    if(Deregress){
      g = c(crossprod(fit,e0[comn])/(crossprod(fit)))
      xb = fit*g
      b0 = e0-xb
      fit = b0+fit*g
    }
    # Output
    hh = e0
    hh[comn] = fit
    hhh = hh[as.character(x)]
    if(anyNA(hhh)) hhh[is.na(hhh)] = 0
    res = e00 - hhh
    OUT = list(g=fit,h=hhh,b=h$b,e=res)
    return(OUT)}
  
  # Fit (loop)
  mu = mean(y,na.rm=T)
  e = yc = y-mu
  B[['Intercept']] = mu
  
  # Fit (loop)
  R2c = 0.5
  if(maxit!=10) pb = txtProgressBar(style = 3)
  for(iter in 1:maxit){
    
    # Fixed effects 
    if(FIX){
      for(i in fxd){
        if(iter>1) e = e+H[[i]]
          go = FCT(e,data[[i]])
          e = go$e
          B[[i]] = go$b
          H[[i]] = go$g
      }
    }
    
    # Random effects
    if(RND){
      # Coefficients
      for(i in rnd){
        # Structured
        if(i%in%ls(X)){
          go = gws(e,i,...)
          e = go$e
          B[[i]] = go$g
          H[[i]] = go$h
          A[[i]] = go$b
        }else{
          # Unstructured 
          if(iter>=1) e+H[[i]]
            go = RCT(e,data[[i]],LMB[i])
            e = go$e
            B[[i]] = go$b
            H[[i]] = go$h
        }
      }
      
      # VarComp & Lambda
      Error = var(e,na.rm=T)
      SSa = sapply(H, function(h) mean(yc*h,na.rm=T) )
      Vg = c(SSa,Error=Error)
      Vg[Vg<(0.001*Vy)] = 0.001*Vy
      Va = Vg[rnd]
      Ve = Vg['Error']
      LMB = Ve/Va; names(LMB) = names(Va)
    }
    
    # Print R2 and check convergence based on Ve
    if(maxit!=10) setTxtProgressBar(pb,iter/maxit)
    R2 = round(1-Ve/Vy,8)
    if(abs(R2c-R2)==0) break()
    R2c = R2
  }
  if(maxit!=10) close(pb)
  
  # Fit model
  fit = rep(B[['Intercept']],length(y))
  for(i in 1:length(H)) fit = fit+H[[i]]
  
  # Wrapping up
  Fitness = list(obs=y,hat=fit,res=e,fits=H)
  OUT = list(Fitness=Fitness,Coefficients=B)
  if(RND){
    wVar = Vg/sum(Vg)
    VC = list( VarComponents = round(Vg,6),
               VarExplained = round(wVar,6) )
    OUT[['VarComp']] = VC;
    if(length(X)>0) OUT[['Structure']] = A
  }  
  
  class(OUT) = 'mixed'
  return(OUT)}


#############################################################################################################                   
                   
mtmixed = function(resp, random=NULL, fixed=NULL, data, X=list(), maxit=10, init=10, regVC=FALSE){
  
  # Get y matrix
  k = length(resp) 
  rownames(data) = paste('Obs',1:nrow(data),sep='')
  y = data[,resp]
  Vy = apply(y,2,var,na.rm=TRUE)
  
  # Remove missing - (if missing for all traits)
  if(anyNA(y)){
    wna = which(rowMeans(is.na(y))==1);
    if(length(wna)>0){ y = y[-wna,]
    data = droplevels.data.frame(data[-wna,]) }}
  
  # Centralized y
  yc = apply(y,2,function(x) x-mean(x,na.rm=T))
  
  # Random effects
  if(!is.null(random)){
    rnd = attr(terms(random),"term.labels")
    LMB = apply(!is.na(y),2, function(q)  sapply(data[rnd], function(x) mean(table(droplevels(x[q])))) )
    LMB = matrix(LMB,nrow=length(rnd),ncol=k,dimnames=list(rnd,resp))
    df0 = apply(!is.na(y),2, function(q) sapply(data[rnd],function(x) ifelse(is.factor(x),length(unique(x[q])),1)) )
    # df0 = apply(!is.na(y),2, function(q) sapply(data[rnd],function(x) ifelse(is.factor(x), mean(table(x[q])) ,1)) )
    RND = TRUE }else{ RND = FALSE }
  
  # Fixed effects or intercept
  if(!is.null(fixed)){
    fxd = attr(terms(fixed),"term.labels")
    df = sum(sapply(data[fxd],function(x) ifelse(is.factor(x),length(unique(x)),1)))
    FIX = TRUE }else{ FIX = FALSE; df=1 }
  n = colSums(!is.na(y))
  B = H = list() # Store coefficients
  
  # Fixed Categorical & Continuous Term
  FCT = function(y,x){
    if(is.factor(x)){
      b = tapply(y,x,mean,na.rm=T);
      b[is.na(b)] = 0
      g = b[x]
      e = y-g
    }else{
      b = cov(y,x,use='pair')/var(x,na.rm=T);
      b = as.vector(b);
      g = x*b
      e = y-g }
    return(list(g=g,b=b,e=e))}
  
  # Random Categorical & Continuous Term
  RCT = function(yy,x,lmb){
    update = function(y0,x,lmb0){
      y00 = y0
      x00 = x
      if( anyNA(y0) ){
        x = x[!is.na(y0)]
        y0 = y0[!is.na(y0)] }
      Mean = function(x) sum(x)/(length(x)+lmb0)
      if(is.factor(x)){
        b = tapply(y0,x,Mean);
        if(anyNA(b)) b[is.na(b)] = 0
        g = b[x00]
      }else{
        b = cov(y0,x,use='pair')/(var(x,na.rm=T)+lmb0/length(y0));
        b = as.vector(b)
        g = b*x00 }
      return(list(h=g,b=b,e=y00-g))
    }
    outs = list()
    for(q in resp) outs[[q]] = update( y0=yy[,q],x,lmb0 = lmb[q])
    return(outs)
  }
  
  # MV parameters for structured terms
  MV = lapply(X, function(X){
    LIST = list()
    LIST[['b']] = matrix(0,ncol(X),k)
    rownames(LIST[['b']]) = colnames(X)
    colnames(LIST[['b']]) = resp
    kk = matrix(0,k,k,dimnames = list(resp))
    LIST[['iG']] = LIST[['vb']] = kk
    diag(LIST[['vb']]) = Vy/ncol(X)
    LIST[['iG']] = solve(LIST[['vb']])
    LIST[['ve']] = 0.5*Vy
    names(LIST[['ve']]) = resp
    return(LIST) } )
  
  # Structured random
  gws = function(e,i){
    x = data[[i]]
    if(iter==1){ e00 = e }else{ e00 = e + H[[i]][as.character(x),] }
    e0 = apply(e00,2,function(e00) tapply(e00,x,mean,na.rm=TRUE) )
    # Add for debugging
    if(any(!rownames(X[[i]])%in%rownames(e0))){
      ww = which(!rownames(X[[i]])%in%rownames(e0))
      ww = rownames(X[[i]])[ww]
      mm = matrix(NA,length(ww),k,dimnames=list(ww,resp))
      e0 = rbind(e0,mm)}
    e0 = e0[rownames(X[[i]]),]
    e0 = apply(e0,2,function(x) x-mean(x,na.rm=T))
    h = mtgsru(Y=e0,X=X[[i]],b=MV[[i]]$b,vb=MV[[i]]$vb,ve=MV[[i]]$ve,iG=MV[[i]]$iG,maxit=init)
    unlist(lapply(h,anyNA)); head(h$b)
    fit = h$hat
    rownames(fit) = rownames(X[[i]])
    comp = list(b=h$b,vb=h$vb,ve=h$ve,iG=h$iG,h2=h$h2,MSx=h$MSx)
    h = fit[as.character(x),]
    # dereg
    dereg = function(x,y){
      b = c(cov(x,y,use='p')/var(x,na.rm = T))
      xb = x*b
      mu = mean(y-xb,na.rm=T)
      ft = mu+xb
      return(ft)}
    for(drg in 1:k) h[,drg] = dereg(h[,drg],e00[,drg])               
    res = e00 - h
    colnames(fit) = resp
    OUT = list(g=fit,h=h,e=res,comp=comp)
    return(OUT)}
  
  # Get ready for loop
  e = y; R2c = 0.5;
  pb = txtProgressBar(style = 3)
  
  # Fit (loop)
  for(iter in 1:maxit){
    
    # Intercept
    if(iter==1){
      mu = apply(e,2,mean,na.rm=T)
      for(q in resp) e[,q] = e[,q]-mu[q]
      B[['Intercept']] = mu
    }else{
      mu = apply(e,2,mean,na.rm=T)
      for(q in resp) e[,q] = e[,q]-mu[q]
      B[['Intercept']] = mu+B[['Intercept']]
    }
    
    # Fixed effects 
    if(FIX){
      for(i in fxd){
        if(iter==1){
          go = apply(e,2,FCT,x=data[[i]])
          e = sapply(go, function(x) x$e )
          B[[i]] = sapply(go, function(x) x$b )
          H[[i]] = sapply(go, function(x) x$g )
        }else{
          E = e+H[[i]]
          go = apply(E,2,FCT,x=data[[i]])
          e = sapply(go, function(x) x$e )
          B[[i]] = sapply(go, function(x) x$b )
          H[[i]] = sapply(go, function(x) x$g )
        }}
    }
    
    # Random effects
    if(RND){
      
      # Coefficients
      for(i in rnd){
        
        # Structured
        if( i%in%ls(X) ){
          go = gws(e,i)
          e = go$e
          B[[i]] = go$g
          H[[i]] = go$h
          MV[[i]]$b = go$comp$b
          MV[[i]]$vb = go$comp$vb
          MV[[i]]$iG = go$comp$iG
          MV[[i]]$ve = go$comp$ve
          MV[[i]]$h2 = go$comp$h2
          MV[[i]]$MSx = go$comp$MSx
          rownames(MV[[i]]$iG) = resp
          colnames(MV[[i]]$vb) = resp
          
        }else{
          
          # Unstructured 
          if(iter==1){
            go = RCT( yy = e,x = data[[i]],lmb = LMB[i,])
            e = sapply(go, function(x) x$e )
            B[[i]] = sapply(go, function(x) x$b )
            H[[i]] = sapply(go, function(x) x$h )
          }else{
            E = e+H[[i]]
            go = RCT( yy = E,x = data[[i]],lmb = LMB[i,])
            e = sapply(go, function(x) x$e )
            B[[i]] = sapply(go, function(x) x$b )
            H[[i]] = sapply(go, function(x) x$h )
          }
        }
      }
      
      # VarComp & Lambda
      Error = colSums(e*e,na.rm=T)
      SSa = sapply(H, function(h) colSums(h*h,na.rm=T) )
      SS = cbind(SSa,Error)
      SS[SS<0] = 0.01
      if( regVC ) SS = sqrt(SS)
      wVar = t(apply(SS,1,function(x) x/sum(x) ))
      Vg = wVar*Vy
      Va = Vg[,rnd]/t(df0)*n
      Ve = Vg[,'Error'] * (n-df)/n
      LMB = t(Ve/Va)
      LMB = matrix(LMB,nrow=length(rnd),ncol=k,dimnames=list(rnd,resp))
      
      
    }
    
    setTxtProgressBar(pb,iter/maxit)
    # Print R2 and check convergence
    R2 = round(1-Ve/Vy,6)
    # cat('Iter ',iter,':',sep='')
    # cat(' R2 =',round(R2,3),'\n')
    if(sum(abs(R2c-R2))==0 & iter>=3){break(); setTxtProgressBar(pb, 1 )} 
    R2c = R2
    
    
  }
  close(pb)
  
  # Fit model
  fit = sapply(B[['Intercept']],function(x) rep(x,nrow(y)) )
  for(i in 1:length(H)) fit = fit+H[[i]]
  rownames(fit) = rownames(data)
  
  # Wrapping up
  Fitness = list(obs=y,hat=fit,res=e,fits=H)
  OUT = list(Fitness=Fitness,Coefficients=B)
  if(RND){
    VC = list( VarComponents = round(Vg,4),
               VarExplained = round(wVar,3) )
    OUT[['VarComp']] = VC;
    if(length(X)>0) OUT[['Structure']] = MV
  }
  
  # Return
  class(OUT) = 'mixed'
  return(OUT)}

#############################################################################################################                   

SimY = function(Z, k=5, h2=0.5, GC=0.5,  seed=123, 
                unbalanced=FALSE, PercMiss=0, BlkMiss=FALSE){
  
  # Store inputs
  trueVal = list(h2=h2,GC=GC,seed=seed)
  n = nrow(Z)
  p = ncol(Z)
  
  # Pick a dataset
  set.seed(seed)
  
  # GC
  if(is.matrix(GC)){
    k = ncol(GC)
    G0 = GC
  }else{
    numGC = (k*(k-1))/2
    GC = rep(GC,numGC)
    G0 = diag(k)
    G0[lower.tri(G0)] = GC
    G0[upper.tri(G0)] = t(G0)[upper.tri(t(G0))]
    GC = mean(GC)
  }
  
  # Genetic parameters
  if(length(h2)==k){
    h20 = h2
  }else{
    h20 = h2
    h2 = rep(h20,k)
  }
  
  # Sample effects
  alpha = 1/sum(apply(Z,2,var))
  trueVal['scaleG'] = alpha
  Vb = G0*alpha
  ev = eigen(Vb, symmetric = TRUE)
  ev$values = ifelse(ev$values<0.0001,0.0001,ev$values)
  UD = ev$vectors %*% diag(sqrt(ev$values))
  beta = matrix(rnorm(p * k), nrow = p)
  trueBeta = UD %*% t(beta)
  
  # True breeding values
  tbv = Z %*% t(trueBeta)
  
  # Residual variances and phenotypes
  ve = (1-h2)/h2;
  E = sapply(ve,function(ve) rnorm(n,0,sqrt(ve)))
  Y = 10 + tbv + E
  colnames(Y) = paste('y',1:k,sep='')
  
  # Unbalance
  if(unbalanced){
    Miss = sample(1:k,n,T)
    for(i in 1:k) Y[which(Miss!=i),i] = NA
  }
  
  # Percentage missing
  if(PercMiss>0){
    if(BlkMiss){
      Ym = matrix(F,10,10)
      col_index = sort(rep(1:10,length.out=ncol(Y)))
      row_index = sort(rep(1:10,length.out=nrow(Y)))
      Ym[sample(100,100*PercMiss)] = T
      for(i in 1:10) for(j in 1:10) if(!Ym[i,j]) Y[which(row_index==i),which(col_index==j)] = NA
    }else{
      Y[sample(length(Y),length(Y)*PercMiss)] = NA
    }
  }
  
  # Output
  return(list(Y=Y,tbv=tbv,settings=trueVal))
  
}

#############################################################################################################                   

SimZ = function(ind = 500, snp = 500, chr=2, F2=TRUE, rec = 0.01){
  # ind = 500; snp = 500; chr=2; F2=TRUE; rec = 0.01
  makeZ = function(...){
    Z = sapply(1:ind,function(a){
      z = c(F,sample(c(T,F),snp-1,T,c(rec,1-rec))); recs = which(z);
      for(i in recs){z[i:snp]=(!z[i-1])};
      z = z*sample(c(T,F),1); if(runif(1)>0.5){ z = !z};
      return(z)}); Z = t(Z); return(Z)}
  makeZF2 = function(...){ makeZ()+makeZ()  }
  if(F2){
    Z = do.call(cbind,lapply(1:chr,makeZF2)) }else{
    Z = do.call(cbind,lapply(1:chr,makeZ)) } 
  return(Z)}
             
#############################################################################################################                   

SimGC = function(k=50, MIX=TRUE,
                 BLK_UNS=FALSE,BLK=FALSE,
                 GRAD_UNS=FALSE,GRAD=FALSE,
                 BLK_GRAD=FALSE){
  kk = round(sqrt(k))
  k0 = ceiling(sqrt(k))
  k00 = floor(sqrt(k))
  if(MIX){
    GC0 = (exp(-as.matrix(dist(1:k0,method = 'man')/(0.5*k0) ))-0.5)*2
    GC = kronecker(GC0,matrix(1,k0,k0))[1:k,1:k]
    GC = GC + 0.5*(1-as.matrix(dist(matrix(runif(k*kk),k,kk ))))
    GC = cov2cor(GC)
    diag(GC)=1; GC = cov2cor(GC)
  }else if(BLK_GRAD){
    GC0 = (exp(-as.matrix(dist(1:k0,method = 'man')/(0.5*k0) ))-0.5)*2
    GC = kronecker(GC0,matrix(1,k0,k0))[1:k,1:k]*0.75
    diag(GC)=1; GC = cov2cor(GC)
  }else if(BLK_UNS){
    GC0 = 1-as.matrix(dist(matrix(runif(k0*k00),k0,k00)))
    GC = kronecker(GC0,matrix(1,k0,k0))[1:k,1:k]
    GC = GC + (1-as.matrix(dist(matrix(runif(k*kk),k,kk ))))
    GC = cov2cor(GC)
  }else if(BLK){
    GC0 = 1-as.matrix(dist(matrix(runif(k0*k00),k0,k00)))
    GC = kronecker(GC0,matrix(1,k0,k0))[1:k,1:k]*0.75
    diag(GC)=1; GC = cov2cor(GC)
  }else if(GRAD_UNS){
    GC = (exp(-as.matrix(dist(1:k,method = 'man')/(0.5*k) ))-0.5)*2
    GC = 0.5*(GC+1-as.matrix(dist(matrix(runif(k*kk),k,kk))))
  }else if(GRAD){
    GC = (exp(-as.matrix(dist(1:k,method = 'man')/(0.5*k) ))-0.5)*2
  }else{
    GC = 1-as.matrix(dist(matrix(runif(k*kk),k,kk)))
  }
  E = eigen(GC,T)
  if(any(E$values<0.01)){
    E$values = ifelse(E$values>0.01,E$values,0.01)
    GC = E$vectors %*% diag(E$values) %*% t(E$vectors)
    GC = cov2cor(GC)}
  GC = round(GC,2)
  return(GC)
}

#############################################################################################################                   

# Main function
mm = function(y,random=NULL,fixed=NULL,data=NULL,
              M=NULL,bin=FALSE,AM=NULL,it=20,verb=TRUE,
              FLM=FALSE,wgtM=FALSE,cntM=FALSE,nPc=3){
  
  # Datasets for testing and debugging
  # require(eMM3); data(cateto);
  # random=~ID+Tester+Env:ID+Tester:ID+Env:Tester;
  # fixed=~Env; data=dt; M=list(ID=Geno);
  # bin=FALSE; AM=NULL; it=10; verb=TRUE
  # FLM=FALSE; wgtM=TRUE; cntM=TRUE; nPc=3
  
  # require(Matrix); require(eMM3)
  # random=~imm1+imm2;
  # fixed=NULL; data=dta; M=list(imm1=Z1,imm2=Z2);
  # bin=FALSE; AM=NULL; it=10; verb=TRUE
  # FLM=FALSE; wgtM=TRUE; cntM=TRUE; nPc=3
  
  # Base algorithm settings
  if(nPc<2) nPc=2
  if(!is.null(data)) data = droplevels.data.frame(data); as = 0
  if(!verb){ cat = function(...) NULL }
  
  # Spline
  if(!is.null(AM)){
    SPACE = function(y,NN) (NN$X%*%y)[,1]/NN$n # Fit
    data = data.frame( data, obs_num_space = 1:nrow(data)) # Track what was removed after QC
    spc = TRUE
  }else{
    spc = FALSE
  } 
  
  ######################
  ### Design Matrices ##
  ######################
  
  # Response variable
  if(!is.null(data)) y = data[[deparse(substitute(y))]]
  # y = data$Yield; # For debugging with example dataset
  # y = dta$NY
  
  Y = y
  sdy = sqrt(var(Y,na.rm = T))
  my = mean(Y,na.rm = T)
  upLim = my+5*sdy
  loLim = my-5*sdy
  w = which(Y<loLim | Y>upLim)
  if(length(w)>0){
    cat("Phenotypes contained",length(w),"outlier(s)\n")
    y[w] = NA
  }
  
  # Fixed effects or intercept
  if(!is.null(fixed)){
    cat("Setting fixed effect design matrices\n")
    fixed = update(fixed,~.-1)
    X = sparse.model.matrix(fixed,data = data)
    #XX = as(crossprod(X),"dgCMatrix")
    XX = as(as(as(crossprod(X), "dMatrix"),"generalMatrix"), "CsparseMatrix")
    
    B = rep(0,ncol(X)); names(B) = colnames(X)
    Lmb_X = (1e-12)*sum(colMeans(X^2)-colMeans(X)^2)
    mu = 0
    FIX = TRUE
  }else{
    mu = 0
    FIX = FALSE
  }
  
  # Random effects
  if(!is.null(random)){
    cat("Setting random effect design matrices\n")
    rnd0 = rnd = attr(terms(random),"term.labels")
    rndInt = grep(':',rnd)
    # Temporarily drop interactions with Markers
    if(length(rndInt)>0){
      tmp = rnd[grep(':',rnd)]
      ms = ifelse(!is.null(M), paste(ls(M),sep='|'), 'NothingReally')
      tmp = grep(ms,tmp,invert=F,value=T)
      rnd = rnd[!rnd%in%tmp]}
    nr0 = length(rnd0)
    nr = length(rnd)
    rndInt = nr0>nr
    Z = list()
    for(i in 1:nr){
      Z[[i]]=sparse.model.matrix(formula(paste('~',rnd[i],'-1')),data=data,drop.unused.levels=TRUE)
      Z[[i]] = as(as(as(Z[[i]], "dMatrix"),"generalMatrix"), "CsparseMatrix")
    } 
    names(Z) = rnd
    U = lapply(Z, function(x){z=rep(0,ncol(x));names(z)=colnames(x);return(z)}  )
    RND = TRUE
  }else{
    nr = nr0 = 0
    rndInt = nr0>nr
    rnd = NULL
    RND = FALSE
  }
  
  # Spatial to Random effects
  if(spc & RND){
    cat("Adding adjacency matrix to random effects\n")
    rnd = c(rnd,'spatial')
    nr = nr+1
    Z[[nr]] = AM$X[data$obs_num_space,data$obs_num_space]
    names(Z)[nr] = 'spatial'
    U[[nr]] = rep(0,nrow(data))
    names(U)[nr] = 'spatial'
  }
  
  # Pre-cook marker terms to redesign matrices for interactions
  if(!is.null(M) ){
    ms = ls(M)
    # Second level logic
    if(cntM | any(sapply(M, anyNA)) ){
      cat('Centralize marker matrices\n')
      # QC'ish
      for(i in ms){
        # Assign correct class
        if(!is.matrix(M[[i]])) M[[i]] = data.matrix(M[[i]])
        # Centralize
        M[[i]] = apply(M[[i]],2,function(x) x-mean(x,na.rm=T))
        # Impute
        if(anyNA(M[[i]])) M[[i]][is.na(M[[i]])] = 0
      }
    }
  }
  
  # Adding back design matrices for interactions
  if(rndInt){
    keyInteractions = rnd0[!rnd0%in%rnd]
    mainTerms = grep(gsub(':','|',paste(keyInteractions,collapse=':')),ms,value=T)
    # Eigen decomposing key main terms
    M_SVD = lapply(M[mainTerms], svd, nu=nPc, nv=nPc)
    Zpc = list()
    for(i in mainTerms){
      cat('Extracting PCs of',i,'\n')
      rownames(M_SVD[[i]]$u) = rownames(M[[i]])
      rownames(M_SVD[[i]]$v) = colnames(M[[i]])
      # Rescale for more meaningful variance components
      k = sqrt(nPc/min(dim(M[[i]])))
      M_SVD[[i]]$d = M_SVD[[i]]$d*k
      M_SVD[[i]]$v = M_SVD[[i]]$v*k
      # Add a value for missing ids
      M_SVD[[i]]$u = rbind(M_SVD[[i]]$u,miss=0)
      # Scale eigenvectors based on their eigenvalues
      M_SVD[[i]]$u = M_SVD[[i]]$u %*% diag(M_SVD[[i]]$d[1:ncol(M_SVD[[i]]$u)])
      tmp_ids = as.character(data[[i]])
      tmp_ids[ ! tmp_ids%in%rownames(M[[i]]) ] = 'miss'
      Zpc[[i]] = M_SVD[[i]]$u[tmp_ids,]
      if(anyNA(Zpc[[i]])) Zpc[[i]][is.na(Zpc[[i]])]=0
      colnames(Zpc[[i]]) = paste('_pc',1:nPc,sep='')}
    # Create interaction matrices
    AxB = function(i){
      tmp = list()
      individual_terms = strsplit(i,':')[[1]]
      for(j in individual_terms){
        if(j%in%ls(Zpc)){ tmp[[j]] = Zpc[[j]] }else{ tmp[[j]] = data[[j]] }}
      AB = model.matrix(as.formula(paste('~',i,'-1')),data=tmp)
      #AB = as(AB,'dgCMatrix')
      #AB = as(as(as(AB, "dMatrix"), "generalMatrix"), "CsparseMatrix")
      AB = as(AB,"CsparseMatrix")
      # AB = sparse.model.matrix(~tmp[[1]]:tmp[[2]],drop.unused.levels=TRUE)
      return(AB)}
    # Add terms back to: nr, rnd, Z, U
    for(i in keyInteractions){
      cat('Creating',i,'matrix\n')
      Z[[i]] = AxB(i)
      U[[i]] = rep(0,ncol(Z[[i]]))
      names(U[[i]]) = colnames(Z[[i]])
      rnd = c(rnd,i)
      nr = nr+1 
    }
  }else{
    keyInteractions = c()
  }
  
  # Missing values
  if(anyNA(y)){
    y0 = y;
    if(FIX) X0 = X
    if(RND) Z0 = Z
    wNA = which(is.na(y))
    y = y[-wNA];
    if(FIX) X = X[-wNA,]
    if(RND) Z = lapply(Z,function(z) z[-wNA,] )
    for(i in 1:nr){
      cm = colSums(Z[[i]]);
      if(any(cm==0)){
        w = which(cm==0);
        Z[[i]] = Z[[i]][,-w]
        Z0[[i]] = Z0[[i]][,-w]
        U[[i]] = U[[i]][-w]
      }
    }
    MIS = TRUE
  }else{MIS = FALSE}
  n0 = length(Y)
  n = length(y)
  
  # Binomial
  if(bin){
    cat("Logit tranformation\n")
    MinY = min(y)-1
    MaxY = max(y)+1
    ry = c(MinY,MaxY)
    y = (y-MinY)/(MaxY-MinY)
    y = log(y/(1-y))
  }
  
  # Random parameters for regularization
  if(RND){
    #ZZ = lapply(Z, lapply(Z, function(z) as(crossprod(z),"dgCMatrix") )
    ZZ = lapply(Z, function(z){as(as(as(crossprod(z), "dMatrix"),"generalMatrix"), "CsparseMatrix") }  )
    Z_cxx = lapply(Z, function(X) sum(colMeans(X^2)-colMeans(X)^2) )
    Lmb_Z = mapply( function(cxx,h2) cxx*((1-h2)/h2), cxx = Z_cxx, h2=0.5)
    trAC22 = sapply(ZZ,function(x) mean(1/(diag(x)+1)))
    fxx = ifelse(FIX,sum(1/diag(XX)),0)
  }
  
  # Starting value for variances
  Vy = var(y)
  Ve = Vy*0.5
  if(RND){ Va = rep(var(y)*0.5,nr); names(Va) = rnd } 
  e = y
  
  # Work on genotypes and structured terms
  if(!is.null(M)){
    
    # Name of structured terms
    ms = ls(M)
    ZXX = lapply(ms, function(z){  ZXX = diag(ZZ[[z]]); names(ZXX)=sub(z,'',colnames(ZZ[[z]])); return(ZXX) })
    names(ZXX) = ms
    
    # Empty lists to store structure stuff
    gFit = gHat = M_xx = M_cxx = Beta = UM = LMB = Gh2 = nrep = list()
    
    # Loop across matrices of M to collect parameters and indeces
    cat("Computing parameters of M terms \n")
    for(i in ms){
      
      if(wgtM){
        
        # Get weights for M
        wgts = sqrt(ZXX[[i]][rownames(M[[i]])])
        if(anyNA(wgts)) wgts[is.na(wgts)] = 0
        
        # Compute cross-products with weights
        M_xx[[i]] = apply(M[[i]],2,crossprod)
        M_cxx[[i]] = sum(M_xx[[i]])/nrow(M[[i]])
        M_xx[[i]] = apply(M[[i]],2,function(x) crossprod(x*wgts) )
        M_xx[[i]] = ifelse(M_xx[[i]]!=0,M_xx[[i]],0.0001)
        M[[i]] = apply(M[[i]],2,function(x){x*wgts})
        
      }else{
        
        # Compute cross-products without weights
        M_xx[[i]] = apply(M[[i]],2,crossprod)
        M_cxx[[i]] = sum(M_xx[[i]])/nrow(M[[i]])
        M_xx[[i]] = ifelse(M_xx[[i]]!=0,M_xx[[i]],0.0001)
        
      }
      
      # Controls
      if( is.null(rownames( M[[i]] )) ) stop(paste("Matrix M",i,"does not have row names to indentify levels"))
      if( !i%in%rnd ) stop(paste("No random effect called",i,"was declared"))
      
      # Rename and match levels
      colnames(Z[[i]]) = gsub(i,"",colnames(Z[[i]]))
      names(U[[i]]) = gsub(i,"",names(U[[i]]))
      proportion_genotyped = round(mean(colnames(Z[[i]])%in%rownames( M[[i]] ))*100,2)
      cat(' (',proportion_genotyped," percent of ",i," have markers)\n",sep='')
      if(proportion_genotyped==0) stop(paste(i,"levels from data and M do not match"))
      
      # Match M to levels
      mm = mean(rownames(M[[i]])%in%colnames(Z[[i]]))*100
      if(mm<100){
        cat(' Note: ',round(100-mm)," percent of ",i," markers have no observations\n",sep='')
        UM[[i]] = M[[i]][ which( ! rownames(M[[i]]) %in% colnames(Z[[i]]) ) , ]
        M[[i]] = M[[i]][ which( rownames(M[[i]]) %in% colnames(Z[[i]]) ) , ]
      }  
      
      # Mapping matrix
      wGen = which(!rownames(M[[i]])%in%colnames(Z[[i]]))
      if(length(wGen)!=0) M[[i]] = M[[i]][-wGen,]
      gFit[[i]] = rep(0,ncol(Z[[i]])); names(gFit[[i]]) = colnames(Z[[i]])
      gHat[[i]] = rep(0,n)
      
      # Other fitting parameters
      Beta[[i]] = rep(0,ncol(M[[i]]))
      LMB[[i]] = rep(M_cxx[[i]],ncol(M[[i]]))
      nrep[[i]] = median(diag(ZZ[[i]]))
      
    }
    
  }else{ms=c();UM=list()}
  
  # Remove overall mean
  m0 = mean(e)
  mu = mu+m0
  e = e-m0
  
  ###########
  ### Loop ##
  ###########
  
  if(verb) pb = txtProgressBar(style = 3)
  resids0 = e
  conv = 0
  
  for(i in 1:it){
    
    # Fixed coefficient update
    if(FIX) upFix = GS2EIGEN(e,X,B,XX,Lmb_X)
    
    # Randomc oefficient update
    if(RND){
      for(j in 1:nr){
        
        # UPDATE TERMS WITH MARKERS
        if(rnd[j]%in%ms){
          
          ## Collect information
          w = rnd[j]
          e = e + gHat[[w]]
          
          # Mapping function
          #Gmap0 = c(crossprod(Z[[w]],e)[,1])/diag(ZZ[[w]])
          Gmap0 = c(crossprod(Z[[w]],Matrix(e) )[,1])/diag(ZZ[[w]])
          Gmap = Gmap0[rownames(M[[w]])]
          
          # Whole-genome regression
          tmpY = c(Gmap)
          Gmap = c(Gmap-M[[w]]%*%Beta[[w]])
          if(FLM){
            G_update = GSFLM(tmpY,Gmap,M[[w]],Beta[[w]],LMB[[w]],M_xx[[w]],M_cxx[[w]],maxit=20)
          }else{
            G_update = GSRR(tmpY,Gmap,M[[w]],Beta[[w]],LMB[[w]],M_xx[[w]],M_cxx[[w]],maxit=20)
          }
          Gh2[[w]] = G_update$h2
          
          # Update coefficients of non-M
          e_tmp = e * 1.0
          u_tmp = U[[w]] * 1.0
          GS2EIGEN(e_tmp,Z[[w]],u_tmp,ZZ[[w]],Lmb_Z[w])
          gFit[[w]] = u_tmp # - mean(u_tmp)
          
          # Update coefficients
          tmp = c(M[[w]]%*%Beta[[w]])
          gFit[[w]][rownames(M[[w]])] = tmp # - mean(tmp)
          U[[w]] = gFit[[w]]
          
          # Update fitted values and residuals
          gHat[[w]] = (Z[[w]]%*%U[[w]])[,1]
          e = e - gHat[[w]]
          
        }else{
          
          # UPDATE TERMS WITHOUT MARKERS
          
          # Update coefficients
          GS2EIGEN(e,Z[[j]],U[[j]],ZZ[[j]],Lmb_Z[j])
          
        }
        
      }
    }
    
    # Update overall intercept
    m0 = mean(e)
    mu = mu+m0
    e = e-m0
    
    # Variance components update
    Ve = abs(crossprod(e,y)[1,1]/(n-ifelse(FIX,ncol(X),1)))
    
    if(RND){
      trAC22 = (unlist(mapply(function(ZZ,lmb) sum(1/(diag(ZZ)+lmb)),ZZ=ZZ,lmb=as.list(Lmb_Z))))
      Va = sapply(U,crossprod)/(sapply(U,length)-trAC22*Lmb_Z)
      Lmb_Z = abs(Ve/Va)
    }
    
    # Convergence and update progress bar
    conv = log(sum((e-resids0)^2))/log(1e-4)
    resids0 = e
    if(verb) setTxtProgressBar(pb, max(i/it,conv))
    if( conv>=1 & ( !is.null(M) | spc ) ) break()
    
  }
  
  # Close pregression bar 
  if(verb) close(pb)
  
  ##################
  ### Wrapping up ##
  ##################
  
  # Empty list of outputs
  fx = coef = list()
  class(fx) = "FLMSS"
  
  # Add fixed effects
  coef[['mu']] = mu
  if(FIX){ coef[['Fxd']] = B } 
  
  # Add random effects
  if(RND){
    coef = c(coef,U)
    fx[['VC']] = round(c(Va,Ve=Ve),4)
  }
  
  # Add marker effects
  if(!is.null(M)){
    names(Beta) = ms
    for(i in ms) names(Beta[[i]]) = colnames(M[[i]])
    fx[['Mrk']] = Beta
  }
  
  # Marker interactions
  if(rndInt){
    # Add PCs to output
    fx[['PCs']] = M_SVD
    # Loop across interactions to best allocate output
    for(i in keyInteractions){
      tmp_terms = strsplit(i,':')[[1]]
      num_terms = length(tmp_terms)
      countMs = sum(tmp_terms%in%mainTerms)
      # Fetch markers effects under simple compound symmetry
      if( countMs==1 & num_terms==2 ){
        
        is_M = which(tmp_terms%in%mainTerms)
        if(is_M==1){  tmp0 = t(matrix(coef[[i]],nrow=nPc))
        rownames(tmp0) = unique(gsub('^.+:','',names(coef[[i]])))} 
        if(is_M==2){  tmp0 = matrix(coef[[i]],ncol=nPc)
        rownames(tmp0) = unique(gsub(':.+','',names(coef[[i]])))}
        
        SnpEff = M_SVD[[ tmp_terms[is_M] ]]$v %*% t(tmp0)
        fx[['Mrk']][[i]] = data.frame(SnpEff)
        
      }
    }
  }
  
  # Add coefficients to output
  fx[['Coef']] = coef
  
  # Fitting model
  GOF = list()
  
  ######## WITH MISSING
  if(MIS){
    fit = rep(0,length(y0))
    if(FIX){
      fit=fit+X0%*%B+mu
    }else{
      fit=fit+mu
    }
    GOF[['Fixed']]=as.vector(fit)
    if(RND){
      for( i in rnd ){
        # Unobserved with markers
        if( (!is.null(M)) & (i%in%ls(UM)) ){
          # Fit unobserved
          ufit = (UM[[i]]%*%Beta[[i]])[,1]
          U[[i]] = c(U[[i]],ufit)
          # Check if unobserved was in the data
          all_obs = as.character(data[[i]])
          all_fits = U[[i]][all_obs]
          if(anyNA(all_fits)) all_fits[is.na(all_fits)] = 0
          # Fit
          GOF[[i]] = as.vector(all_fits)
          fit = fit + all_fits
        }else{
          tmp = Z0[[i]]%*%U[[i]]
          GOF[[i]] = as.vector(tmp)
          zu = tmp
          fit = fit + zu
        }
      } 
    }
    Obs = y0
    
    ######## NO MISSING
  }else{
    fit = rep(0,n)
    if(FIX){
      tmp = fit+X%*%B+mu
    }else{
      tmp = fit+mu  
    }
    fit=fit+mu
    GOF[['Fixed']]=as.vector(fit)
    if(RND){
      for(i in rnd){
        tmp = Z[[i]]%*%U[[i]]
        GOF[[i]]=as.vector(tmp)
        zu = tmp
        fit = fit+zu
      }}
    Obs = y
  }
  
  # Store fitted values
  GOF = cbind(data.frame(Observed=as.vector(Obs),
                       Predicted=as.vector(fit),
                       Residuals=as.vector(Obs-as.vector(fit))),GOF)
  fx[['GOF']] = GOF
  
  # Store spatial info
  if(spc){
    tmp = cbind(AM$dt,sp=U$spatial)
    sp_tmp = by(tmp,tmp$blk,function(Q){sparseMatrix(i=Q$row,j=Q$col,x=Q$sp)})
    fx[['SpVar']] = sp_tmp
  }
  
  # If binary
  if(bin){
    b_fit = exp(fit)/(1+exp(fit));
    b_fit = b_fit*(ry[2]-ry[1])+ry[1]
    fx[['Bin']] = round(as.numeric(b_fit))
  } 
  
  # Log-likelihood, Adj R2, BIC
  LogLik = sum(dnorm(as.numeric(y),as.numeric(fit),sqrt(as.numeric(Ve)),T))
  px = 1+ifelse(FIX,ncol(X),0)+ifelse(RND,length(Z),0)
  Adj_R2 = 1 - Ve/var(Y,na.rm=T)
  BIC = log(n)*px + log(Ve)*n
  Stat = c(LogLik=LogLik,Adj_R2=Adj_R2,BIC=BIC)
  if(!is.null(M)){
    H2 = unlist(Gh2)
    names(H2) = paste('h2',names(H2),sep='.')
    Stat = c(Stat,H2)
  }
  fx[['Stat']] = round(Stat,3)
  
  # Return
  return(fx)
}

# Spline function
NNS = function(blk,row,col,rN=2,cN=2){
  if(is.factor(blk)|is.character(blk)) blk = as.numeric(as.character(blk))
  if(is.factor(row)|is.character(row)) row = as.numeric(as.character(row))
  if(is.factor(col)|is.character(col)) col = as.numeric(as.character(col))
  e = NNSEARCH(blk,row,col,rN,cN)
  e[e==0]=NA
  E = apply(e,1,function(x) x[!is.na(x)] )
  n = sapply(E, length)
  J = unlist(E)
  e = cbind(1:nrow(e),e)
  E = apply(e,1,function(x) x[!is.na(x)] )
  I = unlist(sapply(E,function(x) rep(x[1],length(x)-1)  ))
  X = sparseMatrix(i=I,j=J,x=1)
  colnames(X) = 1:ncol(X)
  rownames(X) = 1:ncol(X)
  dt = data.frame(blk=blk,row=row,col=col)
  OUT = list(X=X,n=n,dt=dt)
  class(OUT) = "NNS"
  return(OUT)}

# Prediction function
predict_FLMSS = function(x, data=NULL, M=NULL){
  
  # Nothing provided
  if(is.null(data)&is.null(M)) return(x$GOF$Predicted)
  
  #########################
  # If data M is provided #
  #########################
  
  if(!is.null(data)){
    # Design matrices
    cat('Getting design matrices\n')
    Xs0 = lapply(names(data), function(aaa) sparse.model.matrix(as.formula(paste('~',aaa,'-1')),data=data) )
    Xs = Xs0[[1]]
    if(length(Xs0)>1){
      for(i in 1:length(Xs0)){
        Xs = cbind(Xs,Xs0[[i]])
      }
    } 
    # Coefficients
    bs = unlist(x$Coef)
    # Collect some more coefficients
    if( (!is.null(M)) & ("Mrk"%in%ls(x)) ){
      for(i in names(M)){
        if(i %in% ls(x$Mrk)){
          cat('Computing marker coefficients',i,'\n')
          tmp = c( M[[i]] %*% x$Mrk[[i]] )
          names(tmp) = paste(i,rownames(M[[i]]),sep='.')
          bs = c(tmp,bs)
          bs = bs[!duplicated(names(bs))]
        }
      }
    }
    # Match overlapping
    #names(bs) = gsub('^.+\\.','',names(bs))
    names(bs) = gsub('-|:|^Fxd|\\.| ','',names(bs))
    #colnames(Xs) = gsub('^.+\\.','',colnames(Xs))
    colnames(Xs) = gsub('-|:|\\.| ','',colnames(Xs))
    i = intersect(names(bs),colnames(Xs))
    # Prediction
    if(length(i)>0){
      cat('Fitting',length(i),'coefficients\n')
      bs = bs[i]
      Xs = data.matrix(Xs[,i])
      prd = x$Coef$mu + c(Xs%*%bs)
    }else{
      prd = rep(x$Coef$mu,nrow(data))
    }
  }
  
  ######################
  # Only M is provided #
  ######################
  
  if( (!is.null(M)) & (is.null(data)) & ("Mrk"%in%ls(x)) ){
    cat('Predicting marker effects\n')
    prd = list()
    # Additive
    for(i in names(M)){
      if(i %in% ls(x$Mrk)){
        cat('Computing',i,'\n')
        prd[[i]] = c( M[[i]] %*% x$Mrk[[i]] ) + x$Coef$mu
        names(prd[[i]]) = rownames(x$Mrk[[i]])
      }
    }
    # Compound symmetry interactions
    if(any( sapply(x$Mrk,class)=="data.frame" )){
      interactions = ls(x$Mrk)[which(sapply(x$Mrk,class)=="data.frame")]
      for(i in interactions){
        Terms = strsplit(i,':')[[1]]
        if(any(Terms%in%ls(M))){
          cat('Computing',i,'\n')
          MainTerm = Terms[which(Terms%in%ls(M))]
          prd[[i]] = apply(x$Mrk[[i]],2,function(x) M[[MainTerm]]%*%x )
          rownames(prd[[i]]) = rownames(M[[MainTerm]])
        }
      }
    }
    # If there was only one possible solution
    if(length(prd)==1) prd=prd[[1]]
  }
  # Return
  return(prd)
}

#############################################################################################################                   

# MegaSEM
SEM = function(Y,Z,PCs=ifelse(RndLatSp,min(30,ncol(Y)),3),
               TOI=NULL,Beta0=NULL,RndLatSp=TRUE){
  cat('Using',PCs,'latent spaces\n')
  k = ncol(Y)
  Mu = colMeans(Y,na.rm=T)
  NamesY0 = colnames(Y)
  toinames = NamesY0[TOI]
  cat('Step 1 - UV\n')
  if(is.null(Beta0)){
    cat('  |')
    Beta0 = sapply(1:k,function(i){
      y = Y[,i]
      w = which(!is.na(y))
      yy = y[w]
      xx = Z[w,]
      beta = MRR3F(matrix(yy),xx,NonLinearFactor=2,weight_prior_h2=0.1)
      if((100*i/k)%%10==0) cat('===')
      return( c(beta$b) )})
    cat('|\n')
  }else{ cat('Skip step1\n') }
  rownames(Beta0) = colnames(Z)
  colnames(Beta0) = NamesY0
  cat('Step 2 - SVD G\n')
  G = Z %*% Beta0; 
  E = (EigenBDCSVD( G )$V)[,1:min(PCs,ncol(Y))]
  cat('Step 3 - SEM\n')
  k = ncol(Y)
  if(is.null(TOI)){ toi = 1:k  }else{ toi = TOI } 
  cat('  |')
  MvBeta = sapply(toi,function(i){
    if( (100*which(toi==i)/length(toi))  %% 10  ==0 ) cat('===')
    y = Y[,i]
    w = which(!is.na(y))
    yy = y[w]
    xx = Z[w,]
    if(RndLatSp){
      X = cbind(G[w,] %*% E)
      b0 = MRR3F(matrix(yy),X,weight_prior_h2=0.1)
      b1 = MRR3F(matrix(yy-b0$hat),xx,NonLinearFactor=2,weight_prior_h2=0.1)
      betaf = c( Beta0 %*% E %*% b0$b ) + c(b1$b)
    }else{
      X = cbind(1,G[w,] %*% E)
      beta = MLM(matrix(yy),X,xx)
      betaf = c( Beta0 %*% E %*% beta$b[-1,] ) + c(beta$u)
    }
    return(betaf)})
  cat('|\n')
  # VC
  G = Z %*% MvBeta
  if(!is.null(TOI)){
    Yc = apply(Y[,toi],2,function(x)x-mean(x,na.rm=T))
    pa = cor(G,Y[,toi],use='p'); GC = 0.5*(pa+t(pa))
    for(i in 1:ncol(G)) G[,i]=G[,i]+Mu[i]
    h2 = 1-c(colMeans((Y[,toi]-G)*Yc,na.rm=T))/apply(Y[,toi],2,var,na.rm=T)
  }else{
    Yc = apply(Y,2,function(x)x-mean(x,na.rm=T))
    pa = cor(G,Y,use='p'); GC = 0.5*(pa+t(pa))
    for(i in 1:ncol(G)) G[,i]=G[,i]+Mu[i]
    h2 = 1-c(colMeans((Y-G)*Yc,na.rm=T))/apply(Y,2,var,na.rm=T)}
  rownames(MvBeta) = colnames(Z)
  if(is.null(TOI)){ colnames(MvBeta) = colnames(Y)  }else{ colnames(MvBeta) = colnames(Y)[TOI] } 
  out = list(Univ=Beta0,SEM=MvBeta,Mu=Mu)
  out = list(mu=Mu,b=MvBeta,GC=GC,hat=G,h2=h2)
  return(out)}

#############################################################################################################                   

mrr = function(Y,X,...) MRR3(Y,X,...)

mrr_float = function(Y,X,...) MRR3F(Y,X,...)

mkr = function(Y,K,...) MRR3(Y,K2X(K),...)

mkr2X = function(Y,K1,K2,...) mrr2X(Y,K2X(K1,...),K2X(K2,...))
