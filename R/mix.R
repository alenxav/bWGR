mixed = function(y,random=NULL,fixed=NULL,data=NULL,X=list(),alg=emML,maxit=10,Deregress=TRUE,...){
  
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
  U = list()
  G = list()
  
  # Fixed effects or intercept
  if(!is.null(fixed)){
    fxd = attr(terms(fixed),"term.labels")
    df = sum(sapply(data[fxd],function(x) ifelse(is.factor(x),length(unique(x)),1)))
    FIX = TRUE }else{ FIX = FALSE; df=1 }
  B = list()
  n = length(y)
  
  # Fixed Categorical & Continuous Term
  FCT = function(y,x){
    if(is.factor(x)){ b = tapply(y,x,mean); e = y-b[x]
    }else{ b = crossprod(y,x)/crossprod(x); b=as.vector(b); e = y-x*b }
    return(list(b=b,e=e))}
  FCT_update = function(y,x,b){
    if(is.factor(x)){ m = tapply(y,x,mean); e = y-m[x]
    }else{ m = crossprod(y,x)/crossprod(x); m=as.vector(m); e = y-x*m }
    b=m+b; return(list(b=m,e=e))}
  
  # Random Categorical & Continuous Term
  RCT = function(y,x,lmb){
    Mean = function(x,lmb) sum(x)/(length(x)+lmb)
    if(is.factor(x)){ b = tapply(y,x,Mean,lmb=lmb); e = y-b[x]
    }else{ b = crossprod(y,x)/(crossprod(x)+lmb); b=as.vector(b); e = y-x*b }
    return(list(b=b,e=e))}
  RCT_update = function(y,x,b,lmb){
    Mean = function(x,lmb) sum(x)/(length(x)+lmb)
    if(is.factor(x)){ m = tapply(y,x,Mean,lmb=lmb); e = y-m[x]
    }else{ m = crossprod(y,x)/(crossprod(x)+lmb); m=as.vector(m); e = y-x*m }
    b=m+b; return(list(b=m,e=e))}
  
  # Structured random
  gws = function(e,i,...){
    x = data[[i]]
    if(iter==1){ e00 = e }else{ e00 = e + U[[i]][as.character(x)]}
    e0 = tapply(e00,x,mean,na.rm=TRUE)[rownames(X[[i]])]
    e0 = as.vector(ifelse(is.na(e0),0,e0))
    # Fit WGR
    h = alg(e0,X[[i]],...)
    fit = as.vector(h$hat)
    names(fit) = rownames(X[[i]])
    # Deregress
    if(Deregress){
     g = c(crossprod(fit,e0)/(crossprod(fit)))
     xb = fit*g
     b0 = e0-xb
     fit = b0+fit*g
    }    
    # Output
    res = e00 - fit[as.character(x)]
    OUT = list(b=fit,g=h$b,e=res)
    return(OUT)}
  
  # Fit (loop)
  e = y; R2c = 0.5
  for(iter in 1:maxit){
    cat('Iter ',iter,':',sep='')
    
    # Fixed effects 
    if(FIX){
      for(i in fxd){
        if(iter==1){
          go = FCT(e,data[[i]])
          e = go$e
          B[[i]] = go$b
        }else{
          go = FCT_update(e,data[[i]],B[[i]])
          e = go$e
          B[[i]] = go$b
        }}
    }else{
      # Intercept only 
      if(iter==1){
        mu = mean(e)
        e = e-mu
        B[['B0']] = mu
      }else{
        mu = mean(e)+B[['B0']]
        e = e-mu
        B[['B0']] = mu 
      }}
    
    # Random effects
    if(RND){
      # Coefficients
      for(i in rnd){
        # Structured
        if(i%in%ls(X)){
          go = gws(e,i,...)
          e = go$e
          U[[i]] = go$b
          G[[i]] = go$g
        }else{
          # Unstructured 
          if(iter==1){
            go = RCT(e,data[[i]],LMB[i])
            e = go$e
            U[[i]] = c(go$b)
          }else{
            go = RCT_update(e,data[[i]],U[[i]],LMB[i])
            e = go$e
            U[[i]] = c(go$b)
          }}}
      # VarComp & Lambda
      Ve = crossprod(y,e)[1,1]/(n-df)
      Va0 = sqrt(sapply(U,sd))
      Va = (Vy-Ve)*(Va0/sum(Va0))      
      LMB = Ve/Va}
    
    # Print R2 and check convergence based on Ve
    R2 = round(1-Ve/Vy,4)
    cat(' R2 =',R2,'\n')
    if(abs(R2c-R2)==0) break()
    R2c = R2}
  
  # Wrapping up
  fitted = y-e
  OUT = list(obs=y,hat=fitted,coef=list(fixed=B))
  if(RND){
    OUT$coef[['random']]=U
    if(length(X)>0) OUT$coef[['markers']]=G
    VC = c(Va,Residual=Ve)
    OUT[['var']]=VC}
  class(OUT) = 'mixed'
  return(OUT)}

mtmixed = function(resp, random=NULL, fixed=NULL, data, X=list(), maxit=5, init=15, regVC=FALSE){
  
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
    e0 = e0[rownames(X[[i]]),]
    # Add for debugging
    if(any(!rownames(X[[i]])%in%rownames(e0))){
      ww = which(!rownames(X[[i]])%in%rownames(e0))
      ww = rownames(X[[i]])[ww]
      mm = matrix(NA,length(ww),k,dimnames=list(ww,resp))
      e0 = rbind(e0,mm)}
    e0 = e0[rownames(X[[i]]),]
    h = mtgsru(Y=e0,X=X[[i]],b=MV[[i]]$b,vb=MV[[i]]$vb,ve=MV[[i]]$ve,iG=MV[[i]]$iG,maxit=init)
    unlist(lapply(h,anyNA)); head(h$b)
    fit = h$hat
    rownames(fit) = rownames(X[[i]])
    comp = list(b=h$b,vb=h$vb,ve=h$ve,iG=h$iG,h2=h$h2,MSx=h$MSx)
    # Output
    h = fit[as.character(x),]
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
      # Error = colSums(yc*e,na.rm=T)
      # SSa = sapply(H, function(h) colSums(yc*h,na.rm=T) )
      SS = cbind(SSa,Error)
      SS[SS<0] = 0.01
      wVar = t(apply(SS,1,function(x) x/sum(x) ))
      if( regVC ) wVar = sqrt(wVar)
      Vg = wVar*Vy
      Va = Vg[,rnd]/t(df0)*n
      Ve = Vg[,'Error']
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
