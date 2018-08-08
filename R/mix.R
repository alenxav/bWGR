mixed = function(y,random=NULL,fixed=NULL,data=NULL,X=list(),alg=emML,...){
  
  # Get y vector
  if(!is.null(data)) y = data[[deparse(substitute(y))]]
  
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
    g = c(crossprod(fit,e0)/(crossprod(fit)))
    xb = fit*g
    b0 = e0-xb
    fit = b0+fit*g
    # Output
    res = e00 - fit[as.character(x)]
    OUT = list(b=fit,g=h$b,e=res)
    # cat(' (WGR on ',i,') ',sep='')
    return(OUT)}
  
  # Fit (loop)
  e = y; R2c = 0.5
  for(iter in 1:10){
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
      Va = (sapply(U,crossprod)+Ve)/(df0-1)
      LMB = sqrt(LMB0*Ve/Va)}
    
    # Print R2 and check convergence based on Ve
    R2 = round(1-Ve/var(y),3)
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