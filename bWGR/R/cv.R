
emCV = function (y, gen, k=5, n=5,
                 Pi=0.75, alpha=0.02, df=10, R2=0.5, avg=TRUE, llo=NULL, tbv=NULL,
                 ReturnGebv = FALSE){
  B0 = list()
  folds = function(Seed, y, gen, k) {
    N = nrow(gen)
    set.seed(Seed)
    Nk = round(N/k)
    w = sample(1:N, Nk)
    Y = y
    y[w] = NA
    f1 = emRR(y[-w], gen[-w, ], R2 = R2, df = df)
    f2 = emEN(y[-w], gen[-w, ], alpha = alpha, R2 = R2)
    f3 = emBL(y[-w], gen[-w, ], alpha = alpha, R2 = R2)
    f4 = emDE(y[-w], gen[-w, ], R2 = R2)
    f5 = emBA(y[-w], gen[-w, ], R2 = R2, df = df)
    f6 = emBB(y[-w], gen[-w, ], Pi = Pi, R2 = R2, df = df)
    f7 = emBC(y[-w], gen[-w, ], Pi = Pi, R2 = R2, df = df)
    f8 = emML(y[-w], gen[-w, ])
    f9 = emBCpi(y[-w], gen[-w, ])
    f10 = lasso(y[-w], gen[-w, ])
    cat("DONE WITH CROSS-VALIDATION CYCLE", Seed, "\n")
    NamesMod = c("emRR", "emEN", "emBL", "emDE", "emBA", 
                 "emBB", "emBC", "emML", "emBCpi","lasso", 
                 "OBSERVATION")
    M = matrix(NA, Nk, length(NamesMod))
    B = matrix(NA, ncol(gen), length(NamesMod)-1)
    colnames(M) = NamesMod
    for (i in 1:(length(NamesMod)-1)){
      M[, i] = gen[w, ] %*% get(paste("f", i,sep = ""))$b
      B[, i] = get(paste("f", i,sep = ""))$b
    } 
    colnames(B) = NamesMod[1:(length(NamesMod)-1)]
    if(is.null(tbv)){
      M[,length(NamesMod)] = Y[w]
    }else{
      M[,length(NamesMod)] = tbv[w]
    }
    B0[[length(B0)+1]] <<- B
    return(M)
  }
  llo_folds = function(lev, y, gen) {
    w = which(llo==lev)
    Nk = length(w)
    Y = y
    y[w] = NA
    f1 = emRR(y[-w], gen[-w, ], R2 = R2, df = df)
    f2 = emEN(y[-w], gen[-w, ], alpha = alpha, R2 = R2)
    f3 = emBL(y[-w], gen[-w, ], alpha = alpha, R2 = R2)
    f4 = emDE(y[-w], gen[-w, ], R2 = R2)
    f5 = emBA(y[-w], gen[-w, ], R2 = R2, df = df)
    f6 = emBB(y[-w], gen[-w, ], Pi = Pi, R2 = R2, df = df)
    f7 = emBC(y[-w], gen[-w, ], Pi = Pi, R2 = R2, df = df)
    f8 = emML(y[-w], gen[-w, ])
    f9 = emBCpi(y[-w], gen[-w, ])
    f10 = lasso(y[-w], gen[-w, ])
    cat("DONE WITH CROSS-VALIDATION CYCLE", lev, "\n")
    NamesMod = c("emRR", "emEN", "emBL", "emDE", "emBA", 
                 "emBB", "emBC", "emML","emBCpi","lasso", 
		 "OBSERVATION")
    M = matrix(NA, Nk, length(NamesMod))
    B = matrix(NA, ncol(gen), length(NamesMod)-1)
    colnames(M) = NamesMod
    for (i in 1:(length(NamesMod)-1)){
      M[, i] = gen[w, ] %*% get(paste("f", i,sep = ""))$b
      B[, i] = get(paste("f", i,sep = ""))$b
    } 
    colnames(B) = NamesMod[1:(length(NamesMod)-1)]
    if(is.null(tbv)){
      M[,length(NamesMod)] = Y[w]
    }else{
      M[,length(NamesMod)] = tbv[w]
    }
    B0[[length(B0)+1]] <<- B
    return(M)
  }
  if(is.null(llo)){
    Seeds = 1:n
    b = lapply(Seeds, FUN = folds, y = y, gen = gen, k = k)
  }else{
    lev = unique(as.character(llo))
    b = lapply(lev, FUN = llo_folds, y = y, gen = gen)
  }
  names(b) = paste("CV_", 1:length(b), sep = "")
  sCV = function(cv) {
    n = length(cv)
    m = ncol(cv$CV_1)
    if(avg){
      dta = matrix(0, 0, m)
      for (i in 1:n) dta = rbind(dta, cv[[i]])
      PA = sort(cor(dta,use='p')[-m, m], decreasing = TRUE)
      return(round(PA, digits = 4))
    }else{
      dta = c()
      for (i in 1:n) dta = rbind(dta, cor(cv[[i]])[-m, m])
      rownames(dta) = paste('CV',1:n,sep='_')
      return(round(dta, digits = 4))
    }
  }
  if(ReturnGebv){
    beta = B0[[1]]
    for(i in 2:length(B0)) beta=beta+B0[[i]]
    beta = beta/length(B0)
    hat = gen%*%beta + mean(y,na.rm = T)
    OUT = list(cv=sCV(b),hat=hat,beta=beta)
  }else{
    OUT = sCV(b)
  }
  return(OUT)
}

mcmcCV = function (y, gen, k = 5, n = 5,
                   it=1500, bi=500, pi=0.95, df=5, R2=0.5, avg=TRUE, llo=NULL, tbv=NULL,
                   ReturnGebv = FALSE){
  B0 = list()
  folds = function(Seed, y, gen, k) {
    N = nrow(gen)
    set.seed(Seed)
    Nk = round(N/k)
    w = sample(1:N, Nk)
    Y = y
    y[w] = NA
    f1 = BayesA(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi)
    f2 = BayesB(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi, pi=pi)
    f3 = BayesC(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi, pi=pi)
    f4 = BayesL(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi)
    f5 = BayesRR(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi)
    f6 = BayesCpi(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi)
    f7 = BayesDpi(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi)
    cat("DONE WITH CROSS-VALIDATION CYCLE", Seed, "\n")
    NamesMod = c("BayesA", "BayesB", "BayesC", "BayesL",
                 "BayesCpi", "BayesDpi", "BayesRR", "OBSERVATION")
    M = matrix(NA, Nk, length(NamesMod))
    B = matrix(NA, ncol(gen), length(NamesMod)-1)
    colnames(M) = NamesMod
    for (i in 1:(length(NamesMod)-1)){
      M[, i] = gen[w, ] %*% get(paste("f", i,sep = ""))$b
      B[, i] = get(paste("f", i,sep = ""))$b
    } 
    colnames(B) = NamesMod[1:(length(NamesMod)-1)]
    if(is.null(tbv)){
      M[,length(NamesMod)] = Y[w]
    }else{
      M[,length(NamesMod)] = tbv[w]
    }
    B0[[length(B0)+1]] <<- B
    return(M)
  }
  llo_folds = function(lev, y, gen) {
    w = which(llo==lev)
    Nk = length(w)
    Y = y
    y[w] = NA
    f1 = BayesA(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi)
    f2 = BayesB(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi, pi=pi)
    f3 = BayesC(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi, pi=pi)
    f4 = BayesL(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi)
    f5 = BayesRR(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi)
    f6 = BayesCpi(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi)
    f7 = BayesDpi(y[-w], gen[-w, ], R2 = R2, df = df, it=it, bi=bi)
    cat("DONE WITH CROSS-VALIDATION CYCLE", lev, "\n")
    NamesMod = c("BayesA", "BayesB", "BayesC", "BayesL",
                 "BayesCpi", "BayesDpi", "BayesRR",
		 "OBSERVATION")
    M = matrix(NA, Nk, length(NamesMod))
    B = matrix(NA, ncol(gen), length(NamesMod)-1)
    colnames(M) = NamesMod
    for (i in 1:(length(NamesMod)-1)){
      M[, i] = gen[w, ] %*% get(paste("f", i,sep = ""))$b
      B[, i] = get(paste("f", i,sep = ""))$b
    } 
    colnames(B) = NamesMod[1:(length(NamesMod)-1)]
    if(is.null(tbv)){
      M[,length(NamesMod)] = Y[w]
    }else{
      M[,length(NamesMod)] = tbv[w]
    }
    B0[[length(B0)+1]] <<- B
    return(M)
  }
  if(is.null(llo)){
    Seeds = 1:n
    b = lapply(Seeds, FUN = folds, y = y, gen = gen, k = k)
  }else{
    lev = unique(as.character(llo))
    b = lapply(lev, FUN = llo_folds, y = y, gen = gen)
  }
  
  names(b) = paste("CV_", 1:length(b), sep = "")
  sCV = function(cv) {
    n = length(cv)
    m = ncol(cv$CV_1)
    if(avg){
      dta = matrix(0, 0, m)
      for (i in 1:n) dta = rbind(dta, cv[[i]])
      PA = sort(cor(dta,use='p')[-m, m], decreasing = TRUE)
      return(round(PA, digits = 4))
    }else{
      dta = c()
      for (i in 1:n) dta = rbind(dta, cor(cv[[i]])[-m, m])
      rownames(dta) = paste('CV',1:n,sep='_')
      return(round(dta, digits = 4))
    }
  }
  if(ReturnGebv){
    beta = B0[[1]]
    for(i in 2:length(B0)) beta=beta+B0[[i]]
    beta = beta/length(B0)
    hat = gen%*%beta + mean(y,na.rm = T)
    OUT = list(cv=sCV(b),hat=hat,beta=beta)
  }else{
    OUT = sCV(b)
  }
  return(OUT)
}

AccByC = function(X1,X2,h2=0.5){
  alpha = 1/sqrt(mean(apply(X1,1,crossprod)))
  V = svd(crossprod(X1),T)$v *alpha
  Z1 = X1 %*% V 
  Z2 = X2 %*% V 
  D = apply(Z1,2,crossprod)
  ve = ((1-h2)/h2)
  VarBhat = 1-1/(D/ve+1)
  sqrtVarBhat = sqrt(VarBhat)
  Acc = apply(Z2,1,function(z) sqrt( crossprod(z*sqrtVarBhat)/(crossprod(z)) )  )
  return(Acc)}