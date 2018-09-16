emCV = function (y, gen, k=5, n=5, Pi=0.75, alpha=0.02, df=10, R2=0.5) {
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
    cat("DONE WITH CROSS-VALIDATION CYCLE", Seed, "\n")
    NamesMod = c("emRR", "emEN", "emBL", "emDE", "emBA", 
                 "emBB", "emBC", "emML", "OBSERVATION")
    M = matrix(NA, Nk, length(NamesMod))
    colnames(M) = NamesMod
    for (i in 1:(length(NamesMod)-1)) M[, i] = gen[w, ] %*% get(paste("f", i,sep = ""))$b
    M[,length(NamesMod)] = Y[w]
    return(M)
  }
  Seeds = 1:n
  b = lapply(Seeds, FUN = folds, y = y, gen = gen, k = k)
  names(b) = paste("CV_", 1:length(b), sep = "")
  sCV = function(cv) {
    n = length(cv)
    m = ncol(cv$CV_1)
    dta = matrix(0, 0, m)
    for (i in 1:n) dta = rbind(dta, cv[[i]])
    PA = sort(cor(dta)[-m, m], decreasing = TRUE)
    return(round(PA, digits = 4))
  }
  return(sCV(b))
}

mcmcCV = function (y, gen, k = 5, n = 5, it=1500, bi=500, pi=0.95, df=5, R2=0.5) {
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
    colnames(M) = NamesMod
    for (i in 1:(length(NamesMod)-1)) M[, i] = gen[w, ] %*% get(paste("f", i,sep = ""))$b
    M[, length(NamesMod)] = Y[w]
    return(M)
  }
  Seeds = 1:n
  b = lapply(Seeds, FUN = folds, y = y, gen = gen, k = k)
  names(b) = paste("CV_", 1:length(b), sep = "")
  sCV = function(cv) {
    n = length(cv)
    m = ncol(cv$CV_1)
    dta = matrix(0, 0, m)
    for (i in 1:n) dta = rbind(dta, cv[[i]])
    PA = sort(cor(dta)[-m, m], decreasing = TRUE)
    return(round(PA, digits = 4))
  }
  return(sCV(b))
}
