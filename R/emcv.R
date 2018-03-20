emCV=function(y, gen, k = 5, n = 5){
  # Cross-validation function
  folds = function(Seed,y,gen,k,it,bi){
    N = nrow(gen)
    # Begin folds
    set.seed(Seed)
    Nk = round(N/k)
    w=sample(1:N,Nk)
    Y=y; y[w]=NA
    f1=emRR(y[-w],gen[-w,])
    f2=emEN(y[-w],gen[-w,])
    f3=emBL(y[-w],gen[-w,])
    f4=emDE(y[-w],gen[-w,])
    f5=emBA(y[-w],gen[-w,])
    f6=emBB(y[-w],gen[-w,])
    f7=emBC(y[-w],gen[-w,])
    cat('DONE WITH CROSS-VALIDATION CYCLE',Seed,'\n')
    NamesMod = c('emRR','emEN','emBL','emDE',
                 'emBA','emBB','emBC','OBSERVATION')
    M = matrix(NA,Nk,length(NamesMod))
    colnames(M) = NamesMod
    for(i in 1:7) M[,i]=gen[w,]%*%get(paste('f',i,sep=''))$b
    M[,8] = Y[w]
    return(M)
  }
  # Running cross-validations
  Seeds=1:n
  b = lapply(Seeds,FUN=folds,y=y,gen=gen,k=k)
  names(b) = paste('CV_',1:length(b),sep='')
  # Summary function
  sCV=function(cv){
    n = length(cv)
    m = ncol(cv$CV_1)
    dta = matrix(0,0,m)
    for(i in 1:n) dta = rbind(dta,cv[[i]])
    # Summary
    PA = sort(cor(dta)[-m,m],decreasing = TRUE)
    return(round(PA,digits=4))
  }
  return(sCV(b))
}