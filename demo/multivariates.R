require(bWGR)
Z = SimZ(ind = 1000)
GC = SimGC(MIX = F)
S = SimY(Z = Z, GC = GC)
Y = S$Y
tbv = S$tbv
fit_mega = MEGAF(Y,Z)$gebv
fit_gsem = GSEMF(Y,Z)$hat

fit_uv = Z %*% FUVBETA(Y,Z)
fit_mrr = MRR3F(Y,Z)$hat
fit_xfa = MRR3(Y,Z,XFA=T)$hat
fit_sem = SEM(Y,Z)$hat
fit_xsem = XSEMF(Y,Z)$hat
fit_zsem = ZSEMF(Y,Z)$hat
fit_gsem = GSEMF(Y,Z)$hat
fit_mega = MEGAF(Y,Z)$gebv

models = grep('fit',ls(),value=T) 
sort(sapply(models, function(x) mean(diag(cor(get(x),tbv)))))
