# Simulate data
Z = SimZ(ind = 1250)
S = SimY(Z,GC = 0.75,k = 10, h2=0.25,PercMiss=0.5)
Y = S$Y
colnames(Z) = paste0("SNP", 1:ncol(Z))
rownames(Z) = rownames(Y) = paste0("IND", 1:nrow(Z))
colnames(Y) = paste0("ENV", 1:ncol(Y))
dta = data.frame(GEN = paste0("IND", 1:nrow(Z)),Y,row.names = NULL)

# Reshape dta from wide to long data frame using native R function reshape
dta = reshape(dta,varying = list(colnames(Y)),direction = 'long',v.names = 'YIELD',timevar = 'ENV',times = colnames(Y),idvar = 'GEN')
dta$ENV = factor(dta$ENV)
dta$GEN = factor(dta$GEN)

# Compare BLUP vs gBLUP
test1 = mm(y=YIELD,random=~GEN+GEN*ENV,fixed=~ENV,data=dta,M=list(GEN=Z))
test2 = mm(y=YIELD,random=~GEN,fixed=~ENV,data=dta)
cor(test1$Coef$GEN,test2$Coef$GEN)

# Plot
par(mfrow=c(1,3))
plot(test1$GOF$Predicted,test1$GOF$Observed,main='Goodness of fit',xlab='Fitted',ylab='Observed',pch=20)
abline(lm(test1$GOF$Observed~test1$GOF$Predicted),col=2)
legend('bottom',paste('Cor =',round(cor(test1$GOF$Observed,test1$GOF$Predicted,use='p'),2)),bty='n')

plot(test1$Coef$GEN,test2$Coef$GEN,xlab='GEBV',ylab='BLUP',pch=20,main='Genomic fit')
abline(lm(test2$Coef$GEN~test1$Coef$GEN),col=2)
legend('bottom',paste('Cor =',round(cor(test1$Coef$GEN,test2$Coef$GEN),2)),bty='n')

plot(test1$Mrk$GEN,main='Marker Effects',ylab='Effect',xlab='SNP',pch=20,type='h')
abline(h=0,col=2,lwd=3)
