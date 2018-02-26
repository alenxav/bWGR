# Load data
data(met,package='NAM')
Obs$Block = factor(Obs$Block)

# Subset family 05
fam = substr(rownames(Gen),6,7)
Gen = Gen[which(fam%in%c('05')),]
Gen = data.matrix(Gen-1)
fam = substr(as.character(Obs$ID),6,7)
Obs = data.frame(Obs,SP=NAM::SPC(Obs$YLD,Obs$Block,Obs$Row,Obs$Col))
Obs = droplevels.data.frame(Obs[which(fam=='05'),])
rm(fam)

# Compare BLUP vs gBLUP
test1 = mixed(y=YLD,random=~ID+Block,fixed=~Year+Row+Col,data=Obs,X=list(ID=Gen))
test2 = mixed(y=YLD,random=~ID+Block,fixed=~Year+Row+Col,data=Obs)
cor(test1$coef$random$ID,test2$coef$random$ID)

# Plot
par(mfrow=c(1,3))
plot(test1$hat,Obs$YLD,main='Goodness of fit',xlab='Fitted',ylab='Observed',pch=20)
plot(test1$coef$random$ID,main='Genomic fit',test2$coef$random$ID,xlab='GEBV',ylab='BLUP',pch=20)
plot(test1$coef$markers$ID,main='Marker Effects',ylab='Effect',xlab='SNP',pch=20)