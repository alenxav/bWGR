Hmat = function (ped,gen=NULL,Diss=FALSE,inbred=FALSE){
  
  # ped - Pedigree. Numeric matrix n by 3.
  # gen - Numeric matrix where rownames regards the ped ID
  # Diss - Logical. If true, it ignores inbreeding.
  # inbred - Logical. Use genomic for FIS or set main diagonal as 2.
  
  n = nrow(ped)
  A = diag(0, n)
  
  if(is.null(gen)){
    
    # WITHOUT GENOTYPES
    for (i in 1:n) {
      for (j in 1:n) {
        if (i > j) {
          A[i, j] = A[j, i]
        } else {
          d = ped[j, 2]
          s = ped[j, 3]
          if (d == 0) Aid = 0 else Aid = A[i, d]
          if (s == 0) Ais = 0 else Ais = A[i, s]
          if (d == 0 | s == 0) Asd = 0 else Asd = A[d, s]
          if (i == j) Aij = 1 + 0.5 * Asd else Aij = 0.5 * (Aid + Ais)
          A[i, j] = Aij
        }}}
    
  }else{
    
    # WITH GENOTYPES 
    G = as.numeric(rownames(gen))
    if(any(is.na(gen))){
      ND = function(g1,g2){
        X = abs(g1-g2)
        L = 2*sum(!is.na(X))
        X = sum(X,na.rm=TRUE)/L
        return(X)} 
    }else{
      ND = function(g1,g2) sum(abs(g1-g2))/(2*length(g1))
    }
    Inbr = function(g) 2-mean(g==1)
    
    # LOOP     
    for (i in 1:n) {
      for (j in 1:n) {
        
        #######################
        if (i > j) {
          A[i, j] = A[j, i]
        } else {
          
          d = ped[j, 2]
          s = ped[j, 3]
          
          if(j%in%G){
            
            if(d!=0&d%in%G){ Aid=2*ND(gen[paste(j),],gen[paste(d),]) }else{if(d==0) Aid=0 else Aid=A[i,d]}
            if(s!=0&s%in%G){ Ais=2*ND(gen[paste(j),],gen[paste(s),]) }else{if(s==0) Ais=0 else Ais=A[i,s]}
            Asd = Inbr(gen[paste(j),])
            
          }else{
            if (d == 0) Aid = 0 else Aid = A[i, d]
            if (s == 0) Ais = 0 else Ais = A[i, s]
            if (d == 0 | s == 0) Asd = 0 else Asd = A[d, s]
          }
          
          if (i == j) Aij = ifelse(Diss,1,ifelse(inbred,ifelse(i%in%G,Asd,2),1+0.5*Asd))
          else Aij = 0.5*(Aid+Ais)
          A[i, j] = round(Aij,6)
        }
        #######################
        
      }}
  }
  return(A)
}

SibZ = function(id,p1,p2){
  lvl = unique(c(unique(as.character(p1)),
                 unique(as.character(p2)),
                 unique(as.character(id))))
  id = factor(as.character(id),levels=lvl)
  p1 = factor(as.character(p1),levels=lvl)
  p2 = factor(as.character(p2),levels=lvl)
  x = p1; if(any(is.na(x))) x[is.na(x)]=0; Z = model.matrix(~x-1)*0.5
  x = p2; if(any(is.na(x))) x[is.na(x)]=0; Z = Z + model.matrix(~x-1)*0.5
  Z = Z[,colnames(Z)!='x0']
  x = sqrt(1-rowSums(Z*Z))
  xx = paste('x',id,sep='')
  for(i in 1:length(x)) Z[i,xx[i]]=x[i]
  return(Z)}
      
