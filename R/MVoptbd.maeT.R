#Subsection 2.2: Function for search of MV-optimal or near-optimal block designs
#SubSubsection 2.2.1 (Function for construction of MV-optimal block designs using treatment exchange algorithm) 
MVoptbd.maeT<-function(trt.N,blk.N,theta,nrep,itr.cvrgval) {
  ii=2
  trco=cbind(matrix(1,trt.N-1),-diag(1,trt.N-1,trt.N-1))
  while(ii<=trt.N-1){
    if (ii==trt.N-1){
      trco1=cbind(matrix(0,1,trt.N-2),matrix(1,trt.N-ii),-diag(1,trt.N-ii,trt.N-ii))}
    else
    {trco1=cbind(matrix(0,trt.N-ii,trt.N-(trt.N-ii+1)),matrix(1,trt.N-ii),-diag(1,trt.N-ii,trt.N-ii))}
    trco=rbind(trco,trco1)
    ii=ii+1
  }
  del.1<-matrix(10^20,trt.N,3)
  desbest.1<-matrix(0,nrep*2,blk.N)
  MVoptbest.1<-matrix(0,nrep,2)
  for(irep in 1:nrep){
    des<-intcbd.mae(trt.N, blk.N)
    if(trt.N==blk.N&trt.N>3&irep<(trt.N-1)) {in.desns=matrix(0,(trt.N-3)*2,blk.N)
    in.desns0=rbind(seq(1,trt.N),c(seq(1,trt.N)[2:trt.N],1))
    for(i in 1:(trt.N-3)) {in.desns01=cbind(rbind(seq(1,(trt.N-i)),c(seq(1,(trt.N-i))[2:(trt.N-i)],1)), rbind(rep(1,i),((trt.N-i+1):trt.N))); in.desns[c((i-1)*2+1,i*2),]=in.desns01}
    in.desns=rbind(rbind(seq(1,trt.N),c(seq(1,trt.N)[2:trt.N],1)),in.desns)
    des=in.desns[c((irep-1)*2+1,irep*2),]}
    cmat<-cmatbd.mae(trt.N,blk.N,theta,des)
    invc=ginv(cmat)
    invcp=trco%*%invc%*%t(trco);
    MVopt =max(diag(invcp));
    MVcold=MVopt
    descold=t(des)
    cdel=100
    ivalMVcold={}
    for (i in 1:blk.N){
      m=1;
      for (m in 1:2){
        j=1;
        for (j in 1:trt.N){
          temp=descold[i,]
          
        if(m==1) {
            if(j==descold[i,1]|j==descold[i,2]) {MVopt=MVcold; del.1[j,]<-c(descold[i,1],(MVcold-MVopt),MVopt); next} else { descold[i,]=c(j,descold[i,2])}}
          if(m==2) {
            if(descold[i,2]==j|j==descold[i,1]) {MVopt=MVcold; del.1[j,]<-c(descold[i,2],(MVcold-MVopt),MVopt); next} else { descold[i,]=c(descold[i,1],j)}}
          trtin<-contrasts(as.factor(t(descold)),contrasts=FALSE)[as.factor(t(descold)),]
          R.trt<-t(trtin)%*%trtin
          if (rankMatrix(R.trt)[1]<trt.N)  {MVopt=MVcold; descold[i,]=temp; if(m==1) {del.1[j,]<-c(descold[i,1],(MVcold-MVopt),MVopt)} else {
           del.1[j,]<-c(descold[i,2],(MVcold-MVopt),MVopt)}; next}
          cmato=cmatbd.mae(trt.N,blk.N, 0,t(descold))
          egv<-sort(eigen(cmato)$values)
          if(egv[2]<0.000001) {MVopt=MVcold; descold[i,]=temp; if(m==1) {del.1[j,]<-c(descold[i,1],(MVcold-MVopt),MVopt)} else {
           del.1[j,]<-c(descold[i,2],(MVcold-MVopt),MVopt)}; next}
          cmat=cmatbd.mae(trt.N,blk.N,theta,t(descold))
          invc=ginv(cmat)
          invcp=trco%*%invc%*%t(trco);
          MVopt =max(diag(invcp));
          del.n<-del.1[j,]<-c(j,(MVcold-MVopt),MVopt)
          descold[i,]=temp
        }
        del.1<-del.1[order(del.1[,3]),]
        delbest=t(del.1[1,])
        if (m==1) {
          if (delbest[1]==descold[i,2]) { descold[i,]= descold[i,]}  else { descold[i,]=c(delbest[1],descold[i,2]); cdel=delbest[2];  MVcold=delbest[3] }} else {
            if (descold[i,1]==delbest[1]) {descold[i,]= descold[i,]} else { descold[i,]=c(descold[i,1],delbest[1]); cdel=delbest[2];  MVcold=delbest[3] }}
      }
      ivalMVcold=rbind(ivalMVcold, c(i,MVcold))
      if(i>itr.cvrgval) if(all(ivalMVcold[c(i-(itr.cvrgval-2),i),2]==ivalMVcold[i-(itr.cvrgval-1),2])) break
    }
    if (irep==1) {desbest.1=t(descold)} else {desbest.1=rbind(desbest.1,t(descold))}
    MVoptbest.1[irep,]=c(irep,MVcold)
  }
  best=MVoptbest.1[order(MVoptbest.1[,2]),]
  nb=best[1,1]
  MVscore<-best[1,2]
  MVoptde<- desbest.1[c((nb-1)*2+1,nb*2),]
  tkmessageBox(title="Search completed",message=paste("Search completed",sep=""))
  cnames=paste0("Ary",1:blk.N)
  dimnames(MVoptde)=list(NULL,cnames)
  MVopt_sum2<-list("v"=trt.N,"b"=blk.N,theta=theta,nrep=nrep,itr.cvrgval=itr.cvrgval, "OptdesF"=MVoptde,"Optcrtsv" =MVscore)
  return(MVopt_sum2)
}#End of SubSubsection 2.2.1 (MVoptbd.maeT function) construction of MV-optimal block design using treatment exchange algorithm

