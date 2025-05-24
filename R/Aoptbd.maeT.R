#Subsection 2.1: Function for search of A-optimal or near-optimal block designs
# SubSubsection 2.1.1 (Function for construction of A-optimal block designs using treatment exchange algorithm) 
Aoptbd.maeT<-function(trt.N,blk.N,theta,nrep,itr.cvrgval) {
  #House keeping
  del.1<-matrix(1000,trt.N,3)
  desbest.1<-matrix(0,nrep*2,blk.N)
  aoptbest.1<-matrix(0,nrep,2)
  #Start iteration
  for(irep in 1:nrep){
    #Initial design with its corresponding Ascore value
    des<-intcbd.mae(trt.N,blk.N)
    if(trt.N==blk.N&trt.N>3&irep<(trt.N-1)) {in.desns=matrix(0,(trt.N-3)*2,blk.N)
    in.desns0=rbind(seq(1,trt.N),c(seq(1,trt.N)[2:trt.N],1))
    for(i in 1:(trt.N-3)) {in.desns01=cbind(rbind(seq(1,(trt.N-i)),c(seq(1,(trt.N-i))[2:(trt.N-i)],1)), rbind(rep(1,i),((trt.N-i+1):trt.N))); in.desns[c((i-1)*2+1,i*2),]=in.desns01}
    in.desns=rbind(rbind(seq(1,trt.N),c(seq(1,trt.N)[2:trt.N],1)),in.desns)
    des=in.desns[c((irep-1)*2+1,irep*2),]}
    cmat<-cmatbd.mae(trt.N,blk.N,theta,des)
    aopt=sum(diag(ginv(cmat)))
    acold=aopt
    descold=t(des)
    #deletion difference
    cdel=100
    i=1;
    ivalacold={}
    for(i in 1:blk.N){
      for (m in 1:2){
        j=1;
        for(j in 1:trt.N){
          temp=descold[i,]
          if(m==1) {
            if(j==descold[i,1]|j==descold[i,2]) {aopt=acold; del.1[j,]<-c(descold[i,1],(acold-aopt),aopt); next} else { descold[i,]=c(j,descold[i,2])}}
          if(m==2) {
            if(descold[i,2]==j|j==descold[i,1]) {aopt=acold; del.1[j,]<-c(descold[i,2],(acold-aopt),aopt); next} else { descold[i,]=c(descold[i,1],j)}}
          
          trtin<-contrasts(as.factor(t(descold)),contrasts=FALSE)[as.factor(t(descold)),]
          R.trt<-t(trtin)%*%trtin
          if (rankMatrix(R.trt)[1]<trt.N)  {aopt=acold; descold[i,]=temp; if(m==1) {del.1[j,]<-c(descold[i,1],(acold-aopt),aopt)} else {
           del.1[j,]<-c(descold[i,2],(acold-aopt),aopt)}; next}
          cmato=cmatbd.mae(trt.N,blk.N, 0,t(descold))
          egv<-sort(eigen(cmato)$values)
          if(egv[2]<0.000001) {aopt=acold; descold[i,]=temp; if(m==1){del.1[j,]<-c(descold[i,1],(acold-aopt),aopt)} else {
           del.1[j,]<-c(descold[i,2],(acold-aopt),aopt)}; next}
          cmat=cmatbd.mae(trt.N,blk.N,theta,t(descold))
          aopt=sum(diag(ginv(cmat)))
          del.n<-del.1[j,]<-c(j,(acold-aopt),aopt)
          descold[i,]=temp
        }
        del.1<-del.1[order(del.1[,3]),]
        delbest=t(del.1[1,])
        if (m==1) {
          if (delbest[1]==descold[i,2]) {descold[i,]=descold[i,]}  else 
          {descold[i,]=c(delbest[1],descold[i,2]); cdel=delbest[2]; acold=delbest[3]}} else {
            if (descold[i,1]==delbest[1]) {descold[i,]= descold[i,]} else 
            {descold[i,]=c(descold[i,1],delbest[1]); cdel=delbest[2]; acold=delbest[3]}}
      }
      ivalacold=rbind(ivalacold, c(i,acold))
      if(i>itr.cvrgval) if(all(ivalacold[c(i-(itr.cvrgval-2),i),2]==ivalacold[i-(itr.cvrgval-1),2])) break
    }
    next.it<- if (irep==1) {desbest.1=t(descold)} else {desbest.1=rbind(desbest.1,t(descold))}
    aoptbest.1[irep,]=c(irep,acold)
  }
  best=aoptbest.1[order(aoptbest.1[,2]),]
  nb=best[1,1]
  Ascore<-best[1,2]
  Aoptde<- desbest.1[c((nb-1)*2+1,nb*2),]
  if(trt.N!=3) {tkmessageBox(title="Search completed",message=paste("Search completed",sep=""))}
  cnames=paste0("Ary",1:blk.N)
  dimnames(Aoptde)=list(NULL,cnames)
  Aopt_sum2<-list("v"=trt.N,"b"=blk.N,theta=theta,nrep=nrep,itr.cvrgval=itr.cvrgval, "OptdesF"=Aoptde,"Optcrtsv" =Ascore)
  return(Aopt_sum2)
}#End of SubSubsection 2.1.1 (Aoptbd.maeT function) construction of A-optimal block design using treatment exchange algorithm
