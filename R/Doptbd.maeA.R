# SubSubsection 2.3.2 (Function for construction of D-optimal block design using array exchange algorithm) 
Doptbd.maeA<-function(trt.N,blk.N,theta,nrep,itr.cvrgval)  {
  #House keeping
  arrays=t(combn(trt.N,2))
  na=dim(arrays)[1]
  del.1<-matrix(10^20,na,3)
  desbest.1<-matrix(0,nrep*2,blk.N)
  doptbest.1<-matrix(0,nrep,2)
  for(irep in 1:nrep){
    des<-intcbd.mae(trt.N,blk.N)
    if(trt.N==blk.N&trt.N>3&irep<(trt.N-1)) {in.desns=matrix(0,(trt.N-3)*2,blk.N)
    in.desns0=rbind(seq(1,trt.N),c(seq(1,trt.N)[2:trt.N],1))
    for(i in 1:(trt.N-3)) {in.desns01=cbind(rbind(seq(1,(trt.N-i)),c(seq(1,(trt.N-i))[2:(trt.N-i)],1)),rbind(rep(1,i),((trt.N-i+1):trt.N))); in.desns[c((i-1)*2+1,i*2),]=in.desns01}
    in.desns=rbind(rbind(seq(1,trt.N),c(seq(1,trt.N)[2:trt.N],1)),in.desns)
    des=in.desns[c((irep-1)*2+1,irep*2),]}
    cmat<-cmatbd.mae(trt.N,blk.N,theta,des)
    degv<-sort(eigen(cmat)$values)
    degvp<-degv[2:length(degv)]
    dopt<-prod(1/degvp) 
    dcold=dopt
    descold=t(des)
    cdel=100
    i=1;
    ivaldcold={}
    for (i in 1:blk.N){
      j=1;
      for (j in 1:na){
        temp=descold[i,]
        if(all(descold[i,]==arrays[j,]))  {dopt=dcold; del.1[j,]<-c(j,(dcold-dopt),dopt); next}
        descold[i,]=arrays[j,]
        trtin<-contrasts(as.factor(t(descold)),contrasts=FALSE)[as.factor(t(descold)),]
        R.trt<-t(trtin)%*%trtin
        if (rankMatrix(R.trt)[1]<trt.N)  {dopt=10^20; del.1[j,]<-c(j,(dcold-dopt),dopt); next}
        cmato=cmatbd.mae(trt.N,blk.N, 0,t(descold))
        egv<-sort(eigen(cmato)$values)
        if(egv[2]<0.000001) {dopt=10^20; del.1[j,]<-c(j,(dcold-dopt),dopt); next}
        cmat=cmatbd.mae(trt.N,blk.N,theta,t(descold))
        degv<-sort(eigen(cmat)$values)
        degvp<-degv[2:length(degv)]
        dopt<-prod(1/degvp) 
        del.n<-del.1[j,]<-c(j,(dcold-dopt),dopt)
        descold[i,]=temp
      }
      del.1<-del.1[order(del.1[,3]),]
      delbest=t(del.1[1,])
      descold[i,]=arrays[delbest[1],]
      dcold=delbest[3]
      cdel=delbest[2]
      ivaldcold=rbind(ivaldcold, c(i,dcold))
      if(i>itr.cvrgval) if(all(ivaldcold[c(i-(itr.cvrgval-2),i),2]==ivaldcold[i-(itr.cvrgval-1),2])) break
    }
    if (irep==1) {desbest.1=t(descold)} else {desbest.1=rbind(desbest.1,t(descold))}
    doptbest.1[irep,]=c(irep,dcold)
  }
  best=doptbest.1[order(doptbest.1[,2]),]
  nb=best[1,1]
  Dscore<-best[1,2]
  Doptde<- desbest.1[c((nb-1)*2+1,nb*2),]
  tkmessageBox(title="Search completed",message=paste("Search completed",sep=""))
  cnames=paste0("Ary",1:blk.N)
  dimnames(Doptde)=list(NULL,cnames)
  Dopt_sum2<-list("v"=trt.N,"b"=blk.N,theta=theta,nrep=nrep,itr.cvrgval=itr.cvrgval, "OptdesF"=Doptde,"Optcrtsv" =Dscore)
  return(Dopt_sum2)
}#End of SubSubsection 2.3.2 (Doptbd.maeA function) construction of D-optimal block design using array exchange algorithm
#End of Subsection 2.3 (Function for construction of D-optimal block design)    
