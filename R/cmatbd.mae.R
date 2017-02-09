#Section 3 computation of the information matrix (C-matrix) for a given block design (des) 
cmatbd.mae<-function(trt.N,blk.N,theta,des){
  k=2
  trtin<-contrasts(as.factor(des),contrasts=FALSE)[as.factor(des),]
  blk.1<-rep(1:blk.N,each=2)
  blkin<-contrasts(as.factor(blk.1),contrasts=FALSE)[as.factor(blk.1),]
  vec.1<-rep(1,blk.N*2)
  R.trt<-t(trtin)%*%trtin
  N.tb<-t(trtin)%*%blkin
  r.trt<-t(trtin)%*%vec.1
  cmat<-R.trt-(1/k)*(N.tb%*%t(N.tb))+theta*((1/k)*(N.tb%*%t(N.tb))-(1/(blk.N*k))*(r.trt%*%t(r.trt)))
  cmat
}#End of Section 3 (computation of C-matrix, cmatbd.mae)
