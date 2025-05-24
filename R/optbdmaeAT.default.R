optbdmaeAT.default<-function(trt.N,blk.N,theta,nrep,itr.cvrgval,Optcrit="",Alg="",...)
{
trt.N=as.numeric(trt.N)
blk.N=as.numeric(blk.N)
theta=as.numeric(theta)
nrep=as.numeric(nrep)
itr.cvrgval=as.numeric(itr.cvrgval)
#"===================================================================================================="
  if(is.na(theta)|theta<0|theta>1){
    tkmessageBox(title="Error",message=paste("Please insert correct value of theta, it should be between 0 and 1 inclusive of 0 and 1, click OK to reset.",sep=""),icon = "error"); 
    stop("Wrong value of 'theta', it should be between 0 and 1, inclusive of 0 and 1")
  }#end of if
  if(is.na(trt.N)|is.na(blk.N)|trt.N!=round(trt.N)|blk.N!=round(blk.N)) {
    tkmessageBox(title="Error",message=paste("Please insert correct format of the number of treatments and arrays. The number of treatments and arrays should be an integer, click OK to reset the values.",sep=""),icon = "error"); 
    stop("Wrong format of 'trt.N' and/or 'blk.N', both should be an integer")
  }#end of if
  if(trt.N<3|blk.N<3){ 
    tkmessageBox(title="Error",message=paste("The number of blocks and treatments should be greater than or equal to 3, click Ok to reset.",sep=""),icon = "error"); 
    stop("Very small value of number of treatments and/or blocks, minimum value of the two is 3")
  }#end of if
  if(trt.N-blk.N>1){ 
    tkmessageBox(title="Error",message=paste("The number of arrays should be greater than or equal to the number of treatments minus one, click Ok to reset.",sep=""),icon = "error"); 
    stop("The number of treatments are larger than the number of arrays minus one (trt.N>blk.N-1)")
  }#end of if
  if(trt.N>10|blk.N>10){ 
    tkmessageBox(title="Information",message=paste("This might take some minutes, please be patient...",sep=""))
  }#end of if
  if(is.na(itr.cvrgval)|itr.cvrgval<2|itr.cvrgval!=round(itr.cvrgval)){
    tkmessageBox(title="Error",message=paste("The number of iteration for the exchange procedure should be a positive integer greater than or equal to two, click OK to reset.",sep=""),icon = "error"); 
    stop("Wrong value of 'itr.cvrgval', it should be greater than or equal to two (only positive integer values)")
  }#end of if
  if(is.na(nrep)|nrep<2|nrep!=round(nrep)){
    tkmessageBox(title="Error",message=paste("The number of replications should be a positive integer greater than or equal to two, click OK to reset.",sep=""),icon = "error"); 
    stop("Wrong value of 'nrep', it should be greater than or equal to two (only positive integer values)")
  }#end of if
#"===================================================================================================="
if(!is.element(Optcrit,c("A","MV","D","E"))){stop("The optimality criterion 'Optcrit' is not correctly sepcified")}
if(!is.element(Alg,c("trtE","arrayE"))){stop("The algorithm 'Alg' is not correctly sepcified")}
if(itr.cvrgval>blk.N) itr.cvrgval<-blk.N
if(Alg=="trtE") {
if(Optcrit=="A") {
optbd_mae<-Aoptbd.maeT(trt.N,blk.N,theta,nrep,itr.cvrgval)} else if(
Optcrit=="MV") {
optbd_mae<-MVoptbd.maeT(trt.N,blk.N,theta,nrep,itr.cvrgval)} else if(
Optcrit=="D") {
optbd_mae<-Doptbd.maeT(trt.N,blk.N,theta,nrep,itr.cvrgval)} else if(
 Optcrit=="E") {
optbd_mae<-Eoptbd.maeT(trt.N,blk.N,theta,nrep,itr.cvrgval)} else{
stop("The optimality criterion is not sepcified")}
optbd_mae$Alg="Treatment exchange"} else if(Alg=="arrayE") {
if(Optcrit=="A") {
optbd_mae<-Aoptbd.maeA(trt.N,blk.N,theta,nrep,itr.cvrgval)} else if(
Optcrit=="MV") {
optbd_mae<-MVoptbd.maeA(trt.N,blk.N,theta,nrep,itr.cvrgval)} else if(
Optcrit=="D") {
optbd_mae<-Doptbd.maeA(trt.N,blk.N,theta,nrep,itr.cvrgval)} else if(
 Optcrit=="E") {
optbd_mae<-Eoptbd.maeA(trt.N,blk.N,theta,nrep,itr.cvrgval)} else{
stop("The optimality criterion is not sepcified")}
optbd_mae$Alg="Array exchange"
} else {stop("The algorithm is not sepcified")}#end of if
optbd_mae$call<-match.call()
optbd_mae$Optcrit<-Optcrit
optbd_mae$Cmat<-cmatbd.mae(optbd_mae$v,optbd_mae$b,optbd_mae$theta,optbd_mae$OptdesF)
trtin <- contrasts(as.factor(optbd_mae$OptdesF), contrasts = FALSE)[as.factor(optbd_mae$OptdesF), ]
vec1 <- rep(1, optbd_mae$b * 2)
vec_trtr <- t(trtin) %*% vec1
optbd_mae$equireplicate<-all(vec_trtr==vec_trtr[1])
optbd_mae$vtrtrep<-t(vec_trtr)
#"======================================================================================"
titleoptbd<-list(c("      --------------------------------------- ",paste("Title: ",Optcrit,"-optimal or near-optimal block design          Date:", format(Sys.time(), "%a %b %d %Y %H:%M:%S"),sep=""),
"      --------------------------------------- "))
write.table(titleoptbd,  file =  file.path(tempdir(), paste(Optcrit,"optbd_",Alg,"_summary.csv",sep = "")),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
parcomb<-list(c("     Parametric combination:", "Number of treatments:", "Number of blocks:", 
"Theta value:", "Number of replications:","Number of exchange iteration:","Algorithm used:", "OPtimality criterion used:"," ","Design obtained:"),
c(" ",optbd_mae$v,optbd_mae$b,optbd_mae$theta,optbd_mae$nrep,optbd_mae$itr.cvrgval,optbd_mae$Alg,paste(Optcrit,"-optimlity criterion",sep="")," "," "))
write.table(parcomb,  file =  file.path(tempdir(), paste(Optcrit,"optbd_",Alg,"_summary.csv",sep = "")),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)

optde<-list("",rbind(paste0("Ary",1:optbd_mae$b),optbd_mae$OptdesF))
write.table(optde,  file =  file.path(tempdir(), paste(Optcrit,"optbd_",Alg,"_summary.csv",sep = "")),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
write.table(list(c("",paste(Optcrit,"-score value:",sep=""), "Equreplicate:",""),c("",optbd_mae$Optcrtsv,optbd_mae$equireplicate,"")), file =    file.path(tempdir(), paste(Optcrit,"optbd_",Alg,"_summary.csv",sep = "")),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)
write.table(list(c("Treatment:", "Treatment replication:"),rbind(1:optbd_mae$v,optbd_mae$vtrtrep)), file =   file.path(tempdir(), paste(Optcrit,"optbd_",Alg,"_summary.csv",sep = "")),append=T ,sep = ",", row.names=FALSE, col.names=FALSE)

optbd_mae$file_loc<-file.path(tempdir(),  paste(Optcrit,"optbd_",Alg,"_summary.csv",sep = ""))
optbd_mae$file_loc2<-paste("Summary of obtained ",Optcrit,"-optimal or near-optimal block design is also saved at:",sep="")
#"======================================================================================"
class(optbd_mae)<-"optbdmaeAT"
optbd_mae
}