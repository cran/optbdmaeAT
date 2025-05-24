#Section 5: Function for the plot of the graphical layout of resultant optimal designs (graphoptbd.mae)
graphoptbd.mae<-function(trt.N, blk.N,theta,OptdesF,Optcrit,cbVal2) {
  #cbValue2.I<-tclVar("0")
  #cbVal2=as.numeric(tclvalue(cbValue2.I))
  trtblkthetano<-paste("(",paste(trt.N, blk.N, theta,sep=", "),")",sep="")
  
  if(cbVal2==0) {
  Alg="trtE";optcN="treatment"} else if(cbVal2==1) {
  Alg="arrayE";optcN="array"} else {stop("The algorithm is not sepcified")
  } #End if (identify the treatment and array exchange algorithm using checkbox value)
  NOptcrtr<-paste(Optcrit,"-optimal",sep="")#name of optimality criteria
  NOptcrtrG<-paste("Graph_layout_",Optcrit,"optbd_",Alg,sep="")#name of folder where the graphical layout will be saved
  NOptcrtrG2<-paste("_Gout",Optcrit,"optbd",Alg,".pdf",sep="")
  NgoutT=paste(NOptcrtr, "block", "design", "for", paste("(v, b, theta) =",trtblkthetano,sep=" "))
  NgoutST=paste("using",optcN,"exchange","algorithim",sep=" ")
  graph.des <- make_graph(as.numeric(as.factor(OptdesF)), directed = FALSE)
  graph.desid <- tkplot(graph.des, canvas.width=515, canvas.height=500,layout=layout.kamada.kawai,vertex.color="cyan",edge.color="black")
  canvas <- tk_canvas(graph.desid)
  padding <- 100
  coords <- norm_coords(layout=layout.kamada.kawai(graph.des), 0+padding, 450-padding,
                        50+padding, 500-padding)
  tk_set_coords(graph.desid, coords)
  width <- as.numeric(tkcget(canvas, "-width"))
  height <- as.numeric(tkcget(canvas, "-height"))
  tkcreate(canvas, "text", width/2, 25, text=NgoutT,
           justify="center", font=tcltk::tkfont.create(family="helvetica",size=15,weight="bold"))
  tkcreate(canvas, "text", width/2, 45, text=NgoutST,
           justify="center", font=tcltk::tkfont.create(family="helvetica",size=15,weight="bold"))
  graph.OutlayoptBlk<-tempdir()
    #file.path(tempdir(),  paste(NOptcrtrG,sep=""))
  if(!file.exists(graph.OutlayoptBlk)) dir.create(graph.OutlayoptBlk)
  obtdes.goutloptBlk<-paste(graph.OutlayoptBlk,paste(trtblkthetano,NOptcrtrG2,sep=""),sep="/")
  pdf(file=obtdes.goutloptBlk)
  plot(graph.des,edge.arrow.size=5, vertex.size=20, margin=0.5,
       layout=layout.kamada.kawai,vertex.color="cyan",edge.color="black")
title(paste("Graphical layout of ", Optcrit,"-optimal or near-optimal block design",sep=""), 
      sub = NULL,cex.main = 1,   font.main= 1, col.main= "black")
mtext(paste(NgoutST," for:",sep=""), line = 0.5, col = "black", font = 1)
mtext(paste("(v, b, theta) =", " (",paste(trt.N, blk.N, theta,sep=", "),")",sep=""), line = -0.50, col = "blue", font = 1)
  dev.off()  
  file_loc<-obtdes.goutloptBlk
  file_loc2<-paste("Graphical layout of obtained", NOptcrtr, "or near-optimal block design is also saved in .pdf at:",sep=" ")
  cat(file_loc2,"\n",file_loc,"\n\n")
}#End Section 5 (plot of the graphical layout of resultant optimal design, graphoptbd.mae)

