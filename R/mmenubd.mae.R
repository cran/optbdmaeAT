#Function for the window of main menu of TCL/TK 
mmenubd.mae<-function(){
  fontHeading<- tkfont.create(family="times",size=40,weight="bold")#,slant="italic")
  fontHeading3<-tkfont.create(family="times",size=10,weight="bold")
  AMVDEopt.top<-tktoplevel()
  tkwm.title(AMVDEopt.top,"Optimal Block Designs For Microarray Experiments")
  tkgrid(tklabel(AMVDEopt.top,text="     optbdmaeAT 1.0.1     ",font=fontHeading))
  tkgrid(tklabel(AMVDEopt.top,text="",font=fontHeading))
  Fixp.butA<-tkbutton(AMVDEopt.top,text="A-Optimal Block Design",font=fontHeading3,command=function()fixparbd.mae("A"),width=30)
  Fixp.butMV<-tkbutton(AMVDEopt.top,text="MV-Optimal Block Design",font=fontHeading3,command=function() fixparbd.mae("MV"),width=30)
  Fixp.butD<-tkbutton(AMVDEopt.top,text="D-Optimal Block Design",font=fontHeading3,command=function() fixparbd.mae("D"),width=30)
  Fixp.butE<-tkbutton(AMVDEopt.top,text="E-Optimal Block Design",font=fontHeading3,command=function() fixparbd.mae("E"),width=30)
  ExitWind.1<-function(){ tkmessageBox(title="Bye...", message=paste("Bye..., Enjoy your optimal design")) 
    tkdestroy(AMVDEopt.top)
  }#End of ExitWind.1
  ExitWin.but<-tkbutton(AMVDEopt.top,text="Exit",font=fontHeading3,command=ExitWind.1,width=15)
  tkgrid(Fixp.butA)
  tkgrid(Fixp.butMV)
  tkgrid(Fixp.butD)
  tkgrid(Fixp.butE)
  tkgrid(ExitWin.but)
  tkgrid(tklabel(AMVDEopt.top,text="",font=fontHeading))   }#end of mmenu.mae function
mmenubd.mae()

