# Viz learning
# c(1,"Gates","Output")
Vizlearning<-function(GatesViz=T, OutputViz=T, NETall=NETall, Step=NULL, Epoch=4,y=y,
                      learning_rate=4, NODESFIX=NULL , NETallActivity=NULL,COST=COST){

  ##############
  # representation of the gates for monitoring? Here or at the end? at the end maybe
  if(GatesViz){

    NamesGrads<-c(grep("dW",colnames(NETall),value = T),grep("db",colnames(NETall),value = T))

    Gradients<-as.matrix(NETall[!NODESFIX, NamesGrads])

    #colnames(Gradients)[unique(which(is.na(Gradients),arr.ind = T)[,2])]
    #NETall$dWeights

    Gradients<-Gradients+1e-8
    ColFonction <- colorRampPalette(c("white","blue","yellow","red"))
    logGradients<-log(abs(Gradients))

    logGradients<-logGradients-min(logGradients)

    par(mfrow=c(2,2))
    par(mar=c(5, 2, 4, 2) + 0.1)
    image(t(logGradients),col=ColFonction(1500),axes=F,
          #ylab=paste("Nodes (b) & Edges (W) \n range =", round(range(logGradients),2) ),
          main="Gradients in logscale")
    #image(t(logGradients),col=grey.colors(500),axes=F,ylab="", main="Gradients in logscale")
    mtext(text = paste("Nodes (b) & Edges (W) \n logrange =", paste(round(range(logGradients),2),collapse = ";")),
          side = 2,at = 0.5, cex=0.7,line = 0.3)
    mtext(text = NamesGrads,side = 1,at = seq(0,1,length.out = length(NamesGrads) ),las=2,cex=0.5,line = 0.3)
    grid(nx = length(NamesGrads), ny = NA, col = 1)


    ###### for paramters
    PARAM<-NETall[,c("Weights","bterm")]
    ColFonction <- colorRampPalette(c("white","blue","yellow","red"))
    image(t(PARAM),col=ColFonction(10),axes=F,
          #ylab=paste("Nodes, range =", round(range(abs(Gates)),2) ),
          main="Parameters")
    mtext(text = paste("Nodes, range = \n", paste(round(range(abs(PARAM)),2),collapse = ";") ),
          side = 2,at = 0.5,cex=0.7,line = 0.3)
    mtext(text = colnames(PARAM),side = 1,at = seq(0,1,length.out = ncol(PARAM)),las=2,cex=0.7,line = 0.3)
    grid(nx = ncol(PARAM), ny = NA, col = 1)


    par(mfrow=c(2,1))
    par(new=T)
    par(mar=c(2, 4, 1, 2))
    matplot(COST,type = "l")
    mtext(text ="Leanring curve for back simul", side = 3, cex=0.8,line = 0.3)
    legend("topright",legend = paste("Learning rate =", round(learning_rate,7)),cex = 0.4)

  }

  if(OutputViz){
    # fancy visualization of the net while learning
    par(mfrow=c(2,1))
    PlotOptNet(NETall = NETall, PDF=F, Optimized = T, NameProj = "During learning", PrintOptNET = F, LEGEND = F)

    if(length(grep("Status",rownames(y)))==1|nrow(y)==1){

      VIZ=ifelse(length(grep("Status",rownames(y)))==1,"Output",rownames(y))

      par(mar=c(2, 4, 1, 2) + 0.1)
      #NETallActivity[[1]]
      if(is.null(Step)){
        PRE<-NETallActivity[[length(NETallActivity)-1]][,VIZ]
        NOW<-NETallActivity[[length(NETallActivity)]][,VIZ]
        DIFF<-(PRE-NOW)!=0
      }else {
        PRE<-NETallActivity[[Step-1]][,VIZ]
        NOW<-NETallActivity[[Step]][,VIZ]
        DIFF<-(PRE-NOW)!=0
      }

      if(length(VIZ)==1){

        plot(as.numeric(NOW),as.numeric(y[VIZ,]),pch=21,
             xlim=c(0,1), ylim=c(0,1), cex.main=0.8,
             xlab="",ylab="True", main=paste( "Patients",paste(range(Epoch),collapse = "->"),"updates for :",VIZ), bg=adjustcolor(ifelse(DIFF,4,3), 0.3),
             col=ifelse(DIFF,4,3))
        mtext(text = "Pred", side = 1, line = 0.3)
        legend("bottomright",legend =  c("Updated",paste("Cost=", round(tail(COST,1),3))),
               pch=16, col=c(4,0), cex=0.5)

        } else{

        # check
        GT<-c(y[VIZ,])

        plot((NOW),(GT),pch=21,
             xlim=c(0,1), ylim=c(0,1), cex.main=0.8,
             xlab="",ylab="True", main=paste( "Patients",range(epoch),"updates for :",VIZ), bg=adjustcolor(ifelse(DIFF,4,3), 0.3),
             col=ifelse(DIFF,4,3))
        mtext(text = "Pred", side = 1, line = 0.3)
        legend("bottomright",legend = c("Updated",paste("Cost=", round(tail(COST,1),3))),
               pch=16, col=c(4,0), cex=0.5)

      }

    } else if(nrow(y)>1){
      par(mar=c(5, 4, 4, 2) + 0.1)

      PlotMat<-cbind(t(y),
              NETallActivity[[length(NETallActivity)-1]][colnames(y),grep("Output",rownames(y),value = T)] )

      boxplot(PlotMat,las=2, main="True    Pred")
      abline(v=ncol(PlotMat)/2+0.5,lty=2)
    }

  }
}
