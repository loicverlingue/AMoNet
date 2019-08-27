#' Plot *AMoNet* objects
#'
#' Can either plot the network with layers and reclassify Boolean interactions based on learned parameters (if \code{network=TRUE})
#' or plot simulations performed (if \code{simulation=TRUE})
#'
#' @param net *AMoNet* object, S3 class
#' @param LEGEND Boolean. If \code{TRUE}, add a legend to the plot
#' @param PDF Boolean. If \code{TRUE}, save a pdf version to the DIR argument
#' @param Npat numeric vector. The sample(s) for which simulations are ploted.
#' @param Species character. The name of gene(s) for which simulations are ploted.
#' @param network boolean. Whether to plot the *AMoNet* network.
#' @param simulations boolean. Whether to plot the results of *AMoNet* simulations.
#' @param DIR path. If \code{PDF=TRUE}, set the path to save pdf. Default is the \code{/tmp/} directory of AMoNet package.
#' @return In case the model is optimized it allows returning the updated Boolean relations.
#' @export
plot.AMoNet<- function(x=net, network=F, simulations=F, LEGEND=F, PDF=F, Npat=NULL,
                       DIR=file.path(getwd(),"tmp"), Species=NULL, ... ){

  if(!network&!simulations&"TotAttractors"%in%names(net)){
    print("Plot simulations")
    simulations=T
  } else {
    print("Plot network")
    network=T
  }

  if(!simulations&!is.null(Npat)){
    print("Npat selection is used to plot results of simulations. Set simultation = TRUE")
    stop()
  }
  if(simulations&!"TotAttractors"%in%names(net) ){
    print("Simulate AMoNet first, otherwise set simulation=FALSE")
    stop()
  }

  #
  list2env(net,envir = globalenv())
  NameProj<-net$call$NameProj
  if(is.null(Species)){
    Species<-union(NETall$source_hgnc,NETall$target_hgnc)
  }

  if(PDF){

    Latt<-length(list.files(DIR,pattern =  NameProj))
    if(Latt>0){
      Latt<-max(as.numeric(gsub(".pdf","", gsub(NameProj,"", list.files(DIR,pattern =  NameProj)))))
    }

    NAME<-paste(DIR,NameProj, Latt+1,".pdf",sep = "") #ifelse(Optimized,"Opt","Init"),
    pdf(NAME)
  }
  if(network){
    ##################
    # plot the network

    # coordinates
    Species<-union(NETall$source_hgnc,NETall$target_hgnc)

    x<-rnorm(n = length(Species),mean = 0.5,sd = 0.2)
    names(x)<-Species

    x[NETall$source_hgnc]<-NETall$Layer
    x[NETall[NETall$Output,"target_hgnc"]]<-max(x)+1

    y<-seq(-0.9,0.9,length.out = length(Species))
    names(y)<-Species

    for(h in unique(NETall$Layer)){
      SPy<-unique(NETall$source_hgnc[NETall$Layer%in%h])
      y[SPy]<-seq(min(y),max(y),length.out = length(SPy))
    }
    SPy<-unique(NETall$target_hgnc[NETall$Output])
    y[SPy]<-seq(min(y),max(y),length.out = length(SPy))

    if(length(grep( "Output",NETall$target_hgnc))>0){
      SPy<-unique(NETall$target_hgnc[NETall$Output])
      y[SPy]<-seq(min(y),max(y),length.out = length(SPy))

    }
    # size
    SIZEnode<-normalized(abs(NETall$bterm[which(NETall$source_hgnc%in%names(x))]))*4

    par(mar=c(0,0,2,2)+0.1)
    if(!"xlim|ylim"%in%names(match.call())){
      plot(1,1,type="n",axes=F, xlim=c(-0.1,max(x)+2), ylim=c(-1,1), cex=SIZEnode, ...)
    } else {
      plot(1,1,type="n",axes=F, cex=SIZEnode, ...)
    }
    mtext(text = paste(NameProj,"\nNumber of species =", length(Species) ),side = 3,line = 0.3 )

    # remove warnings in arrows
    options(warn=-1)

    ACTIV<-NETall$interaction_directed_signed%in%"ACTIVATE"
    arrows(x0 = x[NETall$source_hgnc[ACTIV]],x1 = x[NETall$target_hgnc[ACTIV]],
           y0 = y[NETall$source_hgnc[ACTIV]],y1 = y[NETall$target_hgnc[ACTIV]],
           length = 0.1, col=3, lwd = normalized(abs(NETall$Weights[ACTIV]))*4)

    INHIB<-NETall$interaction_directed_signed%in%"INHIBIT"
    arrows(x0 = x[NETall$source_hgnc[INHIB]],x1 = x[NETall$target_hgnc[INHIB]],
           y0 = y[NETall$source_hgnc[INHIB]],y1 = y[NETall$target_hgnc[INHIB]],
           length = 0.1,angle = 90,col=2, lwd = normalized(abs(NETall$Weights[INHIB]))*4)

    REACTIV<-NETall$Weights>0&INHIB
    if(any(REACTIV)){
      arrows(x0 = x[NETall$source_hgnc[REACTIV]],x1 = x[NETall$target_hgnc[REACTIV]],
             y0 = y[NETall$source_hgnc[REACTIV]],y1 = y[NETall$target_hgnc[REACTIV]],
             length = 0.1, col=5, lwd = normalized(abs(NETall$Weights[REACTIV]))*4)
    }
    REINHIB<-NETall$Weights<0&ACTIV
    if(any(REINHIB)){
      arrows(x0 = x[NETall$source_hgnc[REINHIB]],x1 = x[NETall$target_hgnc[REINHIB]],
             y0 = y[NETall$source_hgnc[REINHIB]],y1 = y[NETall$target_hgnc[REINHIB]],
             length = 0.1,angle = 90,col=6, lwd = normalized(abs(NETall$Weights[REINHIB]))*4)
    }

    options(warn=0)

    points(x,y,cex=SIZEnode)
    text(x,y,union(NETall$source_hgnc,NETall$target_hgnc),cex=0.5,adj = 0)

    if(LEGEND){
      legend("topright",legend = c("activators","inhibitors"), title = "Reclassified",fill = c(5,6), cex=0.8 )
    }


    # reclassify interactions

    NETall$interaction_directed_signed[NETall$Weights<0&ACTIV]<-"INHIBIT"
    NETall$interaction_directed_signed[NETall$Weights>0&INHIB]<-"ACTIVATE"

    net$NETall<<-NETall
    #      class(net)<-"AMoNet"
    #      return(net)
  }

#  if(!is.null(Npat)&"TotAttractors"%in%names(net)){
  if(simulations){
    ############
    # plot simulation for 1 patient
    if("col"%in%names(match.call())){
      matplot(rbind(t(TotAttractors[Npat,"iStates",Species,1]),t(TotAttractors[Npat,"A",Species,])),type='l',
              ylab="Species activities",xlab = "Time steps", main=paste("Simulation of the network \n for sample :", paste(Npat,collapse = ", ") ),
              lty=1, ...)
    } else {
      COL=rainbow(dim(TotAttractors)[3])
      matplot(rbind(t(TotAttractors[Npat,"iStates",Species,1]),t(TotAttractors[Npat,"A",Species,])),type='l',
              ylab="Species activities",xlab = "Time steps",  main=paste("Simulation of the network \n for sample :", paste(Npat,collapse = ", ")),
              lty=1, col=COL, ...)
    }

    if(LEGEND){

      if(length(Npat)<=1){
        legend("topleft",legend = Species,col=COL, lty=1)
      } #else {
        #legend("topleft",legend = paste("Patient",Npat), col=COL, lty=1)
      #}

    }
    #    matplot(as.matrix(c(TotAttractors[Npat,"iStates","Output",1],TotAttractors[Npat,"A","Output",])), type='l', col=1,lwd=3, add=T)
}
  #}else{
    if(PDF){
      dev.off()
    }
}

if(FALSE){
plot.build.AMoNet<- function(object=net$NETall, LEGEND=F, PDF=F,
                       DIR=file.path(getwd(),"tmp"),
                       type="n",axes=F, xlim=c(-0.1,3), ylim=c(-1,1), cex=0.5,
                       ... ){

  list2env(net,envir = globalenv())
  NameProj<-net$call$NameProj
  Species<-union(object$source_hgnc,object$target_hgnc)

  if(PDF){

    Latt<-length(list.files(DIR,pattern =  NameProj))
    if(Latt>0){
      Latt<-max(as.numeric(gsub(".pdf","", gsub(NameProj,"", list.files(DIR,pattern =  NameProj)))))
    }

    NAME<-paste(DIR,NameProj, Latt+1,".pdf",sep = "") #ifelse(Optimized,"Opt","Init"),
    pdf(NAME)
  }

  ##################
    # plot the network

    # coordinates
    x<-rnorm(n = length(Species),mean = 0.5,sd = 0.2)
    names(x)<-Species

    x[object$source_hgnc]<-object$Layer
    x[object[object$Output,"target_hgnc"]]<-max(x)+1

    y<-seq(-0.9,0.9,length.out = length(Species))
    names(y)<-Species

    for(h in unique(object$Layer)){
      SPy<-unique(object$source_hgnc[object$Layer%in%h])
      y[SPy]<-seq(min(y),max(y),length.out = length(SPy))
    }
    SPy<-unique(object$target_hgnc[object$Output])
    y[SPy]<-seq(min(y),max(y),length.out = length(SPy))

    if(length(grep( "Output",object$target_hgnc))>0){
      SPy<-unique(object$target_hgnc[object$Output])
      y[SPy]<-seq(min(y),max(y),length.out = length(SPy))

    }
    # size
    SIZEnode<-normalized(abs(object$bterm[which(object$source_hgnc%in%names(x))]))*4

    par(mar=c(0,0,2,2)+0.1)

    plot(1,1,type="n",axes=F, cex=SIZEnode, xlim=c(-0.1,max(x)+2), ylim=ylim, ...)

    mtext(text = paste(NameProj,"\nNumber of species =", length(Species) ),side = 3,line = 0.3 )

    ACTIV<-object$interaction_directed_signed%in%"ACTIVATE"
    arrows(x0 = x[object$source_hgnc[ACTIV]],x1 = x[object$target_hgnc[ACTIV]],
           y0 = y[object$source_hgnc[ACTIV]],y1 = y[object$target_hgnc[ACTIV]],
           length = 0.1, col=3, lwd = normalized(abs(object$Weights[ACTIV]))*4)

    INHIB<-object$interaction_directed_signed%in%"INHIBIT"
    arrows(x0 = x[object$source_hgnc[INHIB]],x1 = x[object$target_hgnc[INHIB]],
           y0 = y[object$source_hgnc[INHIB]],y1 = y[object$target_hgnc[INHIB]],
           length = 0.1,angle = 90,col=2, lwd = normalized(abs(object$Weights[INHIB]))*4)

    REACTIV<-object$Weights>0&INHIB
    if(any(REACTIV)){
      arrows(x0 = x[object$source_hgnc[REACTIV]],x1 = x[object$target_hgnc[REACTIV]],
             y0 = y[object$source_hgnc[REACTIV]],y1 = y[object$target_hgnc[REACTIV]],
             length = 0.1, col=5, lwd = normalized(abs(object$Weights[REACTIV]))*4)
    }
    REINHIB<-object$Weights<0&ACTIV
    if(any(REINHIB)){
      arrows(x0 = x[object$source_hgnc[REINHIB]],x1 = x[object$target_hgnc[REINHIB]],
             y0 = y[object$source_hgnc[REINHIB]],y1 = y[object$target_hgnc[REINHIB]],
             length = 0.1,angle = 90,col=6, lwd = normalized(abs(object$Weights[REINHIB]))*4)
    }
    points(x,y,cex=SIZEnode)
    text(x,y,union(object$source_hgnc,object$target_hgnc),cex=0.5,adj = 0)

    if(LEGEND){
      legend("topright",legend = c("activators","inhibitors"), title = "Reclassified",fill = c(5,6), cex=0.8 )
    }


    # reclassify interactions

    object$interaction_directed_signed[object$Weights<0&ACTIV]<-"INHIBIT"
    object$interaction_directed_signed[object$Weights>0&INHIB]<-"ACTIVATE"

    net$NETall<<-object
    #      class(net)<-"AMoNet"
    #      return(net)


  if(PDF){
    dev.off()
  }

}
}
