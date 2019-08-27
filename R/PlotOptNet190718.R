
# plot the net with layers and reclassify interactions based on learned parameters

PlotOptNet<-function(NETall=NETall, PDF=T, Optimized=T, NameProj="", PrintOptNET=T, LEGEND=F){
  
  #setwd(gsub("/tmp","",getwd()))
#  getwd()
  if(PDF){
    NameProj=paste(NameProj,"NETallopt",sep = "_")
    
    Latt<-length(list.files(paste(getwd(),"/model/",sep = ""),pattern =  NameProj))
    if(Latt>0){
      Latt<-max(as.numeric(gsub(".pdf","", gsub(NameProj,"", list.files(paste(getwd(),"/model/",sep = ""),pattern =  NameProj)))))
    } 
    
    NAME<-paste("/model/",NameProj, Latt+1,sep = "") #ifelse(Optimized,"Opt","Init"),
    pdf(paste(getwd(),"/",NAME,".pdf",sep = ""))
  }
  # plot
  
  normalized = function(x) {(x-min(x))/(max(x)-min(x))}
  
  # coordinates
  Species<-union(NETall$source_hgnc,NETall$target_hgnc)
  x<-rnorm(n = length(Species),mean = 0.5,sd = 0.2)
  names(x)<-Species
  
  x[NETall$source_hgnc]<-NETall$Layer
  x[NETall[NETall$Output,"target_hgnc"]]<-max(x)+1
  
#  y<-rnorm(n = length(Species),mean = 0.5,sd = 0.2)
  y<-seq(-0.9,0.9,length.out = length(Species))
  names(y)<-Species
  
  for(h in unique(NETall$Layer)){
    SPy<-unique(NETall$source_hgnc[NETall$Layer%in%h])
    y[SPy]<-seq(min(y),max(y),length.out = length(SPy))
  }
  SPy<-unique(NETall$target_hgnc[NETall$Output])
  y[SPy]<-seq(min(y),max(y),length.out = length(SPy))
  
  if(length(grep( "Output",NETall$target_hgnc))>0){
    #SPy<-unique(NETall$target_hgnc[NETall$Layer%in%max(NETall$Layer)])
    #SPo<-unique(NETall$source_hgnc[NETall$Layer%in%max(NETall$Layer)])
    # y[SPy]<-y[SPo]
    SPy<-unique(NETall$target_hgnc[NETall$Output])
    y[SPy]<-seq(min(y),max(y),length.out = length(SPy))
    
  }
  # size
  SIZEnode<-normalized(abs(NETall$bterm[which(NETall$source_hgnc%in%names(x))]))*4
#  if(Optimized){
#    SIZEnode<-normalized(abs(NETall$bterm[which(NETall$source_hgnc%in%names(x))]))*4
#  } else {
#    SIZEnode<-rep(1,length(x))
#    NETall$Weights<-rnorm(nrow(NETall),mean = 0.01,0.01)
#  }
  # 1 or 2 plots
 # ITER<-ifelse(Optimized,c(1:2),1)
  
#  for(iter in ITER){
    
    par(mar=c(0,0,2,2)+0.1)
    
#    plot(1,1,type="n",xlim=c(-0.1,max(x)+2),ylim=c(-0.1,1.1),axes=F,xlab="",ylab="",cex=SIZEnode)
    #  text(x,y,union(NETall$source_hgnc,NETall$target_hgnc),cex=0.5,adj = 0)
    if(!Optimized){
    plot(1,1,type="n",xlim=c(-0.1,max(x)+2),ylim=c(-1,1),axes=F,xlab="",ylab="",cex=SIZEnode)
      mtext(text = paste(NameProj,"\nNumber of species =", length(Species) ),side = 3,line = 0.3 )
      
      ACTIV<-NETall$interaction_directed_signed%in%"ACTIVATE"
      arrows(x0 = x[NETall$source_hgnc[ACTIV]],x1 = x[NETall$target_hgnc[ACTIV]],
             y0 = y[NETall$source_hgnc[ACTIV]],y1 = y[NETall$target_hgnc[ACTIV]],
             length = 0.1, col=3, lwd = normalized(abs(NETall$Weights[ACTIV]))*4)
      
      INHIB<-NETall$interaction_directed_signed%in%"INHIBIT"
      arrows(x0 = x[NETall$source_hgnc[INHIB]],x1 = x[NETall$target_hgnc[INHIB]],
             y0 = y[NETall$source_hgnc[INHIB]],y1 = y[NETall$target_hgnc[INHIB]],
             length = 0.1,angle = 90,col=2, lwd = normalized(abs(NETall$Weights[INHIB]))*4)
      points(x,y,cex=SIZEnode)
      text(x,y,union(NETall$source_hgnc,NETall$target_hgnc),cex=0.5,adj = 0)
      
    } else if(Optimized){
      #plot(1,1,type="n",xlim=c(-0.1,max(x)+2),ylim=c(-0.1,1.1),axes=F,xlab="",ylab="",cex=SIZEnode)
        plot(1,1,type="n",xlim=c(-0.1,max(x)+2),ylim=c(-1,1),axes=F,xlab="",ylab="",cex=SIZEnode)
      mtext(text = paste(NameProj,"\nNumber of species =", length(Species) ),side = 3,line = 0.3 )
      
      ACTIV<-NETall$interaction_directed_signed%in%"ACTIVATE"
      arrows(x0 = x[NETall$source_hgnc[ACTIV]],x1 = x[NETall$target_hgnc[ACTIV]],
             y0 = y[NETall$source_hgnc[ACTIV]],y1 = y[NETall$target_hgnc[ACTIV]],
             length = 0.1, col=3, lwd = normalized(abs(NETall$Weights[ACTIV]))*4)
      
      INHIB<-NETall$interaction_directed_signed%in%"INHIBIT"
      arrows(x0 = x[NETall$source_hgnc[INHIB]],x1 = x[NETall$target_hgnc[INHIB]],
             y0 = y[NETall$source_hgnc[INHIB]],y1 = y[NETall$target_hgnc[INHIB]],
             length = 0.1,angle = 90,col=2, lwd = normalized(abs(NETall$Weights[INHIB]))*4)
      
      REACTIV<-NETall$Weights>0&INHIB
      arrows(x0 = x[NETall$source_hgnc[REACTIV]],x1 = x[NETall$target_hgnc[REACTIV]],
             y0 = y[NETall$source_hgnc[REACTIV]],y1 = y[NETall$target_hgnc[REACTIV]],
             length = 0.1, col=5, lwd = normalized(abs(NETall$Weights[REACTIV]))*4)
      
      REINHIB<-NETall$Weights<0&ACTIV
      arrows(x0 = x[NETall$source_hgnc[REINHIB]],x1 = x[NETall$target_hgnc[REINHIB]],
             y0 = y[NETall$source_hgnc[REINHIB]],y1 = y[NETall$target_hgnc[REINHIB]],
             length = 0.1,angle = 90,col=6, lwd = normalized(abs(NETall$Weights[REINHIB]))*4)
      
      points(x,y,cex=SIZEnode)
      text(x,y,union(NETall$source_hgnc,NETall$target_hgnc),cex=0.5,adj = 0)
      
      if(LEGEND){
        legend("topright",legend = c("activators","inhibitors"), title = "Reclassified",fill = c(5,6), cex=0.8 )
      }
    }
    
  # reclassify interactions
  
  NETall$interaction_directed_signed[NETall$Weights<0&ACTIV]<-"INHIBIT"
  NETall$interaction_directed_signed[NETall$Weights>0&INHIB]<-"ACTIVATE"
  
  if(FALSE){
    # distribution of learned parmaters
    par(mfrow=c(1,2))
    par(mar=c(5, 2, 4, 2) + 0.1)
    plot(density(NETall$Weights), main="Weights distribution")
    plot(density( NETall$bterm[!duplicated(NETall$source_hgnc)] ),main="B terms distribution")
    par(mfrow=c(1,1))
  }
  
  if(PDF){
    dev.off()
  }
  if(PrintOptNET){
    return(NETall)
  }
}
