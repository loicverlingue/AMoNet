
# NetaAnalyser search path into the net and find stuck path
# the nodes involved in stuck path can be considered as inputs
# it returns plot of the paths that can be turned into gif
# and retruns "input nodes" and circuits

# could be important to define a path for back propagation
# to do that, reverse the interactions, and propagate
# then identify where it is stucked
# and tag each species to a unique layer if possible

NetAnalyser<-function(NETall,Propagation=c("foward","back"), PDF=F, PLOT=T){
  ########
  
  ###################
  #1) plot the network
  #NETall<-NETall[order(NETall$source_hgnc),]
  NodesOut<-vector()
  NodesOutpre<-vector()
  NodesOutfreq<-vector() # check at the end the frequencies
  
  if(PLOT){
    # coordinate for the plots
    Species<-union(NETall$source_hgnc,NETall$target_hgnc)
    x<-rnorm(n = length(Species),mean = 0.5,sd = 0.2)
    names(x)<-Species
    x[NETall[NETall$Output,"target_hgnc"]]<-max(x)+0.1
    y<-rnorm(n = length(Species),mean = 0.5,sd = 0.2)
    names(y)<-Species
    #  y[unique(NETall[NETall$Output,"target_hgnc"])]<-seq(min(y),max(y),length.out = length(unique(NETall[NETall$Output,"target_hgnc"])))
    
    #par(mfrow=c(1,1))
    #par(mar=c(0,0,0,0)+0.1)
    #plot(x,y,xlim=c(-0.2,max(x)+0.3),ylim=c(-0.2,1.2), cex=0.1, axes = F)
    #plot(x,y,xlim=c(0.4,0.6),ylim=c(0.4,0.6), cex=0.1, axes = F)
    # text(x,y,Species,cex=0.5,adj = 0)
    #arrows(x0 = x[NETall$source_hgnc],x1 = x[NETall$target_hgnc],y0 = y[NETall$source_hgnc],y1 = y[NETall$target_hgnc],length = 0.05)
    #arrows(x0 = x[NETall$source_hgnc],x1 = x[NETall$target_hgnc],y0 = y[NETall$source_hgnc],y1 = y[NETall$target_hgnc],length = 0.05,col=2)
    #dev.off()
    #####
  }
  
  j=0
  while(length(NodesOutpre)!=length(unique(c(names(NodesOut),names(NodesOutpre))))|j<1 ){ #|j<2 :to put back?
    j=j+1
    print(j)
    if(length(NodesOut)>0){
      NodesOutpre<-c(NodesOut,NodesOutpre)
      NodesOutfreq<-c(NodesOut,NodesOutfreq)
      NodesOutpre<-NodesOutpre[!duplicated(names(NodesOutpre))]
      
      PosNodesOut<-which(NETall$source_hgnc%in%names(NodesOutpre))
      if("foward"%in%Propagation){
      # to order the table with inputs at the begining... not necessary because now in "layers"
        NETall<-NETall[c(PosNodesOut,setdiff(seq(nrow(NETall)),PosNodesOut)),] # maybe not the best way
        NETall[NETall$source_hgnc%in%names(NodesOutpre),"Inputs"]<-1
      }
      
      if(PDF){
        # remove previous file
        fn<-paste(getwd(),"/",list.files(getwd(),pattern = paste("Diffusion",N,".pdf",sep = "")),sep = "")
        if (file.exists(fn)) file.remove(fn)
      }
    }
    
    # explore the network
    # choose the direction of propagation
#    if("back"%in%Propagation&length(grep("Output",colnames(NETall)))>0){
    if("back"%in%Propagation&"Output"%in%colnames(NETall)){
      # reverse the interactions
      NETexplore<-NETall[,c(3,2,1,seq(ncol(NETall))[-c(1:3)])]
      colnames(NETexplore)<-colnames(NETall)
      NewSpecies<-NETexplore[NETexplore$Output,1]
      
    } else if("foward"%in%Propagation&"Output"%in%colnames(NETall)){
      NETexplore<-NETall[!NETall$Output,]
      # NETexplore<-NETexplore[!NETexplore$source_hgnc%in%names(NodesOutpre),] # otherwise paths impossible
      NewSpecies<-Species[1]
    } else {
      NETexplore<-NETall
    }
    
    Species<-union(NETexplore$source_hgnc,NETexplore$target_hgnc)
    #Species<-unique(NETexplore$source_hgnc)
  #  NewSpecies<-Species[1]
  # or select the outputs
  #  NewSpecies<-NETexplore[NETexplore$Output,1]
    NETall$Layer<-0
    MonitorSpecies<-list()
    i=1
    MonitorSpecies[[i]]<-NewSpecies
    L<-1
    Change<-vector()
    
    if(PDF){
      N<-length(list.files(getwd(),pattern = "Diffusion"))+1
      pdf(paste(getwd(),"/tmp/NetDiffusion",N,".pdf",sep = ""))
    }
    while(L[i]!=length(Species)){
      i=i+1
      
      # or
      SOU<-unique(NETexplore[NETexplore$source_hgnc%in%MonitorSpecies[[i-1]],1])
      TAR<-unique(NETexplore[NETexplore$source_hgnc%in%MonitorSpecies[[i-1]],3])
      
      MonitorSpecies[[i]]<-TAR
      L[i]<-length(unique(unlist(MonitorSpecies)))
      
      if(PLOT){
        par(mar=c(5, 2, 2, 2) + 0.1)
        plot(x,y,xlim=c(-0.1,max(x)+0.5),ylim=c(-0.1,1.1),axes=F,xlab="",ylab="")
        text(x,y,union(NETall$source_hgnc,NETall$target_hgnc),cex=0.5,adj = 0)
        arrows(x0 = x[NETall$source_hgnc],x1 = x[NETall$target_hgnc],
               y0 = y[NETall$source_hgnc],y1 = y[NETall$target_hgnc],length = 0.05)
        
        # has to replicate the source node to corresponds to the number of targets
        #  arrows(x0 =x[NewSpecies[,1]],x1 = x[NewSpecies[,2]],
        #         y0 = y[NewSpecies[,1]],y1 = y[NewSpecies[,2]],
        #         length = 0.1,col = 2)
        
        arrows(x0 =x[SOU],x1 = x[TAR],
               y0 = y[SOU],y1 = y[TAR],
               length = 0.1,col = 2)
      }
      
      NewSpecies<-unique(TAR)
      
      # tag new species as belonging to a layer in the NETall table
      # careful, it is the reverse order!! to put back at the end
      # and check of some species are in different layers. cf: NodesOut
      
#      NETall[NETall$source_hgnc%in%NewSpecies[!NewSpecies%in%unlist(MonitorSpecies[1:(i-1)])],"Layer"]<-i-1
      NETall[NETall$source_hgnc%in%NewSpecies[!NewSpecies%in%unlist(MonitorSpecies[1:(i-1)])],"Layer"]<-i-1
#      NETexplore[NETexplore$source_hgnc%in%NewSpecies[!NewSpecies%in%unlist(MonitorSpecies[1:(i-1)])],"Layers" ]<-i
      
      if(L[i]==L[i-1]){
        if(PLOT){
          mtext(text = paste(paste(NewSpecies,collapse = " "),"\nalready explored"))
        }
        Change<-c(Change,rep(i,length(NewSpecies)))
        names(Change)[(length(Change)-length(NewSpecies)+1) :length(Change)]<-NewSpecies
        
        # load new species
        NewSpecies<-sample(Species[!Species%in%names(Change)],1)
        MonitorSpecies[[i]]<-NewSpecies # reload prior selection to unstuck
      } 
    }
    
    
    # put back
    if("back"%in%Propagation){
      # reverse the order of the layers
      NETall$Layer<-abs(NETall$Layer-(max(NETall$Layer)+1))
    }
    
    # add a layer for the output node
 #     NETall[NETall$Output,"Layer"]<-max(NETall$Layer)+1
    
    if(PLOT){
      #plot the layers
      x[NETall$source_hgnc]<-NETall$Layer
      x[NETall[NETall$Output,"target_hgnc"]]<-max(x)+1
      
      for(h in unique(NETall$Layer)){
        SPy<-unique(NETall$source_hgnc[NETall$Layer%in%h])
        y[SPy]<-seq(min(y),max(y),length.out = length(SPy))
      }
      SPy<-unique(NETall$target_hgnc[NETall$Output])
      y[SPy]<-seq(min(y),max(y),length.out = length(SPy))
      
      par(mar=c(5, 2, 2, 2) + 0.1)
      plot(x,y,xlim=c(-0.1,max(x)+2),ylim=c(-0.1,1.1),axes=F,xlab="",ylab="")
      text(x,y,union(NETall$source_hgnc,NETall$target_hgnc),cex=0.5,adj = 0)
      arrows(x0 = x[NETall$source_hgnc],x1 = x[NETall$target_hgnc],
             y0 = y[NETall$source_hgnc],y1 = y[NETall$target_hgnc],length = 0.1)
      
      # plot diffusion
      #    pdf(paste(getwd(),"/NetDiffusionAnalyser",N,".pdf",sep = ""))
      par(mar=c(5, 4, 4, 2) + 0.1)
      plot(seq(i),L,type="l",ylim=c(-5,length(Species)+2),xlim = c(0,i+5),
           ylab="Cumulative N species explored", xlab = "Steps", main="Network analyser")
      if(length(Change)>0){
        arrows(x0 = unique(Change),x1 = unique(Change),y0 = L[unique(Change)]+2,y1 = L[unique(Change)],
               length = 0.05)
      }
    }
    
    
    # which are the nodes not found by propagation?
    NodesOut<-vector()
    for(h in unique(Change) ){
      if(length(Change[Change==h-1])==0){
        NodesOut<-c(NodesOut,Change[Change==h])
      }
    }
    
    if(PLOT){
      text(x = unique(as.numeric(NodesOut)),y = unique(L[as.numeric(NodesOut)]),
           labels = sapply(unique(NodesOut),function(x)paste(names(NodesOut[NodesOut==x]),collapse = " ")),
           srt=90,cex=0.5,adj = 1)
    }
     
#    dev.off()
    
    if(PDF){
      dev.off()
    }
    
    # put the species that stuck the propagation in the net, at the begining and tag : inputs
    Species<-c(Species[Species%in%names(NodesOut)],Species[!Species%in%names(NodesOut)])
    
    # redo
    
  }
  
  ##########
  # put layers in sequences
  ORD<-sort(unique(NETall$Layer))
  LAY<-seq(length(unique(ORD)))
  #names(LAY)<-ORD
  NETall$Layer<-LAY[match(NETall$Layer,ORD)]
  
  ##################
  # Identify circuits
  Circuits<-(lapply(MonitorSpecies,function(x){
    as.numeric(table(MonitorSpecies[duplicated(MonitorSpecies)]%in%x)["TRUE"])
  }))
  names(Circuits)<-lapply(MonitorSpecies, function(x) (paste(x,collapse = " ")))

#  Circuits<-unique(names(Circuits[!is.na(Circuits)])) # test
  
  return(list(Layers=NETall$Layer,Inputs=names(NodesOutpre),InputsTest=NodesOutfreq,Circuits=Circuits))
}
