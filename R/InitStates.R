################
# initial states

InitStates<- function(NETall,istates=istates, Logic=c("Boolean","Sigmoid","tanh"),NumiStates=20){
  Species<-union(NETall$source_hgnc,NETall$target_hgnc)
  Nspecies<-length(Species)
  
  if(Nspecies<20){
    if("Boolean"%in%Logic){
      PossiStates<-expand.grid(rep(list(c(T,F) ),Nspecies))
    } else if("Sigmoid"%in%Logic){
      PossiStates<-expand.grid(rep(list(c(0,1)), Nspecies))
    } else if("tanh"%in%Logic){
      PossiStates<-expand.grid(rep(list(c(-1,1)), Nspecies))
    }
  } else{
    # other way, by dividing the matrix
    
    # for sigmoid and tanh: do random continuous values
    
    if("Boolean"%in%Logic){
      SubPossiStates<-matrix(sample(rep(c(T,F),Nspecies*NumiStates*10),Nspecies*NumiStates),ncol = Nspecies)
    } else if("Sigmoid"%in%Logic){
      SubPossiStates<-matrix(sample(rep(c(0,1),Nspecies*NumiStates*10),Nspecies*NumiStates),ncol = Nspecies)
    # that's a challenge        SubPossiStates<-matrix(  abs(rnorm(Nspecies*NumiStates*10,mean = 0.5,sd=0.4)) ,ncol = Nspecies)
    }else if("tanh"%in%Logic){
      #      SubPossiStates<-matrix(sample(rep(c(-1,0,1),Nspecies*NumiStates*10),Nspecies*NumiStates),ncol = Nspecies)
  #          SubPossiStates<-matrix(  rnorm(Nspecies*NumiStates*10,mean = 0,sd=0.5) ,ncol = Nspecies)
      SubPossiStates<-matrix(sample(rep(c(-1,1),Nspecies*NumiStates*10),Nspecies*NumiStates),ncol = Nspecies)
    }
    
    SubPossiStates<-SubPossiStates[!duplicated(apply(SubPossiStates,1,function(x)paste(x,collapse = ""))),]
    PossiStates<-SubPossiStates
  }
  colnames(PossiStates)<-Species
  
  # that is to force particular distribution for initial states
  if(is.null(istates)){
    istates=rep(0.5,length(Species))
    names(istates)<-Species
  }
  N1<-istates*nrow(PossiStates)
  x<-N1[1]
  if("Boolean"%in%Logic){
    PossiStates<-sapply(N1,function(x){
      sample(c(rep(TRUE,round(x)),rep(FALSE,nrow(PossiStates)-round(x))))
    })
  } else {
    PossiStates<-sapply(N1,function(x){
      sample(c(rep(max(PossiStates),round(x)),rep(min(PossiStates),nrow(PossiStates)-round(x))))
    })
  }
  
  return(PossiStates)
}
