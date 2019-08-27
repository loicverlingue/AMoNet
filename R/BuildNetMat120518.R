#buid netMat
#Layer=4
BuildNetMat<-function(NETall=NETall,Layer=NULL,Type="Weights"){
  
  #COLO<-paste(NETall$source_hgnc,NETall$interaction_directed_signed,NETall$target_hgnc,sep = "_")
 # COLO<-unique(NETall$target_hgnc[NETall$Layer==Layer])
  
  if(is.null(Layer)){
    COLO<-union(NETall$source_hgnc,NETall$target_hgnc)
    ROW<-unique(NETall$source_hgnc)
    
  } else{
    # because there are recurrent relations:
    #COLO<-union(NETall$source_hgnc[NETall$Layer==Layer+1],NETall$target_hgnc[NETall$Layer==Layer])
    COLO<-unique(NETall$target_hgnc[NETall$Layer==Layer])
    
    ROW<-union(NETall$source_hgnc[NETall$Layer==Layer],NETall$target_hgnc[NETall$Layer==Layer-1])
    #ROW<-unique(NETall$source_hgnc[NETall$Layer==Layer])
  }
  
  # build mat
  NetMat<-matrix(data = 0,nrow = length(ROW),ncol = length(COLO))
  rownames(NetMat)<-ROW;colnames(NetMat)<-COLO
  #Weigths<-NETall$Weights
  
 # NETall[NETall$source_hgnc%in%ROW&NETall$target_hgnc%in%COLO,"Weights"]
  
  for(i in which(NETall$source_hgnc%in%ROW&NETall$target_hgnc%in%COLO)){
    # N<-runif(1,0,10)
    #  N<-ifelse(NETall[i,"interaction_directed_signed"]%in%"ACTIVATE",N,-N)
    N<-NETall[i,Type]
    NetMat[NETall[i,"source_hgnc"],NETall[i,"target_hgnc"]]<-N
  }
  return(NetMat)
}
