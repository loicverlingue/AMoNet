###########
# identify outputs
# with immune gene sets
Outputs<-function(NETall,FamilyGene=NULL, FinalOutput=NULL, FilterFamily=1){
  if(!is.null(FamilyGene)){

#    if(file.exists(paste(getwd(),"/data/ImmunGenes.csv",sep = ""))){
#      ImmunGenes<-read.csv2(paste(getwd(),"/data/ImmunGenes.csv",sep = ""), stringsAsFactors = F)
#    } else {
#      ImmunGenes<-read.csv2(paste(gsub("/tmp","", getwd()),"/data/ImmunGenes.csv",sep = ""), stringsAsFactors = F)
#    }
    
#    FamilyGene<-ImmunGenes[ImmunGenes$gene%in%NETall$source_hgnc,c("gene","cell")]
 #   table(FamilyGene$source_hgnc%in%union(NETall$source_hgnc,NETall$target_hgnc))
    
    FamilyGene<-FamilyGene[FamilyGene$source_hgnc%in%union(NETall$source_hgnc,NETall$target_hgnc),]
    # discard output with degree < FilterFamily (<2 connection only for example)
    #ToKeep<-duplicated(FamilyGene$target_hgnc)
    if( !all(table(FamilyGene$target_hgnc)<=FilterFamily) ){
      ToKeep<-FamilyGene$target_hgnc%in%names(table(FamilyGene$target_hgnc))[table(FamilyGene$target_hgnc)>FilterFamily]
      FamilyGene<-FamilyGene[FamilyGene$target_hgnc%in%FamilyGene$target_hgnc[ToKeep],]
    } else {
      print("Careful, few connections with outputs")
    }
    
    if(!is.null(FinalOutput)){
      # add a final node that integrates the outputs (better for BackProp.R), or for outcome
      AddOutput<-data.frame(source_hgnc=rep(unique(FamilyGene$target_hgnc),FinalOutput),
                            interaction_directed_signed=rep("ACTIVATE",FinalOutput),
                            target_hgnc= paste( "Output", if(FinalOutput==1){""}else{
                              #seq(FinalOutput)
                              sort(rep(seq(FinalOutput),length(unique(FamilyGene$target_hgnc))))
                              },sep = "") )
      #sort(rep(seq(FinalOutput),length(unique(FamilyGene$target_hgnc))))
      #AddOutput[AddOutput$target_hgnc=="Output3",c(1,3)]

      FamilyGene<-rbind(FamilyGene,AddOutput)
    }
    
    # randomize interaction with outputs
    FamilyGene$interaction_directed_signed<-sample(rep(c("ACTIVATE","INHIBIT"),nrow(FamilyGene)*100),nrow(FamilyGene))
    
    # put identifiers
    FamilyGene[,setdiff(colnames(NETall), colnames(FamilyGene))]<-NA
    FamilyGene<-FamilyGene[,colnames(NETall)]
    NamesInteraction<-paste(FamilyGene[,1],FamilyGene[,3],sep = ".")
#    table(duplicated(NamesInteraction))
    FamilyGene<-FamilyGene[!duplicated(NamesInteraction),]
    
    # fusion
    #NETall<-rbind(NETall[c(1:4,grep("type",colnames(NETall)))],FamilyGene)
    NETall<-rbind(NETall,FamilyGene)
    
    # add a final node that integrates the outputs (better for BackProp.R)
    #AddOutput<-data.frame(source_hgnc=unique(NETall[NETall$Output,3]),interaction_directed_signed="ACTIVATE",target_hgnc="Output")
    #NETall<-rbind(NETall,AddOutput)
    # tag it
    NETall[,"Output"]<-F
    if(!is.null(FinalOutput)){
    #  NETall[,"Output"]<-F
#      NETall[NETall$source_hgnc%in%AddOutput$source_hgnc,"Output"]<-T  # output tag go to the source
      NETall[grep("Output", NETall$target_hgnc),"Output"]<-T  # output tag go to the source
      #      NETall[NETall$source_hgnc%in%AddOutput$source_hgnc,"Layer"]<-max(NETall$Layer,na.rm = T)+1  
    } else {
      NETall[NETall$target_hgnc%in%FamilyGene[,3],"Output"]<-T
    }
    
    return(NETall)
  } else {
    print("No table with genes to output provided")
  }
 
}


