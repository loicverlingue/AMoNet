#' Computes a list of mutations per sample
#'
#' @param MUTa matrix. Mutation matrix with samples in rows, genes in columns, and numeric values between 0 and 1 for functional impact
#' @return A list of altereed genes, per sample, with a vector numeric values for each patient/gene with signs corresponding to either gain or loss of function
#'
MutMatToList<-function(MUTa=MUTa){
  MUTl <- lapply(seq(nrow(MUTa)),function(b){
    a=MUTa[b,]

    if(any(is.na(a))){
      NAME<-names(a)[!is.na(a)]
      a<-a[!is.na(a)]
      names(a)<-NAME
    }
    if("NaN"%in%a){
      NAME<-names(a)[a!="NaN"]
      a<-a[a!="NaN"]
      a<-rep(1,length(a))
      names(a)<-NAME
    }
    if(any(a==0)){
      NAME<-names(a)[a!=0]
      a<-a[a!=0]
      names(a)<-NAME
    }

    if(!is.numeric(a)){
      NAME<-names(a)
      Annot<-CGS[CGS$Gene.Symbol%in%names(a),c("Gene.Symbol","Role.in.Cancer")]
      a<-rep(1,length(NAME))
      a[Annot$Gene.Symbol[grep("TSG",Annot$Role.in.Cancer)]]<-(-1)
      names(a)<-NAME
    }
    return(a)
  })
  names(MUTl)<-rownames(MUTa)
  return(MUTl)
}
