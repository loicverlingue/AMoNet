#' Loads cBioportal data
#'
#' cBioportal stores publically available TCGA data.
#'
#' @param object AMoNet object.
#' @param Species character vectors. Gene names in Hugo nomenclature.
#' @param Param character vectors. Parameters to randomly search for.
#' @param RestrictUnique boolean. In case of multiple patients with same molecular features, these are removed.
#' @param PLOT boolean. Should the \code{LoadCleanTCGA} function return a waterfall plot of the mutational landscape?
#' Default is \code{FALSE}
#' @param FunctionalAnnot boolean. Should the mutations be binary (1 or 0, if \code{FunctionalAnnot=FALSE}) or quantitative (if \code{FunctionalAnnot=TRUE}).
#' Default is \code{TRUE}.
#' @param organ regexp or unique character. Data will be loaded for the type of organs selected. \code{NULL} is accepted and will select default value, \code{Default$organ="lung"}.
#'
#' @details
#' requires cgdsr package
#' requires internet connection
#'
#' @return A list with data downloaded from cBioportal, and corresponding to Default parameters.
#' y is a descrete time survival matrix whose number of intervals depend on the number of output of the model.
#' SurvData is the classical survival data of 2 dimensions with rowns: Time, Status and patients in colmuns.
#' MUTa is a mutation data matrix, either binary if FunctionalAnnot=F or numeric is FunctionalAnnot=T
#' MUT is a list format of MUTa generated with the \code{MutMatToList} function
#' Init is a normalized expression matrix
#'
#' @export
#'
LoadCleanTCGA<-function(object, Species=NULL,
                        Param="", RestrictUnique=T,
                        PLOT=F, FunctionalAnnot=T, organ=Default$organ){
  if(is.null(Species)){
    Species<-union(object$NETall$source_hgnc,object$NETall$target_hgnc)
  }

  if(is.null(organ)){
    organ<-object$Parameters$Default$organ
  } else {
    object$Parameters$Default$organ<-organ
  }

  if(length(Species)>800){ # todo improve

    TCGAdata<-LoadcBioportal(Genes = Species[Species%in%CGS$Gene.Symbol[CGS$Hallmark=="Yes"]],
                             Organ = organ, ClinicNeeded = T,
                             RNANeeded = object$Parameters$Default$EXPinit, MutNeeded = object$Parameters$Default$MUTinit,
                             FunctionalAnnot=FunctionalAnnot,
                             NormalizeRNA = T, PDF = F, Tests=T)
  } else {
    TCGAdata<-LoadcBioportal(Genes = Species, Organ = organ, ClinicNeeded = T,
                             RNANeeded = object$Parameters$Default$EXPinit, MutNeeded = object$Parameters$Default$MUTinit,
                             FunctionalAnnot=FunctionalAnnot,
                             NormalizeRNA = object$Parameters$Default$EXPinit, PDF = F, Tests=T)
  }

  ### adjust max boundary minibatch
  if("MiniBatch"%in%Param){
    object$Parameters$Boundaries$MiniBatch[2] <-  ceiling(log2(nrow(TCGAdata$CLINIC)))
    object$Parameters$Default <- HyperP(C="MiniBatch",Default = object$Parameters$Default,Boundaries = object$Parameters$Boundaries)
  }

  # remove follow up = 0
  REMOVE<-is.na(TCGAdata$CLINIC$OS_MONTHS)
  REMOVE<-REMOVE|TCGAdata$CLINIC$OS_MONTHS==0

  # if want to remove early censored patients
  REMOVE<-REMOVE|TCGAdata$CLINIC$OS_MONTHS==0|(TCGAdata$CLINIC$OS_MONTHS<quantile(TCGAdata$CLINIC$OS_MONTHS,0.4,na.rm=T)&TCGAdata$CLINIC$OS_STATUS=="LIVING")
  TCGAdata$CLINIC<-TCGAdata$CLINIC[!REMOVE,]
  TCGAdata$MUT<-TCGAdata$MUT[!REMOVE,]
  TCGAdata$EXP<-TCGAdata$EXP[!REMOVE,]

  # remove genes that are not in species - wrongly selected in LoadTCGAdata
  TCGAdata$MUT<-TCGAdata$MUT[,Species[Species%in%colnames(TCGAdata$MUT)]]
  TCGAdata$EXP<-TCGAdata$EXP[,Species[Species%in%colnames(TCGAdata$EXP)]]

  # censor long survivors
  TCGAdata$CLINIC[TCGAdata$CLINIC$OS_MONTHS>quantile(TCGAdata$CLINIC$OS_MONTHS,0.95),"death"]<-0
  TCGAdata$CLINIC[TCGAdata$CLINIC$OS_MONTHS>quantile(TCGAdata$CLINIC$OS_MONTHS,0.95),"OS"]<-as.numeric(quantile(TCGAdata$CLINIC$OS_MONTHS,0.95))

  # if Npat in the parameters to randomize: learning curve
  if(object$Parameters$Default$Npat<nrow(TCGAdata$CLINIC)){

    Selec<-sample(seq(nrow(TCGAdata$CLINIC)),object$Parameters$Default$Npat)
    Clin=TCGAdata$CLINIC[Selec,]
    MUTa=TCGAdata$MUT[Selec,]
    EXP=TCGAdata$EXP[Selec,]
    dim(Clin);dim(MUTa);dim(EXP)
  } else {
    object$Parameters$Default$Npat<-nrow(TCGAdata$CLINIC)
    Clin=TCGAdata$CLINIC
    MUTa=TCGAdata$MUT
    EXP=TCGAdata$EXP
    dim(Clin);dim(MUTa);dim(EXP)
  }

  print(paste("Data:",dim(MUTa)[1],"patients &",dim(MUTa)[2],"genes"))


  #MUT<-MutMatToList(MUTa = MUTa) # list

  # can be usefull # image at the end with results
    MUTmat<-apply(MUTa,1, function(x){
      ifelse(!is.na(x)&!x%in%"NaN",1,0)
    })

    FreqAssoc<-sort(table(apply(MUTmat,2,function(x){
      paste(rownames(MUTmat)[x==1],collapse = "_")
    })),decreasing = T)
    if(PLOT){
      par(mfrow=c(2,1))
      par(mar=c(2,4,4,2))
      ORD1<-order(rowSums(MUTmat))
      MUTmat<-MUTmat[ORD1,]
      image(t(MUTmat),col=c(0,1), xlab = "Patients", axes=F, main="Mutations")
      mtext(colnames(MUTa)[ORD1],2,las=2,cex = 0.4, at =  normalized( seq(ncol(MUTa))) )
      barplot(FreqAssoc*100/ncol(MUTmat),las=2,cex.names = 0.4,ylab = "Mutation associations %",space = 0)
      if(length(table(FreqAssoc==1))==2){
        legend("topright",legend = paste("Unique profiles =", round(table(FreqAssoc==1)["TRUE"]*100/ncol(MUTmat)),"%"),
               bty = "n")
      } else {
        legend("topright",legend = "Unique profiles = 100%", bty = "n")
      }
    }

  if(RestrictUnique){
    # remove big groups
    FreqAssoc<-FreqAssoc/nrow(MUTmat)
    PatientAssoc<-apply(MUTmat,2,function(x){
      paste(rownames(MUTmat)[x==1],collapse = "_")
    })
    #table(FreqAssoc<0.05)
    if(!all(FreqAssoc<0.05)){
      Select<-which(PatientAssoc%in%names(FreqAssoc[FreqAssoc<0.05]))

      Clin=Clin[Select,,drop=F]
      MUTa=MUTa[Select,]
      EXP=EXP[Select,]
    }
    dim(Clin);dim(MUTa);dim(EXP)
    print(paste("Restricted to:",dim(MUTa)[1],"patients &",dim(MUTa)[2],"genes"))
  }

    # mut in lists
    MUTl<-MutMatToList(MUTa = MUTa)
  ###############
  # for survival
  Y_base_raw<-Clin[,c("OS_MONTHS","OS_STATUS" )]
  Y_base_raw[,2]<-ifelse(Y_base_raw[,2]=="LIVING",0,1)
  Y_base_raw<-t(Y_base_raw)
  rownames(Y_base_raw)<-c("Output","Status")

  # 2. with descrete time estimation
  NYsurv<-SurvEstimator(Interval = ifelse(object$Parameters$Default$Interval>1,
                                          object$Parameters$Default$Interval,10), time = Y_base_raw["Output",],
                        status = Y_base_raw["Status",], PDF = F)

  if(object$Parameters$Default$Interval>1){
    Y_base<-NYsurv$ProbMat
    rownames(Y_base)<-paste("Output",seq(object$Parameters$Default$Interval),sep = "")
  } else {
    Y_base<-NYsurv$MeanProbMat
    Y_base<-t(NYsurv$MeanProbMat)
    rownames(Y_base)<-"Output"
  }
  colnames(Y_base)<-colnames(Y_base_raw)
  Clin<-Y_base

# list2env(list(y=Clin,SurvData=Y_base_raw,MUTa=MUTa,MUT=MUTl,Init=EXP), envir = .GlobalEnv)
  DATA<-list(y=Clin,SurvData=Y_base_raw,MUTa=MUTa,MUT=MUTl,Init=EXP)
  object$Data<-DATA

  # update Default and Boundaries in .GlobalEnv
  assign("Default",object$Parameters$Default,envir = .GlobalEnv)
  assign("Boundaries",object$Parameters$Boundaries,envir = .GlobalEnv)

  return(object)
}
