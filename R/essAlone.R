####
# ess
# data are already loaded
#dim(Clin);dim(MUTa);dim(EXP)

####
#DIR<-"C:/Users/L_VERLINGUE/Desktop/ModelK/Rpack/torque/results/GridSearch10/"
#FILES<-list.files(DIR, pattern = NameProjbase)
# for first run of TCGA grid search when using pre-optimized nets from GTEx : match when includes "MeanWinit"
#NETdata<-SelectGTExOptNet(NameProjbase = "LungRegular", # ifelse(Survival&"MeanWinit"%in%C,gsub("TCGA","",NameProjbase),NameProjbase),
#                          dir=DIR, Survival = T, Default=NULL,
#                          addOuputs =T,
#                          MiniBatch=NULL)
#NETall1<-NETdata$NETall # the best net
#Default<-NETdata$Default # the best default parameters from previous optimisations

ess<-function(Quant=30, EXP, MUTa, Clin, CGS, NETall1, Default, FixNodes, Plot=F){

  # without EXP
  no_cores <- detectCores()
  print(paste("run ess with",no_cores,"cores"))
  cl <- makeCluster(no_cores)

  EXPORT<-c("InitStates","BoolSimul", "MUTa","EXP", "Quant","FixNodes",
            "Clin","PlotOptNet","addWeights", #"Outputs","NetAnalyser"
            "BuildNetMat","CGS","survConcordance","Surv","NETall1","Default")
  clusterExport(cl,EXPORT, envir = environment())

  ScreenEss<-parLapply(cl,seq(Quant),function(Q){


    Species<-union(NETall1$source_hgnc,NETall1$target_hgnc)

    # add weigths to new nodes or remove new nodes
    #NETall1[!NETall1$target_hgnc%in%FixNodes,3]
    NETallWeights<-addWeights(NETall = NETall1[!NETall1$target_hgnc%in%FixNodes,],
                              SdWinit = Default$SdWinit,
                              MeanWinit = Default$MeanWinit ,Scaling = F,
                              Adam = T,LSTM = T)
    NETall1[!NETall1$target_hgnc%in%FixNodes,]<-NETallWeights

    print("simulation from random initial states to compute future initial stable states")
    PossiStates<-InitStates(NETall = NETall1, istates = NULL, Logic = "Sigmoid", NumiStates = 30)

    TotAttractors<-BoolSimul(NETall=NETall1, Logic = "Sigmoid", Mode = "LAYER",
                             iStates=PossiStates, MinSteps = 8,
                             LSTM = F, GOF = NULL, LOF =NULL,
                             Parallel = F, Inputs = NULL)

    iStates<-t(colMeans(TotAttractors[,"A",,dim(TotAttractors)[4]]))
    colnames(iStates)<-Species

    ######
    # muts
    MUTl <- lapply(seq(nrow(MUTa)),function(b){
      a=MUTa[b,]
      if(any(is.na(a))){
        a<-a[!is.na(a)]
      }
      NAME<-names(a)[a!="NaN"]
      a<-a[a!="NaN"]
      names(a)<-NAME
      Annot<-CGS[CGS$Gene.Symbol%in%names(a),c("Gene.Symbol","Role.in.Cancer")]

      LOF<-a[Annot$Gene.Symbol[grep("TSG",Annot$Role.in.Cancer)]]
      if(length(LOF)==0){LOF<-NULL}

      GOF<-a[!names(a)%in%names(LOF)]
      if(length(GOF)==0){GOF<-NULL}
      return(list(LOF=LOF,GOF=GOF))
    })
    names(MUTl)<-rownames(MUTa)
    LOFl<-lapply(MUTl,function(a){names(a$LOF)})
    GOFl<-lapply(MUTl,function(a){names(a$GOF)})

    #######
    # initial states

    iStates<-t(replicate(nrow(MUTa), iStates , simplify = TRUE))
    colnames(iStates)<-Species
    rownames(iStates)<-rownames(MUTa)

    #######
    # if EXP
    if(!is.null(EXP)){
      for(i in colnames(Clin)){
        iStates[i,colnames(EXP)]<-as.numeric(EXP[i,]) # match the col order
      }
    }

    # simulate
    TotAttractors<-BoolSimul(NETall = NETall1,
                             Logic = "Sigmoid", Mode = "LAYER",
                             iStates = iStates, Steepness = 1,
                             MinSteps = 8,#Default$iteration
                             Inputs=NULL, LSTM = T,
                             GOF = GOFl,
                             LOF = LOFl, Parallel = F)

    #TotAttractors[1,"A","Output",]
    #par(mar=c(5,4,4,2)+0.1)
    #  matplot(t(TotAttractors[1,"A",,]),type='l',ylim=c(0,1))
    #hist(TotAttractors[,"A","Output",dim(TotAttractors)[4]],breaks=100)

    Pred<-TotAttractors[,"A","Output",dim(TotAttractors)[4]]

    if(ncol(Clin)>5){
      #CONC=survConcordance(Surv(Y_base_raw[1,],Y_base_raw[2,])~Pred)
      CONC=survConcordance(Surv(Clin[1,],rep(1,ncol(Clin)))~Pred)
      CONC=CONC$concordance
    } else{
      CONC=NA
    }
    #      COSTS=survConcordance(Surv(Y_base_raw[1,Train],Y_base_raw[2,Train])~Pred[Train])
    #      COSTS= COSTS$concordance

    COSTS<-mean((TotAttractors[,"A","Output",dim(TotAttractors)[4]]-t(Clin))^2)

    return(list(NETall=NETall1, Cost=COSTS, Conc=CONC,Pred=Pred))

  })

  stopCluster(cl)


  Res<-t(unlist(sapply(ScreenEss,function(x){
    unlist(x[2:3])
  })))

  Cost<-min(Res[,2])[1]

  if(Plot){

    #length(ScreenEss);ScreenEss[[1]]$Conc;ScreenEss[[1]]$Cost
    par(mfrow=c(1,1))
    par(mar=c(5,4,4,2)+0.1)
    #Res<-t(unlist(sapply(ScreenEss,function(x){
    #  unlist(x[2:3])
    #})))
    plot(Res,ylab="Conc",xlab="Cost",xlim=c(0,1),ylim=c(0,1))

    #Res[which(Res[,1]<0.1),]
    #Res[which(Res[,2]<0.5),1]
    #which(Res[,2]<0.5)
    #Res[Res[,2]<0.5,2]

    par(mfrow=c(3,3))
    par(mar=c(2,2,2,2)+0.1)
    lapply(ScreenEss[which(Res[,2]<0.5)],function(x){
      plot(x$Pred,Clin,xlim = c(0,1),ylim=c(0,1))
      legend('bottomright',legend =  round(unlist(x[2:3]),3),cex=0.7, bty = "n" )
    })

    # or with cost only
    par(mfrow=c(3,3))
    par(mar=c(2,2,2,2)+0.1)
    lapply(ScreenEss[which(Res[,1]<0.1)],function(x){
      plot(x$Pred,Clin,xlim = c(0,1),ylim=c(0,1))
      legend('bottomright',legend =  round(unlist(x[2]),3),cex=0.7, bty = "n" )
    })

    #  par(mfrow=c(1,1))
    #  par(mar=c(5,4,4,2)+0.1)
    #  plot(ScreenEss[[which(Res[,2]==min(Res[,2]))[1] ]]$Pred,Clin,xlim = c(0,1),ylim=c(0,1))
#  } else {

  }

  NETall1<-ScreenEss[[which(Res[,2]==Cost)[1] ]]$NETall

  return(list(Best=NETall1,ScreenEss=ScreenEss,Cost=Cost ))
}

####
# ess
# data are already loaded
#dim(Clin);dim(MUTa);dim(EXP)

####
#DIR<-"C:/Users/L_VERLINGUE/Desktop/ModelK/Rpack/torque/results/GridSearch10/"
#FILES<-list.files(DIR, pattern = NameProjbase)
# for first run of TCGA grid search when using pre-optimized nets from GTEx : match when includes "MeanWinit"
#NETdata<-SelectGTExOptNet(NameProjbase = "LungRegular", # ifelse(Survival&"MeanWinit"%in%C,gsub("TCGA","",NameProjbase),NameProjbase),
#                          dir=DIR, Survival = T, Default=NULL,
#                          addOuputs =T,
#                          MiniBatch=NULL)
#NETall1<-NETdata$NETall # the best net
#Default<-NETdata$Default # the best default parameters from previous optimisations

ess<-function(Quant=30, EXP, MUTa, Clin, CGS, NETall1, Default, FixNodes, Plot=F){

  # without EXP
  no_cores <- detectCores()
  print(paste("run ess with",no_cores,"cores"))
  cl <- makeCluster(no_cores)

  EXPORT<-c("InitStates","BoolSimul", "MUTa","EXP", "Quant","FixNodes",
            "Clin","PlotOptNet","addWeights", #"Outputs","NetAnalyser"
            "BuildNetMat","CGS","survConcordance","Surv","NETall1","Default")
  clusterExport(cl,EXPORT, envir = environment())

  ScreenEss<-parLapply(cl,seq(Quant),function(Q){


    Species<-union(NETall1$source_hgnc,NETall1$target_hgnc)

    # add weigths to new nodes or remove new nodes
    #NETall1[!NETall1$target_hgnc%in%FixNodes,3]
    NETallWeights<-addWeights(NETall = NETall1[!NETall1$target_hgnc%in%FixNodes,],
                              SdWinit = Default$SdWinit,
                              MeanWinit = Default$MeanWinit ,Scaling = F,
                              Adam = T,LSTM = T)
    NETall1[!NETall1$target_hgnc%in%FixNodes,]<-NETallWeights

    print("simulation from random initial states to compute future initial stable states")
    PossiStates<-InitStates(NETall = NETall1, istates = NULL, Logic = "Sigmoid", NumiStates = 30)

    TotAttractors<-BoolSimul(NETall=NETall1, Logic = "Sigmoid", Mode = "LAYER",
                             iStates=PossiStates, MinSteps = 8,
                             LSTM = F, GOF = NULL, LOF =NULL,
                             Parallel = F, Inputs = NULL)

    iStates<-t(colMeans(TotAttractors[,"A",,dim(TotAttractors)[4]]))
    colnames(iStates)<-Species

    ######
    # muts
    MUTl <- lapply(seq(nrow(MUTa)),function(b){
      a=MUTa[b,]
      if(any(is.na(a))){
        a<-a[!is.na(a)]
      }
      NAME<-names(a)[a!="NaN"]
      a<-a[a!="NaN"]
      names(a)<-NAME
      Annot<-CGS[CGS$Gene.Symbol%in%names(a),c("Gene.Symbol","Role.in.Cancer")]

      LOF<-a[Annot$Gene.Symbol[grep("TSG",Annot$Role.in.Cancer)]]
      if(length(LOF)==0){LOF<-NULL}

      GOF<-a[!names(a)%in%names(LOF)]
      if(length(GOF)==0){GOF<-NULL}
      return(list(LOF=LOF,GOF=GOF))
    })
    names(MUTl)<-rownames(MUTa)
    LOFl<-lapply(MUTl,function(a){names(a$LOF)})
    GOFl<-lapply(MUTl,function(a){names(a$GOF)})

    #######
    # initial states

    iStates<-t(replicate(nrow(MUTa), iStates , simplify = TRUE))
    colnames(iStates)<-Species
    rownames(iStates)<-rownames(MUTa)

    #######
    # if EXP
    if(!is.null(EXP)){
      for(i in colnames(Clin)){
        iStates[i,colnames(EXP)]<-as.numeric(EXP[i,]) # match the col order
      }
    }

    # simulate
    TotAttractors<-BoolSimul(NETall = NETall1,
                             Logic = "Sigmoid", Mode = "LAYER",
                             iStates = iStates, Steepness = 1,
                             MinSteps = 8,#Default$iteration
                             Inputs=NULL, LSTM = T,
                             GOF = GOFl,
                             LOF = LOFl, Parallel = F)

    #TotAttractors[1,"A","Output",]
    #par(mar=c(5,4,4,2)+0.1)
    #  matplot(t(TotAttractors[1,"A",,]),type='l',ylim=c(0,1))
    #hist(TotAttractors[,"A","Output",dim(TotAttractors)[4]],breaks=100)

    Pred<-TotAttractors[,"A","Output",dim(TotAttractors)[4]]

    if(ncol(Clin)>5){
      #CONC=survConcordance(Surv(Y_base_raw[1,],Y_base_raw[2,])~Pred)
      CONC=survConcordance(Surv(Clin[1,],rep(1,ncol(Clin)))~Pred)
      CONC=CONC$concordance
    } else{
      CONC=NA
    }
    #      COSTS=survConcordance(Surv(Y_base_raw[1,Train],Y_base_raw[2,Train])~Pred[Train])
    #      COSTS= COSTS$concordance

    COSTS<-mean((TotAttractors[,"A","Output",dim(TotAttractors)[4]]-t(Clin))^2)

    return(list(NETall=NETall1, Cost=COSTS, Conc=CONC,Pred=Pred))

  })

  stopCluster(cl)


  Res<-t(unlist(sapply(ScreenEss,function(x){
    unlist(x[2:3])
  })))

  Cost<-min(Res[,2])[1]

  if(Plot){

    #length(ScreenEss);ScreenEss[[1]]$Conc;ScreenEss[[1]]$Cost
    par(mfrow=c(1,1))
    par(mar=c(5,4,4,2)+0.1)
    #Res<-t(unlist(sapply(ScreenEss,function(x){
    #  unlist(x[2:3])
    #})))
    plot(Res,ylab="Conc",xlab="Cost",xlim=c(0,1),ylim=c(0,1))

    #Res[which(Res[,1]<0.1),]
    #Res[which(Res[,2]<0.5),1]
    #which(Res[,2]<0.5)
    #Res[Res[,2]<0.5,2]

    par(mfrow=c(3,3))
    par(mar=c(2,2,2,2)+0.1)
    lapply(ScreenEss[which(Res[,2]<0.5)],function(x){
      plot(x$Pred,Clin,xlim = c(0,1),ylim=c(0,1))
      legend('bottomright',legend =  round(unlist(x[2:3]),3),cex=0.7, bty = "n" )
    })

    # or with cost only
    par(mfrow=c(3,3))
    par(mar=c(2,2,2,2)+0.1)
    lapply(ScreenEss[which(Res[,1]<0.1)],function(x){
      plot(x$Pred,Clin,xlim = c(0,1),ylim=c(0,1))
      legend('bottomright',legend =  round(unlist(x[2]),3),cex=0.7, bty = "n" )
    })

    #  par(mfrow=c(1,1))
    #  par(mar=c(5,4,4,2)+0.1)
    #  plot(ScreenEss[[which(Res[,2]==min(Res[,2]))[1] ]]$Pred,Clin,xlim = c(0,1),ylim=c(0,1))
#  } else {

  }

  NETall1<-ScreenEss[[which(Res[,2]==Cost)[1] ]]$NETall

  return(list(Best=NETall1,ScreenEss=ScreenEss,Cost=Cost ))
}

