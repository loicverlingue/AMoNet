#' Wrapper function to run the AMonet workflow, including hyper-parameters grid search, on TCGA data
#'
#' @description AMonet workflow comprises building, training and simulating a network model to predict survival of patients from TCGA genomics WES data.
#'
#' @param Param character vector: the hyper-parameter(s) to randomly test whithin this run.
#' @param DIR directory: directory of the models storing
#' @param NameProj character: name of your project. If at least one model with the name of your project NameProj is in the DIR directory, the script will load it and use it.
#' @param NewNet boolean: when set to FALSE, if at least one model with the name of your project NameProj is in the DIR directory, the script will load it and use it. If set to TRUE, a new molecular network will be built using arguments GENESman and default parameters (Interval, MinConnect, nblayers).
#' @param GENESman character vector: gene(s)' name(s) (in Hugo nomenclature) used to build the AMoNET network. If "all", all genes related to SelectMECA (gene sets) will be selected. If "frequent", 95\% most frequently altered genes from the cancer gene census (CGS) list in the cohort will be selected.
#' @param treatmt character vector: gene(s)' name(s) targeted by simulated treatment(s). Added to GENEman toused to build the AMoNET network.
#' @param SelectMECA character in regexp format: selection of one or several gene sets related to biological mecanisms. Names gene sets are available with the command \code{print(names_MECA)}.
#' @param organ character in regexp format: selection of one or several organs or cancer types. Names of organs are available with the command \code{print(names_MECA)}
#' @param eSS boolean to perform multistart or not. Default is set to \code{FALSE}. Recommandation is to use \code{eSS=TRUE} within a hyper-parameters search and in case of \code{Param = "MeanWinit"} and/or \code{Param = "SdWinit"} eSS \code{TRUE}.
#' @param ... any argument(s) present in Default can be set here by the user. \code{print(Default)} for information on Default hyper-parameters.
#'
#' @details
#' This function can be run:
#' 1) Alone to train an AMoNet model with fixed hyper-parameters are either user-defined or the default one are used. For default hyper-parameters values call \code{print(Default)}.
#' 2) In array and iteratively from command line together with the \code{PlotAndPredict} function to perform hyper-parameters search. Hyperparameter selected in the \code{Param} argument are thus randomly selected in batches and best results selected.
#'
#' Requirements:
#' Uses the pre-configured /model and /tmp paths in AMoNet package directories (checks todo)
#'
#'
#' @return todo AMoNet object
#'
RunTCGAopt<-function(Param=c("nblayers", "MinConnect"), DIR=file.path(getwd(),"model"),
                     NameProj="HallmarksLungHN", GENESman=c("EGFR", "KRAS", "TP53", "MTOR"),
                     treatmt="", SelectMECA="HALLMARK", organ="luad",
                     eSS=F, NewNet=T, no_cores=3, KeepData=T, PartitionSplit=0.7, ...){

  set.seed(NULL)

  NameProjbase<-NameProj

  ###############
  # initiate AMoNet
  if(is.null(GENESman)&is.null(treatmt)){
    print("You should select a list of genes in initial function to build the net")
    stop()
  }

  net<-AMoNet(GENESman = GENESman, treatmt = treatmt)

  # package should be loaded with AMoNet package
  # todo
if(FALSE){
  HyperP<-function(C=Param, Default=Default, Boundaries=Boundaries){

    # define distribution you want
    LogDistrib<-function(min,max,N){
      sample(exp(seq(log(abs(min)),log(abs(max)),length.out = 100)),N)
    }

    NormDistrib<-function(min,max,N){
      sample(seq(min,max,length.out = 100),N)
    }

    for(Cp in C){
      if(is.numeric(Boundaries[[Cp]])){
        if(is.integer(Boundaries[[Cp]])){
          Default[[Cp]]<- sample(seq(Boundaries[[Cp]][1],Boundaries[[Cp]][2]),1)
        } else if(all(Boundaries[[Cp]]>=0.1)){
          Default[[Cp]]<- NormDistrib(min = Boundaries[[Cp]][1],max = Boundaries[[Cp]][2], N = 1)
        }else{
          Default[[Cp]]<- LogDistrib(min = Boundaries[[Cp]][1],max = Boundaries[[Cp]][2], N = 1)
        }
      } else if(is.character(Boundaries[[Cp]])|is.logical(Boundaries[[Cp]])){
        Default[[Cp]]<- sample(Boundaries[[Cp]],1)
      }
    }

    if("MiniBatch"%in%C){
      Default$MiniBatch<-sample(2^seq(Boundaries$MiniBatch[1],Boundaries$MiniBatch[2]),1)
    }

    list2env(Default, envir = .GlobalEnv)
    return(Default)
  }
}

  # update Default parameters
  CALL<-mget(names(formals()))
  CNames<-intersect(names(Default),names(CALL))
  Default[CNames] <- CALL[CNames]

  if(FALSE){
  #### control optional arguments
  # totest
  #CLnames<-list(...)
  CLnames<-match.call()
  CLnames<-CLnames[names(CLnames)%in%names(Default)]
  if(length(CLnames)>0){
    Default[names(CLnames)]<-CLnames
    print("updated:")
    print(Default[names(CLnames)])
  }
  }

  ########
  # Genes sets for phenotypes: take all here
  GenesSelec<-rbind(AMoNet::GenesSelectImmuno, AMoNet::GenesSelecHall)
  MECA<-unique(GenesSelec$target_hgnc) # do an interactcive function to choose MECA
  if(!is.null("SelectMECA")){
    MECA<-grep(SelectMECA,MECA,value = T)
  }
  if("all"%in%tolower(GENESman)){
    GENESman<-unique(GenesSelec$source_hgnc[GenesSelec$target_hgnc%in%MECA])
  }


  ###############
  # load networks or build
  if(!dir.exists(DIR)){
    dir.create(DIR,showWarnings = T)
    print("New direction -> new nets")
  }

  FILES<-list.files(DIR,pattern = NameProjbase)

  if(NewNet){
    FILES<-NULL
  }

  if(length(FILES)>1){
    # to change++++++
    # for first run of TCGA grid search when using pre-optimized nets from GTEx : match when includes "MeanWinit"
    NETdata<-SelectGTExOptNet(NameProjbase = NameProjbase, NETall = NULL,
                              DIR=DIR, Default=Default,
                              addOutputs = NULL, MiniBatch=Default$MiniBatch,
                              ValSelect=T)

    NETall1<-NETdata$NETall # the best net
    Default<-NETdata$Default # the best default parameters from previous optimisations
    Boundaries<-NETdata$Boundaries
    Train<-NETdata$NETallProp$TrainSplit$Train
    Val<-NETdata$NETallProp$TrainSplit$Val
    organ<-NETdata$NETallProp$organ

    # update Hyperparameters with new ones
    Default<-HyperP(C=Param, Default = Default, Boundaries = Boundaries)


  } else if(length(FILES)==1){
 # change with preBuild function

    # update Hyperparameters with new ones
    Default<-HyperP(C=Param,Default = Default,Boundaries = Boundaries)

    if(length(grep("csv",FILES))>0){
      NETall1<-read.csv2(file = paste(DIR,FILES[1],sep = ""), row.names = 1, header = T,stringsAsFactors = F)
      NETall1$X<-NULL

      net<-AMoNet(GENESman=unique(NETall1$source_hgnc),treatmt=treatmt)
      net<-build.AMoNet(net,
                        InteractionBase = NETall1,
                        nblayers = Default$nblayers, MinConnect = Default$MinConnect,
                        RestrictedBuilding = T, RestrictedBase = F, FilterCGS = F,
                        MeanWinit = Default$MeanWinit, SdWinit = Default$SdWinit,
                        MECA=NULL, Phenotypes=NETall1[NETall1$target_hgnc%in%NETall1$source_hgnc[NETall1$Output],1:3],
                        Interval = Default$Interval, LSTM = Default$LSTM, Optimizer = Default$Optimizer,
                        KeepPhenotypes=T, WRITE = F, no_cores = no_cores, NameProj = NameProjbase)

    } else {
      # to check
      x<-load(file.path(DIR,FILES))
      net<-get(x)
      #    NETall1<-NETallProp$NETallList[[length(NETallProp$NETallList)]]
      #    Default<-NETallProp$Parameters$Default # the best default parameters from previous optimisations
      #    Boundaries<-NETallProp$Parameters$Boundaries
      #    Train<-NETallProp$TrainSplit$Train
      #    Val<-NETallProp$TrainSplit$Val
      #    organ<-NETallProp$organ
    }

  } else { # this is used for the first net building

    # update Hyperparameters with new ones
    Default<-HyperP(C=Param,Default = Default, Boundaries = Boundaries)

    net<-build.AMoNet(object = net,
                      InteractionBase = OMNI,
                      nblayers = Default$nblayers, MinConnect = Default$MinConnect,
                      RestrictedBuilding = T, RestrictedBase = T, FilterCGS = F,
                      MeanWinit = Default$MeanWinit, SdWinit = Default$SdWinit,
                      Phenotypes=GenesSelec, MECA=MECA,
                      Interval = Default$Interval, LSTM = Default$LSTM, Optimizer = Default$Optimizer,
                     KeepPhenotypes=T,no_cores = no_cores, NameProj = NameProjbase)
  }

  Species<-union(net$NETall$source_hgnc,net$NETall$target_hgnc)

  ###############
  print("load TCGA data")
  # run the function to load data

  DATA<-LoadCleanTCGA(net, Species = Species,
                      Param = Param, RestrictUnique = F)

  #DATA$Perturb<-MutMatToList(MUTa = DATA$Perturb) # to list
  #lapply(DATA, dim)
  list2env(DATA,envir = globalenv())

  # rename project
  QUANTS<-unlist(sapply(Param,function(PR){
    PARAM<-eval(parse(text=PR) ,envir = Default)
    NAME<-paste(PR,ifelse(is.na(as.numeric(PARAM)),PARAM,as.numeric(PARAM)),sep = "_") #round(as.numeric(PARAM),4)
    return(NAME)
  }))

  NameProj<-net$call$NameProj<-paste(NameProjbase, paste(QUANTS,collapse = "_"),"_", sep = "")
  print(paste("Welcome to project", NameProj))

  ###############
  # compute initial states
  RandomiStates(net)

  #iStates<-matrix(runif(ncol(y)*length(Species)),nrow = ncol(y), ncol = length(Species))
  #rownames(iStates)<-colnames(y)
  #colnames(iStates)<-Species

  ##### EXPi # to improve by using it as constrained activity value
  if(Default$EXPinit){
    for(i in colnames(y)){
      iStates[i,colnames(Init)]<-as.numeric(Init[i,]) # match the col order
    }
  }

  if(Default$LSTM){ # to check if doesn impare LSTM
    Ct<-t(replicate(ncol(y), rep(1,length(Species)) , simplify = TRUE))
    colnames(Ct)<-Species
    rownames(Ct)<-colnames(y)
  } else{
    Ct=NULL
  }

  ####################
  print("run learning")

  if(Default$FixWeights){
    FixNodes<-Species[!Species%in%c(net$NETall$source_hgnc[net$NETall$target_hgnc%in%"Output"],"Output")]
  } else {
    FixNodes<-NULL
  }

  if(!"TrainSplit"%in%names(net)){
    net<-split(net,PartitionSplit)
    list2env(net$TrainSplit, envir = globalenv())
  } else {
    list2env(net$TrainSplit, envir = globalenv())

    Train<-Train[Train%in%colnames(y)]
    Val<-Val[Val%in%colnames(y)]

    MorePat<-colnames(y)[!colnames(y)%in%c(Train,Val)]
    if(length(MorePat)>0){

      MPTrain<-sample(MorePat,size = round(length(MorePat)*0.7) )
      MPVal<-MorePat[!MorePat%in%MPTrain]

      Train<-c(Train,MPTrain)
      Val<-c(Val,MPVal)
      net$TrainSplit$Train<<-Train
      net$TrainSplit$Val<<-Val
    }
  }

if(FALSE){
  if(!exists("Train")|!exists("Val")){

    Train<-sample(colnames(y),ncol(y)*0.7)
    Val<-colnames(y)[!colnames(y)%in%Train]

  } else if(is.null(Train)|is.null(Val)){
    Train<-sample(colnames(y),ncol(y)*0.7)
    Val<-colnames(y)[!colnames(y)%in%Train]

  } else {

    Train<-Train[Train%in%colnames(y)]
    Val<-Val[Val%in%colnames(y)]

    MorePat<-colnames(y)[!colnames(y)%in%c(Train,Val)]
    if(length(MorePat)>0){

      MPTrain<-sample(MorePat,size = round(length(MorePat)*0.7) )
      MPVal<-MorePat[!MorePat%in%MPTrain]

      Train<-c(Train,MPTrain)
      Val<-c(Val,MPVal)

    }
  }
}
  # cap minibatches
  Default$MiniBatch<-min(Default$MiniBatch,length(Train))
  if(Default$MiniBatch==length(Train)){Default$Optimizer=NULL}

  print(paste("Training on",length(Train),"; validation on",length(Val), "patients"))

  if(!eSS){

    tic=Sys.time()

    net<-train.AMoNet(net,y=y[,Train,drop=F],
           MUT=MUT[Train], treatmt=NULL,
           Init = Init,
           ValMut = Default$ValMut,
           iStates=iStates[Train,],
           Ct=Ct[Train,],
           MiniBatch=Default$MiniBatch,
           alpha = Default$alpha, lambda = Default$lambda,
           iteration=Default$iteration,
           beta1=Default$beta1, beta2=Default$beta2,
           Parallelize= Default$Parallelize, no_cores = no_cores,
           learning_rate=Default$learningrate,
           LearningRateDecay = Default$LearningRateDecay,
           adaptive_iStates= Default$adaptive_iStates,
           Optimizer=Default$Optimizer,
           Logic = Default$Logic, Mode = Default$Mode, ModeBack = Default$Mode,
           MinStepsForward = Default$MinStepsForward,
           MinStepsBackward = Default$MinStepsBackward,
           gradClipping=Default$gradClipping,
           LSTM=Default$LSTM,
           FixNodes = FixNodes, NameProj = net$call$NameProj,
           PDF=F, GIF=F, Visualize=NULL, KeepData=F,  KeepTraining=T)

    if(FALSE){
      NETallProp<-RunBackSimul(NETall=NETall1, y=y[,Train,drop=F],
                               MUT=MUT[Train], ValMut = Default$ValMut,
                               Init=Init[Train,], treatmt = NULL,
                               iStates=iStates[Train,],
                               Ct=Ct[Train,],
                               MiniBatch=Default$MiniBatch,
                               iteration=round(Default$iteration),
                               beta1=Default$beta1, beta2=Default$beta2,
                               Parallelize= Default$Parallelize, no_cores=no_cores,
                               learning_rate=Default$learningrate,
                               LearningRateDecay = Default$LearningRateDecay,
                               adaptive_iStates=Default$adaptive_iStates,
                               Optimizer=Default$Optimizer,
                               Logic = Default$Logic,
                               Mode = Default$Mode, ModeBack = Default$Mode,
                               MinStepsForward = Default$MinStepsForward,
                               MinStepsBackward = Default$MinStepsBackward,
                               gradClipping=Default$gradClipping,
                               LSTM=Default$LSTM,
                               FixNodes = FixNodes, NameProj = NameProj,
                               PDF=F, GIF=F, Visualize="Output",
                               alpha=Default$alpha, lambda=Default$lambda)
    }

    toc=Sys.time()-tic
    print(toc)

  } else {
    # to check in to S3
    ResEss<-ess(Quant=3, Init, MUT, y, CGS, NETall1, Default, FixNodes, Plot=T)

    net[["NETall"]]<-ResEss$Best
    net$history[["NETallList"]]<-list(ResEss$Best)
    net$history[["NETallActivity"]]<-NULL # to check
    net$history[["Cost"]]<-ResEss$Cost

  }

  ############ save
  print("update parameters")

  net[["Parameters"]]$Default<-Default
  net[["Parameters"]]$Boundaries<-Boundaries
  net[["TrainSplit"]]$Train<-Train
  net[["TrainSplit"]]$Val<-Val
  net$Data$Surv<-DATA$SurvData


  Latt<-length(list.files(DIR, pattern =  NameProj))
  if(Latt>0){
    CHAR<-as.numeric(gsub(".*_","", gsub(".Rdata","", gsub(NameProj ,"", list.files(DIR, pattern =  NameProj)))))
    Latt<-max(CHAR[nchar(CHAR)<=2],na.rm = T)
  }

  print("save")

  DIRSave<-paste(DIR,"/", NameProj, Latt+1,".Rdata",sep = "") # keep it for saving validation
  save(net, file = DIRSave)

  return(net)
}

#net<-RunTCGAopt(Param="", DIR=file.path(getwd(),"model"),
#                       NameProj="HNSCC_AMoNet", GENESman="MTOR",
#                       treatmt="", SelectMECA="HALLMARK", organ="hnsc",KeepData = T,
#                       eSS=F, NewNet=T, MinConnect=1, nblayers=4, alpha=0)


#' Plots and predictions within AMoNet grid search pipeline
#'
#' @param net *AMoNet* object.
#' @param DIR path to file to plot. Default in tmp/
#'
#' @return
#' Metriccs for training and validation
#'
PlotAndPredict<-function(net, DIR=file.path(getwd(),"tmp")){

  pdf(paste(DIR, NameProj,Latt+1,".pdf",sep = ""))

  ## visualize learning phase
  plot(net$history)

  ## relation of phenotype weights to outputs

  NETall1<-net$history$NETallList[[length(net$history$NETallList)]]
  par(mar=c(5, 4, 4, 2) + 0.1)
  #par(mar=c(15, 4, 4, 2) + 0.1)
  barplot(NETall1[NETall1$Output,"Weights"],names.arg =NETall1[NETall1$Output,1], las=2, cex.names = 0.5)

  # predict

  #predict(net)

  net<-predict(net, newiStates = net$iStates,
          newInit = NULL, newMUT = net$Data$MUT,
          newtreatmt = NULL, newy = net$Data$y)

  plot(predsAMoNet, xlim=c(0,1))
  #net$Predict$metrics$CindexTrain
  #plot(net,Npat = c(1,2),PDF = F,LEGEND = T, ylim=c(0,1), Species = "Output", col=c(2,3))
  if(is.null(names(net$TotAttractors[,1,1,1]))){
    names(net$TotAttractors[,1,1,1])<-rownames(net$iStates)
  }

  # Train & Val split
  Train<-net$TrainSplit$Train; Val<-net$TrainSplit$Val

  predsTrain<-predict(net, newiStates = net$iStates[Train,],
                       newInit = NULL, newMUT = net$Data$MUT[Train],
                       newtreatmt = NULL, newy = net$Data$y[,Train],
                       RETURN = T)

  plot(predsTrain, xlim=c(0,1))

  predsVal<-predict(net, newiStates = net$iStates[Val,],
                      newInit = NULL, newMUT = net$Data$MUT[Val],
                      newtreatmt = NULL, newy = net$Data$y[,Val],
                      RETURN = T)

  plot(predsVal, xlim=c(0,1))
  ###

  return(list(TrainMetrics=predsTrain$metrics,ValMetrics=predsVal$metrics))
}
