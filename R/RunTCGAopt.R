#' Wrapper function to run the AMonet workflow, including hyper-parameters grid search, on TCGA data
#'
#' @description AMonet workflow comprises building, training and simulating a network model to predict survival of patients from TCGA genomics WES data.
#'
#' @param Param character vector: the hyper-parameter(s) to randomly test whithin this run.
#' @param DIR directory: directory where the models are stored
#' @param NameProj character: name of your project. If at least one model with the name of your project NameProj is in the DIR directory, the script will load it and use it.
#' @param NewNet boolean: when set to FALSE, if at least one model with the name of your project NameProj is in the DIR directory, the script will load it and use it. If set to TRUE, a new molecular network will be built using arguments GENESman and default parameters (Interval, MinConnect, nblayers).
#' @param GENESman character vector: gene(s)' name(s) (in Hugo nomenclature) used to build the AMoNET network. If "all", all genes related to SelectMECA (gene sets) will be selected. If "frequent", 95\% most frequently altered genes from the cancer gene census (CGS) list in the cohort will be selected.
#' @param treatmt character vector: gene(s)' name(s) targeted by simulated treatment(s). Added to GENEman toused to build the AMoNET network.
#' @param SelectMECA character in regexp format: selection of one or several gene sets related to biological mecanisms. Names gene sets are available with the command \code{print(names_MECA)}.
#' @param organ character in regexp format: selection of one or several organs or cancer types. Names of organs are available with the command \code{print(names_MECA)}
#' @param eSS boolean. If \code{TRUE}, performs multistart without training. Default is set to \code{FALSE}. Recommandation is to use \code{eSS=TRUE} within a hyper-parameters search and in case of hyper-parameters for building the net such as: \code{Param = "MeanWinit"}, \code{Param = "SdWinit"}, \code{Param = "nblayers"} or \code{Param = "MinConnect"}.
#' @param Default list. Default hyper-parameters used for the workflow (building and training).
#' @param Boundaries list. Boundaries for Default hyper-parameters for random selection.
#'
#' @details
#' This function can be run:
#' 1) Alone to train an AMoNet model with fixed hyper-parameters are either user-defined or the default one are used. For default hyper-parameters values call \code{print(Default)}.
#' 2) In array and iteratively from command line together with the \code{PlotAndPredict} function to perform hyper-parameters search. Hyperparameter selected in the \code{Param} argument are thus randomly selected in batches and best results selected.
#'
#' Requirements:
#' Uses the pre-configured /model and /tmp paths in AMoNet package directories (checks todo)
#'
#' To set manually any hyper-paramters, change it into the \code{print(Default)} object.
#'
#' Setting \code{eSS=TRUE} is a convinient way to select best networks initialisations parameters. It is not recommended for learning parameters.
#'
#' @return todo AMoNet object
#' @export
RunTCGAopt<-function(Param=c("nblayers", "MinConnect"), DIR=getwd(),
                     NameProj="LUNG_AMoNet_ess", GENESman=c("EGFR","MTOR"),
                     treatmt=NULL, SelectMECA="HALLMARK", organ="lung",
                     eSS=F, NewNet=T, KeepData=T, PartitionSplit=0.7,
                     Default=AMoNet::Default, Boundaries=AMoNet::Boundaries){

  set.seed(NULL)

  NameProjbase<-NameProj

  ###############
  # initiate AMoNet
  if(is.null(GENESman)&is.null(treatmt)){
    print("You should select a list of genes in initial function to build the net")
    stop()
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
  # initialize the net
  net<-AMoNet(GENESman = GENESman, treatmt = treatmt)
  net$Parameters$Default<-Default
  net$Parameters$Boundaries<-Boundaries

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

  # check and set directories
  if(length(grep("/model",DIR))==0){
    DIR<-file.path(DIR,"model/")
  }

  if(!dir.exists(DIR)){
    dir.create(DIR,showWarnings = T)
    print(paste("New direction to store networks :", DIR))
  }

  ###############
  # load networks or build

  FILES<-list.files(DIR, pattern = NameProjbase)

  if(NewNet){
    FILES<-NULL
  }

  if(length(FILES)>1){
    # ok

    net<-SelectGridSearch(NameProjbase = NameProjbase,
                              DIR=DIR, Default=net$Parameters$Default,
                              ValSelect=T)

    # update Hyperparameters with new ones
    net$Parameters$Default<-HyperP(C=Param, Default = net$Parameters$Default,
                                   Boundaries = net$Parameters$Boundaries)


  } else if(length(FILES)==1){
    # change with preBuild function
    # update Hyperparameters with new ones
    x<-load(file.path(DIR,FILES))
    net<-get(x)
    net$Parameters$Default<-HyperP(C=Param, Default = net$Parameters$Default,
                                   Boundaries = net$Parameters$Boundaries)

  } else { # this is used for the first net building

    # update Hyperparameters with new ones
    net$Parameters$Default<-HyperP(C=Param,Default = net$Parameters$Default,
                                   Boundaries = net$Parameters$Boundaries)

    net<-build.AMoNet(object = net,
                      InteractionBase = OMNI,
                      nblayers = net$Parameters$Default$nblayers,
                      MinConnect = net$Parameters$Default$MinConnect,
                      RestrictedBuilding = T, RestrictedBase = T, FilterCGS = T,
                      MeanWinit = net$Parameters$Default$MeanWinit,
                      SdWinit = net$Parameters$Default$SdWinit,
                      Phenotypes=GenesSelec, MECA=MECA,
                      Interval = net$Parameters$Default$Interval,
                      LSTM = net$Parameters$Default$LSTM,
                      Optimizer = net$Parameters$Default$Optimizer,
                     KeepPhenotypes=T, no_cores = net$Parameters$Default$no_cores,
                     NameProj = NameProjbase)
  }

  Species<-union(net$NETall$source_hgnc,net$NETall$target_hgnc)

  ###############
  print("load TCGA data")
  # run the function to load data

  net<-LoadCleanTCGA(net, Species = Species,
                      Param = Param, RestrictUnique = F, organ = organ)

  #DATA$Perturb<-MutMatToList(MUTa = DATA$Perturb) # to list
  #lapply(DATA, dim)
  #list2env(DATA,envir = .GlobalEnv)

  # rename project
  QUANTS<-unlist(sapply(Param,function(PR){
    PARAM<-eval(parse(text=PR) ,envir = Default)
    NAME<-paste(PR,ifelse(is.na(as.numeric(PARAM)),PARAM,as.numeric(PARAM)),sep = "_") #round(as.numeric(PARAM),4)
    return(NAME)
  }))

  net$call$NameProj<-paste(NameProjbase, paste(QUANTS,collapse = "_"),"_", sep = "")
  NameProj<-net$call$NameProj
  print(paste("Welcome to project", NameProj))

  ###############
  # compute initial states
  net<-RandomiStates(net)

  #iStates<-matrix(runif(ncol(y)*length(Species)),nrow = ncol(y), ncol = length(Species))
  #rownames(iStates)<-colnames(y)
  #colnames(iStates)<-Species

  ##### EXPi # to improve by using it as constrained activity value
  if(net$Parameters$Default$EXPinit){
    for(i in colnames(y)){
      net$iStates[i,colnames(net$Data$Init)]<-as.numeric(net$Data$Init[i,]) # match the col order
    }
  }

  if(net$Parameters$Default$LSTM){ # to check if doesnt impare LSTM
    Ct<-t(replicate(ncol(y), rep(1,length(Species)) , simplify = TRUE))
    colnames(Ct)<-Species
    rownames(Ct)<-colnames(y)
  } else{
    Ct=NULL
  }

  ####################
  print("run learning")

  if(net$Parameters$Default$FixWeights){
    FixNodes<-Species[!Species%in%c(net$NETall$source_hgnc[net$NETall$target_hgnc%in%"Output"],"Output")]
  } else {
    FixNodes<-NULL
  }

  if(!"TrainSplit"%in%names(net)){
    net<-split(net,PartitionSplit)
    #list2env(net$TrainSplit, envir = globalenv())
  } else {
    #list2env(net$TrainSplit, envir = globalenv())

    net$TrainSplit$Train<-net$TrainSplit$Train[net$TrainSplit$Train%in%colnames(net$Data$y)]
    net$TrainSplit$Val<-net$TrainSplit$Val[net$TrainSplit$Val%in%colnames(net$Data$y)]

    MorePat<-colnames(net$Data$y)[!colnames(net$Data$y)%in%c(net$TrainSplit$Train,net$TrainSplit$Val)]
    if(length(MorePat)>0){

      MPTrain<-sample(MorePat,size = round(length(MorePat)*PartitionSplit))
      MPVal<-MorePat[!MorePat%in%MPTrain]

      net$TrainSplit$Train<-c(net$TrainSplit$Train,MPTrain)
      net$TrainSplit$Val<-c(net$TrainSplit$Val,MPVal)
#      net$TrainSplit$Train<<-Train
#      net$TrainSplit$Val<<-Val
    }
  }

  # cap minibatches
  net$Parameters$Default$MiniBatch<-min(net$Parameters$Default$MiniBatch,length(net$TrainSplit$Train))

  if(net$Parameters$Default$MiniBatch==length(net$TrainSplit$Train)){net$Parameters$Default$Optimizer=NULL}

  print(paste("Training on", length(net$TrainSplit$Train),"; validation on",length(net$TrainSplit$Val), "patients"))

  if(!eSS){

    tic=Sys.time()

    net<-train.AMoNet(net, y= net$Data$y[,net$TrainSplit$Train,drop=F],
           MUT= net$Data$MUT[net$TrainSplit$Train], treatmt=NULL,
           Init = net$Data$Init[,net$TrainSplit$Train],
           ValMut = net$Parameters$Default$ValMut,
           iStates= net$iStates[net$TrainSplit$Train,],
           Ct=Ct[net$TrainSplit$Train,],
           MiniBatch= net$Parameters$Default$MiniBatch,
           alpha = net$Parameters$Default$alpha, lambda = net$Parameters$Default$lambda,
           iteration= net$Parameters$Default$iteration,
           beta1= net$Parameters$Default$beta1, beta2= net$Parameters$Default$beta2,
           Parallelize= net$Parameters$Default$Parallelize, no_cores = net$Parameters$Default$no_cores,
           learning_rate= net$Parameters$Default$learningrate,
           LearningRateDecay = net$Parameters$Default$LearningRateDecay,
           adaptive_iStates= net$Parameters$Default$adaptive_iStates,
           Optimizer= net$Parameters$Default$Optimizer,
           Logic = net$Parameters$Default$Logic,
           Mode = net$Parameters$Default$Mode, ModeBack = net$Parameters$Default$Mode,
           MinStepsForward = net$Parameters$Default$MinStepsForward,
           MinStepsBackward = net$Parameters$Default$MinStepsBackward,
           gradClipping= net$Parameters$Default$gradClipping,
           LSTM= net$Parameters$Default$LSTM,
           FixNodes = FixNodes, NameProj = net$call$NameProj,
           PDF=F, GIF=F, Visualize=NULL, KeepData=F,  KeepTraining=T)

    toc=Sys.time()-tic
    print(toc)

  } else {
    # to check in to S3
    net<-ess(net, y= net$Data$y[,net$TrainSplit$Train,drop=F],
                MUT= net$Data$MUT[net$TrainSplit$Train],
                Init = net$Data$Init[,net$TrainSplit$Train],
                iStates= net$iStates[net$TrainSplit$Train,],
                treatmt=NULL) # CGS, NETall1, Default,ValMut = net$Parameters$Default$ValMut,FixNodes=FixNodes, Plot=F
      # Quant=3, no_cores=net$Parameters$Default$no_cores
  }

  ############ save
  print("update parameters")

  Latt<-length(list.files(DIR, pattern =  NameProj))
  if(Latt>0){
    CHAR<-as.numeric(gsub(".*_","", gsub(".Rdata","", gsub(NameProj ,"", list.files(DIR, pattern =  NameProj)))))
    Latt<-max(CHAR[nchar(CHAR)<=2],na.rm = T)
  }

  print("save")

  DIRSave<-paste(DIR,"/", NameProj, Latt+1,".Rdata",sep = "") # keep it for saving validation
  net$DIRSave<-DIRSave
  save(net, file = DIRSave)

  return(net)
}

#' Plots and predictions within AMoNet grid search pipeline
#'
#' @param net *AMoNet* object, after training and with split defined (net$TrainSplit)
#'
#' @return
#' Metrics for training and validation in \code{$predict_[]} from *AMoNet* object
#' Plots of training procedures and predictions in the tmp/ directory
#' @export
PlotAndPredict<-function(net){

  DIR<-gsub("Rdata","pdf", gsub("/model/","/tmp/", net$DIRSave))
  if(!dir.exists(gsub("tmp/.*","tmp/", DIR))){
    dir.create(gsub("tmp/.*","tmp/", DIR),showWarnings = T)
    print(paste("Learning plots are in",DIR))
  }

#  if( dir.exists(gsub("tmp/.*","tmp/", DIR)) ){


#    Latt<-length(list.files(DIR,pattern = net$call$NameProj))
#    if(Latt>0){
#      Latt<-max(as.numeric(gsub(".pdf","", gsub(net$call$NameProj,"", list.files(DIR,pattern = net$call$NameProj)))))
#    }

#    NAME<-paste(DIR,net$call$NameProj, Latt+1,".pdf",sep = "") #ifelse(Optimized,"Opt","Init"),
#    pdf(NAME)

    pdf(DIR)
#    }
  ## visualize learning phase
    if(net$call$train_call$iteration>1){
      plot(net$history)
    }

  ## relation of phenotype weights to outputs

  NETall1<-net$history$NETallList[[length(net$history$NETallList)]]
  par(mar=c(5, 4, 4, 2) + 0.1)
  #par(mar=c(15, 4, 4, 2) + 0.1)
  barplot(NETall1[NETall1$Output,"Weights"],names.arg =NETall1[NETall1$Output,1], las=2, cex.names = 0.5,
          main="Relation of the phenotype layer with outputs")

  # predict whole base with split
    net<-predict(net, SplitType = "Train")

    plot(net$Predict_Train, xlim=c(0,1), main = paste("Training predictions for \n",net$call$NameProj))

    net<-predict(net, SplitType = "Val")

    plot(net$Predict_Val, xlim=c(0,1), main = paste("Validation predictions for \n",net$call$NameProj))

    if(!is.null(DIR)){
      dev.off()
    }

    if(FALSE){

    if(is.null(names(predsAll$TotAttractors[,1,1,1]))){
      names(predsAll$TotAttractors[,1,1,1])<-rownames(predsAll$iStates)
    }
    predsAll$Predict_Train$metrics$Cindex


  #net$Predict$metrics$CindexTrain
  #plot(net,Npat = c(1,2),PDF = F,LEGEND = T, ylim=c(0,1), Species = "Output", col=c(2,3))

  # Train & Val split
  Train<-net$TrainSplit$Train; Val<-net$TrainSplit$Val

  predsTrain<-predict(net, newiStates = net$iStates[net$TrainSplit$Train,],
                       newInit = NULL, newMUT = net$Data$MUT[net$TrainSplit$Train],
                       newtreatmt = NULL, newy = net$Data$y[,net$TrainSplit$Train],
                      Logic = net$Parameters$Default$Logic, Mode = net$Parameters$Default$Mode,
                      MinSteps = net$Parameters$Default$MinStepsForward,LSTM = net$Parameters$Default$LSTM,
                      Parallelize = net$Parameters$Default$Parallelize, no_cores = net$Parameters$Default$no_cores,
                      ValMut = net$Parameters$Default$ValMut)
  plot(predsTrain$Predict)

  predsVal<-predict(net, newiStates = net$iStates[net$TrainSplit$Val,],
                    newInit = NULL, newMUT = net$Data$MUT[net$TrainSplit$Val],
                    newtreatmt = NULL, newy = net$Data$y[,net$TrainSplit$Val],
                    Logic = net$Parameters$Default$Logic, Mode = net$Parameters$Default$Mode,
                    MinSteps = net$Parameters$Default$MinStepsForward,LSTM = net$Parameters$Default$LSTM,
                    Parallelize = net$Parameters$Default$Parallelize, no_cores = net$Parameters$Default$no_cores,
                    ValMut = net$Parameters$Default$ValMut)

  plot(predsVal$Predict, xlim=c(0,1))
  predsVal$Predict$metrics
  ###

  return(list(TrainMetrics=predsTrain$metrics,ValMetrics=predsVal$metrics))
    }

    print("save")
    save(net, file = net$DIRSave)

  return(net)
}
