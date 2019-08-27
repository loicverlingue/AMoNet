
#' Wrapper function to run the AMonet workflow on TCGA data
#'
#' AMonet workflow comprise building, training and simulating a network model to predict survival of patients from TCGA genomics WES data.
#'
#' @param Param character vector: the hyperparameters to randomly test
#' @param DIR directory: directory of the models storing
#' @param NameProj character: name of your project. If at least one model with the name of your project NameProj is in the DIR directory, the script will load it and use it.
#' @param NewNet boolean: when set to FALSE, if at least one model with the name of your project NameProj is in the DIR directory, the script will load it and use it. If set to TRUE, a new molecular network will be built using arguments GENESman and default parameters (Interval, MinConnect, nblayers).
#' @param GENESman character vector: gene(s)' name(s) (in Hugo nomenclature) used to build the AMoNET network. If "all", all genes related to SelectMECA (gene sets) will be selected. If "frequent", 95\% most frequently altered genes from the cancer gene census (CGS) list in the cohort will be selected.
#' @param treatmt character vector: gene(s)' name(s) targeted by simulated treatment(s). Added to GENEman toused to build the AMoNET network.
#' @param SelectMECA character in regexp format: selection of one or several gene sets related to biological mecanisms. Names gene sets are available with the command \code{print(names_MECA)}.
#' @param organ character in regexp format: selection of one or several organs or cancer types. Names of organs are available with the command \code{print(names_MECA)}
#' @param eSS boolean to perform multistart or not. Default is set to \code{FALSE}. Recommandation is to use \code{eSS=TRUE} within a hyperparameter search and in case of \code{Param = "MeanWinit"} and/or \code{Param = "SdWinit"} eSS \code{TRUE}.
#' @param ... any argument(s) present in Default can be set here. \code{print(Default)} for information on Default hyper-parameters
#'
#' @examples
#' NETallProp<-RunTCGAopt(Param="", DIR=file.path(getwd(),"model"),
#'                       NameProj="HallmarksLungHN", GENESman=c("EGFR", "KRAS", "TP53", "MTOR"),
#'                       treatmt="", SelectMECA="HALLMARK", organ="luad",
#'                       eSS=F, NewNet=T, MinConnect=6, alpha=0)
#'
#' @export
#'
RunTCGAopt<-function(Param=c("nblayers", "MinConnect"), DIR=file.path(getwd(),"model"),
         NameProj="HallmarksLungHN", GENESman=c("EGFR", "KRAS", "TP53", "MTOR"),
         treatmt="", SelectMECA="HALLMARK", organ="luad",
         eSS=F, NewNet=T, ...){

  set.seed(NULL)

  if(FALSE){
    # if command line use:  CommandArgs()
    ######################
    #args<-"--Param nblayers MinConnect --DIR /data/tmp/lverling/PossiLandscape/ --NameProj HallmarksLungHN --GENESman EGFR KRAS TP53 MTOR --treatmt --Interval 10 --SelectMECA HALLMARK --organ luad"


    # retreive the hyperparameters to search from the commandline
    args <- commandArgs(trailingOnly = TRUE)

    #args <- gsub("\r","",args)

    hh <- paste(unlist(args),collapse=' ')
    listoptions <- unlist(strsplit(hh,'--'))[-1]
    options.args <- sapply(listoptions,function(x){
      unlist(strsplit(x, ' '))[-1]
    })
    options.names <- sapply(listoptions,function(x){
      option <-  unlist(strsplit(x, ' '))[1]
    })
    names(options.args) <- unlist(options.names)

    options.args<-lapply(options.args,function(Correct){
      if(length(Correct)==0){
        return("")
      } else if(all(!is.na(as.numeric(Correct)))){
        return(as.numeric(Correct))
      } else if(exists(Correct)){
        return(eval(parse(text = Correct)))
      } else {
        return(Correct)
      }
    })

    print("options are:")
    print(options.args)

    ##### other param
    # eSS : dont run the learning to directly test the multi starting points
    eSS<-ifelse("MeanWinit"%in%unlist(options.args),T,F)
    LSTM<-F
    Logic<-"Sigmoid"
    Optimizer="Momentum"
    #Adam<-F
    NameProj="HallmarksLonger";NameProjbase<-NameProj
    organ="lung"
    organ="lusc|luad|hnsc"
    SelectMECA<-""
    SelectMECA<-"HALLMARK"
    GENESman="all"
    GENESman=c("KRAS", "PDCD1","CTLA4","CD274","EGFR")
    GENESman=c("KRAS", "EGFR","TP53")
    treatmt=""
    Interval=10
    NewNet=F
    no_cores=8

    # if from command line:
    if(exists("options.args")){
      list2env(options.args,globalenv())
    }
  }
if(FALSE){

  Param=c("nblayers", "MinConnect")
  DIR=file.path(getwd(),"model")
  NewNet=T
  NameProj="HallmarksLungHN"
  GENESman=c("EGFR", "KRAS", "TP53", "MTOR")
  treatmt=""
  Interval=10
  SelectMECA="HALLMARK"
  organ="luad"
  no_cores=4
  eSS=F
  LSTM=F
  Logic="Sigmoid"
  Optimizer="Momentum"
  Interval=10
  NewNet=F

}

  NameProjbase<-NameProj
  # package should be loaded with AMoNet package
  # todo
  if(FALSE){
    print("Loading packages")
    library(reshape2)
    library(parallel)
    library(stringr)
    library(splines)
    library(cgdsr)
    library(survival)
    library(magick)
    library(matrixStats)
    library(survivalROC)
    library(pec)
  }

if(FALSE){
#setwd(paste(getwd(),"/R/",sep = ""))
print("source functions")

source("R/NetBuilding.R")
source("R/addWeights.R")
source("R/NetAnalyser170418.R")

#source("R/PlotAttractors240718.R")
source("R/plotProbaState141118.R")
source("R/BuildNetMat120518.R")
source("R/InitStates.R")
source("R/Outputs_function_130718.R")
source("R/SelectGTExOptNet.R")

source("R/PlotOptNet190718.R")
#source("DownloadGTEx.R")
#source("ProcessGTEx180718.R")
#source("LoadTCGAdata040918.R")
source("R/LoadcBioportal150119.R")
#source("LoadRecountClinicTCGA.R")
source("R/NETforOpt.R")

#source("R/BoolSimul170319.R")
source("R/BoolSimul200519.R")
source("R/MutMatToList.R")
source("R/Back_BoolSimul290119.R")
#source("R/RunBackSimul_141118.R")
source("R/RunBackSimul_200519.R")
source("R/SurvEstimator261018.R")
source("R/Vizlearning.R")

source("R/essAlone.R")

#setwd("/data/tmp/lverling/")
#setwd("C:/Users/L_VERLINGUE/Desktop/ModelK/Rpack/ArtMolNet/")
}

######
# default hyperparameters
if(FALSE){
  Default<-list(learningrate=0.003, MeanWinit=0.001, SdWinit = 0.001,
              MiniBatch=2^11, SimAnnealMaxSd=0.001,
              beta1=0.9,beta2=0.999,gradClipping=NULL,
              iteration=5, LearningRateDecay=NULL,
              Parallelize=T, MinStepsBackward=1,
              MinStepsForward=5, adaptive_iStates=T,
              MUTinit=T, EXPinit=F, ValMut=50,
              Logic="Sigmoid", Mode="LAYER",FixWeights=F,
              Interval=10, Npat=1e4, Optimizer="Momentum",
              MinConnect=6, nblayers=2, alpha=1, lambda=1e-2,
              LSTM=F, no_cores=4, organ="lung")

Bondaries<-list(learningrate=c(1e-5,0.05), MeanWinit=c(10^-4,0.99), SdWinit = c(10^-4,0.99),
                MiniBatch=as.integer(c(3,11)), SimAnnealMaxSd=c(10e-4,0.1),
                beta1=c(0.1,0.99), beta2=c(0.99,0.9999999),gradClipping=as.integer(c(10,100)),
                iteration=as.integer(c(30,80)), LearningRateDecay=c(NULL,"linear","exponential"),
                Parallelize=c(TRUE,FALSE), MinStepsBackward=as.integer(c(1,10)),
                MinStepsForward=as.integer(c(2,25)), adaptive_iStates=c(TRUE,FALSE),
                MUTinit=c(TRUE,FALSE),EXPinit=c(TRUE,FALSE),
                ValMut=as.integer(c(2,500)), Mode=c("ASYNC","LAYER"),
                FixWeights=c(TRUE,FALSE), Interval=as.integer(c(5,50)),
                Npat=as.integer(c(10,1e4)), MinConnect=as.integer(c(1,5)),
                nblayers=as.integer(c(1,10)),
                alpha=as.integer(c(1,2)), lambda=c(1e-5,0.9),
                Optimizer=c(NULL,"Adam","Momentum","RMSprop"),
                Logic=c("Boolean","Sigmoid","tanh", "ReLU"),
                Mode = c("LAYER","ASYNC","SYNC"),LSTM=c(T,F),
                no_cores=as.integer(c(2,20)) )

#save(Default,file = file.path(getwd(),"data/Default.RData"))
}

HyperP<-function(C=Param, Default=Default, Bondaries=Bondaries){

  # define distribution you want
  LogDistrib<-function(min,max,N){
    sample(exp(seq(log(abs(min)),log(abs(max)),length.out = 100)),N)
  }

  NormDistrib<-function(min,max,N){
    sample(seq(min,max,length.out = 100),N)
  }

  for(Cp in C){
    if(is.numeric(Bondaries[[Cp]])){
      if(is.integer(Bondaries[[Cp]])){
        Default[[Cp]]<- sample(seq(Bondaries[[Cp]][1],Bondaries[[Cp]][2]),1)
      } else if(all(Bondaries[[Cp]]>=0.1)){
        Default[[Cp]]<- NormDistrib(min = Bondaries[[Cp]][1],max = Bondaries[[Cp]][2], N = 1)
      }else{
        Default[[Cp]]<- LogDistrib(min = Bondaries[[Cp]][1],max = Bondaries[[Cp]][2], N = 1)
      }
    } else if(is.character(Bondaries[[Cp]])|is.logical(Bondaries[[Cp]])){
      Default[[Cp]]<- sample(Bondaries[[Cp]],1)
    }
  }

  if("MiniBatch"%in%C){
    Default$MiniBatch<-sample(2^seq(Bondaries$MiniBatch[1],Bondaries$MiniBatch[2]),1)
  }

  list2env(Default, envir = .GlobalEnv)
  return(Default)
}

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

########
#
GenesSelec<-rbind(GenesSelecImuno,GenesSelecHall)
MECA<-unique(GenesSelec$target_hgnc) # do an interactcive function to choose MECA
if(!is.null("SelectMECA")){
  MECA<-grep(SelectMECA,MECA,value = T)
}
if("all"%in%tolower(GENESman)){
  GENESman<-unique(GenesSelec$source_hgnc[GenesSelec$target_hgnc%in%MECA])
}


###############
# load networks

if(!dir.exists(DIR)){
  dir.create(DIR,showWarnings = T)
  print("New direction -> new nets")
}

FILES<-list.files(DIR,pattern = NameProjbase)

if(NewNet){
  FILES<-NULL
}

if(length(FILES)>1){

  # for first run of TCGA grid search when using pre-optimized nets from GTEx : match when includes "MeanWinit"
  NETdata<-SelectGTExOptNet(NameProjbase = NameProjbase, NETall = NULL,
                            DIR=DIR, Default=Default, Bondaries=Bondaries,
                            addOutputs = NULL, MiniBatch=Default$MiniBatch,
                            ValSelect=T)

  NETall1<-NETdata$NETall # the best net
  Default<-NETdata$Default # the best default parameters from previous optimisations
  Bondaries<-NETdata$Bondaries
  Train<-NETdata$NETallProp$TrainSplit$Train
  Val<-NETdata$NETallProp$TrainSplit$Val
  organ<-NETdata$NETallProp$organ

  # update Hyperparameters with new ones
  Default<-HyperP(C=Param, Default = Default, Bondaries = Bondaries)

  # check outputs are in the good format
  if(Default$Interval>1&FALSE){
    if(!all(as.numeric(unique(gsub("Output","",NETdata$NETall[NETdata$NETall$Output,3])))==seq(Default$Interval))){
      NETdata<-SelectGTExOptNet(NameProjbase = NameProjbase, NETall = NULL,
                                dir=DIR, Default=Default, Bondaries=Bondaries,
                                addOutputs = Default$Interval, MiniBatch=Default$MiniBatch)
    }
  }

  # initialize colors of the net
   NETall1$interaction_directed_signed<-ifelse(NETall1$Weights>0,"ACTIVATE","INHIBIT")


} else if(length(FILES)==1){

  # update Hyperparameters with new ones
  Default<-HyperP(C=Param,Default = Default,Bondaries = Bondaries)

  if(length(grep("csv",FILES))>0){
    NETall1<-read.csv2(file = paste(DIR,FILES[1],sep = ""), row.names = 1, header = T,stringsAsFactors = F)
    NETall1$X<-NULL
    NETalldata<-NETforOpt(WRITE = F,NETall = NULL, NameProj=NameProjbase,
                        GENESman=unique(NETall1$source_hgnc),
                        treatmt="", InteractionBase = NETall1, RestrictedBuilding = F,
                        MeanWinit = Default$MeanWinit, SdWinit = Default$SdWinit,
                        addOutputs = Default$Interval, FilterCGS = F,
                        LSTM = Default$LSTM, Adam = !is.null(Default$Optimizer),
                        Activity = F, KeepPhenotypes = T,
                        Phenotypes =  NETall1[NETall1$target_hgnc%in%NETall1$source_hgnc[NETall1$Output],1:3],
                        MECA = NULL,
                        nblayers = Default$nblayers, MinConnect = Default$MinConnect,
                        RestrictedBase = F)

  NETall1<-NETalldata$NETall
  } else {
    x<-load(file.path(DIR,FILES))
    NETallProp<-get(x)
    NETall1<-NETallProp$NETallList[[length(NETallProp$NETallList)]]
    Default<-NETallProp$Parameters$Default # the best default parameters from previous optimisations
    Bondaries<-NETallProp$Parameters$Bondaries
    Train<-NETallProp$TrainSplit$Train
    Val<-NETallProp$TrainSplit$Val
    organ<-NETallProp$organ
  }

} else { # this is used for the first net building

  if(is.null(GENESman)&is.null(treatmt)){
    print("You should select a list of genes in initial function to build the net")
    stop()
  }

  # update Hyperparameters with new ones
  Default<-HyperP(C=Param,Default = Default,Bondaries = Bondaries)

  #Default$MinConnect=6
  NETalldata<-NETforOpt(WRITE = F, NameProj=NameProjbase,
                        GENESman=GENESman,treatmt=treatmt,
                        InteractionBase = OMNI,
                        MinConnect = Default$MinConnect, RestrictedBuilding = T,
                        MeanWinit = Default$MeanWinit, SdWinit = Default$SdWinit,
                        nblayers = Default$nblayers, addOutputs = Default$Interval,
                      FilterCGS = F, LSTM = Default$LSTM, Adam = !is.null(Default$Optimizer),
                        Activity = F,  KeepPhenotypes=T,
                      Phenotypes=GenesSelec, #[GenesSelec$target_hgnc%in%MECA,],
                      MECA=MECA, RestrictedBase = T, no_cores=no_cores)

   NETall1<-NETalldata$NETall
}

Species<-union(NETall1$source_hgnc,NETall1$target_hgnc)


###############
print("load TCGA data")
# run the function to load data
DATA<-LoadCleanTCGA(NETall = NETall1, Default = Default, Bondaries = Bondaries,
                Param = Param, RestrictUnique = T)
#lapply(DATA, dim)
list2env(DATA,envir = globalenv())
MUT<-MutMatToList(MUTa = Perturb) # list

# rename project
QUANTS<-unlist(sapply(Param,function(PR){
  PARAM<-eval(parse(text=PR) ,envir = .GlobalEnv)
  NAME<-paste(PR,ifelse(is.na(as.numeric(PARAM)),PARAM,as.numeric(PARAM)),sep = "_") #round(as.numeric(PARAM),4)
  return(NAME)
}))

NameProj<-paste(NameProjbase, paste(QUANTS,collapse = "_"),"_", sep = "")
print(paste("Welcome to project", NameProj))

###############
# compute initial states
iStates<-matrix(0.5,nrow = ncol(Target), ncol = length(Species))
rownames(iStates)<-colnames(Target)
colnames(iStates)<-Species

##### EXPi # to improve by using it as contrained activity value
if(Default$EXPinit){
  for(i in colnames(Target)){
    iStates[i,colnames(Init)]<-as.numeric(Init[i,]) # match the col order
  }
}

if(Default$LSTM){ # to check if doesn impare LSTM
  Ct<-t(replicate(ncol(Target), rep(1,length(Species)) , simplify = TRUE))
  colnames(Ct)<-Species
  rownames(Ct)<-colnames(Target)
} else{
  Ct=NULL
}

####################
print("run learning")

if(Default$FixWeights){
  FixNodes<-Species[!Species%in%c(NETall1$source_hgnc[NETall1$target_hgnc%in%"Output"],"Output")]
} else {
  FixNodes<-NULL
}

if(!exists("Train")|!exists("Val")){

  Train<-sample(colnames(Target),ncol(Target)*0.7)
  Val<-colnames(Target)[!colnames(Target)%in%Train]

} else if(is.null(Train)|is.null(Val)){

  Train<-sample(colnames(Target),ncol(Target)*0.7)
  Val<-colnames(Target)[!colnames(Target)%in%Train]

} else {

  Train<-Train[Train%in%colnames(Target)]
  Val<-Val[Val%in%colnames(Target)]

  MorePat<-colnames(Target)[!colnames(Target)%in%c(Train,Val)]
  if(length(MorePat)>0){

    MPTrain<-sample(MorePat,size = round(length(MorePat)*0.7) )
    MPVal<-MorePat[!MorePat%in%MPTrain]

    Train<-c(Train,MPTrain)
    Val<-c(Val,MPVal)

  }
}

# cap minibatches
Default$MiniBatch<-min(Default$MiniBatch,length(Train))
if(Default$MiniBatch==length(Train)){Default$Optimizer=NULL}

print(paste("Training on",length(Train),"; validation on",length(Val), "patients"))

if(!eSS){

    tic=Sys.time()
    NETallProp<-RunBackSimul(NETall=NETall1, y=Target[,Train,drop=F],
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
                             Discretize = F, LSTM=Default$LSTM,
                             FixNodes = FixNodes, NameProj = NameProj,
                             PDF=F, GIF=F, Visualize="Output",
                             alpha=Default$alpha, lambda=Default$lambda)

    toc=Sys.time()-tic
    print(toc)

} else {

  ResEss<-ess(Quant=3, Init, Perturb, Target, CGS, NETall1, Default, FixNodes, Plot=T)

  NETallProp<-list()
  NETallProp[["NETall"]]<-ResEss$Best
  NETallProp[["NETallList"]]<-list(ResEss$Best)
  NETallProp[["NETallActivity"]]<-NULL # ot check
  NETallProp[["Cost"]]<-ResEss$Cost

}

############ save
print("save")

NETallProp[["Parameters"]]$Default<-Default
NETallProp[["Parameters"]]$Bondaries<-Bondaries
NETallProp[["TrainSplit"]]$Train<-Train
NETallProp[["TrainSplit"]]$Val<-Val
NETallProp[["organ"]]<-organ

Latt<-length(list.files(DIR, pattern =  NameProj))
if(Latt>0){
  CHAR<-as.numeric(gsub(".*_","", gsub(".Rdata","", gsub(NameProj ,"", list.files(DIR, pattern =  NameProj)))))
  Latt<-max(CHAR[nchar(CHAR)<=2],na.rm = T)
}

DIRSave<-paste(DIR, NameProj, Latt+1,".Rdata",sep = "") # keep it for saving validation
save(NETallProp, file = DIRSave)

return(NETallProp)
}


#NETallProp<-RunTCGAopt(Param="", DIR=file.path(getwd(),"model"),
#                       NameProj="BRCA_AMoNet", GENESman="MTOR",
#                       treatmt="", SelectMECA="HALLMARK", organ="brca",
#                       eSS=F, NewNet=T, MinConnect=6, alpha=0)



#####################
PlotAndPredict<-function(NETallProp=NETallProp, Data=TCGAdata, DIR=file.path(getwd(),"tmp")){

  MUTmat<-apply(Data$Perturb,1, function(x){
    ifelse(!is.na(x)&!x%in%"NaN",1,0)
  })

if(!eSS){
  print("do plots")

  # retrieve data
#  iStates<-NETallProp$NETallActivity[[length(NETallProp$NETallActivity)]]
  Default<-NETallProp$Parameters$Default
  list2env(NETallProp$TrainSplit,envir = globalenv())

  # compute initial states
  iStates<-matrix(0.5,nrow = ncol(Target), ncol = length(Species))
  rownames(iStates)<-colnames(Target)
  colnames(iStates)<-Species

  ##### EXPi # to improve by using it as contrained activity value
  if(Default$EXPinit){
    for(i in colnames(Target)){
      iStates[i,colnames(Init)]<-as.numeric(Init[i,]) # match the col order
    }
  }

  if(Default$LSTM){ # to check if doesn impare LSTM
    Ct<-t(replicate(ncol(Target), rep(1,length(Species)) , simplify = TRUE))
    colnames(Ct)<-Species
    rownames(Ct)<-colnames(Target)
  } else{
    Ct=NULL
  }


  # to change
 # DIR<-paste(getwd(),"/tmp/",sep = "")

  pdf(paste(DIR, NameProj,Latt+1,".pdf",sep = ""))

  # load new net
  N<-length(NETallProp$NETallActivity)
  NETall1<-NETallProp$NETallList[[N]]

  PlotOptNet(NETall1,PDF = F,Optimized = T,PrintOptNET = F,LEGEND = F, NameProj = NameProj)

  # mut
  par(mfrow=c(2,1))
  par(mar=c(4,4,2,2)+0.1)
  image(t(MUTmat),col=c(0,1), xlab = "Patients", axes=F, main="Mutations")
  mtext(colnames(Perturb)[ORD1],2,las=2,cex = 0.4, at =  normalized( seq(ncol(Perturb))) )
  barplot(FreqAssoc*100/ncol(MUTmat),las=2,cex.names = 0.4,ylab = "Mutation associations %",space = 0)
  if(length(table(FreqAssoc==1))==2){
    legend("topright",legend = paste("Unique profiles =", round(table(FreqAssoc==1)["TRUE"]*100/ncol(MUTmat)),"%"),
           bty = "n")
  } else {
    legend("topright",legend = "Unique profiles = 100%", bty = "n")
  }

  # learning
  # plot cost
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)+0.1)
  matplot(NETallProp$Cost,type = "l", main='Cost')

  # plot weights
  # iStates<-NETallProp$NETallActivity[[11]]
  LWeights<-lapply(NETallProp$NETallList,function(x){
    x$Weights
  })
  LWeights<-(do.call('rbind',LWeights))
  matplot(LWeights,type='l',add=F, main="Evolution of weights through learning")

  VAR<-as.matrix(apply(LWeights[,],2,sd))
  par(mar=c(3,13,4,2))
  barplot(VAR[VAR>quantile(VAR,0.89,na.rm = T)],names.arg = apply(NETall1[VAR>quantile(VAR,0.89,na.rm = T),c(1,3)],1,function(x){
    paste(x,collapse = "_")}), las=2, cex.names = 0.5, horiz = T,space = 0, cex.axis = 0.5,
    main="Top weights' correction")

  # variations (sd) in outputs across patients
  if(Default$Interval>1){
    OutputsEvol<-lapply(NETallProp$NETallActivity,function(x){
      mean(apply(x[,grep("Output",colnames(x))],1,sd))
      #sd(x[,grep("Output",colnames(x))])
    })
  } else{
    OutputsEvol<-lapply(NETallProp$NETallActivity,function(x){
      sd(x[,"Output"])
    })
  }
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)+0.1)
  matplot(unlist(OutputsEvol)[-1],type='l',add=F, main="Evolution of std deviations \n in outputs through learning")

   if(Default$Interval==1){
    Pred_train_mat<-NETallProp$NETallActivity[[N]][,"Output",drop=F]
    Ytrain<-as.data.frame(t(Y_base_raw)) # do that for IPCW
    plot(Pred_train_mat,Ytrain[Train,1],pch=ifelse(Ytrain[Train,2]==0,2,1),main="Train")
    CindexTrain<-survConcordance(Surv(Ytrain[Train,1],Ytrain[Train,2])~Pred_train_mat)
  } else {
    #  N=1
    par(mfrow=c(2,1))
    par(mar=c(4,4,2,2)+0.1)
    Pred_train_mat<-NETallProp$NETallActivity[[N]][,paste("Output",seq(Default$Interval),sep = "")]

    Pred_train_mat<-Pred_train_mat[colnames(Target[,Train]),]

    Target_train<-Target[,Train]
    ORD<-colnames(Target_train)[order(colMeans(Target_train))]
    image(t(Target_train[,ORD]),main="Train");image(Pred_train_mat[ORD,])

    par(mfrow=c(1,1))
    par(mar=c(4,4,2,2)+0.1)
    Pred_train<-rowMeans(Pred_train_mat)

    Ytrain<-as.data.frame(t(Y_base_raw)) # do that for IPCW

    if(!all(names(Pred_train)==rownames(Ytrain[Train,]))){
      Pred_train<-Pred_train[rownames(Ytrain[Train,])]
    }

    plot(Pred_train,Ytrain[Train,1],pch=ifelse(Ytrain[Train,2]==0,2,1),main="Train")
    CindexTrain<-survConcordance(Surv(Ytrain[Train,1],Ytrain[Train,2])~Pred_train)
    #  IPCW <- cindex(as.matrix(Pred_train), cens.model="marginal",
    #                 formula = Surv(Output,Status)~., data = Ytrain[Train,])
  }

  legend("topleft",legend = c(paste("C-index =", 1-round(CindexTrain$concordance,3),
                                    "; se =", round(CindexTrain$std.err,3)),
                              paste("MSE =",round(median(tail(NETallProp$Cost[,1],Default$MiniBatch)),3))#,
                              #paste("IPCW =",round(as.numeric(IPCW$AppCindex),3))
  ),cex=0.5)


  ### val
  # compute list of mut per patient from mat
  #MUTl<-MutMatToList(MUTa = Perturb)
  #length(MUTl$GOFl)

  ### simulate

  # load new net
  # NETall1<-NETallProp$NETallList[[N]]

  # simulate for every one

  TotAttractors<-BoolSimul(NETall=NETall1, Logic = Default$Logic, Mode = Default$Mode,
                           iStates=iStates, #[Val,],
                           Parallel = Default$Parallelize, no_cores=no_cores,
                           MinSteps = Default$MinStepsForward,
                           Discretize = F,
                           LSTM = Default$LSTM, Ct=Ct,#[Val,],
                           MUT = MUT, treatmt=NULL,
                           ValMut = Default$ValMut)

  if(is.null(names(TotAttractors[,1,1,1]))){
    names(TotAttractors[,1,1,1])<-rownames(iStates)
  }
  # sapply(1:15,function(x){
  #    cor(TotAttractors[x,"A",,5],TotAttractors1[x,"A",,5])
  #  })

  matplot(rbind(t(TotAttractors[1,"iStates",,1]),t(TotAttractors[1,"A",,])),
          type='l',main="1 patient", ylim=c(0,1),
          ylab="Nodes' activities", xlab = "Simulation steps")

  #TotAttractors[,"A",1:3,1]
  if(Default$Interval>1){
    Pred_val_mat<-TotAttractors[match(Val,rownames(iStates)),"A",paste("Output",seq(Default$Interval),sep = ""),dim(TotAttractors)[4]]
    #rownames(iStates[match(Val,rownames(iStates)),])==Val
    rownames(Pred_val_mat)<-Val
  } else {
    Pred_val_mat<-t(TotAttractors[rownames(iStates)%in%Val,"A","Output",dim(TotAttractors)[4]])
    colnames(Pred_val_mat)<-Val
  }


  #colnames(Y_base[,Val])==rownames(Pred_val_mat)
  #Pred_val_mat<-Pred_val_mat[colnames(Y_base[,Val]),]

  Cost <- mean( (t(Target[,Val])-Pred_val_mat)^2 )
  Cost

  # save cost val to perform learning curve
  NETallProp[["CostVal"]]<-Cost

  # re-save
  save(NETallProp, file = DIRSave )

  if(Default$Interval>1){
    #Pred_train_mat<-predict(base_model,as.array(as.matrix(XtrainN[Train,])) )
    par(mfrow=c(2,1))

    Target_val<-Target[,Val]
    ORD<-colnames(Target_val)[order(colMeans(Target_val))]
    image(t(Target_val[,ORD]),main="Val");image(Pred_val_mat[ORD,])

    # image(t(Target[,Val]),main="Val");image(Pred_val_mat)

    par(mfrow=c(1,1))
    par(mar=c(4,4,2,2)+0.1)
    Pred_val<-rowMeans(Pred_val_mat)

    Yval<-as.data.frame(t(Y_base_raw)) # do that for IPCW

    plot(Pred_val,Yval[Val,1],pch=ifelse(Yval[Val,2]==0,2,1),main="Val")
    CindexVal<-survConcordance(Surv(Yval[Val,1],Yval[Val,2])~Pred_val)
    # IPCW <- cindex(as.matrix(Pred_val), cens.model="marginal",
    #                 formula = Surv(Output,Status)~., data = Yval[Val,])

    legend("topleft",legend = c(paste("C-index =", 1-round(CindexVal$concordance,3),
                                      "; se =", round(CindexVal$std.err,3)),
                                paste("MSE =",round(as.numeric(Cost),3))#,
                                #paste("IPCW =",round(as.numeric(IPCW$AppCindex),3))
    ),cex=0.5)
    #  TrainIPCW<-as.numeric(IPCW$AppCindex)
    #  TrainIPCW;1-Cindex$concordance

  } else {
    plot(Pred_val_mat,Y_base_raw[1,Val],pch=ifelse(Y_base_raw[2,Val]==0,2,1),main="Val")
    CindexVal<-survConcordance(Surv(Y_base_raw[1,Val],Y_base_raw[2,Val])~as.numeric(Pred_val_mat))
    # IPCW <- cindex(as.matrix(Pred_val), cens.model="marginal",
    #                 formula = Surv(Output,Status)~., data = Yval[Val,])

    legend("topleft",legend = c(paste("C-index =", 1-round(CindexVal$concordance,3),
                                      "; se =", round(CindexVal$std.err,3)),
                                paste("MSE =",round(as.numeric(Cost),3))#,
                                #paste("IPCW =",round(as.numeric(IPCW$AppCindex),3))
    ),cex=0.5)
    #  TrainIPCW<-as.numeric(IPCW$AppCindex)
    #  TrainIPCW;1-Cindex$concordance
  }

  Cindex<-c(CindexTrain$concordance , CindexVal$concordance)
  names(Cindex)<-c("Train","Val")
  NETallProp[["Cindex"]]<-Cindex

  # re-re-save
  save(NETallProp, file = DIRSave )

  dev.off()


}

print("end")
}
