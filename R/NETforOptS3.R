#' Build the AMoNet network
#'
#' @param object *AMoNet* object, S3 class. if NULL (default) \code{build.AMonet} function will create a new *AMoNet* object from scratch.
#'   If you already have a AMoNet model, \code{build.AMoNet} can use it and run following argument on it (eg Activity, Interval, ...).
#' @param InteractionBase data frame. The iteraction base used to build the network. Should be a data.frame in a sif-like format, with at least 3 colmuns: source_hgnc, interaction_directed_signed, target_hgnc
#' @param MeanWinit numeric. The mean value of the normal distribution used to initiate the weights.
#' @param SdWinit numeric. The standard deviation of the normal distribution used to initiate the weights.
#' @param RestrictedBuilding boolean. During the network buidling phase, restrict new species only to feedback connections within the net.
#' @param nblayers integer. Number of layers from the gene queries to explore within the InteractionBase (number of search steps). Should be >=1.
#' @param MinConnect integer. During the filtering phase, what is the minimal value of connections to keep a node. If MinConnect=1, you'll discard all nodes that are only "passengers".
#' @param FilterCGS boolean. Do you want to keep cancer genes only
#' @param LSTM boolean. Will you use LSTM unit update in the learning phase?
#' @param Optimizer boolean. Will you use another optimizer than GD (ie Adam, Momentum, RMSprop), to initiate related weights.
#' @param Interval integer or NULL. Number of outputs to the network. To simulate descrete time survival, use Interval to set the number of time intervals. For label dimensions of 1, set \code{Interval=1}. Outputs will be connected to phenotypes. Default is 10.
#' @param Phenotypes data frame of the same format as InteractionBase but with species as source_hgnc and names of phenotypes as target_hgnc. Can be a table of gene sets for example. Default loaded are Hallmarks and immune gene sets.
#' @param KeepPhenotypes boolean. Phenotypes are required to connect output(s)
#' @param MECA character vector or NULL. the phenotypes you want to keep in the network.
#' @param RestrictedBase boolean. Restrict the building of the nets to genes related to the phenotypes selected with MECA.
#' @param NameProj character. The name of your project.
#' @param WRITE boolean. Should the network be saved as csv?
#'
#' @details
#' Run \code{print(Default)} to check and optionnaly set default values.
#'
#' @return An object of class *AMoNet*. Stores the name of the project, the network, the initial states and initial Ct if LSTM.
#' @export
build.AMoNet<-function(object, # GENESman=c("PDCD1", "CD274", "CTLA4","KRAS", "EGFR","TP53"),
                    #treatmt="",
                    MeanWinit=Default$MeanWinit, SdWinit=Default$SdWinit,
                    RestrictedBuilding=F,
                    nblayers = Default$nblayers, MinConnect=Default$MinConnect,
                    FilterCGS=F, LSTM=Default$LSTM, Optimizer = Default$Optimizer,
                    Interval=Default$Interval,  InteractionBase = OMNI,
                    Phenotypes=NULL,  KeepPhenotypes=F, MECA=names_MECA,
                    RestrictedBase=F, NameProj="My AMoNet",
                    no_cores=Default$no_cores, WRITE=F ){

  # update Default parameters
  CALL<-mget(names(formals()))
  CNames<-intersect(names(object$Parameters$Default),names(CALL))
  object$Parameters$Default[CNames] <- CALL[CNames]

#  if(!is.null(addOutputs)){
#    object$Parameters$Default$Interval <- addOutputs
#  }

  if("NETall"%in%names(object)){
    NETall<-object$NETall
  } else if("GENESman"%in%names(object)){
    GENESman<-object$GENESman
    NETall<-NULL
    treatmt<-object$treatmt
  } else {
    object<-AMoNet()
    GENESman<-object$GENESman
    treatmt<-object$treatmt
    NETall<-NULL
  }
  #
#  if(!is.null(object)){
#    NETall<-object$NETall
#  } else {
#    NETall<-NULL
#  }
  ########
  # load and define the FamilyGenes
  # Gene sets to add, and menu to select your mechanism of interest
  if(is.null(Phenotypes)){
    #GenesSelecImuno<-read.csv2(paste(getwd(),"/data/ImmunGenes2.csv",sep = ""), stringsAsFactors = F, row.names = 1) # generic 13/7/18
    #GenesSelecHall<-read.csv2(paste(getwd(),"/data/HallmarksGenes.csv",sep = ""), stringsAsFactors = F, row.names = 1) # generic 13/7/18
    Phenotypes<-rbind(GenesSelecImuno,GenesSelecHall)
  }

  #  if(is.null(Phenotypes)&KeepPhenotypes){
  #    print("A phenotype table is needed when KeepPhenotypes is TRUE")
  #    KeepPhenotypes=F
  #  }

  if(!is.null(Interval)&!KeepPhenotypes){
    print("Penotypes are needed to build outputs")
    KeepPhenotypes=T
  }
  ########
  # load data net
  if(length(MECA)==0){
    #  print("Keep all phenotypes")
    MECA<-unique(Phenotypes$target_hgnc)

  }

  if(FilterCGS){
    Phenotypes<-Phenotypes[Phenotypes$source_hgnc%in%CGS$Gene.Symbol,]
    InteractionBase<-InteractionBase[InteractionBase$source_hgnc%in%CGS$Gene.Symbol,]
    InteractionBase<-InteractionBase[InteractionBase$target_hgnc%in%CGS$Gene.Symbol,]
  }

  if(is.null(NETall)){
    Species<-union(InteractionBase$source_hgnc,InteractionBase$target_hgnc)
    if(RestrictedBase){
      LENG<-length(intersect(Phenotypes$source_hgnc[Phenotypes$target_hgnc%in%MECA],Species))
    } else {
      LENG<-length(intersect(Phenotypes$source_hgnc, Species))
    }
    print(paste("Building for",paste(head(GENESman),collapse = ","),"..., from",length(Species),"species, from which" , LENG, "are connected to", length(unique(MECA)), "biological functions"))
    #print("remenber, the size of the final network depend on the nb of species and sometimes the nb of layers (nblayers) you allow between species")
  }

  ####
  # build NET
  if(is.null(NETall)){
    print("build net")
    #  NET<-InteractionBase[InteractionBase$source_hgnc%in%union(GENESman,treatmt)|InteractionBase$target_hgnc%in%union(GENESman,gsub("_TTT","",treatmt)),]
    #nblayers=2
    if(RestrictedBase){
      GENES<-c(GENESman,Phenotypes$source_hgnc[Phenotypes$target_hgnc%in%MECA])
    } else {
      GENES<-c(GENESman,Phenotypes$source_hgnc)
    }

    NETall<-NetBuilding(GENESman=GENESman, treatmt=treatmt,
                        OMNI = InteractionBase, MinConnect = MinConnect,
                        VeryRestricted = RestrictedBuilding,
                        nblayers = nblayers,
                        FamilyGene = GENES,no_cores=no_cores ) #

    if(is.null(NETall)){
      return(list(NETall=NETall))
    }

    if(!all(GENESman%in%union(NETall$source_hgnc,NETall$target_hgnc))){
      print(paste("Gene(s) not retained after cleaning :", paste(GENESman[!GENESman%in%union(NETall$source_hgnc,NETall$target_hgnc)],collapse = ", ")))
    }
  }

  # but necessary to peform the layer based fonctions
  if(!is.null(Interval)){
    # addings to NETall
    # outputs

    # check if want to have only selected phenotypes or all
    NETall<-Outputs(NETall, FamilyGene = Phenotypes[Phenotypes$target_hgnc%in%MECA,],
                    FinalOutput = Interval, FilterFamily = 0) # Interval=1 for adding a single output


    # how much species are connected to outputs
    CONout<-table(union(NETall$source_hgnc,NETall$target_hgnc)%in%NETall$source_hgnc[NETall$Output])
    print(paste("On",sum(CONout), "species in the net,", as.numeric(CONout["TRUE"]),
                "species are connected to",length(unique(NETall[NETall$Output,3])), "outputs"))

    print("net analysis into layers, warnings in the plots but no pb")

    Analysis<-NetAnalyser(NETall, Propagation = "back",PDF=F, PLOT=F) # warnings due to autoregulated nodes (FOXP3): no danger

    NETall$Layer<-Analysis$Layers

    # add weigths to new nodes or remove new nodes
    NETall<-addWeights(NETall = NETall, SdWinit = object$Parameters$Default$SdWinit,
                       MeanWinit = object$Parameters$Default$MeanWinit , Scaling = T,
                       Adam = !is.null(object$Parameters$Default$Optimizer), LSTM = LSTM)

    #  PlotOptNet(NETall,PDF = F,Optimized = F,PrintOptNET = F,LEGEND = F)

  } else {

    # just adding phenotypes to perform the NetAnalysis, then remove it
    if(!"Output"%in%colnames(NETall)){
      NETall<-Outputs(NETall, FamilyGene = Phenotypes[Phenotypes$target_hgnc%in%MECA,], FinalOutput = NULL) # FinalOutput=T for adding a unique output
    }

    print("net analysis into layers, warnings in the plots but no pb")

    Analysis<-NetAnalyser(NETall, Propagation = "back",PDF=F, PLOT = F) # warnings due to autoregulated nodes (FOXP3): no danger
    NETall$Layer<-Analysis$Layers

    print("add weights")
    NETall<-addWeights(NETall, MeanWinit = abs(MeanWinit), SdWinit = abs(SdWinit),
                       Scaling=T, Adam = !is.null(object$Parameters$Default$Optimizer), LSTM = LSTM)

    if(!KeepPhenotypes){
      print("removing phenotypes")
      NETall<-NETall[!NETall$Output,]
    }
  }


  if(WRITE){
    Latt<-length(list.files(paste(getwd(),"/model/",sep = ""),pattern =  NameProj))
    if(Latt>0){
      Latt<-max(as.numeric(gsub(".csv","", gsub(NameProj,"", list.files(paste(getwd(),"/model/",sep = ""),pattern =  NameProj)))))
    }
    write.csv2(NETall,paste(getwd(),"/model/",NameProj,Latt+1,".csv",sep = ""))
  }

  # output
  #object<<-list()
  netLIST<-list(call=list(build_call=CALL), NETall=NETall,iStates=NULL,Ct=NULL) # Parameters=list(object$Parameters$Default=object$Parameters$Default, Boundaries=Boundaries)
  DIFF<-setdiff(names(object),names(netLIST))

  object<-c(object[DIFF],netLIST)
  class(object)<-"AMoNet"

  # update Default in .GlobalEnv
  #  Default<<-object$Parameters$Default
  assign("Default",object$Parameters$Default,envir = .GlobalEnv)
  # net<<-net
 # class(net$NETall) <- "build.AMoNet"
  return(object)
}
