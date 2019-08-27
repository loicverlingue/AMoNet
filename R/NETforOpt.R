
#' This function do the first steps of building of the net (for future optimization)
#' @param NETall: data frame or NULL. to build a netwrok or your network table.
#' @param GENESman: Character vector. The genes query, in symbols corresponding to the InteractionBase.
#' @param treatmt: Character vector. The genes targets of your treatment, in symbols corresponding to the InteractionBase.
#' @param InteractionBase: data frame. The iteraction base used to build the network. Should be a data.frame in a sif-like format, with at least 3 colmuns: source_hgnc, interaction_directed_signed, target_hgnc
#' @param MeanWinit: Numeric. The mean value of the normal distribution used to initiate the weights.
#' @param SdWinit: Numeric. The standard deviation of the normal distribution used to initiate the weights.
#' @param RestrictedBuilding: Boolean. During the network buidling phase, restrict new species only to feedback connections within the net.
#' @param nblayers: Integer. Number of layers from the gene queries to explore within the InteractionBase (number of search steps). Should be >=1.
#' @param MinConnect: Integer. During the filtering phase, what is the minimal value of connections to keep a node. If MinConnect=1, you'll discard all nodes that are only "passengers".
#' @param FilterCGS: Boolean. Do you want to keep cancer genes only
#' @param LSTM: Boolean. Will you use LSTM unit update in the learning phase?
#' @param Adam: Boolean. Will you use another optimizer than GD (ie Adam, Momentum, RMSprop), to initiate related weights.
#' @param Activity:  Boolean. Run simulationn of the network to return activity values for each species.
#' @param addOutputs: Integer or NULL. Do you need outputs to the net, for example to simulate survival per interval. Specify the number of outputs. Outputs will be connected to phenotypes.
#' @param Phenotypes: data frame of the same format as InteractionBase but with species as source_hgnc and names of phenotypes as target_hgnc. Can be a table of gene sets for example. Default loaded are Hallmarks and immune gene sets.
#' @param KeepPhenotypes: Boolean. Phenotypes are required to connect output(s)
#' @param MECA : Character vector or NULL. the phenotypes you want to keep in the network.
#' @param RestrictedBase: Boolean. Restrict the building of the nets to genes related to the phenotypes selected with MECA
#' @param NameProj: Character. The name of your project
NETforOpt<-function(WRITE=F, NETall=NULL, GENESman=c("PDCD1", "CD274", "CTLA4","KRAS", "EGFR","TP53"),
                    treatmt=c("PDCD1", "CD274", "CTLA4"), InteractionBase = OMNI,
                    MeanWinit=0, SdWinit=0.1, RestrictedBuilding=F,
                    nblayers = 1, MinConnect=1, FilterCGS=F, LSTM=F, Adam = F,
                    Activity=T, addOutputs=NULL,
                    Phenotypes=NULL,  KeepPhenotypes=F, MECA=NULL,
                    RestrictedBase=F, NameProj="Model",no_cores=no_cores ){


  ########
  # load and define the FamilyGenes
  # Gene sets to add, and menu to select your mechanism of interest
  if(is.null(Phenotypes)){
    GenesSelecImuno<-read.csv2(paste(getwd(),"/data/ImmunGenes2.csv",sep = ""), stringsAsFactors = F, row.names = 1) # generic 13/7/18
    GenesSelecHall<-read.csv2(paste(getwd(),"/data/HallmarksGenes.csv",sep = ""), stringsAsFactors = F, row.names = 1) # generic 13/7/18
    Phenotypes<-rbind(GenesSelecImuno,GenesSelecHall)
  }

  #  if(is.null(Phenotypes)&KeepPhenotypes){
  #    print("A phenotype table is needed when KeepPhenotypes is TRUE")
  #    KeepPhenotypes=F
  #  }

  if(!is.null(addOutputs)&!KeepPhenotypes){
    print("Penotypes are needed to build outputs")
    KeepPhenotypes=T
  }
  ########
  # load data net
  # print("load files")

  #  OMNI<-read.csv2(paste(getwd(),"/data/OMNI_AngioMAPKImmuno.csv",sep = ""), stringsAsFactors = F, row.names = 1)

  #union(OMNI$source_hgnc,OMNI$target_hgnc)[!union(OMNI$source_hgnc,OMNI$target_hgnc)%in%rownames(DATA$expMat)]

  # if you want to reduce to specific functions:
  #MECA<-unique(GenesSelec$target_hgnc)[c(97,98,88,89)]
  #MECA<-c(MECA,unique(GenesSelecImuno$target_hgnc))
  if(length(MECA)==0){
    #  print("Keep all phenotypes")
    MECA<-unique(Phenotypes$target_hgnc)

  }

  if(FilterCGS){
    CGS<-read.csv(paste(getwd(),"/data/CancerGeneCensusCOSMIC.csv",sep = ""))
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
  # Cleaning the net
  # add an interactive function where you can choose the species you want to keep and discard

  if(FALSE){
    print("clean net")

    Unwanted<-c("PLCG2","PLCG1","CTNND1","RARA","PKN1","VAV2","FGR","LRIG1","RASA1","PTK2B",
                "CAV1","PRKACA","CSNK2A1","PAK1", "PLD1","ITGB4","EZR","PIK3R3","TNK2","BCAR1",
                "AMHR2","WASF1","ATF2","BMPR1B","ARHGDIA","NCK1","ACVR1","GNA13",
                "NUMB","PLCG2","PRKCA","PLCG1","RIN1","RIPK1","RARA","ICAM1","PKN1","VAV2","FGR","RASA1",
                "PTK2B","SMURF2","STK4","PRKCD", "BLM","IER3","MAPK9","RFWD2","HSP90AA1","HIPK2","PRKDC",
                "NPM1","NR3C1","KAT2B","EIF2AK2","PLK1","PIN1","HNRNPK","BCL2L1","DAPK1","PPP1R13L","LATS2",
                "RCHY1","RRM2B","NR4A1","HSPB1","KAT5","TRIM24","DAG","CDC42","SIAH1","RELA","HSPA8",
                "AURKA","CAV1","PRKACA","YWHAZ","VDR","PAK1","PYCARD","FYN","PLCB1","CBLB","PTK2", "CRKL","GRAP",
                "VAV1","PAK2","YWHAH","PAK3","HIC1","YWHAQ","PPM1D","DDB2","SHC1","NFATC4","TEK","SGK1","PML","BID",
                "PLK2","EPHA2","APAF1","HTT","RPS6KA2","DUSP1","HSP90AB1","RALGDS","PIDD1","GADD45A",
                "BNIP3L","ITGB4","NOXA1","PMAIP1","EZR","HGF","CFLAR","ITGB2","ETS1","TBET","GNA13","MAD2L1",
                "TBK1","SFN","ACVR1","HGS","NCK1","BMPR1B","BHLHE41","PTGS2","RGS16","AMHR2","TP53INP1","KLF4")
    InteractionBase<-InteractionBase[!InteractionBase$source_hgnc%in%Unwanted,]
    InteractionBase<-InteractionBase[!InteractionBase$target_hgnc%in%Unwanted,]
    Phenotypes<-Phenotypes[!Phenotypes$source_hgnc%in%Unwanted,]
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
      #NETall<-OMNI[OMNI$source_hgnc%in%union(GENESman,treatmt)|OMNI$target_hgnc%in%union(GENESman,gsub("_TTT","",treatmt)),]
      return(list(NETall=NETall))
    }

    if(!all(GENESman%in%union(NETall$source_hgnc,NETall$target_hgnc))){
      print(paste("Gene(s) not retained after cleaning :", paste(GENESman[!GENESman%in%union(NETall$source_hgnc,NETall$target_hgnc)],collapse = ", ")))
    }
  }

  # but necessary to peform the layer based fonctions
  if(!is.null(addOutputs)){
    # addings to NETall
    # outputs

    if(FALSE){
      MECA<-unique(Phenotypes$target_hgnc)[c(97,98,88,89)]
      MECA<-c(MECA,unique(PhenotypesImuno$target_hgnc))
    }

    #NETall<-Outputs(NETall,FamilyGene = Phenotypes, FinalOutput = NULL)

    # check if want to have only selected phenotypes or all
    NETall<-Outputs(NETall, FamilyGene = Phenotypes[Phenotypes$target_hgnc%in%MECA,],
                    FinalOutput = addOutputs, FilterFamily = 0) # addOutputs=1 for adding a single output


    # how much species are connected to outputs
    CONout<-table(union(NETall$source_hgnc,NETall$target_hgnc)%in%NETall$source_hgnc[NETall$Output])
    print(paste("On",sum(CONout), "species in the net,", as.numeric(CONout["TRUE"]),
                "species are connected to",length(unique(NETall[NETall$Output,3])), "outputs"))

    print("net analysis into layers, warnings in the plots but no pb")

    Analysis<-NetAnalyser(NETall, Propagation = "back",PDF=F, PLOT=F) # warnings due to autoregulated nodes (FOXP3): no danger

    NETall$Layer<-Analysis$Layers

    # add weigths to new nodes or remove new nodes

    NETall<-addWeights(NETall = NETall,SdWinit = Default$SdWinit,
                       MeanWinit = Default$MeanWinit ,Scaling = F,
                       Adam = Adam, LSTM = LSTM)

    #  PlotOptNet(NETall,PDF = F,Optimized = F,PrintOptNET = F,LEGEND = F)

  } else {

    # just adding phenotypes to perform the NetAnalysis, then remove it
    #  print("put outputs")
    #  setwd("/data/tmp/lverling/")
    if(!"Output"%in%colnames(NETall)){
      NETall<-Outputs(NETall, FamilyGene = Phenotypes[Phenotypes$target_hgnc%in%MECA,], FinalOutput = NULL) # FinalOutput=T for adding a unique output
    }

    # how much species are connected to outputs
    # CONout<-table(union(NETall$source_hgnc,NETall$target_hgnc)%in%NETall$source_hgnc[NETall$Output])
    #  print(paste("On",sum(CONout), "species in the net,", as.numeric(CONout["TRUE"]),
    #              "species are connected to",length(unique(NETall[NETall$Output,3])), "outputs"))

    #
    print("net analysis into layers, warnings in the plots but no pb")
    #  setwd("/data/tmp/lverling/")

    Analysis<-NetAnalyser(NETall, Propagation = "back",PDF=F, PLOT = F) # warnings due to autoregulated nodes (FOXP3): no danger
    NETall$Layer<-Analysis$Layers
    #table(NETall$Layer)
    # add weights
    print("add weights")
    NETall<-addWeights(NETall, MeanWinit = abs(MeanWinit), SdWinit = abs(SdWinit),
                       Scaling=T, Adam = Adam, LSTM = LSTM)

    if(!KeepPhenotypes){
      print("removing phenotypes")
      NETall<-NETall[!NETall$Output,]
    }
  }

  if(Activity){

    # do a simulation
    print("a simulation")

    PossiStates<-InitStates(NETall = NETall, istates = NULL, Logic = "Sigmoid", NumiStates = 10)
    if(dim(PossiStates)[1]!=10){
      PossiStates = NULL
    }

    #    TotAttractors<-BoolSimul(NETall = NETall, Logic = "Sigmoid", Mode = "ASYNC", iStates = PossiStates,
    #                             ComputeZ = T, MinSteps = 10, Discretize = F, LSTM = T, GOF = NULL, LOF = NULL,
    #                             Inputs = NULL)
    TotAttractors<-BoolSimul(NETall = NETall,
                             Logic = "Sigmoid", Mode = "LAYER",
                             iStates =PossiStates,
                             MinSteps =4,
                             Inputs=NULL, LSTM = LSTM, Discretize = F,
                             GOF = NULL, LOF = NULL, Parallel = F)

    # representation
    #    SpeciesActivity<-plotProbaState(Species = union(NETall$source_hgnc,NETall$target_hgnc),
    #                                    TotAttractors = TotAttractors, Level = "A", PDF = F, NameProj = NameProj)

    #  SpeciesActivity<-plotProbaState(Level = "A",Species = union(NETall$source_hgnc,NETall$target_hgnc),
    #                 TotAttractors = TotAttractors,Plot = T, PDF = F, Legend = F )

    #matplot(rbind(t(TotAttractors[2,"iStates",,1]),t(TotAttractors[2,"A",,])),
    #        type='l',main="1 Val patient",ylim=c(0,1))

    iStates<-t(colMeans(TotAttractors[,"A",,dim(TotAttractors)[4]]))
    if(LSTM){
      Ct<-t(colMeans(TotAttractors[,"Ct",,dim(TotAttractors)[4]]))
    }
    # iStates<-SpeciesActivity[nrow(SpeciesActivity),] # compute initial states from attractors: speed up learning
    #  iStates[which(iStates==0)]<-iStates[which(iStates==0)]+1e-8

    # put activity corresponding to target sepcies
    for(i in colnames(SpeciesActivity)){
      NETall[NETall$target_hgnc%in%i,"Activity"]<-SpeciesActivity[nrow(SpeciesActivity),i]
    }
    NETall$Activity[which(NETall$Activity==0)]<-NETall$Activity[which(NETall$Activity==0)]+1e-8
  }

  if(WRITE){
    Latt<-length(list.files(paste(getwd(),"/model/",sep = ""),pattern =  NameProj))
    if(Latt>0){
      Latt<-max(as.numeric(gsub(".csv","", gsub(NameProj,"", list.files(paste(getwd(),"/model/",sep = ""),pattern =  NameProj)))))
    }
    write.csv2(NETall,paste(getwd(),"/model/",NameProj,Latt+1,".csv",sep = ""))
  }
  if(Activity&LSTM){
    return(list(NETall=NETall,iStates=iStates,Ct=Ct))
  }else if(Activity&!LSTM){
    return(list(NETall=NETall,iStates=iStates))
  } else {
    return(list(NETall=NETall))
  }
}

