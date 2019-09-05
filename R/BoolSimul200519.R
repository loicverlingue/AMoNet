#' Simulates the network model
#'
#' Called within the \code{simulate()} *AMoNet* objects function
#'
#' @param NETall data frame. The network in format (similar to .sif with) columns "source", "interaction", "target", additionally storing the architecture (column "Layer"), the parameters and optionnaly Boolean logic.
#' @param Logic character. The update function used for the simulations. Choices are: "Sigmoid" (default), "Boolean","tanh" or "ReLU".
#' @param Mode character. The type of simulation to perform. Choice are "LAYER" (layer based, default), "ASYNC" (assynchronous update), "SYNC" (synchronous update)
#' @param iStates matrix. Initial states of the whole network model. Either set from data or randomized.
#' @param MUT list. Corresponds to mutations with names and corresponding values (list of vectors). If MUT=NULL, a wild type *AMoNet* is simulated.. Patients ordered the same as y.
#' @param treatmt list. The same for than newMUT, with treatments' targets and corresponding values (list of vectors).
#' @param ValMut numeric. Multicplicative factor for the effect of mutations. Default to 50.
#' @param MinSteps numeric. Number of simulation steps. Default to 5. Ideally set the number of simulation steps to achieve stable states of *AMoNet* activities.
#' @param LSTM boolean. Should a Long-Short Term Memory (LSTM) unit architecture be used?
#' @param Ct matrix. If LSTM=TRUE, Ct is used to set the initial states of the correponding values in LSTM units.
#' @param Parallel boolean. Should the simulation run in parallel (TRUE by default)? Not recommended for simulations of less than 10 examples.
#' @param no_cores numeric. If Parallel=TRUE, set the number of cores to parallelize on. Default is 4. Can detect and set to available no_cores if inferior to user defined no_cores. to
#' @param Steepness numeric. Can be used to modulate the stepness of the sigmoid curve (default to 1).
#'
#' @details
#' The ASYNC mode perform a fully assynchroneous simulation that may be relatively long. The LAYER mode is an innovation in *AMoNet* to further speed up the simulations and nevertheless optain similar activities to the ASYNC mode. SYNC mode is not recommended.
#'
#' @return TotAttractors, a 4-dimensions array with: D1 the number of examples, D2 the results of the linear function (z=ax+b), D3 the results of the activator function (defined by the \code{Logic} argument) and D4 the number of iterations.
#'
BoolSimul<-function(NETall=NETall, Logic="Sigmoid", Mode="LAYER",
                    iStates=NULL, MUT=NULL, treatmt=NULL, ValMut=100, MinSteps=5,
                   LSTM=F, Ct=NULL, Parallel=TRUE, no_cores=4, Steepness=1 ){ # ComputeZ=TRUE,
#  Tic<-Sys.time()

  # integrate treatm to LOF and GOF
  if(!is.null(treatmt)){
    if(!is.list(treatmt)|length(treatmt)!=nrow(iStates)){ # to check
      print("treatm should be a list and length(treatmt)==nrow(iStates)")
    } else {
      MUT<-lapply(seq(length(treatmt)),function(TT){
        c(MUT[[TT]][!names(MUT[[TT]])%in%names(treatmt[[TT]])],treatmt[[TT]])
      })
    }
  }

  ##
  # load functions used
  Sigmoid<-function(z){1/(1+exp(-z*Steepness))}
  tanh<-function(z){(exp(z)-exp(-z))/(exp(z)+exp(-z))}
  ReLU<-function(z){max(0,z)}
  mReLU<-function(z){ifelse(z<=0,z*0.2,z ) }

  if( !(is.list(MUT)|is.null(MUT)) ) {
    stop(print("MUT should be either lists of alterations or NULL"))
  }
 # if(!Boolean){Sigmoid=T}else{Sigmoid=F}
  Species<-union(NETall$target_hgnc,NETall$source_hgnc)

  # do an array
  NAMESLSTM<-c("iStates","InputGate","ZInputGate",
               "ForgetGate","ZForgetGate",
               "OutputGate","ZOutputGate",
               "Ct","Gt","ZGt","A")
  NAMESregular<-c("iStates","A","Z")

  if(!"Boolean"%in%Logic){
    TotiState<-array(1, dim = c(nrow(iStates),ifelse(LSTM,length(NAMESLSTM),length(NAMESregular)),length(Species), MinSteps),
                     dimnames = list(NULL,if(LSTM){NAMESLSTM}else{NAMESregular},Species,NULL) )
  } else {
    TotiState<-array(T, dim = c(nrow(iStates),2,length(Species), MinSteps),
                     dimnames = list(NULL,c("iStates","A"),Species,NULL) )
  }

  TotiState[,"iStates",,1]<-iStates

  if(!is.null(Ct)&LSTM){
    TotiState[,"Ct",,1]<-iStates
  }

 #Ospecies=c("EGFR","TP53")
  #Ospecies="TP53"

  ####### functions

  Update<-function(NETall=NETall, Ospecies="ERBB2", Round=Round,
                   TotiState=TotiState,NameWeigts="Weights",
                   NameB="bterm", FUN=Logic){ #GOF=NULL, LOF=NULL,

    W<-BuildNetMat(NETall = NETall, Layer = NULL,Type = NameWeigts) # sparse matrix of weights
    # to be sure it is in the rigth order, but maybe not that efficient...
    #Bterms<-apply(B[,Ospecies],2,function(x){sum(unique(x))})

    if(length(Ospecies)==1){
      Bterms<-as.matrix(rep(unique(NETall[NETall$target_hgnc%in%Ospecies,NameB]),dim(TotiState)[1]))
    } else {
      B<-BuildNetMat(NETall = NETall, Layer = NULL,Type = NameB) # weights for this layer
      Bterms<-t(replicate(as.numeric(apply(B[,Ospecies,drop=F],2,function(x){sum(unique(x))})),n = dim(TotiState)[1],simplify = T))
    }

    # because iState is updated with previous A in Simulation fonction
    TOT<-TotiState[,"iStates",rownames(W),Round]%*%W[,Ospecies,drop=F]+Bterms

    A=match.fun(Logic)(TOT)

    return(list(A=A,Z=TOT))
  }

  ############
  # other updates formalisms

  if("Boolean"%in%Logic){

    Update<-function(NETall,Ospecies,TotiState,Round,
                     NameWeigts=NULL, NameB = NULL,
                     FUN=Logic){

      #SP="Output2"
      ResSP<-sapply(Ospecies,function(SP){


        TabAND<-NETall[NETall$target_hgnc%in%SP,]

        GREP<-TabAND$LogicRule%in%"OR"

        LOGOR<-as.matrix(TotiState[,"iStates",TabAND[GREP,1],Round])
        if(length(LOGOR)>0){

          TabANACT<-t(replicate(dim(LOGOR)[1],TabAND[GREP,"interaction_directed_signed"]%in%"ACTIVATE",simplify = T))
          if(dim(TabANACT)[1]!=dim(LOGOR)[1]){TabANACT<-t(TabANACT)}
          dim(TabANACT);dim(LOGOR)
          LOGOR<-!xor(TabANACT,LOGOR)
        }
        # the species with whom is the AND connection is the other species with a AND annotation!
        GREP<-TabAND$LogicRule%in%"AND"

        LOGAND<-as.matrix(TotiState[,"iStates",TabAND[GREP,1],Round])

        if(length(LOGAND)>0){
          TabANACT<-t(replicate(dim(LOGAND)[1],TabAND[GREP,"interaction_directed_signed"]%in%"ACTIVATE",simplify = T))
         # dim(TabANACT);dim(LOGAND)
          if(dim(TabANACT)[1]!=dim(LOGAND)[1]){TabANACT<-t(TabANACT)}

          LOGAND<-!xor(TabANACT,LOGAND)
        }

        Res<- sapply(seq(dim(TotiState)[1]),function(LOGX){
          any(LOGOR[LOGX,])&all(LOGAND[LOGX,])
        })

        return(Res)
      })
      return(list(A=ResSP))
    }
  }  # add others: Hill...


  #Ospecies<-'FOXP3'
  if(LSTM){
    # to change by using Logic functions...

    ModeUpdate<-function(NETall=NETall,Ospecies=Ospecies,MUT=MUT,
                         Logic=Logic, Round = Round, ValMut=ValMut,
                         TotiState=TotiState, Update.=Update, Parallel=Parallel){

      if(Parallel){
        # parallelise
        FunDf<-list(First=c("Whi","bi", Logic),
                    Sec=c("Whf","bf", Logic),
                    Third=c("Who","bo", Logic),
                    Fourth=c("Weights","bterm", "tanh"))

        EXPORT<-c("FunDf","NETall","Update.","Ospecies", "TotiState", "Round",
                  "Sigmoid", "tanh", "ReLU", "BuildNetMat")
        parallel::clusterExport(cl,EXPORT, envir = environment())

        Lgates<-parallel::parLapply(cl,FunDf,function(FunDfx){
          Update.(NETall=NETall,Ospecies=Ospecies,TotiState=TotiState,
                  NameWeigts=FunDfx[[1]], NameB = FunDfx[[2]],
                  FUN=FunDfx[[3]], Round = Round)
        })

        names(Lgates)<-c("InputGate","ForgetGate","OutputGate","Gt")
        list2env(Lgates,envir = .GlobalEnv)

      } else {
        InputGate<-Update.(NETall=NETall,Ospecies=Ospecies,TotiState=TotiState,
                           NameWeigts="Whi", NameB = "bi", FUN=Logic, Round = Round)
        ForgetGate<-Update.(NETall=NETall,Ospecies=Ospecies,TotiState=TotiState,
                            NameWeigts="Whf", NameB = "bf", FUN=Logic, Round = Round)
        OutputGate<-Update.(NETall=NETall,Ospecies=Ospecies,TotiState=TotiState,
                            NameWeigts="Who", NameB = "bo", FUN=Logic, Round = Round)
        Gt<-Update.(NETall=NETall,Ospecies=Ospecies,TotiState=TotiState,
                    NameWeigts="Weights", NameB = "bterm", FUN="tanh", Round = Round)
      }
      #  Ctpre<-as.matrix(Ct[,NETall[NETall$target_hgnc%in%Ospecies,"source_hgnc"],drop=F])
      # Round?
#      names(TotiState[1,,1,1])

#      TotiState[,"Ct",Ospecies,Round]<-ForgetGate$A*TotiState[,"Ct",Ospecies,Round]+InputGate$A*Updated1$A
#      TotiState[,"Ht",Ospecies,Round]<-OutputGate$A*match.fun(Logic)(TotiState[,"Ct",Ospecies,Round])
      TotiState[,"InputGate",Ospecies,Round]<-InputGate$A
      TotiState[,"ZInputGate",Ospecies,Round]<-InputGate$Z
      TotiState[,"ForgetGate",Ospecies,Round]<-ForgetGate$A
      TotiState[,"ZForgetGate",Ospecies,Round]<-ForgetGate$Z
      TotiState[,"OutputGate",Ospecies,Round]<-OutputGate$A
      TotiState[,"ZOutputGate",Ospecies,Round]<-OutputGate$Z
      TotiState[,"Gt",Ospecies,Round]<-Gt$A
      TotiState[,"ZGt",Ospecies,Round]<-Gt$Z

     # last steps:
      TotiState[,"Ct",Ospecies,Round]<-ForgetGate$A*TotiState[,"Ct",Ospecies,Round]+InputGate$A*Gt$A
      TotiState[,"A",Ospecies,Round]<-OutputGate$A*match.fun(Logic)(TotiState[,"Ct",Ospecies,Round])

      # in case of mutation, A is altered
      if(!is.null(MUT)){
        Trash<-sapply(seq(dim(TotiState)[1]),function(NPat){
          if(any(Species%in%names(MUT[[NPat]]))){
            TotiState[NPat,"A",names(MUT[[NPat]]),Round]<<- MUT[[NPat]]*ValMut #match.fun(Logic)(1)
          #  TotiState[NPat,"A",Species%in%GOF[[NPat]],Round]<<- 1
          }
        })
      }

      return(TotiState)
    }

  } else {
    ModeUpdate<-function(NETall=NETall,Ospecies=Ospecies, MUT=MUT,
                         Logic=Logic, Round = Round,
                         TotiState=TotiState, Update.=Update, Parallel=Parallel){

      Updated1<-Update.(NETall=NETall,Ospecies=Ospecies,TotiState=TotiState,
                        NameWeigts="Weights", NameB = "bterm",  #GOF=GOF, LOF=LOF,
                        FUN=Logic, Round = Round)

      for(Name in names(Updated1)){
        TotiState[,Name,Ospecies,Round]<-Updated1[[Name]]
      }
      #TotiState[,"A",Ospecies,Round]<-Updated1$A
      #TotiState[,"Z",Ospecies,Round]<-Updated1$Z
      # in case of mutation, A is altered
      if(!is.null(MUT)){

        Trash<-sapply(seq(dim(TotiState)[1]),function(NPat){
          if(any(Species%in%names(MUT[[NPat]]))){
            if("Boolean"%in%Logic){
              TotiState[NPat,"A",names(MUT[[NPat]]),Round] <<- MUT[[NPat]]>0
            } else {
              TotiState[NPat,"A",names(MUT[[NPat]]),Round] <<- MUT[[NPat]]*ValMut
            }
          }
        })
      }

      return(TotiState)

    }
  }

  # put connectors
#  if(length(grep("LogicRule",colnames(NETall)))==0){
  if(!"LogicRule"%in%colnames(NETall)&"Boolean"%in%Logic){
    NETall[,"LogicRule"]<-"OR" # begin with OR connectors only
    NETall[NETall$interaction_directed_signed%in%"INHIBIT","LogicRule"]<-"AND"
  }

  ################
  # initial states
  # initial states changed the 02/05/18
  # turn to matrices the 14/11/18
  if(is.null(iStates)){
    iStates<-t(rnorm(length(Species),0,5))
    colnames(iStates)<-Species
    ListAttractors<-list(iStates)
  } else if(is.vector(iStates)){
    ListAttractors<-list(iStates)
    iStates<-t(iStates)
    colnames(iStates)<-Species
  } else if(is.matrix(iStates)){
    ListAttractors<-as.list(data.frame(t(iStates)))
    ListAttractors<-lapply(ListAttractors,function(x){
      names(x)<-Species
      return(x)
    })
  } else if(is.data.frame(iStates)){
    ListAttractors<-as.list(data.frame(t(iStates)))# as.list(iStates)
    ListAttractors<-lapply(ListAttractors,function(x){
      names(x)<-Species
      return(x)
    })
    iStates<-as.matrix(iStates)

  } else if(is.list(iStates)){
    ListAttractors<-iStates
    iStates<-do.call("rbind",iStates)

  } else {
    print("iStates should be vector, data.frame, matrix or list of vectors,
          with name of all species as columns and example as rows")
  }

  if(length(ListAttractors[[1]])!=length(Species)){
    print("iState should be vector, data.frame, matrix or list of vectors,
          with name of all species as columns and example as rows")
  }


  #iStates<-ListAttractors[[9]]

  ##########
  # simulations

  #TotiState=TotiState[1,,,,drop=F]
    Simulation<-function(TotiState=TotiState, NETall=NETall, Mode=Mode, MinSteps=MinSteps,
                         LSTM=LSTM, MUT=MUT, Parallel=Parallel){ #,ComputeZ=ComputeZ

      Species<-union(NETall$target_hgnc,NETall$source_hgnc)
      Nspecies<-length(Species)
      if(is.null(MinSteps)){MinSteps<-Nspecies}

      # table to store the updates
      # tables should be dimensions: [N, Type, Species, iterations], with Type = A, Z and/or LSTM values ...
      ########  ###############
    # random selection of the first one / select how many? tot, 1/2, 1? -> now with layers

    if("LAYER"%in%Mode){
     # CHECK<-TRUE
      NLayer<-1
      ULayer<-sort(unique(NETall$Layer))
      Round<-1

      while(Round<=MinSteps){
          # Nlayer+1 updated at the end of the while loop
        if(NLayer>length(ULayer)){
          NLayer<-1
          TotiState[,"iStates",,Round]<-TotiState[,"A",,Round-1]
          if(LSTM){
            #udpate Ct to next round
            TotiState[,"Ct",,Round]<-TotiState[,"Ct",,Round-1]
          }
        }
        # select new species
          NewSpecies<-unique(NETall$target_hgnc[NETall$Layer==ULayer[NLayer]]) # target

        #  if(LSTM){
            TotiState<-ModeUpdate(NETall = NETall, Ospecies = NewSpecies,
                                  TotiState=TotiState,MUT=MUT,Logic = Logic,
                                   Update.=Update, Round = Round, Parallel = Parallel) # Ct = TotiState[,"Ct",,],

            #Ospecies<-NewSpecies[3]
#            TotiState[,"A",NewSpecies[3],]
            #### to check
          #  iStates[,Ospecies,drop=F]<<-Updated1$Ht
          #  Ct[,Ospecies,drop=F]<<-Updated1$Ct

          #  TotUpdates[[length(TotUpdates)+1]] <- Updated
          #  TotUpdates[,colnames(Updated)]<-as.data.frame(Updated)
          #  rownames(TotUpdates)<-rownames(Updated)

   #         #  return(UpdatedLSTM)
  #        } else{
  #          Updated1<-Update(NETall = NETall, Ospecies=Ospecies,iStates=iStates, NameWeigts = "Weights", NameB = "bterm")
          #  iStates[,Ospecies,drop=F]<<-Updated1$A
  #        }

         # names(Updated1)
        #  x<-"A"

  #        sapply(names(Updated1),function(x){
  #          TotiState[,x,Ospecies,Round]<<-Updated1[[x]]
  #        })



   #       return(Updated1)})

       # species that changes:
        #ChangedSpecies<-names(Updated)[Updated!=iStates[NewSpecies]]
                # how to keep the name of the species in transitions? trying with ChangedSpecies
      #  if(ChangedSpecies!=""){
      #    rownames(TotiState)[nrow(TotiState)]<-paste(paste(ChangedSpecies,collapse = "_"),length(grep(paste(ChangedSpecies,collapse = "_"),rownames(TotiState)))+1,sep = "_")
      #  } else {
      #    rownames(TotiState)[nrow(TotiState)]<-nrow(TotiState)
      #  }
      # NewSpecies<-NETall$target_hgnc[NETall$Layer==NLayer+1]

      # NewSpecies<-unique(NETall[NETall$source_hgnc%in%NewSpecies,"target_hgnc"])

        # stoping rule
        # do next layer
        NLayer<-NLayer+1
        if(NLayer>length(ULayer)){
          Round<-Round+1
        }
        if(Round>MinSteps){
          Round<-Round-1
          break
          }
      } # end of while loop

    #  TotiState[,"A",,10]

    } else if("ASYNC"%in%Mode){

      SpeciesOrd<-unique(c(NETall$source_hgnc[order(NETall$Layer)],NETall$target_hgnc[order(NETall$Layer)]))
      Round<-1
      while(Round<=MinSteps){

        # update iState withe previous A
        if(Round>1){
          TotiState[,"iStates",,Round]<-TotiState[,"A",,Round-1]
          if(LSTM){
            #udpate Ct to next round
            TotiState[,"Ct",,Round]<-TotiState[,"Ct",,Round-1]
          }
        }

        for(NewSpecies in SpeciesOrd){
         # TotiState<-ModeUpdate(NETall = NETall, Ospecies = NewSpecies,
        #                        TotiState=TotiState, Update.=Update, Round = Round)
          TotiState<-ModeUpdate(NETall = NETall, Ospecies = NewSpecies,
                                TotiState=TotiState,MUT=MUT,Logic = Logic,
                                Update.=Update, Round = Round, Parallel = Parallel) # Ct = TotiState[,"Ct",,],

        }
        Round<-Round+1

      } #end of while loop


      } else if("SYNC"%in%Mode){

        Round<-1
        while(Round<=MinSteps){

          if(Round>1){
            TotiState[,"iStates",,Round]<-TotiState[,"A",,Round-1]
            if(LSTM){
              #udpate Ct to next round
              TotiState[,"Ct",,Round]<-TotiState[,"Ct",,Round-1]
            }
          }

            #TotiState<-ModeUpdate(NETall = NETall, Ospecies = Species,
            #                      iStates = iStates, TotiState=TotiState,
          #                        Update.=Update, Round = Round) # Ct = TotiState[,"Ct",,],
            TotiState<-ModeUpdate(NETall = NETall, Ospecies = NewSpecies,
                                  TotiState=TotiState,MUT=MUT, Logic = Logic,
                                  Update.=Update, Round = Round, Parallel = Parallel) # Ct = TotiState[,"Ct",,],

          Round<-Round+1
        } #end of while loop

      } # end of SYNC

    return(TotiState)

  }  # end of simulation function

  #  TotiState[1,,"EGFR",]

  #### run simulations
  # Initiate cluster
  if(Parallel){

      no_cores_detect <- parallel::detectCores()
      if(no_cores_detect<no_cores){
        print(paste("Changed",no_cores,"to",no_cores_detect,"detected cores"))
        no_cores<-no_cores_detect
      }

    if(no_cores>nrow(iStates)/2){
      no_cores1<-no_cores
      no_cores<-floor(nrow(iStates)/2)
      print(paste("Changed",no_cores1,"to",no_cores,"cores"))
    }
  cl <- parallel::makeCluster(no_cores)

  }

  if(Parallel&!LSTM){

  NperB<-ceiling(nrow(iStates)/no_cores)

    EpochMini<-list()
    #NMini=1
    #length(seq((NMini*NperB-NperB)+1,NMini*NperB ))
    #length( seq((NMini*NperB-NperB)+1,ncol(iStates)))
    for(NMini in seq(no_cores)){
      if(NMini == no_cores){ # last one
        EpochMini[[NMini]]<- seq((NMini*NperB-NperB)+1,nrow(iStates))
      } else {
        EpochMini[[NMini]]<- seq((NMini*NperB-NperB)+1,NMini*NperB )
      }
    }


    #Inputs=NULL
    #Mode="LAYER"
    #Logic="tanh"
    #MinSteps=30
    #GOF=NULL
    #LOF=NULL
    #LSTM=T
    #Discretize=F

    EXPORT<-c("Update","Simulation","ListAttractors","NETall", "Mode", "Logic","ModeUpdate",
              "MUT","MinSteps","LSTM","TotiState","EpochMini",
              "BuildNetMat", "Steepness","ValMut")
    if("Sigmoid"%in%Logic){
      EXPORT<-c(EXPORT,"Sigmoid")
    } else if("tanh"%in%Logic){
      EXPORT<-c(EXPORT,"tanh")
    }else if("ReLU"%in%Logic){
      EXPORT<-c(EXPORT,"ReLU")
    }

    parallel::clusterExport(cl,EXPORT, envir = environment())
  # NMini<-EpochMini[[2]]

    LAttractorsAZ<-parallel::parLapply(cl,EpochMini,function(NMini){
      Simulation(TotiState = TotiState[NMini,,,,drop=F] , NETall=NETall,
                 Mode=Mode, MinSteps=MinSteps, LSTM = LSTM,
                 MUT=MUT[NMini,drop=F], Parallel=Parallel)
    })

    Trash<-lapply(seq(length(EpochMini)),function(NumPat){
      TotiState[EpochMini[[NumPat]],,,]<<-LAttractorsAZ[[NumPat]]
    })

    #names(TotiState1[1,1,,dim(TotiState)[4]])
    #dim(TotiState)

    parallel::stopCluster(cl)

  } else {

      TotiState<-Simulation(TotiState = TotiState , NETall=NETall,
                 Mode=Mode, MinSteps=MinSteps,LSTM = LSTM,
                 MUT=MUT, Parallel=Parallel)

      if(Parallel){
        stopCluster(cl)
      }
    }

  return(TotiState)

}
