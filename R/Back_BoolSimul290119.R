#' Back simulation, corresponds to the diffusion of the gradients
#'
#' @description \code{BackBoolSimul} is taking forward simulations to derivate loss, flip the net, and simulate the gradients
#'
#' @param NETall data frame. The network in format (similar to .sif with) columns "source", "interaction", "target", additionally storing the architecture (column "Layer"), the parameters and optionnaly Boolean logic.
#' @param Logic character. The update function used for the simulations. Choices are: "Sigmoid" (default), "Boolean","tanh" or "ReLU".
#' @param Mode character. The type of simulation to perform. Choice are "LAYER" (layer based, default), "ASYNC" (assynchronous update), "SYNC" (synchronous update)
#' @param Y numeric or matrix. labels
#' @param MinSteps numeric. Number of back-simulation steps. Default to 1.
#' @param FixNodes character. Freeze parameters of selected nodes.
#' @param Parallelize boolean. Under devellopment for back-simulation.
#' @param TotiState array. Results of the simulations with the \code{simulate()} function.
#' @param gradClipping nunmeric. Maximum values allowed for gradients otherwise clipped. Gradient scaling in dev.
#' @param LSTM boolean. Should a Long-Short Term Memory (LSTM) unit architecture be used?
#' @param alpha numeric. Regularization term. Set to 1 for L1N (Lasso) regularization, set to 2 for L2N (Ridge) regularization, between 1 and 2 with lambda 1 for L1N + L2N (Elastic Net) regularization (penalization parameter is alpha).
#' @param lambda numeric. the lambda value for Lasso, Ridge or Elastic Net regularization. Set to 0 corresponds to no regularization. When alpha between 1 and 2 (Elastic Net regularization), lambda should be set to 1 and penalization parameter is alpha.
#'
#' @details
#' \code{MinSteps} is corresponding to generic argument \code{nsim} in S3 function \code{simulate()}.
#'
#' @return gradients of parameters
#' c(list(TotdZState=TotdZState,TotdWState=TotdWState,TotdbState=TotdbState),Cache)
#'
BackBoolSimul<-function(NETall=NETall, Logic=c("Sigmoid","tanh"), Mode=c("LAYER","ASYNC","SYNC"),
                    Y=Y, MinSteps=MinSteps, FixNodes=FixNodes, #learning_rate =learning_rate,
                    Parallelize=Parallelize, TotiState=NULL,gradClipping=F, LSTM=LSTM,
                    alpha=0, lambda=0.1){

  # todo: parallellize

  # pool of functions
if("Sigmoid"%in%Logic){
    activefun<-function(z){1/(1+exp(-z))}
    derivfun<-function(z){activefun(z)*(1-activefun(z))}
  } else if("tanh"%in%Logic){
    activefun<-function(z){(exp(z)-exp(-z))/(exp(z)+exp(-z))}
    derivfun<-function(z){1-(activefun(z)^2)}
  }
  activetanh<-function(z){(exp(z)-exp(-z))/(exp(z)+exp(-z))}
  derivtanh<-function(z){1-(activefun(z)^2)}

  # sepcies
  OutputsNETall<-unique(NETall$target_hgnc[NETall$Output])
  InputsNETall<-unique(NETall$source_hgnc[NETall$Layer%in%1])
  Species<-union(NETall$target_hgnc,NETall$source_hgnc)

  ##############
  # flip the Net
  NETexplore<-NETall[,c(3,2,1,seq(ncol(NETall))[-c(1:3)])]
  colnames(NETexplore)<-colnames(NETall)
  # flip the layers
  NETexplore$Layer<-max(NETexplore$Layer)-NETexplore$Layer+1
  NETexplore$Inputs<-NETexplore$source_hgnc%in%OutputsNETall
  NETexplore$Output<-NETexplore$target_hgnc%in%InputsNETall

  # function to update gradients
  UpdateLSTMgrad<-function(dX=dX, NAMESnode=NULL, Layer=Layer, NETexplore=NETexplore,
                           TotiState=TotiState, LSTM=LSTM){

    NAMESGates<-names(TotiState[1,,1,dim(TotiState)[4]])
    if(dim(TotiState)[1]>1){
      ValuesPerNames<-lapply(NAMESGates,function(NG){
        TotiState[,NG,,dim(TotiState)[4]]
      })
    } else {
      ValuesPerNames<-lapply(NAMESGates,function(NG){
        t(TotiState[,NG,,dim(TotiState)[4]])
        #rownames(M)<-Species
      })
    }
    names(ValuesPerNames)<-NAMESGates


    # remove iState (not usefull) and change name of Ct so that it does change parent Ct
    ValuesPerNames$iStates<-NULL
    names(ValuesPerNames)[names(ValuesPerNames)%in%"Ct"]<-"Ctl"

    list2env(ValuesPerNames,.GlobalEnv)
    # activefun.<-function(z){1/(1+exp(-z))}
      #derivfun<-function(z){(1/(1+exp(-z)))*(1-(1/(1+exp(-z))))}
        # dX <-> dht in paper's equation

    ##############
    # Input gate
    # element wise derivative of Ct
    #unlist(LSTMvalues["ZOutputGate",names(Y)])

    # when dc_next is calculated, I've choosed to use it instead of Cactivity: to check
    # compute dCt before the function: to take into account dCt-1!
    # turn LSTMvalues list or dataframe into variables
    if(LSTM){

      # find target nodes to uptdate
      TargetNodes<-unique(NETexplore[NETexplore$source_hgnc%in%NAMESnode,3])
# ????does TargetNodes needs to be sorted????????? in Gt[,TargetNodes] for example
    #  NAMESnode==TargetNodes
    #  NAMESnode==colnames(iStates)
    #  W<-BuildNetMat(NETall = NETexplore,Layer = Layer, Type = "Weights")
    #  TargetNodes<-TargetNodes[match(colnames(W[NAMESnode,]),TargetNodes)] #==colnames(W[NAMESnode,])
      if(is.null(TargetNodes)){ print("TargetNodes problem") }


      dCt <- dX*derivtanh(Ctl[,NAMESnode,drop=F])*OutputGate[,NAMESnode,drop=F]

      # to parallelize?

      ##############
      # input Gate

      Wi<-BuildNetMat(NETall = NETexplore,Layer = Layer, Type = "Whi")

      dIt<-dCt*Gt[,NAMESnode,drop=F]
    dZi<-as.matrix(dIt*derivfun(ZInputGate[,NAMESnode,drop=F]))
    #dWhi<-t(Gt[,unique(NETexplore[NETexplore$source_hgnc%in%NAMESnode,3]),drop=F])%*%dZi
    dWhi<-t(Gt[,TargetNodes,drop=F])%*%dZi
    dbi<-colMeans(dZi) # colMeans added
    dXi <- dZi %*% Wi[NAMESnode,TargetNodes]

    ##############
    # output gate
    Wo<-BuildNetMat(NETall = NETexplore,Layer = Layer, Type = "Who")
    dOt<-dX*derivfun(Ctl[,NAMESnode,drop=F])# dX!!!
    dZo<-as.matrix(dOt*derivfun(ZOutputGate[,NAMESnode,drop=F]))
#    dWho<-t(Gt[,unique(NETexplore[NETexplore$source_hgnc%in%NAMESnode,3]),drop=F])%*%dZo
    dWho<-t(Gt[,TargetNodes,drop=F])%*%dZo
    dbo<-colMeans(dZo)
    dXo <- dZo %*% Wo[NAMESnode,TargetNodes]

    ##############
    # Forget gate
    Wf<-BuildNetMat(NETall = NETexplore,Layer = Layer, Type = "Whf")
    dFt<-dCt*(Ctl[,NAMESnode,drop=F]) ##
    dZF<-dFt*derivfun(ZForgetGate[,NAMESnode,drop=F])
    dWhf<-t(Gt[,TargetNodes,drop=F])%*%dZF
    dbf<-colMeans(dZF)
    dXf <- dZF %*% Wf[NAMESnode,TargetNodes]

    ###########
    # and the nodes activity...

    W<-BuildNetMat(NETall = NETexplore,Layer = Layer, Type = "Weights")
    dGt<-dCt*InputGate[,NAMESnode,drop=F]
    dZg<-dGt*derivfun(ZGt[,NAMESnode,drop=F])
    dWhg<-t(Gt[,TargetNodes,drop=F])%*%dZg
    dbg<-colMeans(dZg)
    dXg <- dZg %*% W[NAMESnode,TargetNodes]

    #colnames(dWhg)==NAMESnode
    #rownames(dWhg)==TargetNodes
    #rownames(dWhg)==colnames(dXg)

    # As X = ht-1 was used in multiple gates, the gradient must be accumulated here
    dXt <- (dXg + dXo + dXi + dXf)
    colnames(dXt)<-TargetNodes

    # compute next dC
    # maybe dCt-1 is computed only from dCt?
    dc_next <- dXt*derivfun(Ctl[,TargetNodes,drop=F])*OutputGate[,TargetNodes,drop=F]
    dc_next <- ForgetGate[,TargetNodes,drop=F] * dc_next

    if(gradClipping){
      Clip<-3
      dWho[dWho>Clip]<-Clip; dbo[dbo>Clip]<-Clip; dWhi[dWhi>Clip]<-Clip; dbi[dbi>Clip]<-Clip
      dWhf[dWhf>Clip]<-Clip; dbf[dbf>Clip]<-Clip; dWhg[dWhg>Clip]<-Clip; dbg[dbg>Clip]<-Clip
      dXt[dXt>Clip]<-Clip
      # and negatives too
      dWho[dWho<(-Clip)]<-(-Clip); dbo[dbo<(-Clip)]<-(-Clip); dWhi[dWhi<(-Clip)]<-(-Clip); dbi[dbi<(-Clip)]<-(-Clip)
      dWhf[dWhf<(-Clip)]<-(-Clip); dbf[dbf<(-Clip)]<-(-Clip); dWhg[dWhg<(-Clip)]<-(-Clip); dbg[dbg<(-Clip)]<-(-Clip)
      dXt[dXt<(-Clip)]<-(-Clip)

      # try z scores, but NAs!
    }

    return(list(dWho=dWho,dbo=dbo,dWhi=dWhi,dbi=dbi,dWhf=dWhf,dbf=dbf,
                dWhg=dWhg,dbg=dbg,dZg=dZg,dXt=dXt,dc_next=dc_next)) #,dc_next=dc_next
    } else {
      TargetNodes<-unique(NETexplore[NETexplore$source_hgnc%in%NAMESnode,3])
      if(is.null(TargetNodes)){ print("TargetNodes problem") }

      ##30/01/19

      dZ<-dX*derivfun(Z[,NAMESnode,drop=F])

      dW<-t(A[,TargetNodes,drop=F])%*%dZ
      db<-colMeans(dZ)

      W<-BuildNetMat(NETall = NETexplore, Layer = Layer, Type = "Weights")

      dXt <- dZ %*% W[NAMESnode,TargetNodes]

      return(list(dW=dW,db=db,dZ=dZ,dXt=dXt)) #,dc_next=dc_next

      if(FALSE){
      ###
      dZ=(dX)*(derivfun(Zactivity[NAMES]))  # check if value too big -> Inf --> NaN --> error

      dZStates[colnames(dZ)]<-dZ
      TotdZState<-rbind(TotdZState,dZStates) #!

      # 2. derivative of the weights
      #  dim(t(dZ))
      #   dim(t(Xactivity[NETexplore[NETexplore$source_hgnc%in%In,3]]))
      #       print('2nd step')
      dW=t(dZ)%*%t(Xactivity[NETexplore[NETexplore$source_hgnc%in%NAMES,3]]) # Xactivity at stable state for target species!

      dW<-reshape2::melt(dW)
      Rnames<-paste(dW[,1],dW[,2],sep = ".")
      dW<-dW[!duplicated(Rnames),]
      rownames(dW)<-paste(dW[,1],dW[,2],sep = ".")

      dW<-dW[rownames(dW)%in%names(dWStates),]
      dWStates[rownames(dW)]<-dW[,3]
      TotdWState<-rbind(TotdWState,dWStates) #!

      #  db=dZ
      dbStates[colnames(dZ)]<-dZ #! for real inputs
      TotdbState<-rbind(TotdbState,dbStates) #!
      }

    }
  } # end of update gradient functcion
  # put into shape compatible with NETall... not simple

  # function to compute and store (to keep track of) the gradients in LSTM
  StoreLSTMgrad<-function(Cache=Cache, NamesInteraction=NamesInteraction, grads=grads){

    # put in right format
   # dWhState<-sapply(grep("dW",names(grads),value = T ),function(Variable){(Variable)})
    DWi<-grep("dW",names(grads),value = T)
    #DW="dWho"
    dWhState<-lapply(DWi,function(DW){
      # reshape it
#      dW<-dWshape(grads[grep(DW,names(grads),value = T)],NamesInteraction)

      dW<-reshape2::melt(grads[[DW]])
      Rnames<-paste(dW[,1],dW[,2],sep = ".")
      dW<-dW[!duplicated(Rnames),]
      rownames(dW)<-paste(dW[,2],dW[,1],sep = ".")
      dW<-dW[rownames(dW)%in%NamesInteraction,]
    #  dW<-dW[NamesInteraction,]

   #   if(gradClipping&nrow(dW)>=1){dW[,3]<-ifelse(dW[,3]<3,dW[,3],3)}

      DWCache<-Cache[[grep(DW,names(Cache),value = T)]]
      DWCache[rownames(dW),]<-dW[,3]

#      dW<-dWshape(grads[[DW]],NamesInteraction)
      # include in cache
      #Cache[[DW]][,rownames(dW)]<-dW[,3]
      return(DWCache)
    })

    names(dWhState)<-DWi

    # dWhState$dWhg

    #Variable=grep("db",names(grads),value = T)[1]
    DBi<-grep("db",names(grads),value = T)
    #DB="dbg"
    #if(!all(DBi%in%names(grads))){print("DBi pb")}
    dbState<-lapply(DBi,function(DB){
      db<-grads[[DB]]
      #class(db)
      #if(class(db)!="matrix"){db<-as.matrix(db)}
   #   if(gradClipping){db<-ifelse(db<3,db,3)}
      # include in cache
     # colnames(Cache[[grep(DB,names(Cache),value = T)]])%in%colnames(db)
     # if(!all(!is.na(db))){print(DB%in%names(grads))}

      DBCache<-Cache[[grep(DB,names(Cache),value = T)]]
      DBCache[if(!is.null(colnames(db))){colnames(db)}else{names(db)},]<-db
      #Cache[[grep(DB,names(Cache),value = T)]][,colnames(db)]<-db

      return(DBCache)
    })
      names(dbState)<-DBi

      # to change when no LSTM
      if(LSTM){
        Addings<-c("dXt","dc_next")
      } else {
        Addings<-c("dXt")
      }
      #DX="dc_next"

      dXState<-lapply(Addings,function(DX){
        dx<-grads[[DX]]
        DXCache<-Cache[[grep(DX,names(Cache),value = T)]]
        DXCache[,if(!is.null(colnames(dx))){colnames(dx)}else{names(dx)}]<-dx
        return(DXCache)
      })
      names(dXState)<-Addings
      # gradient clipping?

      return(c(dWhState,dbState,dXState))
  }

  ##########
  # simulations
    Simulation<-function(NETexplore=NETexplore,  LSTM=LSTM, gradClipping=gradClipping,
                          Mode=Mode, MinSteps=MinSteps,Y=Y, NAMES=rownames(Y),
                         FixNodes=FixNodes, TotiState=TotiState){

      # names species
      Species<-union(NETexplore$target_hgnc,NETexplore$source_hgnc)
      Nspecies<-length(Species)
      OutNet<-unique(NETexplore$target_hgnc[NETexplore$Layer==max(NETexplore$Layer)])
      Out<-gsub("_.*","", OutNet )
      InNet<-unique(NETexplore$source_hgnc[NETexplore$Layer==1])
      In<-unique(NETexplore$source_hgnc[NETexplore$Inputs]) # check


      ##### back prop
      # initiate
      #list(dWho=dWho,dbo=dbo,dWhi=dWhi,dbi=dbi, dWhf=dWhf,dbf=dbf, dWhg=dWhg, dbg=dbg)
      NamesInteraction<-paste(NETexplore[,1],NETexplore[,3],sep = ".")

      if(LSTM){

        Cache<-list(
          dWho=data.frame(A=rep(1,length(NamesInteraction)),row.names = NamesInteraction),
          dbo=data.frame(A=rep(1,length(Species)),row.names = Species),
          dWhi=data.frame(A=rep(1,length(NamesInteraction)),row.names = NamesInteraction),
          dbi=data.frame(A=rep(1,length(Species)),row.names = Species),
          dWhf=data.frame(A=rep(1,length(NamesInteraction)),row.names = NamesInteraction),
          dbf=data.frame(A=rep(1,length(Species)),row.names = Species),
          dWhg=data.frame(A=rep(1,length(NamesInteraction)),row.names = NamesInteraction),
          dbg=data.frame(A=rep(1,length(Species)),row.names = Species),
          #dc_next=t(data.frame(A=derivfun(Cactivity),row.names = Species)),
          dc_next=derivfun(TotiState[,"Ct",,dim(TotiState)[4]]),
          dXt=derivfun(TotiState[,"A",,dim(TotiState)[4]]))

        if(dim(TotiState)[1]==1){
          Cache$dc_next<-t(Cache$dc_next)
          Cache$dXt<-t(Cache$dXt)
        }
             # bdg=  dZgStates=t(data.frame(A=rep(1,length(Species)),row.names = Species))
        #)
        #lapply(Cache,dim)
        #names(Cache)%in%names(grads)
        # this could be the measure for stability detection
        #TotdWState should have ncol=nrow(NETexplore) and nrow= nb of iterations
        TotdWState<-t(Cache$dWhg)
        TotdZState<-t(Cache$dbg) # because dbg=dZg, and need to keep track of Z I think... to check
        TotdbState<-t(Cache$dbg)

        #  Totdc_next<-Cache$dc_next
      #  TotdXt<-Cache$dXt
      #  hist(TotdbState[2,])
      #  TotdbState[2,][TotdbState[2,]==1]
      } else{
        ### 30/1/19
        Cache<-list(
          dW=data.frame(A=rep(1,length(NamesInteraction)),row.names = NamesInteraction),
          db=data.frame(A=rep(1,length(Species)),row.names = Species),
          dXt=derivfun(TotiState[,"A",,dim(TotiState)[4]]))

        TotdWState<-t(Cache$dW)
        TotdZState<-t(Cache$db) # because dbg=dZg, and need to keep track of Z I think
        TotdbState<-t(Cache$db)

        if(FALSE){
        dWStates<-rnorm(nrow(NETexplore),mean = 0,sd = 10)
        dWStates<-rep(0,nrow(NETexplore))
        names(dWStates)<-NamesInteraction

        dbStates<-rep(0,length(Species))
        names(dbStates)<-Species
        dZStates<-derivfun(Zactivity)

        TotdWState<-as.data.frame(t(dWStates),row.names = "")
        TotdbState<-as.data.frame(t(dbStates),row.names = "")
        TotdZState<-as.data.frame(t(dZStates),row.names = "")
        }
      }

      ###############
      # first step:
      # 1. derivative of the cost <-> inputs
      # for 1 patient and each layers

      #dX<-2*(t(Y[NAMES,])-TotiState[,"A",NAMES,dim(TotiState)[4]])#/NPat
      dX<-2*(TotiState[,"A",NAMES,dim(TotiState)[4]]-t(Y) )#/NPat

      # with regularization
      if(1==alpha){
        Wcost<-BuildNetMat(NETall = NETexplore, Layer = 1, Type = "Weights")

      #  dX<-2*(TotiState[,"A",NAMES,dim(TotiState)[4]]-t(Y)) +
      #    lambda*diag(sign(Wcost)) #/dim(TotiState)[1]
        #table(diag(sign(Wcost)))
          dX<-2*(TotiState[,"A",NAMES,dim(TotiState)[4]]-t(Y)) +
            lambda*diag(sign(Wcost)%*%t(Wcost)) #/dim(TotiState)[1]
      } else if (2==alpha){
        Wcost<-BuildNetMat(NETall = NETexplore, Layer = 1, Type = "Weights")
        # pb of dimensions if mat mul between 2Z(A-y), so picewise mult choosen, to check
        #dX<-2*(TotiState[,"Z",NAMES,dim(TotiState)[4]])*(TotiState[,"A",NAMES,dim(TotiState)[4]]-t(Y))+
        #  2*lambda*diag(Wcost) #/dim(TotiState)[1]
          dX<-2*(TotiState[,"A",NAMES,dim(TotiState)[4]]-t(Y)) +
            lambda*diag(Wcost)*2 #/dim(TotiState)[1]
      }

      if("LAYER"%in%Mode){

      #1st step : update measured species
      grads <- UpdateLSTMgrad(dX = dX, NAMESnode = NAMES, NETexplore = NETexplore,
                              Layer=NULL, TotiState = TotiState, LSTM=LSTM)

    #                          Xactivity=Xactivity,Zactivity=Zactivity, Cactivity=Cactivity,
    #                          LSTMvalues = LSTMvalues)
      Cache<-StoreLSTMgrad(Cache=Cache, NamesInteraction=NamesInteraction,grads=grads)

      if(LSTM){
        TotdWState<-rbind(TotdWState,t(Cache$dWhg))
        TotdZState<-rbind(TotdZState,t(Cache$dbg)) # ...
        TotdbState<-rbind(TotdbState,t(Cache$dbg))
      } else {
        TotdWState<-rbind(TotdWState,t(Cache$dW))
        TotdZState<-rbind(TotdZState,t(Cache$db)) # ...
        TotdbState<-rbind(TotdbState,t(Cache$db))
        }
      #dim(TotdWState)
      #2nd step : diffuse the gradients for others, per layers
      CHECK<-TRUE
#     NLayer<-1
      # begin with measured layer
      #NLayer<-min(NETexplore$Layer[NETexplore$source_hgnc%in%NAMES])

      # begin with first layer that contains non-measured species
    NLayer<-min(NETexplore$Layer[which(!NETexplore$source_hgnc%in%NAMES)])
#      ULayer<-sort(unique(NETexplore$Layer))
      ULayer<-c(seq(NLayer,max(NETexplore$Layer)),seq(1,NLayer-1))
      NLayer<-1

      MonitorSpecies<-NAMES
      SetCut<-max(3,round(log(MinSteps)))
      Round=1

  #    print("2nd step")
      # do at least 2 runs ( or 1) of simulationsn of each layers
      while(ifelse(nrow(TotdWState)<MinSteps,TRUE,CHECK)){


        # Nlayer+1 updated at the end of the while loop
        if(NLayer>length(ULayer)){
          NLayer<-1
          Round<-Round+1
        }
        # print(paste("Iter layer",NLayer))

        # if(ULayer[NLayer]==0){
        #  Species<-unique(NETexplore$source_hgnc[NETexplore$Layer==2]) # source first then targets
        #} else {

        # in LAYER mode, no FixeNode in the back simulation strategy
        NewSpecies<-unique(NETexplore$target_hgnc[NETexplore$Layer==ULayer[NLayer]]) # target

        ####
        if(MinSteps==1){
          NewSpecies<-NewSpecies[!NewSpecies%in%MonitorSpecies]
          if(length(NewSpecies)==0){
            break() ########## update only once!!!
          }
        }

        # species selection for first round for the layer based simulation
        if(Round==1){
        # for first round, discard species already updated in 1st step
          NewSpecies<-NewSpecies[!NewSpecies%in%NAMES]
          # if all species already updated,
          if(length(NewSpecies)==0){
            NLayer<-1
            Round<-Round+1
            next
            #break() ########## update only once!!!
            }
        }
        #}

        MonitorSpecies<-c(MonitorSpecies,NewSpecies)

        #length(NewSpecies)
        # adjuction matrix
      W<-BuildNetMat(NETall = NETexplore, Layer = ULayer[NLayer],Type = "Weights") # weights for this layer

      #  NewSpecies<-NewSpecies[match(colnames(W),NewSpecies)] #==colnames(W[NAMESnode,])
      #  if(is.null(NewSpecies)){ print("NewSpecies problem") }

    #    if(LSTM){

          # prepare the next step for grad diffusion by selecting connected nodes
          #dX_next<-t(TotdXt[nrow(TotdXt),rownames(W)])
          #dX_next<-(Cache$dXt[,rownames(W)])
#          dX_next<-(Cache$dXt[,colnames(W)])
          dX_next<-(Cache$dXt[,NewSpecies,drop=F])

          #  dim(dX_next)
          # dc_next<-dX_next*derivfun(Cactivity[colnames(dX_next)])*unlist(LSTMvalues["OutputGate",colnames(dX_next)])
        #  dc_next <- as.vector(unlist(LSTMvalues["ForgetGate",colnames(dX_next)]) * dc_next)
        #  names(dc_next)<-colnames(dX_next)

          # run the gradient computation

          # update NAMES  with dx_next????
          #grads <- UpdateLSTMgrad(dX = dX_next, NAMESnode = NewSpecies, Xactivity=Xactivity,Zactivity=Zactivity, Cactivity=Cactivity,
          #                        LSTMvalues = LSTMvalues, NETexplore = NETexplore, Layer=NULL) #ULayer[NLayer]

          grads <- UpdateLSTMgrad(dX = dX_next, NAMESnode = NewSpecies,
                                  NETexplore = NETexplore, Layer=NULL,
                                  TotiState = TotiState,  LSTM=LSTM) #ULayer[NLayer]

          #lapply(grads,dim)
          Cache<-StoreLSTMgrad(Cache=Cache, NamesInteraction=NamesInteraction,grads = grads)

          #lapply(Cache,dim)
          if(LSTM){

            # store it
            TotdWState<-rbind(TotdWState,t(Cache$dWhg))
            TotdZState<-rbind(TotdZState,t(Cache$dbg))
            TotdbState<-rbind(TotdbState,t(Cache$dbg))
        # matplot(TotdWState,type='l')
          # TotdWState[3,][TotdWState[3,]!=0.1]
         # Totdc_next<-rbind(Totdc_next,Cache$dc_next)
        #  TotdXt<-rbind(TotdXt,Cache$dXt)

        } else {
          TotdWState<-rbind(TotdWState,t(Cache$dW))
          TotdZState<-rbind(TotdZState,t(Cache$db)) # ...
          TotdbState<-rbind(TotdbState,t(Cache$db))

        }

      if(FALSE){
        W<-BuildNetMat(NETall = NETexplore, Layer = ULayer[NLayer],Type = "Weights") # weights for this layer
      #  dim((W))
      #  dim(as.matrix(dZStates[colnames(W)]))
      #  dim(as.matrix(TotdZState[nrow(TotdZState),colnames(W)]))

        print('dx step')
        dim(W);dim(dZStates[,rownames(W)])
        dX<-as.vector(t(W)%*%dZStates[,rownames(W)])  # size [nb genes, nb patients] # how dim fits?
     #   length(dX)
        dX<-t(t(W)%*%as.matrix(TotdZState[nrow(TotdZState),colnames(W)]))  # size [nb genes, nb patients]

#        dX<-t(t(W[NewSpecies,,drop=F])%*%TotdZState[nrow(TotdZState),NewSpecies,drop=F])  # size [nb genes, nb patients]

           print('dz step')

       #   length(dX)
      #    length(Zactivity[colnames(W)])
          dZ<-dX*derivfun(Zactivity[colnames(W)]) # [nb outputs, nb patients] #
          # dZ<-dX*derivfun(Zactivity[,NewSpecies,drop=F]) # [nb outputs, nb patients] #

      #    dim(as.matrix(dZ))
      #    dim(t(Xactivity[rownames(W)]))
          dW<-as.matrix(dZ)%*%t(Xactivity[rownames(W)]) # Xactivity at stable state for target species!
      #    dim(dW)
          dW<-reshape2::melt(dW)
          Rnames<-paste(dW[,2],dW[,1],sep = ".")
          # Rnames[Rnames%in%names(dWStates)]
          dW<-dW[!duplicated(Rnames),]
          rownames(dW)<-paste(dW[,2],dW[,1],sep = ".")
          dW<-dW[rownames(dW)%in%names(dWStates),]
          dWStates[rownames(dW)]<-dW[,3]
          TotdWState<-rbind(TotdWState,dWStates) #!

          print('db step')

          db<-dZ
          dbStates[names(db)]<-dZ #! for real inputs
          TotdbState<-rbind(TotdbState,dbStates) #!
        }

        #MonitorSpecies<-c(MonitorSpecies,NewSpecies)

        # stoping rule
        # it is where you can define the length of the simulation
        if(nrow(TotdWState)>100){break}
        if(MinSteps==1){
          if(NLayer==length(ULayer)){break}
          #  if(all(Species[!Species%in%FixNodes]%in%unique(MonitorSpecies))){break}
          if(all(Species%in%unique(MonitorSpecies))){break}
        }

          # if stabilized (to keep?)
        ROWAttractors<-apply(round(TotdWState,3),1,function(y)paste(y,collapse = ""))

        if(length(unique(tail(ROWAttractors,MinSteps)))<3){
          # if last iState is in the list of the last 3 iStates, returns CHECK=F
          #CHECK<-!ROWAttractors[length(ROWAttractors)]%in%unique(tail(ROWAttractors[-length(ROWAttractors)],3))
          CHECK<-1!=length(unique(tail(ROWAttractors[-length(ROWAttractors)], SetCut)))

        }else {
          CHECK<-TRUE
        }

        # do next layer
        NLayer<-NLayer+1

      } # end of while loop

      #  nrow(TotdWState)
      #matplot(TotdZState,type='l')
#      matplot(TotdWState,type='l',main="dWeights")
  #    matplot(TotdbState,type='l')

    } else if("ASYNC"%in%Mode){
  #    print("1st step")
      #1st step : update measured species
      #grads <- UpdateLSTMgrad(dX = dX, NAMESnode = NAMES,
      #                        Xactivity=Xactivity,Zactivity=Zactivity, Cactivity=Cactivity,
      #                        LSTMvalues = LSTMvalues,NETexplore = NETexplore, Layer=NULL)
      grads <- UpdateLSTMgrad(dX = dX, NAMESnode = NAMES, NETexplore = NETexplore, Layer=NULL,
                              TotiState = TotiState,  LSTM=LSTM)
      Cache<-StoreLSTMgrad(Cache=Cache, NamesInteraction=NamesInteraction,grads=grads)


      if(LSTM){
        TotdWState<-rbind(TotdWState,t(Cache$dWhg))
        TotdZState<-rbind(TotdZState,t(Cache$dbg)) # ...
        TotdbState<-rbind(TotdbState,t(Cache$dbg))
      } else {
        TotdWState<-rbind(TotdWState,t(Cache$dW))
        TotdZState<-rbind(TotdZState,t(Cache$db)) # ...
        TotdbState<-rbind(TotdbState,t(Cache$db))
      }

        #2nd step : diffuse the gradients for others asynchronously
        #    NewSpecies<-sample(Species,1)
    #  NewSpecies<-sample(In,1)

     # NewSpecies<-Species[!Species%in%NAMES]

      ### take the Species following species used to compute dx at 1st step:
      NewSpecies<-unique(NETexplore[NETexplore$source_hgnc%in%NAMES,"target_hgnc"])

    #  # put the species not measured, therefore not updated yet, upfront:
    #  NewSpecies<-NewSpecies[c(which(!NewSpecies%in%NAMES),which(NewSpecies%in%NAMES))]
      # for first step, take only the species not measured, therefore not updated yet
      NewSpecies<-NewSpecies[which(!NewSpecies%in%NAMES)]

      # species selection for first round for the asynchronous based simulation
      # if the species haven't been updated already in 1st steps, then loop:
      if(!all(Species%in%NAMES)|MinSteps>1){

 #       print("2nd step")

      MonitorSpecies<-c(NAMES,NewSpecies)
      SetCut<-max(3,round(log(MinSteps)))
      CHECK<-TRUE

        while(ifelse(nrow(TotdWState)<MinSteps,TRUE,CHECK)){
          #Yspecies<-NewSpecies[1]
          # no loop is ok when only 1 species measured, eg "output" (& is <-> Layer)
          # but not if many species measured (& every species measured <-> sync)
   #       Yspecies="PP2A"
   for(Yspecies in NewSpecies){
            #if(LSTM){

              # prepare the next step for grad diffusion by selecting connected nodes
              #dX_next<-as.matrix(TotdXt[nrow(TotdXt),Yspecies])
              #dX_next<-(Cache$dXt[,NewSpecies,drop=F])
              dX_next<-(Cache$dXt[,Yspecies,drop=F])
              #colnames(dX_next)<-NewSpecies

              #dc_next<-dX_next*derivfun(Cactivity[colnames(dX_next)])*unlist(LSTMvalues["OutputGate",colnames(dX_next)])
              #dc_next <- as.vector(unlist(LSTMvalues["ForgetGate",colnames(dX_next)]) * dc_next)
              #names(dc_next)<-colnames(dX_next)

              grads <- UpdateLSTMgrad(dX = dX_next, NAMES=Yspecies,NETexplore = NETexplore, Layer=NULL,
                                      TotiState = TotiState,  LSTM=LSTM)
              Cache<-StoreLSTMgrad(Cache=Cache, NamesInteraction=NamesInteraction,grads = grads)

              if(LSTM){
                TotdWState<-rbind(TotdWState,t(Cache$dWhg))
                TotdZState<-rbind(TotdZState,t(Cache$dbg)) # ...
                TotdbState<-rbind(TotdbState,t(Cache$dbg))
              } else {
                TotdWState<-rbind(TotdWState,t(Cache$dW))
                TotdZState<-rbind(TotdZState,t(Cache$db)) # ...
                TotdbState<-rbind(TotdbState,t(Cache$db))
              }


          if(FALSE){
              TabAND<-NETexplore[NETexplore$target_hgnc%in%Yspecies,]
            #function (a*dZ)/m  *  g'(Z)
            # length(TabAND[,"Weights"])
            #  length(dZStates[TabAND[,1]])
            # TabAND[,1:3]
            #TotdZState[nrow(TotdZState),TabAND[,1]]
            dX<-as.matrix(TabAND[,"Weights"]%*%TotdZState[nrow(TotdZState),TabAND[,1]])  # size [nb genes, nb patients]
            colnames(dX)<-Yspecies

              dZ<-dX*derivfun(Zactivity[unique(TabAND[,3])]) # [nb outputs, nb patients] #

              dZStates[colnames(dZ)]<-dZ
              TotdZState<-rbind(TotdZState,dZStates) #!

            #  dim(as.matrix(dZ))
            #  dim(t(Xactivity[TabAND[,3]]))

              #dW<-as.matrix(dZ)%*%t(Xactivity[NETall[NETall$source_hgnc%in%NewSpecies,1]]) # Xactivity at stable state for target species!
              dW<-as.matrix(dZ)%*%t(Xactivity[TabAND$source_hgnc]) # Xactivity at stable state for target species!
              # dW
              #  rownames(dW)%in%NewSpecies
              #  colnames(dW)%in%NewSpecies
              dW<-reshape2::melt(dW)
              Rnames<-paste(dW[,2],dW[,1],sep = ".")
              # Rnames[Rnames%in%names(dWStates)]
              dW<-dW[!duplicated(Rnames),]
              rownames(dW)<-paste(dW[,2],dW[,1],sep = ".")

              dW<-dW[rownames(dW)%in%names(dWStates),]
              dWStates[rownames(dW)]<-dW[,3]
              TotdWState<-rbind(TotdWState,dWStates) #!

              db<-dZ
              dbStates[names(db)]<-dZ #! for real inputs
              TotdbState<-rbind(TotdbState,dbStates) #!
            }

     # adding 7/11/18
          MonitorSpecies<-c(MonitorSpecies,Yspecies)
          if(MinSteps==1){
            #  if(all(Species[!Species%in%FixNodes]%in%unique(MonitorSpecies))){break}
            if(all(Species%in%unique(MonitorSpecies))){break}
          }

         } # end of for loop

          #matplot((TotdXt),xlim = c(2,nrow(TotdXt)),type='l')
     #     matplot((TotdWState),type='l')
          #matplot((TotdZState),type='l')
         # matplot((Totdc_next),type='l')

          # if other pass, diffuse to all species if needed
          NewSpecies<-unique(NETexplore[NETexplore$source_hgnc%in%NewSpecies,"target_hgnc"])

          #NewSpecies<-NewSpecies[!NewSpecies%in%FixNodes]
          if(length(NewSpecies)==0){NewSpecies<-sample(Species,1)}
        # if the last of monitorspecies are the same as NewSpecies, change species <-> unstuck
        if(all(NewSpecies%in%tail(MonitorSpecies,length(NewSpecies)))){
#          NewSpecies<-sample(Species,1)
          NewSpecies<-sample(Species,1)
        }

#        MonitorSpecies<-c(MonitorSpecies,NewSpecies)

        # stoping rule
        if(nrow(TotdWState)>100){break}
        if(MinSteps==1){
        #  if(nrow(TotdWState)>=max(NETexplore$Layer)){break}
          if(all(Species%in%unique(MonitorSpecies))){break}
        }

        ROWAttractors<-apply(round(TotdWState,3),1,function(y)paste(y,collapse = ""))
#        if(length(unique(tail(ROWAttractors,Nspecies*2)))==1){

        if(length(unique(tail(ROWAttractors,MinSteps)))<3){
         # CHECK<-!ROWAttractors[length(ROWAttractors)]%in%unique(tail(ROWAttractors[-length(ROWAttractors)],3))
        #  CHECK<-!all(unique(tail(ROWAttractors[-length(ROWAttractors)],round(MinSteps/5) ))%in%ROWAttractors[length(ROWAttractors)])

          CHECK<-1!=length(unique(tail(ROWAttractors[-length(ROWAttractors)], SetCut)))

          }else{
          CHECK<-TRUE
        }

        } # end of while loop

      } # end of if()

      #  nrow(TotdWState)
  #    matplot(TotdZState,type='l')
 #     matplot(TotdWState,type='l')
  #    matplot(TotdbState,type='l')

      } else if("SYNC"%in%Mode){

        ######################
        ###### NO USE ##########
        ######################

        # in SYNC mode, no Fixnode in the updates
        #NewSpecies<-sample(names(iStates),1)
        #MonitorSpecies<-NewSpecies
        CHECK<-TRUE
        SetCut<-max(3,round(log(MinSteps)))

          while(ifelse(nrow(TotiState)<MinSteps,TRUE,CHECK)){

            W<-BuildNetMat(NETall = NETexplore, Layer =NULL,Type = "Weights") # weights for this layer
            dim(t(W))
            dim(as.matrix(dZStates[colnames(W)]))

            dX<-as.vector(t(W)%*%dZStates[rownames(W)])  # size [nb genes, nb patients]
            length(dX)
            length(Zactivity[colnames(W)])
            dZ<-dX*derivfun(Zactivity[colnames(W)]) # [nb outputs, nb patients] #
            dim(as.matrix(dZ))
            dim(t(Xactivity[rownames(W)]))
            dW<-as.matrix(dZ)%*%t(Xactivity[rownames(W)]) # Xactivity at stable state for target species!
            dim(dW)
            dW<-reshape2::melt(dW)
            Rnames<-paste(dW[,2],dW[,1],sep = ".")
            # Rnames[Rnames%in%names(dWStates)]
            dW<-dW[!duplicated(Rnames),]
            rownames(dW)<-paste(dW[,2],dW[,1],sep = ".")
            dW<-dW[rownames(dW)%in%names(dWStates),]
            dWStates[rownames(dW)]<-dW[,3]
            TotdWState<-rbind(TotdWState,dWStates) #!

            db<-dZ
            dbStates[names(db)]<-dZ #! for real inputs
            TotdbState<-rbind(TotdbState,dbStates) #!

            # stoping rule
      ROWAttractors<-apply(TotdWState,1,function(y)paste(y,collapse = ""))

      #unique(tail(ROWAttractors,log2(Nspecies)))

      if(length(unique(tail(ROWAttractors,MinSteps)))<3){
        #CHECK<-!ROWAttractors[length(ROWAttractors)]%in%unique(tail(ROWAttractors[-length(ROWAttractors)],3 ))
        CHECK<-1!=length(unique(tail(ROWAttractors[-length(ROWAttractors)], SetCut)))

      }else {
        CHECK<-TRUE
        }
      if(nrow(TotdWState)>500){break}
        } # #end of while loop
      } # end of SYNC

#      if(LSTM){
        return(c(list(TotdZState=TotdZState,TotdWState=TotdWState,TotdbState=TotdbState),Cache) )
 #     } else{
#        return(list(TotdZState=TotdZState,TotdWState=TotdWState,TotdbState=TotdbState))
#      }

  }  # end of simulation function


    # vectorized matrix calculations
    LAttractorsAZ<-Simulation(NETexplore=NETexplore, TotiState = TotiState,
                 Y = Y, Mode=Mode, MinSteps = MinSteps,
                 NAMES=rownames(Y), FixNodes=FixNodes, #R_batch=R_batch[Npat,,drop=F],
                 LSTM=LSTM)

  #  LAttractorsAZ<-Simulation(NETexplore=NETexplore, Xactivity=Xactivity,Zactivity=Zactivity,
  #             Cactivity = Cactivity, Mode = Mode, FixNodes=FixNodes,
  #             MinSteps = MinSteps, Y = Y, NAMES=rownames(Y), #R_batch=R_batch,
  #             LSTM=LSTM, LSTMvalues=LSTMvalues)
    LAttractorsAZ<-list(LAttractorsAZ)


   return(LAttractorsAZ)
}

