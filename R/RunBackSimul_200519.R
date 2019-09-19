#' Runs the back simulation or back propagation of the gradients
#'
#' \code{RunBackSimul} process minibatches and run :
#' 1. simulations
#' 2. compute error on outputs
#' 3. flip the net
#' 4. simulation of the derivate of the functions
#' 5. correct the parameters with gradients
#'
#' @param NETall data frame. The network in format (similar to .sif with) columns "source", "interaction", "target", additionally storing the architecture (column "Layer"), the parameters and optionnaly Boolean logic.
#' @param y numeric or matrix. labels
#' @param MUT list. Corresponds to mutations with names and corresponding values (list of vectors). If MUT=NULL, a wild type *AMoNet* is simulated.. Patients ordered the same as y.
#' @param treatmt list. The same for than newMUT, with treatments' targets and corresponding values (list of vectors).
#' @param Init matrix. Used to set some of the initial states (iStates). To force some initial states during learning set adaptive_iStates=F. Init=NULL otherwise.
#' @param iStates matrix. Initial states of the whole network model. Either set from data or randomized.
#' @param Optimizer character. Options are "Adam","Momentum","RMSprop" or NULL. NULL corresponds to gradient descent and can be used with batch learning. Use of optimizers are recommended to mini batch learning.
#' @param beta1 numeric. Momentum term. Adam is defined by Momentum / RMSprop.
#' @param beta2 numeric. RMSprop term. Adam is defined by Momentum / RMSprop.
#' @param MiniBatch integer. nb of sample per MiniBatch. For batch learning, use MiniBatch = total nb of examples
#' @param iteration integer. Number of iterations, corresponding to epochs, i.e. one path through the full data.
#' @param Visualize character or integer. Options : \code{c(1,"Gates","Output")} Real time visualization of learning, requires selecting the visualization target: "Gates","Output", or all nodes activities of a patients by selecting an integer (a patient's number) or an ID. Default is no visualization: set to \code{NULL}, .
#' @param alpha numeric. Regularization term. Set to 1 for L1N (Lasso) regularization, set to 2 for L2N (Ridge) regularization, between 1 and 2 with lambda 1 for L1N + L2N (Elastic Net) regularization (penalization parameter is alpha).
#' @param lambda numeric. the lambda value for Lasso, Ridge or Elastic Net regularization. Set to 0 corresponds to no regularization. When alpha between 1 and 2 (Elastic Net regularization), lambda should be set to 1 and penalization parameter is alpha.
#' @param Logic character. The update function used for the simulations. Choices are: "Sigmoid" (default), "Boolean","tanh" or "ReLU".
#' @param Mode character. The type of simulation to perform. Choice are "LAYER" (layer based, default), "ASYNC" (assynchronous update), "SYNC" (synchronous update)
#' @param MinStepsForward numeric. Number of simulation steps. Default to 5. Ideally set the number of simulation steps to achieve stable states of *AMoNet* activities.
#' @param MinStepsBackward numeric. Number of back-simulation steps. Default to 1, as in a standard back-propagation of gradient in neural nets.
#' @param FixNodes Character. Freeze parameters of selected nodes.
#' @param Parallelize boolean. Should the simulation run in parallel (TRUE by default)? Not recommended for simulations of less than 10 examples.
#' @param no_cores numeric. If Parallel=TRUE, set the number of cores to parallelize on. Default is 4. Can detect and set to available no_cores if inferior to user defined no_cores. to
#' @param gradClipping nunmeric. Maximum values allowed for gradients otherwise clipped. Gradient scaling in dev.
#' @param LSTM boolean. Should a Long-Short Term Memory (LSTM) unit architecture be used?
#' @param ValMut numeric. Multicplicative factor for the effect of mutations. Default to 50.
#'
#'
#' @details
#' Preferably use the \code{train()} function from AMoNet package to run the training.
#'
#' Default hyperparameters used for training are stored in \code{Default}.
#'
#' For character objects in hyperparameters (for eg Optimizer can be "Adam", "RMSprop, "Momentum" or NULL), options are stored in \code{Boundaries}
#'
#' \code{MUT} and \code{treatmt} parameters can be generated for a matrix (samples * genes) filled with values and \code{NAs}, using the \code{MutMatToList()} functcion.
#'
#' \code{MinStepsForward} is corresponding to generic argument \code{nsim} in S3 function \code{simulate()}.
#'
#' Setting argument \code{MinStepsBackward} > 1 can result in vanishing gradient. Nevertheless, using \code{LSTM=TRUE} in this context may help.
#'
#' @return a list with the network dataframe NETall, and for each epoch:
#' the training Cost, the network weights (NETallList) and the simulation activities (NETallActivity)
#'
#' @examples
#' \dontrun{
#' #see the \code{?train()} function from AMoNet package
#' }
RunBackSimul<-function(NETall=NETall, y=y, MUT=NULL, treatmt=NULL,Init=NULL, iStates=NULL, Ct=NULL, MiniBatch=Default$MiniBatch,
                       Optimizer=Default$Optimizer, beta1=Default$beta1, beta2=Default$beta2,
                       iteration=Default$iteration, learning_rate=Default$learningrate,
                       adaptive_iStates=Default$adaptive_iStates, FixNodes=NULL,
                       Parallelize=Default$Parallelize, no_cores=Default$no_cores,
                       Logic = Default$Logic, Mode = Default$Mode, ModeBack=Default$Mode,
                       MinStepsForward = Default$MinStepsForward,MinStepsBackward=Default$MinStepsBackward,
                       LSTM=Default$LSTM, gradClipping=Default$gradClipping, LearningRateDecay=Default$LearningRateDecay,
                       ValMut=Default$ValMut, PDF=F, GIF=F, NameProj="AMoNet",
                       Visualize=Default$Visualize, alpha=Default$alpha, lambda=Default$lambda ){

  ### call
  Call<-list(Logic=Logic, Mode=Mode, MUT=!is.null(MUT), Init=!is.null(Init))

  ### checks
  if(!is.null(MUT)&!is.null(Init)){
    if(length(MUT)!=nrow(Init)){stop("MUT list and Init table/matrices must have same number of examples (eq patients)")}
  }
  if(!is.null(MUT)&!is.null(treatmt)){
    if(length(MUT)!=length(treatmt)){stop("MUT and treatmt lists must have same length (eq patients)")}
  }
  if(!is.null(MUT)){
    if(length(MUT)!=ncol(y)){stop("length MUT and ncol y table/matrices should be the same (eq patients)")}
  }
  if(!is.null(Init)){
    if(nrow(Init)!=ncol(y)){stop("nrow Init and ncol y table/matrices should be the same (eq patients)")}
  }
  if(!is.null(treatmt)){
    if(length(treatmt)!=ncol(y)){stop("length MUT and ncol y table/matrices should be the same (eq patients)")}
  }

 # if(file.exists(paste(getwd(),"/data/CancerGeneCensusCOSMIC.csv",sep = ""))){
 #   CGS<-read.csv(paste(getwd(),"/data/CancerGeneCensusCOSMIC.csv",sep = ""),stringsAsFactors = F)
#  } else {
#    print("If you simulate mutations, be sure to have the Cancer Gene Census (from Cosmic) table in /data/ directory")
#  }

  if(is.null(iStates)|is.null(nrow(iStates))){
    print("Be sure to provide a matrix of initial states for each species = iStates
          for eg: simulate once with multiple random initial states to compute the attractors and use it as iStates")
    stop()
  }

  #########
  # settings
  # store initial learning_rate
  learning_rate_init<-learning_rate

  # store species
  OutputsNETall<-unique(NETall$target_hgnc[NETall$Output])
  InputsNETall<-unique(NETall$source_hgnc[NETall$Layer%in%1])
  Species<-union(NETall$target_hgnc,NETall$source_hgnc)

  # compute minibatches

  NMiniBatch<-ceiling(ncol(y)/MiniBatch)
  if(NMiniBatch==1){print(paste("Batch learning with",ncol(y),"examples" ))}

  EpochMiniBatch<-list()

  for(NMiniBatch in seq(ceiling(ncol(y)/MiniBatch))){
    if(NMiniBatch == ceiling(ncol(y)/MiniBatch) ){ # last one
      EpochMiniBatch[[NMiniBatch]]<- seq((NMiniBatch*MiniBatch-MiniBatch)+1,ncol(y))
    } else {
      EpochMiniBatch[[NMiniBatch]]<- seq((NMiniBatch*MiniBatch-MiniBatch)+1,NMiniBatch*MiniBatch)
    }
  }

  # initialise cache
  NETallList<-list()
  NETallList[[1]]<-NETall
  NETallActivity<-list()
  NETallActivity[[1]]<-iStates

  COST<-data.frame()
  StopAdam=NULL

  # Vizualizations

  if(is.null(Visualize)){
    print("no visualization selected")
    PDF=FALSE;GIF=FALSE
  }
  if(PDF&GIF){
    print("Select either PDF or GIF: PDF automatically selected")
    GIF=F
  }

  if(PDF){
    Latt<-length(list.files(paste(getwd(),"/tmp/",sep = ""),pattern =  NameProj))
    if(Latt>0){
      Latt<-max(as.numeric(gsub(".pdf|.gif","", gsub(NameProj,"", list.files(paste(getwd(),"/tmp/",sep = ""),pattern =  NameProj)))),na.rm = T)
    }
    pdf(paste(paste(getwd(),"/tmp/",sep = ""),NameProj,Latt+1,".pdf",sep = ""))
  }

  if(GIF){
    print("Run a gif on the type of vizualization selected")
    img <- magick::image_graph(600, 340, res = 60)
  }

  par(mar=c(5, 4, 4, 2) + 0.1)

  ############
  # run the optimization procedure
  # Loop over the iterations
  #iter=1
  for(iter in seq(iteration)){

    # shuffle minibatches
    if(NMiniBatch>1){
      ORD<-sample(seq(ncol(y)))
      y<-y[,ORD,drop=F]
      iStates<-iStates[ORD,]
      if(LSTM){
        Ct<-Ct[ORD,]
      }
    }

    #epoch=1
    for(epoch in seq(NMiniBatch)){

      if(iter==1&epoch==1){
        print("Waiting for first epoch...")
      }
      Tic<-Sys.time()

      # select subset of data for minibatches
      Y_epoch<-y[,EpochMiniBatch[[epoch]],drop=F]

      iStates_epoch<-iStates[EpochMiniBatch[[epoch]],Species,drop=F]
      if(LSTM){
        Ct_epoch<-Ct[EpochMiniBatch[[epoch]],Species,drop=F]
      } else {
        Ct_epoch=NULL
      }

      if( !is.null(MUT) ){
        MUT_epoch<-MUT[EpochMiniBatch[[epoch]]]
      } else {
        MUT_epoch<-NULL
      }

      if(!is.null(treatmt)){
        treatmt_epoch<-treatmt[EpochMiniBatch[[epoch]]]
      } else {
        treatmt_epoch<-NULL
      }

      ############
      #### do the simulation
      TotAttractors<-BoolSimul(NETall=NETall, Logic = Logic, Mode = Mode,
                               iStates=iStates_epoch,
                               Parallel = Parallelize, no_cores=no_cores,
                               MinSteps = MinStepsForward,
                               LSTM = LSTM, Ct=Ct_epoch,
                               MUT=MUT_epoch,treatmt = treatmt_epoch,
                               ValMut = ValMut)

      if(!all(rownames(Y_epoch)%in%names(TotAttractors[1,"A",,dim(TotAttractors)[4]]))){
        print("Pb with names to compute Costs")
        print(rownames(Y_epoch))
        print(names(TotAttractors[1,"A",,dim(TotAttractors)[4]]))
      }
      ######
      # compute cost with penalization

      W<-BuildNetMat(NETall = NETall[NETall$Output,],Layer = NULL, Type = "Weights")

      Cost <- mean( (t(Y_epoch)-TotAttractors[,"A",rownames(Y_epoch),dim(TotAttractors)[4]])^2 ) +
        (2-alpha)*lambda*norm(W,"1") +
        (alpha-1)*lambda*norm(W,"f")

 #   RealCost <- mean( (t(Y_epoch)-TotAttractors[,"A",rownames(Y_epoch),dim(TotAttractors)[4]])^2 )

      ### cost storing
      COST<-rbind(COST,Cost)

      # manage gradient clipping # to improve
      if(nrow(COST)>2&!is.null(gradClipping)){
        # slop calculated relatively to iteration: at the end of learning can flatten
        ModCost<-lm(as.numeric(t(tail(COST[,1],iter+2)))~seq(length(tail(COST[,1],iter+2)))) #which(COST[,1]%in%tail(COST[,1],iter+2)))#seq(length(tail(COST[,1],iter+1))))
        dCOST[nrow(COST)]<-as.numeric(ModCost$coeff[2]) # derivate of the cost, which allows to determine when to stop the leanring
        # should use it to modulate the level to wich you do gradclipping?!
        CLIP<-F
        # or check if slope to important
        CLIP<-ifelse(as.numeric(ModCost$coefficients[2]*iter*epoch+ModCost$coefficients[1])<=0.00001,T,CLIP)
        # or gradients too high
        CLIP<-ifelse(max(abs(NETall[,grep("Vd|dW|db",colnames(NETall))]))>gradClipping,T,CLIP)
        #  which(NETall[,grep("Vd",colnames(NETall))]>gradClipping,arr.ind = T)
        #  NETall[126,c(1,3)]
      } else{
        dCOST<-1
        CLIP<-F
      }

      ####################
      # do the back propagation <=> 1 step back simulation recommended

      LAttractorsAZ<-BackBoolSimul(NETall = NETall, Logic = Logic, Mode = ModeBack,
                                   MinSteps = MinStepsBackward, FixNodes=FixNodes,
                                   Y=Y_epoch, gradClipping = F, TotiState = TotAttractors,
                                   Parallelize = F, LSTM = LSTM, alpha = alpha, lambda = lambda)
#names(LAttractorsAZ[[1]])

      ###########
      # compute the parameters' update

      States<-lapply(grep("^dW|^db", names(LAttractorsAZ[[1]]),value = T),function(PARAM){
        State<-lapply(LAttractorsAZ,function(p){
          p[[PARAM]]
        })

        if(length(LAttractorsAZ)==ncol(Y_epoch)){
          State<-t(do.call('cbind',State))
          rownames(State)<-colnames(Y_epoch)
          State<-as.data.frame(colMeans(State))
        } else{
          State<-as.data.frame(State)
        }
        return(State)
      })

      names(States)<-grep("^dW|^db", names(LAttractorsAZ[[1]]),value = T)

      if(length(grep("g",names(States)))>0){
        names(States)[grep("g",names(States))]<-c("dWeights","dbterm")
      } else if (length(States)==2){
        names(States)<-c("dWeights","dbterm")
      }

      # nodes that are not updated:
      # either manually selected or that are not in the ground truth table
      # manualy <-> influence toward the node
      if(!is.null(FixNodes) ){
        #NODESFIX<-NETall[,1]%in%FixNodes|NETall[,3]%in%FixNodes
        NODESFIX<-NETall[,3]%in%FixNodes
      } else {
        NODESFIX<-rep(FALSE,nrow(NETall))
      }

      #if not in the ground truth table
      NoUpdate<-unlist(unique(lapply(States,function(S){
        rownames(S[which(S==1),,drop=F]) # initialize 1 in Cache in Back_BoolSimul121018
      })))

      NODESFIX<-NODESFIX|NETall[,3]%in%NoUpdate

      #################
      # run the descent/optimisation algorithms

      if(any(c("Adam","Momentum","RMSprop")%in%Optimizer) ){

        # 1. momentum
        for(VAdam in grep("VdW",colnames(NETall),value = T) ){
          NETall[,gsub("V","",VAdam)]<-States[[gsub("V","",VAdam)]]
          NETall[,VAdam]<-beta1*NETall[,VAdam]+(1-beta1)*NETall[,gsub("V","",VAdam)]
        }

        for(VAdam in grep("Vdb",colnames(NETall),value = T) ){
          NETall[,gsub("V","",VAdam)]<-States[[gsub("V","",VAdam)]][NETall$target_hgnc,]
          NETall[,VAdam]<-beta1*NETall[,VAdam]+(1-beta1)*NETall[,gsub("V","",VAdam)]
        }

        if("Momentum"%in%Optimizer){
          for(XAdam in gsub("Vd","", grep("Vd",colnames(NETall),value = T)) ){
            NETall[!NODESFIX,XAdam]<-NETall[!NODESFIX,XAdam]-learning_rate*NETall[!NODESFIX,paste("Vd",XAdam,sep = "")]
          }
        }

        # 2. RSMprop
        for(VAdam in grep("SdW",colnames(NETall),value = T) ){
          NETall[,gsub("S","",VAdam)]<-States[[gsub("S","",VAdam)]]
          NETall[,VAdam]<-beta2*NETall[,VAdam]+(1-beta2)*(NETall[,gsub("S","",VAdam)]^2)
          #NETall[,VAdam]<-NETall[,VAdam]/(1-beta2)
        }

        for(VAdam in grep("Sdb",colnames(NETall),value = T) ){
          NETall[,gsub("S","",VAdam)]<-States[[gsub("S","",VAdam)]][NETall$target_hgnc,]
          NETall[,VAdam]<-beta2*NETall[,VAdam]+((1-beta2)*(NETall[,gsub("S","",VAdam)]^2))
          #NETall[,VAdam]<-NETall[,VAdam]/(1-beta2)
        }

        if("RMSprop"%in%Optimizer){
          #XAdam = gsub("Sd","", grep("Sd",colnames(NETall),value = T))[1]
          for(XAdam in gsub("Sd","", grep("Sd",colnames(NETall),value = T)) ){
            # NETall[!NODESFIX,XAdam]<-NETall[!NODESFIX,XAdam]-learning_rate*NETall[!NODESFIX,paste("Vd",XAdam,sep = "")]
            NETall[!NODESFIX,XAdam]<-NETall[!NODESFIX,XAdam]-learning_rate*(NETall[!NODESFIX,paste("d",XAdam,sep = "")]/(sqrt(NETall[!NODESFIX,paste("Sd",XAdam,sep = "")])+1e-8))
          }
        }

        ############
        # grad clip Adam here. To improve
         if(CLIP){

          NETall[,grep("Vd",colnames(NETall))]<-1e-8
          NETall[,grep("Sd",colnames(NETall))]<-1e-8

          # with scaling ?
          if(FALSE){
            # SC="VdWeights"
            for(SC in c(grep("Vd",colnames(NETall),value = T),grep("Sd",colnames(NETall),value = T)) ){
              S<-scale.default(NETall[,SC])
              if(!any(is.na(S))){
                NETall[,SC]<-S[,1]
              }
            }
          }

        }

        # 3. Adam

       if("Adam"%in%Optimizer){
          for(XAdam in gsub("Vd","", grep("Vd",colnames(NETall),value = T)) ){
            NETall[!NODESFIX,XAdam]<-NETall[!NODESFIX,XAdam]-learning_rate*(NETall[!NODESFIX,paste("Vd",XAdam,sep = "")]/(sqrt(NETall[!NODESFIX,paste("Sd",XAdam,sep = "")])+1e-8 ) )
          }
        }


      } else { #if ! Adam|Momentum|RMSprop &LSTM

       NamesGrads<-grep("dW",colnames(NETall),value = T)
        NamesGrads<-gsub("d","", NamesGrads[NamesGrads%in%names(States)])
        for(WEI in NamesGrads ){
          NETall[!NODESFIX,paste("d",WEI,sep = "")]<-States[[paste("d",WEI,sep = "")]][!NODESFIX,]
          NETall[!NODESFIX,WEI]<-NETall[!NODESFIX,WEI]-learning_rate*NETall[!NODESFIX,paste("d",WEI,sep = "")]
        }

        NamesGrads<-grep("db",colnames(NETall),value = T)
        NamesGrads<-gsub("d","", NamesGrads[NamesGrads%in%names(States)])
        for(BTR in NamesGrads ){
          NETall[!NODESFIX,paste("d",BTR,sep = "")]<-States[[paste("d",BTR,sep = "")]][NETall$target_hgnc,][!NODESFIX]
          NETall[!NODESFIX,BTR]<-NETall[!NODESFIX,BTR]-learning_rate*NETall[!NODESFIX,paste("d",BTR,sep = "")]
        }

      }


      ##################
      # update and store:
      # net
      NETallList[[length(NETallList)+1]]<-NETall

      ################
      # activities
      # put stable states as initial states for next simulation and update in NETall

      iStates1<-iStates


      TotAttractors<-BoolSimul(NETall = NETall, Logic = Logic, Mode = Mode,
                               iStates = iStates_epoch, MinSteps = MinStepsForward,
                               LSTM = LSTM, Ct=Ct_epoch, MUT = MUT_epoch, treatmt = treatmt_epoch,
                               Parallel = Parallelize, no_cores=no_cores, ValMut = ValMut) #

      #### plot
      if(!any(c("Gates","Output")%in%Visualize)&!is.null(Visualize)){

        par(mfrow=c(1,2))
        MAIN<-paste("Pat", ifelse(is.numeric(Visualize),colnames(Y_epoch)[Visualize],Visualize) , "simulation")
        matplot(rbind(t(TotAttractors[Visualize,"iStates",,1]),t(TotAttractors[Visualize,"A",,])),type='l',
                main=MAIN,
                ylim=c(0,1), ylab="Species activities", xlab="Updates")
        plot(TotAttractors[Visualize,"A",rownames(Y_epoch),dim(TotAttractors)[4]],Y_epoch[,Visualize],
             main=MAIN,ylab="True",xlab="Pred",xlim=c(0,1),ylim=c(0,1))

      }


      if(dim(TotAttractors)[1]>1){
        iStates1[EpochMiniBatch[[epoch]],]<-TotAttractors[,"A",colnames(iStates),dim(TotAttractors)[4]]
      } else {
        iStates1[EpochMiniBatch[[epoch]],]<-t(TotAttractors[1,"A",colnames(iStates),dim(TotAttractors)[4]])
      }

      NETallActivity[[length(NETallList)]]<-iStates1

      if(adaptive_iStates){
        if(iter==1){

          iStates<-t(replicate(ncol(y), colMeans(TotAttractors[,"A",Species,dim(TotAttractors)[4]]) , simplify = TRUE))
          colnames(iStates)<-Species
          rownames(iStates)<-colnames(y)
          if(LSTM){
            Ct<-t(replicate(ncol(y), colMeans(TotAttractors[,"Ct",,dim(TotAttractors)[4]]) , simplify = TRUE))
            colnames(Ct)<-Species
            rownames(Ct)<-colnames(y)
          }
        } else {
          iStates<-iStates1
        }
        if(!is.null(Init)){
          for(i in seq(nrow(Init))){
            iStates[i,colnames(Init)]<-as.numeric(Init[i,]) # match the col order
          }
        }
        # update also the Ct activity
        if(LSTM){
          Ct[EpochMiniBatch[[epoch]],]<-TotAttractors[,"Ct",,dim(TotAttractors)[4]]
        }
      }
      rm(iStates1)

      ##############
      # representation of the gates for monitoring

      if(any(c("Gates","Output")%in%Visualize)){

        Vizlearning(GatesViz = "Gates"%in%Visualize&LSTM,
                    OutputViz = "Output"%in%Visualize, NETall = NETall,
                    Epoch=EpochMiniBatch[[epoch]],
                    learning_rate=learning_rate, NODESFIX = NODESFIX, y=y ,
                    NETallActivity = NETallActivity, COST=COST)

      }

      ###################
      #### learning rate decay
      if(!is.null(LearningRateDecay)){
        if('linear'%in%LearningRateDecay){
          #linear decay
          learning_rate=learning_rate/(1+learning_rate*epoch*iter)
        } else if('exponential'%in%LearningRateDecay){
          #exponential decay
          learning_rate=learning_rate*(0.95^(epoch*iter))
        }
      }

      ###################
      # print statistics for 1 epoch & estimate whole running duration
      if(iter==1&epoch==1){

        Toc<-Sys.time()
        Diff<-Toc-Tic

        total<-as.numeric(Diff)*iteration*NMiniBatch

        print("Statistics for optimization")
        if(Parallelize){

          #no_cores <- detectCores() # - 1
          #no_cores<-6
          print(paste("   Number of parallel cores:",no_cores))
        }
        print("   1 epoch took :"); print(Diff)
        print(paste("   Estimated time for", iteration, "iteration(s) and", NMiniBatch, "mini-batch(es):", total, units(Diff)))

        # create progress bar
        pb <- txtProgressBar(min = 0, max = total, style = 3)
        setTxtProgressBar(pb, as.numeric(Diff), title = "Progress status")

      } else {
        Toc<-Sys.time()
        Diff1<-Toc-Tic
        Diff<-Diff+Diff1
        setTxtProgressBar(pb, as.numeric(Diff))
        print(paste("Cost =",Cost))
        if(CLIP){ print("reinitialize gradients") } #
      }

    } # end of epoch

  } # end for loop
  setTxtProgressBar(pb, total)

  if(PDF){
    dev.off()
  }
  if(GIF){
    dev.off()

    # do gif
    if(iteration<30){FPS<-2}else if(iteration<100){FPS<-4}else{FPS<-10}

    animation <- magick::image_animate(img, fps = FPS)

    #### file names
    Latt<-length(list.files(paste(getwd(),"/tmp/",sep = ""),pattern =  NameProj))
    if(Latt>0){
      Latt<-max(as.numeric(gsub(".gif|.pdf","", gsub(NameProj,"", list.files(paste(getwd(),"/tmp/",sep = ""),pattern =  NameProj)))),na.rm = T)
    }
    Latt<-paste(paste(getwd(),"/tmp/",sep = ""),NameProj,Latt+1,".gif",sep = "")

    #print(animation)
    magick::image_write(animation, Latt)
  }

  names(COST)<-"Cost"
  return(list(NETall=NETall, Cost=COST, NETallList=NETallList, NETallActivity=NETallActivity) )

}
