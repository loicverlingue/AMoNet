###
#' Initiate an AMoNet object with gene selection.
#' @aliases AMoNet
#' @usage AMoNet(GENESman=NULL,treatmt="")
#' @param GENESman character vector. A vector of gene(s) selected to build the *AMoNet* object.
#' @param treatmt character vector. A vector of gene(s) targeted by treatment(s). Used in to build the *AMoNet* object.
#' @details
#' The default \code{AMoNet()} function initiate an *AMoNet* object.
#' It can be invoqued either directly by inputing GENESman, or in an interactive mode by just runing \code{AMoNet()}.
#' Initiation consist in storing in *AMoNet* object the gene(s) and treatment(s) queries, default hyper-parameters and boundaries (and import the latter to \code{globalenv()}).
#' \code{AMoNet()} is also included in the \code{buid.AMoNet()} function: the user will be asked to enter gene(s)' names interactively, if not performed before.
#' @examples
#' \dontrun{
#' net<-AMoNet() # prints an interactive question to enter gene(s)' names.
#' }
#' @export
#' @import grDevices graphics stats utils
AMoNet <- function(GENESman, treatmt, ...) UseMethod("AMoNet")

#' @export
AMoNet.default<-function(GENESman=NULL,treatmt=NULL){
  if(is.null(GENESman)){
    GENES <- readline("Select a list of genes to build AMoNet (Hugo symbols separated with comma):\n")
    GENES <- unlist(strsplit(GENES, " |,|, | ,"))
    GENES <- setdiff(GENES, "")
  } else {
    GENES<-GENESman
  }

  # load Default and Boundaries in .GlobalEnv to be writable and accessible in functions calls
 # data("Default",package = "AMoNet", envir = .GlobalEnv)
 #  data("Boundaries",package = "AMoNet", envir = .GlobalEnv)
  data("CGS",package = "AMoNet", envir = .GlobalEnv)
  data("OMNI",package = "AMoNet", envir = .GlobalEnv)

  net<-list(GENESman=GENES, treatmt=treatmt, Parameters=list(Default=Default, Boundaries=Boundaries))
  class(net)<-"AMoNet"

  return(net)
}

#' Split samples from a *AMoNet* object
#' @aliases split
#' @usage split(x,f=0.7, drop=F, RETURN=F)
#' @param x *AMoNet* object
#' @param f numeric. Paratition split comprised in in [0,1].
#' @param RETURN boolean. Do the function returns values or just store it in the *AMoNet* object?
#' @details Can run when available \code{$Data$y} is available in the *AMoNet* object. Can be generated either by the user or by the \code{LoadCleanTCGA()} function. See vignette
#' @return if \code{RETURN=TRUE} the function returns the Train and Validation split
#' @examples
#'
#' @export
split.AMoNet<-function(x,f=0.7, drop=F){

  sdtout<-try(x$Data$y,silent = T)
  if(class(sdtout)=="try-error" ){
    print("Need $Data$y in your AMoNet objbect")
    stop()
    }

  Train<-sample(colnames(x$Data$y),ncol(x$Data$y)*f)
  Val<-colnames(x$Data$y)[!colnames(x$Data$y)%in%Train]

  x$TrainSplit$Train<-Train
  x$TrainSplit$Val<-Val

    return(x)
}

#' Generates random initial states
#'
#' @description Used prior to \code{simulate()} or \code{AMoNet()} with available \code{$Data$y} in the *AMoNet* object
#' @aliases RandomiStates
#' @usage RandomiStates(object, RETURN=F)
#'
#' @param object *AMoNet* object, S3 class.
#'
#' @return iStates, a dataframe that contains initial states
#' @examples
#' \dontrun{RandomiStates(object)}
#' @export
RandomiStates<-function(object){
  Species<-union(object$NETall$source_hgnc,object$NETall$target_hgnc)
  iStates<-matrix(runif(ncol(object$Data$y)*length(Species)),nrow = ncol(object$Data$y), ncol = length(Species))
  rownames(iStates)<-colnames(object$Data$y)
  colnames(iStates)<-Species
  object$iStates<-iStates
  return(object)
}

# build.AMoNet in NETforOptS3.R
# plot in PlotOptNet

#' Simulates *AMoNet*
#'
#' @aliases simulate
#' @usage
#' simulate(object, nsim=Default$MinStepsForward, seed=NULL, Init=NULL, iStates=NULL, MUT=NULL, treatmt=NULL, Logic=Default$Logic, Mode=Default$Mode, ValMut=Default$ValMut, LSTM=Default$LSTM, Ct=NULL, Parallel=Default$Parallelize, no_cores=Default$no_cores, Steepness=1)
#'
#' @param object *AMoNet* object, S3 class.
#' @param Logic character. The update function used for the simulations. Choices are: "Sigmoid" (default), "Boolean","tanh" or "ReLU".
#' @param nsim numeric. Number of simulation steps. Default to 5. Ideally set the number of simulation steps to achieve stable states of *AMoNet* activities.
#' @param seed generic parameter. Not used in *AMoNet* simulation as it doesn't involve stochastic behavior.
#' @param Mode character. The type of simulation to perform. Choice are "LAYER" (layer based, default), "ASYNC" (assynchronous update), "SYNC" (synchronous update)
#' @param iStates matrix. Initial states can be either a vector, data.frame, matrix or list of vectors, with name of all species as columns and example as rows
#' @param MUT list. Corresponds to mutations with names and corresponding values (list of vectors). If MUT=NULL, a wild type *AMoNet* is simulated.. Patients ordered the same as y.
#' @param treatmt list. The same for than newMUT, with treatments' targets and corresponding values (list of vectors).
#' @param Init matrix. Used to set some of the initial states (iStates). To force some initial states during learning set adaptive_iStates=F. Init=NULL otherwise.
#' @param ValMut numeric. Multicplicative factor for the effect of mutations. Default to 50.
#' @param LSTM boolean. Should a Long-Short Term Memory (LSTM) unit architecture be used?
#' @param Ct matrix. If LSTM=TRUE, Ct is used to set the initial states of the correponding values in LSTM units.
#' @param Parallel boolean. Should the simulation run in parallel (TRUE by default)? Not recommended for simulations of less than 10 examples.
#' @param no_cores numeric. If Parallel=TRUE, set the number of cores to parallelize on. Default is 4. Can detect and set to available no_cores if inferior to user defined no_cores. to
#' @param Steepness numeric. Can be used to modulate the stepness of the sigmoid curve (default to 1).
#'
#' @details
#' The S3 \code{simulate()} function calls the core \code{BoolSimul()} function. Simulations are vectorized and optinally paralellized.
#' The ASYNC mode perform a fully assynchroneous simulation that may be relatively long. The LAYER mode is an innovation in *AMoNet* to further speed up the simulations and nevertheless optain similar activities to the ASYNC mode. SYNC mode is not recommended.
#'
#' @return *AMoNet* object with \code{$TotAttractors}, a 4-dimensions array with: D1 the number of examples, D2 the results of the linear function (z=ax+b), D3 the results of the activator function (defined by the \code{Logic} argument) and D4 the number of iterations.
#' @importFrom stats simulate
#' @export
simulate.AMoNet<-function(object, nsim=Default$MinStepsForward,
                          seed=NULL,
                          Init=NULL, iStates=NULL, MUT=NULL,
                          treatmt=NULL, Logic=Default$Logic,
                          Mode=Default$Mode,
                          ValMut=Default$ValMut,
                          LSTM=Default$LSTM, Ct=NULL, Parallel=Default$Parallelize,
                          no_cores=Default$no_cores, Steepness=1 ){

  # update Default parameters
  CALL<-mget(names(formals()))
  CNames<-intersect(names(object$Parameters$Default),names(CALL))
  object$Parameters$Default[CNames] <- CALL[CNames]

  object$Parameters$Default["MinStepsForward"]<-CALL["nsim"]
  object$Parameters$Default["Parallelize"]<-CALL["Parallel"]

  if(is.matrix(MUT)){
    MUT<-MutMatToList(MUTa = MUT)
  }

 # if(!is.null(iStates)){
 #   object$iStates<-iStates
#  }

  # is it ok in .GlobalEnv?
#  list2env(object, envir = .GlobalEnv)

  Species<-union(object$NETall[,1],object$NETall[,3])
  if(is.null(Init)&is.null(MUT)&is.null(treatmt)&is.null(iStates)){
    print("As no Init, iStates, Mut or treatm defined, simulation will take 10 random initial states")
    MUT<-as.list(rep("",10))
  }

  if(is.null(iStates)){ # or not same dim
    Nb<-max(length(MUT),length(treatmt))
    iStates<-matrix(runif(Nb*length(Species)),nrow = Nb, ncol = length(Species))
    colnames(iStates)<-Species
   # object$iStates<-iStates
  }

  if(!is.null(Init)) {
    for(i in seq(nrow(Init))){
      iStates[i,colnames(Init)]<-as.numeric(Init[i,]) # match the col order
    }
   # object$iStates<-iStates
  }

  if(is.null(Ct)&object$Parameters$Default$LSTM) {
    Ct<-iStates
  }

  TotAttractors<-BoolSimul(NETall=object$NETall, Logic=Logic, Mode=Mode,
                           iStates=iStates, MUT=MUT, treatmt=treatmt,
                           ValMut=ValMut, MinSteps=nsim,
                           LSTM=LSTM, Ct=Ct, Parallel=Parallel,
                           no_cores=no_cores, Steepness=Steepness)

  object[["TotAttractors"]]<-TotAttractors
  object$call$simul_call<-match.call()

  # update Default in .GlobalEnv
#  Default<<-object$Parameters$Default
  assign("Default",object$Parameters$Default,envir = .GlobalEnv)

 return(object)
}

#plot simulate in plotAMoNet function


# to check
if(FALSE){
summary.AMoNet<-function(net,Species=NULL){
  #list2env(net,envir = globalenv())
  if(is.null(Species)){
    Species<-union(net$NETall$source_hgnc,net$NETall$target_hgnc)
  }

  if(exists("TotAttractors")){

      if(length(Species)==1){
        Res<-mean(net$TotAttractors[,"A",Species,dim(net$TotAttractors)[4]])
      } else{
        Res<-apply(net$TotAttractors[,"A",Species,dim(net$TotAttractors)[4]],2,mean)
        names(Res)<-Species
      }
      class(Res) <- "summary.AMoNet"
      return(Res)

  } else {
    print("simulate or predict first")
  }
}

print.summary.AMoNet<-function(x,...){
  print(x)
}
#summary(net, Species = "Output")
#Res<-summary(net, Species = "Output")
#class(Res)
#View(Res)
#as.numeric(Res)
}
########


#' AMoNet training
#'
#' Runs the \code{RunBackSimul} function (see \code{?RunBackSimul} for detailed description)
#'
#' @aliases train
#' @usage
#' train(object, y, MUT=NULL, treatmt=NULL, Init=NULL, iStates=NULL, KeepTraining=T, KeepData=T, Ct=NULL, MiniBatch=Default$MiniBatch, Optimizer=Default$Optimizer, beta1=Default$beta1, beta2=Default$beta2, iteration=Default$iteration, learning_rate=Default$learningrate, adaptive_iStates=Default$adaptive_iStates, FixNodes=NULL, Parallelize=Default$Parallelize, no_cores=Default$no_cores, Logic = Default$Logic, Mode = Default$Mode, ModeBack=Default$Mode, MinStepsForward = Default$MinStepsForward,MinStepsBackward=Default$MinStepsBackward, LSTM=Default$LSTM, gradClipping=Default$gradClipping, LearningRateDecay=Default$LearningRateDecay, ValMut=Default$ValMut, PDF=F, GIF=F, NameProj=net$call$build_call$NameProj, Visualize=Default$Visualize, alpha=Default$alpha, lambda=Default$lambda)
#'
#' @param object *AMoNet* object, S3 class.
#' @param y numeric or matrix. labels
#' @param MUT list. Corresponds to mutations with names and corresponding values (list of vectors). If MUT=NULL, a wild type *AMoNet* is simulated.. Patients ordered the same as y.
#' @param treatmt list. The same for than newMUT, with treatments' targets and corresponding values (list of vectors).
#' @param Init matrix. Used to set some of the initial states (iStates). To force some initial states during learning set adaptive_iStates=F. Init=NULL otherwise.
#' @param iStates matrix. Initial states of the simulations.
#' @param KeepTraining boolean. If \code{TRUE} (default) stores the training in the *AMoNet* object, which is usefull to plot the training step.
#' @param KeepData boolean. If \code{TRUE} (default) stores the data in the *AMoNet* object. Usefull to run prediction on same data as training.
#' @param Optimizer character. Options are "Adam","Momentum","RMSprop" or NULL. NULL corresponds to gradient descent and can be used with batch learning. Use of optimizers are recommended to mini batch learning.
#' @param beta1 numeric. Momentum term. Adam is defined by Momentum / RMSprop.
#' @param beta2 numeric. RMSprop term. Adam is defined by Momentum / RMSprop.
#' @param MiniBatch integer. nb of sample per MiniBatch. For batch learning, use MiniBatch = total nb of examples
#' @param iteration integer. Number of iterations, corresponding to epochs, i.e. one path through the full data.
#' @param FixNodes character. Freeze parameters of selected nodes.
#' @param Visualize character or integer. Options : \code{c(1,"Gates","Output")} Real time visualization of learning, requires selecting the visualization target: "Gates","Output", or all nodes activities of a patients by selecting an integer (a patient's number) or an ID. Default is no visualization: set to \code{NULL}, .
#' @param alpha numeric. Regularization term. Set to 1 for L1N (Lasso) regularization, set to 2 for L2N (Ridge) regularization, between 1 and 2 with lambda 1 for L1N + L2N (Elastic Net) regularization (penalization parameter is alpha).
#' @param lambda numeric. the lambda value for Lasso, Ridge or Elastic Net regularization. Set to 0 corresponds to no regularization. When alpha between 1 and 2 (Elastic Net regularization), lambda should be set to 1 and penalization parameter is alpha.
#' @param Parallelize boolean. Should the simulation run in parallel (TRUE by default)? Not recommended for simulations of less than 10 examples.
#' @param no_cores numeric. If Parallel=TRUE, set the number of cores to parallelize on. Default is 4. Can detect and set to available no_cores if inferior to user defined no_cores.
#' @param gradClipping nunmeric. Maximum values allowed for gradients otherwise clipped. Gradient scaling in dev.
#' @param LSTM boolean. Should a Long-Short Term Memory (LSTM) unit architecture be used?
#' @param MinStepsForward numeric. Number of simulation steps. Default to 5. Ideally set the number of simulation steps to achieve stable states of *AMoNet* activities.
#' @param MinStepsBackward numeric. Number of back-simulation steps. Default to 1, as in a standard back-propagation of gradient in neural nets.
#' @param PDF boolean. Whether to save a pdf from visualizations of the learning phase.
#' @param GIF boolean. Whether to save a gif from visualizations during the learning phase.
#' @param NameProj character. Name of the project.
#' @param Ct matrix. If LSTM=TRUE, Ct is used to set the initial states of the correponding values in LSTM units.
#' @param learning_rate numeric. Learning rate value.
#' @param adaptive_iStates boolean. Do the initial states of the simulations be updated at each epoch? Deafult is TRUE. Set to FALSE, the initial states of simulations will always be the same.
#' @param Logic character. The update function used for the simulations. Choices are: "Sigmoid" (default), "Boolean","tanh" or "ReLU".
#' @param Mode character. The type of simulation to perform. Choice are "LAYER" (layer based, default), "ASYNC" (assynchronous update), "SYNC" (synchronous update)
#' @param ModeBack character. The type of back-simulation to perform. Choice are "LAYER" (layer based, default), "ASYNC" (assynchronous update), "SYNC" (synchronous update)
#' @param LearningRateDecay character. Choice are "linear", "exponential" or NULL.
#' @param ValMut numeric. Multicplicative factor for the effect of mutations. Default to 50.
#'
#' @details
#' The \code{train.AMoNet()} function mainly calls \code{RunBackSimul()} function
#'
#' @return stores in *AMoNet* object a list with the network dataframe NETall, and for each epoch, if \code{KeepTraining=TRUE}:
#' the training Cost, the network weights (NETallList) and the simulation activities (NETallActivity)
#'
#' @export
train.AMoNet <- function(object, y, MUT=NULL, treatmt=NULL,
                           Init=NULL, iStates=NULL,
                           KeepTraining=T, KeepData=T,
                         Ct=NULL, MiniBatch=Default$MiniBatch,
                         Optimizer=Default$Optimizer, beta1=Default$beta1, beta2=Default$beta2,
                         alpha=Default$alpha, lambda=Default$lambda,
                         iteration=Default$iteration, learning_rate=Default$learningrate,
                         adaptive_iStates=Default$adaptive_iStates, FixNodes=NULL,
                         Parallelize=Default$Parallelize, no_cores=Default$no_cores,
                         Logic = Default$Logic, Mode = Default$Mode, ModeBack=Default$Mode,
                         MinStepsForward = Default$MinStepsForward,MinStepsBackward=Default$MinStepsBackward,
                         LSTM=Default$LSTM, gradClipping=Default$gradClipping, LearningRateDecay=Default$LearningRateDecay,
                         ValMut=Default$ValMut, PDF=F, GIF=F, NameProj="My_AMoNet",
                         Visualize=Default$Visualize){

  # update Default parameters
  CALL<-mget(names(formals()))
  CNames<-intersect(names(object$Parameters$Default),names(CALL))
  object$Parameters$Default[CNames] <- CALL[CNames]

  # retrieve the net
  NETall <- object$NETall
  NameProj <- object$call$build_call$NameProj

  # set data
  if(is.matrix(MUT)){
    MUT<-MutMatToList(MUTa = MUT)
  }
  #y <- as.matrix(y)

  if("TotAttractors"%in%names(object)&is.null(iStates)){
    iStates<-object$TotAttractors[,"A",,dim(object$TotAttractors)[4]]
  } else if(is.null(iStates)&"iStates"%in%names(object)){
    iStates<-object$iStates
  } else if(is.null(iStates)){
    print("Either provide iStates or run a simulation of *AMoNet* first")
    stop()
  }

  # train
  NETallProp <- RunBackSimul(NETall = NETall,
                             y = y,
                             MUT = MUT,
                             iStates = iStates,
                             NameProj = NameProj,
                             treatmt=treatmt,
                             Ct=Ct, MiniBatch = MiniBatch,
                             Init=Init, Optimizer=Optimizer, beta1=beta1, beta2=beta2,
                             iteration=iteration, learning_rate=learning_rate,
                             adaptive_iStates=adaptive_iStates, FixNodes=FixNodes,
                             Parallelize=Parallelize, no_cores=no_cores,
                             Logic = Logic, Mode = Mode, ModeBack=ModeBack,
                             MinStepsForward = MinStepsForward,MinStepsBackward=MinStepsBackward,
                             LSTM=LSTM, gradClipping=gradClipping, LearningRateDecay=LearningRateDecay,
                             ValMut=ValMut, PDF=F, GIF=F,
                             Visualize=Visualize, alpha=alpha, lambda=lambda)

  # store
  if(KeepTraining){
    object$history<-NETallProp
    class(object$history)<-"AMoNet.history"
  }
  if(KeepData){
    DATA<-list(y=y,MUT=MUT,treatmt=treatmt,Init=Init)
    NewData<-setdiff(names(DATA),names(object$Data))
    #print(paste("saving new data:",NewData,collapse = ", "))
    object$Data[NewData]<-DATA[NewData]
  }

  object$NETall<-NETallProp$NETallList[[length(NETallProp$NETallList)]]
  object$Cost<-NETallProp$Cost
  object$call$train_call<-match.call()

  # update Default in globalenv()
  #  Default<<-object$Parameters$Default
  assign("Default",object$Parameters$Default,envir = .GlobalEnv)

  iStates<-NETallProp$NETallActivity[[length(NETallProp$NETallActivity)]]

  if(all(dim(iStates)==dim(object$iStates))){
    object$iStates<-iStates
  }
  class(object)<-"AMoNet"
  return(object)
}

# to improve
#' Print training cost
#' @param x of class "AMoNet.history" created in the \code{$history} of the *AMoNet* object.
#' @export
print.AMoNet.history<-function(x=net$history, ...){
  tail(x$Cost,1)
}

#' Plot history of training performed with the \code{AMoNet()} function
#' @aliases plot
#' @usage plot(x, ...)
#' @param x of class *AMoNet.history*. Can be foudn in *AMoNet* object \code{$history}
#' @param ... additional optional arguments.
#' @export
plot.AMoNet.history<-function(x, Interval=Default$Interval, ...){
  if(class(x)!="AMoNet.history"){
    print("Input the $history of trained AMoNet object")
  }
  # cost
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)+0.1)
  matplot(x$Cost,type = "l", main='Cost', xlab='Epochs', ylab = "Cost")

  # outputs
  if(Interval>1){
    OutputsEvol<-lapply(x$NETallActivity,function(x){
      mean(apply(x[,grep("Output",colnames(x))],1,sd))
      #sd(x[,grep("Output",colnames(x))])
    })
  } else{
    OutputsEvol<-lapply(x$NETallActivity,function(x){
      sd(x[,"Output"])
    })
  }

  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)+0.1)
  matplot(unlist(OutputsEvol),type='l',add=F, xlab = "Epochs", ylab = "SD",
          main="Evolution of std deviations \n in outputs through learning")

  # weigths
  LWeights<-lapply(x$NETallList,function(x){
    x$Weights
  })
  LWeights<-(do.call('rbind',LWeights))
  matplot(LWeights,type='l',add=F, main="Evolution of weights through learning")

  # per edges
  NETall<-x$NETallList[[length(x$NETallList)]]

  VAR<-as.matrix(apply(LWeights,2,sd))
  par(mar=c(3,8,4,2))
  barplot(VAR[VAR>=quantile(VAR,0.89)],names.arg = apply(NETall[VAR>=quantile(VAR,0.89),c(1,3)],1,function(z){ paste(z,collapse = "_")}), las=2, cex.names = 0.5, horiz = T,space = 0, cex.axis = 0.5, main="Top weights' correction")
}



#' Perform predictions from a trained *AMoNet* object
#'
#' @aliases predict
#' @usage
#' predict(object, newy=NULL, newInit=NULL, newiStates=NULL, newMUT=NULL, newtreatmt=NULL, Logic = Default$Logic, Mode = Default$Mode, MinSteps = Default$MinStepsForward, LSTM = Default$LSTM, Parallelize = Default$Parallelize, no_cores = Default$no_cores, ValMut = Default$ValMut)
#'
#' @param object *AMoNet* object
#' @param newy numeric or matrix. new labels used to compute metrics
#' @param newInit matrix. Used to set some of the initial states (iStates). Init=NULL otherwise.
#' @param newiStates matrix. Initial states of the simulations.
#' @param newMUT list. Corresponds to mutations with names and corresponding values (list of vectors). If MUT=NULL, wild type *AMoNet* is simulated. Patients ordered the same as newy.
#' @param newtreatmt list. In the same for than newMUT, with treatments' targets and corresponding values (list of vectors).
#' @param ValMut numeric. Multicplicative factor for the effect of mutations. Default to 50.
#' @param SplitType character. With split set to use between: "Train", "Val" or NULL
#'
#' @details
#' The \code{predict()} runs new simulations of the *AMoNet* object with new data.
#' Can be used to evaluate the performance of the model on a validation set.
#'
#' @return *AMoNet* object with \code{$Predict} of class *predict.AMoNet* and results of simulations in \code{$TotAttractors}
#'
#' @export
predict.AMoNet<-function(object, newy=NULL, newInit=NULL, newiStates=NULL, newMUT=NULL,
                         newtreatmt=NULL, SplitType=NULL){

  # todo/tocheck #  not feasible to retrieve ... arguments I think
  # update Default parameters
 # CALL<-mget(names(formals()))
 # CNames<-intersect(names(object$Parameters$Default),names(CALL))
#  object$Parameters$Default[CNames] <- CALL[CNames]

  if(is.null(newInit)&is.null(newiStates)&is.null(newMUT)&is.null(newtreatmt)){
    #print("Use load data to perform predictions")
    list2env(object$Data,envir = .GlobalEnv)

    newtreatmt=NULL
    #list2env(Default,envir = globalenv())
    #print(dim(iStates))
    iStates<-object$iStates
  } else {
    Init=newInit; iStates=newiStates; MUT=newMUT;treatmt=newtreatmt
  }

  if(!is.null(newtreatmt)){
    treatmt=list(rep(newtreatmt,length(MUT)))
  }


  if(is.null(iStates)){
    print("Genrate random initial states as newiStates is NULL")
    object<-RandomiStates(object)
    iStates<-object$iStates
  }

  if(is.null(newy)){
    y<-object$Data$y
  } else{
    y<-newy
  }

  if(!is.null(SplitType)){
    Init = Init[net$TrainSplit[[SplitType]]]
    iStates = iStates[net$TrainSplit[[SplitType]],]
    MUT=MUT[net$TrainSplit[[SplitType]]]
    treatmt = treatmt[net$TrainSplit[[SplitType]]]
    y=y[,net$TrainSplit[[SplitType]]]

  }

  object<-simulate(object, nsim = object$Parameters$Default$MinStepsForward,
                   Init = Init, iStates = iStates, MUT=MUT, treatmt = treatmt,
                   Logic = object$Parameters$Default$Logic, Mode = object$Parameters$Default$Mode,
                   Parallel = object$Parameters$Default$Parallelize, no_cores = object$Parameters$Default$no_cores,
                   LSTM = object$Parameters$Default$LSTM, ValMut = object$Parameters$Default$ValMut)

  ###
  # Predictions
  if(object$Parameters$Default$Interval>1){

    Pred_mat<-object$TotAttractors[,"A",paste("Output",seq(object$Parameters$Default$Interval),sep = ""),dim(object$TotAttractors)[4]]
    rownames(Pred_mat)<-rownames(iStates)

  } else {

    Pred_mat<-as.matrix(object$TotAttractors[,"A","Output",dim(object$TotAttractors)[4]])
    rownames(Pred_mat)<-rownames(iStates)
    colnames(Pred_mat)<-"Output"

  }

  Cost <- mean( (t(y)-Pred_mat)^2 )

  # store predictions
  object$Predict<-list(Pred_mat=Pred_mat)

  # store metrics
  object$Predict$metrics<-list(Cost=Cost)


  if(object$Parameters$Default$Interval>1){

    # go back to censored data

    newy<-apply(y,2,function(Y){
      Time<-max(which(Y==1))
      Status<-ifelse(Y[Time+1]==0,1,0)
      Status<-ifelse(is.na(Status),0,Status)
      return(c(Time,Status))
    })
    newy<-t(newy)

    # train metrics
    Pred<-rowMeans(Pred_mat)

    Cindex<-survival::survConcordance(survival::Surv(newy[,1],newy[,2])~as.numeric(Pred))

    # store
    object$Predict$metrics$Cindex<-c(Cindex$concordance,Cindex$std.err)
  } else {
    object$Predict$metrics$Cindex<-NULL
  }
  object$Predict$newy<-newy

  class(object$Predict)<-"predict.AMoNet"

  #
  names(object)[names(object)=="Predict"]<-paste("Predict",SplitType,sep="_")

  # update Default in globalenv()
  #  Default<<-object$Parameters$Default
  assign("Default",object$Parameters$Default, envir = .GlobalEnv)

  return(object)

}

#' Plot predictions generated from the \code{predicte()} function
#' @aliases plot
#' @usage plot(x=net$Predict, ...)
#' @param x object of class *predict.AMoNet* generated with the \code{predict()} function and accessible with *AMoNet* object subobject \code{object$Predict}.
#' @param ... additional arguments for the \code{plot()} function.
#' @return The relation between true and predicted values.
#' @export
plot.predict.AMoNet<-function(x=net$Predict,
                              main= "predictions",
                              ylab="True survival time intervals",
                              xlab="Predicted mean survival probability",
                              pch=1,
                              ...){

  Pred<-rowMeans(x$Pred_mat)
  newy<-x$newy

  par(mfrow=c(1,1))
  par(mar=c(4,4,2,2)+0.1)

  plot(Pred,newy[,1],main=main,ylab=ylab,xlab=xlab,pch=ifelse(newy[,2]==0,2,1), ...)
  #CindexTrain<-survConcordance(Surv(Ytrain[Train,1],Ytrain[Train,2])~Pred_train)
  #}
  legend("topleft",legend = c(paste("C-index =", 1-round(x$metrics$Cindex["concordant"],3),
                                    "; std =", round(x$metrics$Cindex["std(c-d)"],3)),
                              paste("MSE =",round(x$metrics$Cost,3))#,
                              #paste("IPCW =",round(as.numeric(IPCWtrain$AppCindex),3))
  ),cex=0.5)

}


#' Write the network of the *AMoNet* object
#' @description Used to export a .csv file for visualization and other applications
#' @param object *AMoNet* object, S3 class.
#' @param file a (writable binary-mode) connection or the name of the file where the data will be saved.
#' @export
write.AMoNet.network<-function(object, file=stop("'file' must be specified")){
  write.csv(x = object$NETall, file = file)
}
