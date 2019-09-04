#' Use already built network in AMoNet
#'
#' @param object *AMoNet* object.
#' @param NETall data frame. The newtork you provide to *AMoNet*. See the "details" section for format requirements.
#' @param MeanWinit numeric. The mean value of the normal distribution used to initiate the weights.
#' @param SdWinit numeric. The standard deviation of the normal distribution used to initiate the weights.
#' @param LSTM boolean. Will you use LSTM unit update in the learning phase?
#' @param Optimizer boolean. Will you use another optimizer than GD (ie Adam, Momentum, RMSprop), to initiate related weights.
#' @param Outputs character vector. The name(s) of the outputs of the network.
#' @param NameProj character. The name of your project.
#'
#' @details
#' The network files provided in the NETall argument should have the following format:
#' * It should be a R data.frame object
#' * First column name is "source_hgnc" and stores the names of source species in character format.
#' * Second column name is "interaction_directed_signed" and encode the type of interaction for each row with the following characters: "ACTIVATE",  "INHIBIT".
#' * Third column name is "target_hgnc" and stores the names of target species in character format.
#' * Rownames are not used.
#'
#' Example of the first 3 rows of a data frame supported:
#'
#' source_hgnc interaction_directed_signed target_hgnc
#' SRC                     INHIBIT          PTEN
#' PTEN                    ACTIVATE         E2F1
#' CSNK2A1                 INHIBIT          PTEN
#'
#' Names of genes are in Hugo nomenclature. It is necessary to use data from TCGA with the \code{LoadCleanTCGA()} function, for exeample.
#' For any other purpose, the nomenclaure is at the discretion of the user.
#'
#' Run \code{print(Default)} to check and optionnaly set default values.
#'
#' @return An object of class *AMoNet*. Stores the name of the project, the network, the initial states \code{iStates=NULL} and \code{Ct=NULL} if LSTM.
#' @export
prebuilt.AMoNet<-function(object, NETall=NULL,
                       MeanWinit=Default$MeanWinit, SdWinit=Default$SdWinit,
                       LSTM=Default$LSTM, Optimizer = Default$Optimizer,
                       Outputs="Output", NameProj="My AMoNet"){

  # update Default parameters
  CALL<-mget(names(formals()))
  CNames<-intersect(names(object$Parameters$Default),names(CALL))
  object$Parameters$Default[CNames] <- CALL[CNames]

 if(is.null(NETall)){
    print("Please provide a network file in NETall argument, or use the build.AMoNet() function instead")
    stop()
  }

  # add output column
  NETall[NETall[,3]%in%Outputs,"Output"]<-TRUE

  #Analyse net
  Analysis<-NetAnalyser(NETall, Propagation = "back",PDF=F, PLOT=F) # warnings due to autoregulated nodes (FOXP3): no danger
  NETall$Layer<-Analysis$Layers

  # add weigths to new nodes or remove new nodes
  NETall<-addWeights(NETall = NETall, SdWinit = object$Parameters$Default$SdWinit,
                     MeanWinit = object$Parameters$Default$MeanWinit , Scaling = T,
                     Adam = !is.null(object$Parameters$Default$Optimizer), LSTM = object$Parameters$Default$LSTM)


  # output
  netLIST<-list(call=list(build_call=CALL), NETall=NETall,iStates=NULL,Ct=NULL) # Parameters=list(object$Parameters$Default=object$Parameters$Default, Boundaries=Boundaries)
  DIFF<-setdiff(names(object),names(netLIST))

  object<-c(object[DIFF],netLIST)
  class(object)<-"AMoNet"

  # update Default in .GlobalEnv
  assign("Default",object$Parameters$Default,envir = .GlobalEnv)

  return(object)
}
