#' Create and Initialize the weights of the network
#' \code{addWeights}
#' @param NETall df, the AMoNet network architecture
#' @param SdWinit numeric, standard deviation of the random initialization, using the \code{rnorm} base function
#' @param MeanWinit numeric, mean of the random initialization, using the \code{rnorm} base function
#' @param Scaling bool, scales weights, Xavier initialization
#' @param Adam bool, create weigths for advanced optimizers such as Adam, including and used also for Momentum or RMSprop
#' @param LSTM bool, is the unit LSTM
#' @return an updated AMoNet with weights

addWeights<-function(NETall=NETall, SdWinit=2, MeanWinit=1, Scaling=T, Adam=T, LSTM=T){
  # add Weights for each interaction
   # NETall[,"Weights"]<-runif(nrow(NETall),0,MAX)
    NETall[,"Weights"]<-rnorm(nrow(NETall),MeanWinit,SdWinit)

    # maybe better to do a new column with positive or negative value, to multiply with in the simulations/optimisation
    #  NETall[,"Value_interaction"]<-ifelse(NETall$interaction_directed_signed%in%"ACTIVATE",  1,-1)
    #plot(density(NETall$Weights))
    # or
    # Xavier initialisation
    if(Scaling){
      ORD<-order(NETall$Layer)
      layers_dims<-table(NETall$Layer[!duplicated(NETall$source_hgnc)][ORD])
      layers_dims[paste(length(layers_dims)+1)]<-0

      for(l in sort(unique((NETall$Layer)))){
        NETall[NETall$Layer==l-1,"Weights"]<-NETall[NETall$Layer==l-1,"Weights"] * sqrt(1/layers_dims[l])
      }
  }
  NETall[,"Weights"]<-ifelse(NETall$interaction_directed_signed%in%"ACTIVATE",  abs(NETall[,"Weights"]+10^-8),- abs(NETall[,"Weights"]+10^-8))

  NETall[,"bterm"]<-0 # no symetry problem

  if(Adam){
    NETall[,"dWeights"]<-0
    NETall[,"dbterm"]<-0
    NETall[,"VdWeights"]<-0
    NETall[,"SdWeights"]<-0
    NETall[,"Vdbterm"]<-0
    NETall[,"Sdbterm"]<-0
  } else {
    NETall[,"dWeights"]<-0
    NETall[,"dbterm"]<-0
  }
  if(LSTM){
    #gates: i=input, f=forget, o=output
    NETall[,"Whi"]<-NETall$Weights+rnorm(nrow(NETall),mean = 0,sd = 0.01)
    NETall[,"bi"]<-0.1 # input gate should be should be open (>1) at the begining
    NETall[,"Whf"]<-NETall$Weights+rnorm(nrow(NETall),mean = 0,sd = 0.01)
    NETall[,"bf"]<-0.1  # forget gate should be open (>1)/shut off at the begining
    NETall[,"Who"]<-NETall$Weights +rnorm(nrow(NETall),mean = 0,sd = 0.01)
    NETall[,"bo"]<-0.1 # idem
#    NETall[,"Ct"]<-1 # ?
    if(Adam){
      NETall[,"VdWhi"]<-0;NETall[,"SdWhi"]<-0;NETall[,"Vdbi"]<-0;NETall[,"Sdbi"]<-0
      NETall[,"VdWhf"]<-0;NETall[,"SdWhf"]<-0;NETall[,"Vdbf"]<-0;NETall[,"Sdbf"]<-0
      NETall[,"VdWho"]<-0;NETall[,"SdWho"]<-0;NETall[,"Vdbo"]<-0;NETall[,"Sdbo"]<-0
    } else {
      NETall[,"dWhi"]<-0;NETall[,"dbi"]<-0;NETall[,"dWhf"]<-0;NETall[,"dbf"]<-0
      NETall[,"dWho"]<-0;NETall[,"dbo"]<-0
    }
  }

  return(NETall)
}
