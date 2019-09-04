#' Randomise hyper parameters
#'
#' @param C characeter vector. Name(s) of the hyper-parameters to optimize
#' @param Default list. Default hyper-parameters from \code{Default}
#' @param Bondaries list. Default hyper-parameters bondaries from \code{Bondaries}
HyperP<-function(C=C, Default=Default, Bondaries=Bondaries){

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
