#' Randomise hyper parameters
#'
#' @param C characeter vector. Name(s) of the hyper-parameters to optimize
#' @param Default list. Default hyper-parameters from \code{Default}
#' @param Boundaries list. Default hyper-parameters Boundaries from \code{Boundaries}
#' @export
HyperP<-function(C=C, Default=Default, Boundaries=Boundaries){

  # define distribution you want
  LogDistrib<-function(min,max,N){
    sample(exp(seq(log(abs(min)),log(abs(max)),length.out = 100)),N)
  }

  NormDistrib<-function(min,max,N){
    sample(seq(min,max,length.out = 100),N)
  }

  for(Cp in C){
    if(is.numeric(Boundaries[[Cp]])){
      if(is.integer(Boundaries[[Cp]])){
        Default[[Cp]]<- sample(seq(Boundaries[[Cp]][1],Boundaries[[Cp]][2]),1)
      } else if(all(Boundaries[[Cp]]>=0.1)){
        Default[[Cp]]<- NormDistrib(min = Boundaries[[Cp]][1],max = Boundaries[[Cp]][2], N = 1)
      }else{
        Default[[Cp]]<- LogDistrib(min = Boundaries[[Cp]][1],max = Boundaries[[Cp]][2], N = 1)
      }
    } else if(is.character(Boundaries[[Cp]])|is.logical(Boundaries[[Cp]])){
      Default[[Cp]]<- sample(Boundaries[[Cp]],1)
    }
  }

  if("MiniBatch"%in%C){
    Default$MiniBatch<-sample(2^seq(Boundaries$MiniBatch[1],Boundaries$MiniBatch[2]),1)
  }

  list2env(Default, envir = .GlobalEnv)
  return(Default)
}
