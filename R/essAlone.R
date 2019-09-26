#' ess: multi start simulations
#'
#' @param net *AMoNet* object
#' @param y numeric or matrix. labels
#' @param MUT list. Corresponds to mutations with names and corresponding values (list of vectors). If MUT=NULL, a wild type *AMoNet* is simulated.. Patients ordered the same as y.
#' @param treatmt list. The same for than newMUT, with treatments' targets and corresponding values (list of vectors).
#' @param Init matrix. Used to set some of the initial states (iStates). To force some initial states during learning set adaptive_iStates=F. Init=NULL otherwise.
#' @param iStates matrix. Initial states of the simulations.
#' @param treatmt list. The same for than newMUT, with treatments' targets and corresponding values (list of vectors).
#' @export
ess<-function(net, Init, MUT, iStates, y, treatmt){ # Quant=30, no_cores=3, Default,  CGS, NETall1, ValMut,

  if(FALSE){
  #  print(environment())
  cl <- parallel::makeCluster(no_cores)

  parallel::clusterExport(cl,c("net","Init", "MUT", "iStates", "y","treatmt", "predict.AMoNet","simulate.AMoNet"),
                          envir = environment())

  ScreenEss<-parallel::parLapply(cl,seq(Quant),function(Q){

    Species<-union(net$NETall$source_hgnc,net$NETall$target_hgnc)

    net2<-predict(net, newy = y, newInit = Init, newiStates = iStates,
                  newMUT = MUT, newtreatmt = treatmt)
    return(net2)
  })

  parallel::stopCluster(cl)

  Res<-t(unlist(sapply(ScreenEss,function(net2){
    net2$Predict_$metrics$Cost
  })))

  net<-ScreenEss[[which(Res==min(Res))[1] ]]
  }

  Species<-union(net$NETall$source_hgnc,net$NETall$target_hgnc)

  net<-predict(net, newy = y, newInit = Init, newiStates = iStates,
                newMUT = MUT, newtreatmt = treatmt)

  Res<-net$Predict_$metrics$Cost

  net$history[["NETallList"]]<-list(net$NETall)
  net$history[["NETallActivity"]]<-NULL # to check
  net$history[["Cost"]]<-Res
  net$call$train_call<-list(Optimizer="eSS", iteration=1)

  return(net)
}

