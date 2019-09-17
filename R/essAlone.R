#' ess: multi start simulations
#'
#' @param net *AMoNet* object
#' @param Quant numeric. Number of evaluations
#' @param y numeric or matrix. labels
#' @param MUT list. Corresponds to mutations with names and corresponding values (list of vectors). If MUT=NULL, a wild type *AMoNet* is simulated.. Patients ordered the same as y.
#' @param treatmt list. The same for than newMUT, with treatments' targets and corresponding values (list of vectors).
#' @param Init matrix. Used to set some of the initial states (iStates). To force some initial states during learning set adaptive_iStates=F. Init=NULL otherwise.
#' @param iStates matrix. Initial states of the simulations.
#' @param treatmt list. The same for than newMUT, with treatments' targets and corresponding values (list of vectors).
#' @param no_cores numeric. If Parallel=TRUE, set the number of cores to parallelize on. Default is 4. Can detect and set to available no_cores if inferior to user defined no_cores.
#' @export
ess<-function(net, Quant=30, Init, MUT, iStates, y,
              treatmt, no_cores=3){ #Default,  CGS, NETall1, ValMut,

  cl <- makeCluster(no_cores)

  clusterExport(cl,c("net","Init", "MUT", "iStates", "y","treatmt"))
#  clusterEvalQ(cl, library(AMoNet))
  # clusterExport(cl,c("net","Init", "MUT", "iStates", "y","treatmt"), envir = .GlobalEnv)
#  clusterExport(cl,c("predict","simulate"), envir = as.environment("package:AMoNet"))


  ScreenEss<-parLapply(cl,seq(Quant),function(Q){

    library(AMoNet) # maybe not the best way but... cant find how to use & export AMoNet::predict
    Species<-union(net$NETall$source_hgnc,net$NETall$target_hgnc)

    net2<-predict(net, newy = y, newInit = Init, newiStates = iStates,
                  newMUT = MUT, newtreatmt = treatmt)
    return(net2)
  })

  stopCluster(cl)

  Res<-t(unlist(sapply(ScreenEss,function(net2){
    net2$Predict_$metrics$Cost
  })))

  net<-ScreenEss[[which(Res==min(Res))[1] ]]

  net$history[["NETallList"]]<-net
  net$history[["NETallActivity"]]<-NULL # to check
  net$history[["Cost"]]<-Res

  return(net)
}

