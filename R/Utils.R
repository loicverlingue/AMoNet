
#' normalisation
#' @param x vector or matrix to normalise between [0,1]
normalized<-function(x) {(x-min(x))/(max(x)-min(x))}

#' retreive the hyperparameters to search from the commandline and write it on globalenv()
CommandArgs<-function(){

    args <- commandArgs(trailingOnly = TRUE)
    hh <- paste(unlist(args),collapse=' ')
    listoptions <- unlist(strsplit(hh,'--'))[-1]
    options.args <- sapply(listoptions,function(x){
      unlist(strsplit(x, ' '))[-1]
    })
    options.names <- sapply(listoptions,function(x){
      option <-  unlist(strsplit(x, ' '))[1]
    })
    names(options.args) <- unlist(options.names)


  #args <- gsub("\r","",args)


  options.args<-lapply(options.args,function(Correct){
    if(length(Correct)==0){
      return("")
    } else if(all(!is.na(as.numeric(Correct)))){
      return(as.numeric(Correct))
    } else if(exists(Correct)){
      return(eval(parse(text = Correct)))
    } else {
      return(Correct)
    }
  })

  print("options are:")
  print(options.args)

  # if from command line:
  if(exists("options.args")){
    list2env(options.args,globalenv())
  }

  return(options.args)
}

