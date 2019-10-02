
### optionally install and load AMoNet
# todo: put after to put dir in arguments?
XV<-try(library(AMoNet, lib.loc = "C:/Users/L_VERLINGUE/Documents/R/win-library/3.4/"), silent = T)

if("try-error"%in%class(XV)){
  print("Failed to load from local")

#  print("Installing AMoNet from gitlab.")
#  devtools::install_git("https://lverling:Lololo.1@gitlab.curie.fr/lverling/amonet.git", dependencies = TRUE,
#                        destdir="~/R/x86_64-redhat-linux-gnu-library/3.4")
#  library(AMoNet)
}

# set Default parameters if user defined
#Default$iteration=1
#Default$nblayers=1
#Default$MinConnect=6
#Default$MiniBatch=64
#Default$lambda

########
# remove global environment variables that may interfere with cmd and generation of new arguments
rm(list = intersect(ls(), names(Default)))

set.seed(NULL)
######################
#arguments<-"--Param nblayers MinConnect --NameProj HallmarksLung --GENESman EGFR MTOR --treatmt --Interval 10 --SelectMECA HALLMARK --organ luad"

# retrieve arguments eithers form command lines or from environments
XV<-try(arguments)
if("try-error"%in%class(XV)){
  print("Arguments retreived from command line")
  arguments <- commandArgs(trailingOnly = TRUE)
} else {
  print("Arguments retrieved from R environment")
}

# retreive the hyperparameters to search from the commandline

#arguments <- gsub("\r","",arguments)

hh <- paste(unlist(arguments),collapse=' ')
listoptions <- unlist(strsplit(hh,'--'))[-1]
options.arguments <- sapply(listoptions,function(x){
  unlist(strsplit(x, ' '))[-1]
})
options.names <- sapply(listoptions,function(x){
  option <-  unlist(strsplit(x, ' '))[1]
})
names(options.arguments) <- unlist(options.names)
#Correct<-options.arguments[[1]]
options.arguments<-lapply(options.arguments,function(Correct){
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
print(options.arguments)

# update default values - used whith NewNet=T
InterDef<-intersect(names(Default), names(options.arguments))
Default[InterDef]<-options.arguments[InterDef]

#
arguments.needed<-setdiff(names(formals(RunTCGAopt)),c("Default","Boundaries"))
InterArg<-intersect(names(options.arguments),arguments.needed)

# run workflow
net<-do.call(RunTCGAopt, c(options.arguments[InterArg],
                           list(Default=Default,Boundaries=Boundaries)))

# predict
net<-PlotAndPredict(net)
