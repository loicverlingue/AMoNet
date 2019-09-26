
### optionally install and load AMoNet
XV<-try(library(AMoNet),silent = T)

if("try-error"%in%class(XV)){
  print("Installing AMoNet from gitlab.")
  devtools::install_git("https://lverling:Lololo.1@gitlab.curie.fr/lverling/amonet.git", dependencies = TRUE,
                        destdir="~/R/x86_64-redhat-linux-gnu-library/3.4")
  library(AMoNet)
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
#args<-"--Param nblayers MinConnect --NameProj HallmarksLung --GENESman EGFR MTOR --treatmt --Interval 10 --SelectMECA HALLMARK --organ luad"

# retrieve arguments eithers form command lines or from environments
XV<-try(args)
if("try-error"%in%class(XV)){
  print("Arguments retreived from command line")
  args <- commandArgs(trailingOnly = TRUE)
} else {
  print("Arguments retrieved from R environment")
}

# retreive the hyperparameters to search from the commandline

#args <- gsub("\r","",args)

hh <- paste(unlist(args),collapse=' ')
listoptions <- unlist(strsplit(hh,'--'))[-1]
options.args <- sapply(listoptions,function(x){
  unlist(strsplit(x, ' '))[-1]
})
options.names <- sapply(listoptions,function(x){
  option <-  unlist(strsplit(x, ' '))[1]
})
names(options.args) <- unlist(options.names)
#Correct<-options.args[[1]]
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

# update default values
InterDef<-intersect(names(Default), names(options.args))
Default[InterDef]<-options.args[InterDef]

#print("Arguments with no match:")
#print(options.args[InterArg])

#
args.needed<-setdiff(names(formals(RunTCGAopt)),c("Default","Boundaries"))
InterArg<-intersect(names(options.args),args.needed)

# run workflow
net<-do.call(RunTCGAopt, c(options.args[InterArg],
                           list(Default=Default,Boundaries=Boundaries)))
#x<-load(file.path(getwd(),"model", list.files(file.path(getwd(),"model"))[1]))
#net<-get(x)

# predict
net<-PlotAndPredict(net)
