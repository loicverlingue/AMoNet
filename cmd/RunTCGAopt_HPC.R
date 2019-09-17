
library(AMoNet)

# set Default parameters if user defined
Default$iteration=1
Default$nblayers=1
Default$MinConnect=3
Default$MiniBatch=64
net<-RunTCGAopt(Param=c("nblayers","MinConnect"), DIR=file.path(getwd(),"model"),
                NameProj="LUNG_AMoNet", GENESman=c("KRAS","MTOR"),
                treatmt=NULL, SelectMECA="HALLMARK", organ="lung",KeepData = T,
                eSS=F, NewNet=T, Default=Default, Boundaries = Boundaries)

net<-PlotAndPredict(net, DIR=NULL)

save(net, file = net$DIRSave)

###
net<-RunTCGAopt(Param=c("learningrate","lambda"), DIR=file.path(getwd(),"model"),
                NameProj="LUNG_AMoNet", GENESman=c("KRAS","MTOR"),
                treatmt=NULL, SelectMECA="HALLMARK", organ="lung", KeepData = T,
                eSS=F, NewNet=F, Default=Default, Boundaries = Boundaries)
