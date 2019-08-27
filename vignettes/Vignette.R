## ------------------------------------------------------------------------
#install.packages("AMoNet", dependencies=T)
library(AMoNet)
if(TRUE){
    print("Loading packages")
    library(reshape2)
    library(parallel)
    library(stringr)
    library(splines)
    library(cgdsr)
    library(survival)
    library(magick)
    library(matrixStats)
    library(survivalROC)
    library(pec)
  }

set.seed(1234)

## ------------------------------------------------------------------------
GenesSelec<-rbind(GenesSelecImuno,GenesSelecHall)
MECA<-unique(GenesSelec$target_hgnc)  

## ------------------------------------------------------------------------

NETalldata<-NETforOpt(GENESman=c("EGFR","TP53"),treatmt="", InteractionBase = OMNI, nblayers = 2, MinConnect = 2, RestrictedBuilding = T, MeanWinit = 0.1, SdWinit = 1, Phenotypes=GenesSelec, MECA=MECA, RestrictedBase = T, addOutputs = 1, FilterCGS = F, LSTM = F, Adam = T, Activity = F,  KeepPhenotypes=T, WRITE = F, no_cores = 3)

NETall1<-NETalldata$NETall


## ------------------------------------------------------------------------
 PlotOptNet(NETall1,PDF = F,Optimized = F,PrintOptNET = F,LEGEND = F, NameProj = "MyAMoNet")

# retreive the species in your network
Species<-union(NETall1$source_hgnc,NETall1$target_hgnc)


## ------------------------------------------------------------------------
#Generate random initial states to compute future initial stable states
iStates<-matrix(0.5,nrow = 50, ncol = length(Species))
colnames(iStates)<-Species

# Simulate
TotAttractors<-BoolSimul(NETall=NETall1, Logic = "Sigmoid", Mode = "LAYER", iStates=iStates, MinSteps = 10, Discretize = F, LSTM = F, MUT=NULL, Parallel = T, no_cores = 3)

## ------------------------------------------------------------------------
Npat<-1
matplot(rbind(t(TotAttractors[Npat,"iStates",,1]),t(TotAttractors[Npat,"A",,])),type='l', ylab="Species activities",xlab = "Time steps", main="Simulation of the network")
matplot(as.matrix(c(TotAttractors[Npat,"iStates","Output",1],TotAttractors[Npat,"A","Output",])), type='l', col=1,lwd=3, add=T)


## ------------------------------------------------------------------------
# prepare mutations data
MUTa<-matrix(NA,nrow = 2,ncol = length(Species))
colnames(MUTa)<-Species
rownames(MUTa)<-c("Pat1","Pat2")#,"Pat3")

# mutate it
if(FALSE){
  MUTa[2,"TP53"]<-"mut"
  MUTa[2,"EGFR"]<-"mut"
} else{ # or
  MUTa[2,"TP53"]<-(-1)
  MUTa[2,"EGFR"]<-1
}
# transform matrix into list
MUTl<-MutMatToList(MUTa = MUTa[,Species])

# reduce to 2 initial states
iStates<-iStates[1:2,]


## ------------------------------------------------------------------------
TotAttractors<-BoolSimul(NETall=NETall1, Logic = "Sigmoid", Mode = "LAYER", iStates=iStates, MinSteps = 10, Discretize = F, LSTM = F, MUT=MUTl, ValMut = 2, Parallel = F)

Npat<-2
matplot(rbind(t(TotAttractors[Npat,"iStates",,1]),t(TotAttractors[Npat,"A",,])), type='l',main=paste(names(MUTl[[Npat]]),"mutated"), ylab="Species activities",xlab = "Time steps")

matplot(rbind(t(TotAttractors[,"iStates","Output",1]),t(TotAttractors[,"A","Output",])), type='l', ylab="Output activities",xlab = "Time steps", main="Differences in output porbabilities")

## ------------------------------------------------------------------------
Clin<-matrix(c(0.9,0.2),nrow = 1)
colnames(Clin)<-c("Pat1","Pat2")
rownames(Clin)<-"Output"


