##########################
# select best optimised nets

# parameters:
# todo

#DIR<-"C:/Users/L_VERLINGUE/Desktop/ModelK/Rpack/ArtMolNet/results/Hallmarks/"
#NameProjbase="Hallmarks"

#' Title
#'
#' @param NameProjbase
#' @param DIR
#' @param NETall
#' @param Default
#' @param Adam
#' @param LSTM
#' @param Bondaries
#' @param addOutputs
#' @param MiniBatch
#' @param ValSelect boolean. Whether the metrics for validation should be used for model selection. Default is \code{TRUE}
#'
#' @return
#' @export
#'
#' @examples
SelectGTExOptNet<-function(NameProjbase="LungRegular", DIR=file.path(getwd(),"/model"), NETall=NETall,
                           Default=Default,Adam=!is.null(Default$Optimizer),LSTM=Default$LSTM,
                           Bondaries=Bondaries,
#                           addOutputs=Default$Interval,MiniBatch=Default$MiniBatch,
                           ValSelect=T){
  if(is.null(NETall)){
    COSTS<-list()
    FILES<-list.files(DIR, pattern = NameProjbase)
    print(FILES)

    if(length(FILES)>1){
      #f<-FILES[1]

      ### learning curves, cost first, or cindex if available
      # for validation+++ or train
      # todo: add other metrics: accuracy, F1, etc...

      COST<-lapply(FILES,function(f){
        load(file.path(DIR,f))

        if(ValSelect){
          if("Cindex"%in%names(NETallProp)){
            print("Select on validation C-index")
            COSTsim<-as.numeric(1-NETallProp$Cindex["Val"])
          } else {
            print("Select on validation Cost")
            COSTsim<-NETallProp$CostVal
          }
        } else {
          if("Cindex"%in%names(NETallProp)){
            print("Select on training C-index")
            COSTsim<-as.numeric(1-NETallProp$Cindex["Train"])
            # } else if(!is.null(MiniBatch)){
            #   MNB<-round(length(unlist(NETallProp$TrainSplit$Train))/NETallProp$Parameters$Default$MiniBatch)
            #    if(MNB==0){MNB=10}
            #    COSTsim<-median(as.numeric(tail(unlist(NETallProp$Cost),MNB)))
          } else {
            print("Select on training Cost")
          #  if(!is.null(Default$MiniBatch)){
                 MNB<-round(length(unlist(NETallProp$TrainSplit$Train))/NETallProp$Parameters$Default$MiniBatch)
                  if(MNB==0){MNB=1}
                  COSTsim<-median(as.numeric(tail(unlist(NETallProp$Cost),MNB)))
          #  } else {
              # last train cost
          #    COSTsim<-tail(as.numeric(unlist(NETallProp$Cost)),1)
          #  }
          }
        }
        return(COSTsim)
      })

      names(COST)<-FILES
      COST<-unlist(COST)

      BESTone<-names(which(COST==min(COST)))[1]

      print(paste("Best AMoNet is :", gsub(".Rdata","",BESTone) ))

      load(file.path(DIR,BESTone))
#  net<-prebuilt.AMoNet(net, NETall=NULL, MeanWinit = 0.1, SdWinit = 0.1, LSTM=F, Optimizer = Default$Optimizer, Outputs="Output", NameProj="AMoNet_TCGA")

      # plot it
      NETallpre<-NETallProp$NETallList[[length(NETallProp$NETallList)]]
      NETall<-PlotOptNet(NETall = NETallpre,PDF = F,Optimized = T,
                         NameProj = paste(NameProjbase,"\nmin cost=", round(min(COST),4) ) ,
                         PrintOptNET = T)

      ###################
      # retrieve hyperparamters
      if(!is.null(Default)){

        VARname<-unique(gsub(".Rdata","", gsub(NameProjbase,"", unlist(strsplit(FILES,"_")))))
        VARname<-VARname[is.na(as.numeric(VARname))]
        VARname<-VARname[!VARname%in%""]
        VARname<-VARname[VARname%in%names(Default)]
        #VARname

        RES<-list()
        #i=VARname[1]
        for(i in VARname){
          # value of the variable from the file names
          RES[[i]]<-abs(as.numeric(unlist(lapply(strsplit(FILES,"_"),function(x)x[grep(i,x)+1]))))
          if(any(is.na(RES[[i]]))){
            RES[[i]]<-unlist(lapply(strsplit(FILES,"_"),function(x)x[grep(i,x)+1]))
          }
          # position of the name in the file names
          CHECK<-lapply(strsplit(FILES,"_"),function(x)x[grep(i,x)+1])
          names(RES[[i]])<-which(lapply(CHECK,length)!=0)
        }

        # do a linear regression of the parameters one by one to deduce best shot
        par(mfrow=c(2,2))
        par(mar=c(5,4,4,2)+0.1)
        for(i in VARname){
          VAR<-RES[[i]]
          CO<-COST[as.numeric(names(RES[[i]]))]
          if(is.numeric(RES[[i]])){
            plot(CO, VAR, main=i,ylab=i,xlab=paste( ifelse(ValSelect,"Validation","Training"),"Cost") )
            MODEL<-lm(VAR~CO) # bof
            if(!any(is.na(MODEL$coefficients))){
              abline(MODEL)
              PRED<-MODEL$coefficients[1]+MODEL$coefficients[2]*10^-8
              PVALUE<-summary(MODEL)$coeff[2,4]
              COR<-cor(CO,VAR)
              legend("topright",legend = c(paste("pred best =", round(PRED,3)),paste("R2 =",round(COR,3) ),
                                           paste("p =",round(PVALUE,3) )), cex=0.7)
            } else {
              PVALUE=1
              COR<-0.5
              PRED<-0
            }

            if(!any(is.na(c(PVALUE,COR,PRED)))){
              if((PVALUE<0.05|abs(COR)>0.5)&PRED>0){
                # predicted hyperparameters from linear regression
                Default[[i]]<-as.numeric(PRED)

                Bondaries[[i]]<-as.numeric(confint(MODEL,level = 0.2)[1,])

                #MODEL$df.residual
              } else{
                Default[[i]]<-as.numeric(VAR[CO==min(CO)])

                #Bondaries[[i]]<-as.numeric(VAR[tail(order(CO,decreasing = F), 2 )])
                # add the difference to median of the 3 best results
                MedianBest<-median(VAR[tail(order(CO,decreasing = T), 3)])-VAR[CO==min(CO)]
                Bondaries[[i]]<-sort(VAR[CO==min(CO)]+c(MedianBest,-MedianBest))

              }
            }
          } else {
            CO<-COST[as.numeric(names(RES[[i]]))]
            bp<-boxplot(CO~RES[[i]], ylab=paste( ifelse(ValSelect,"Validation","Training"),"Cost"), main = i)
            PRED<-bp$names[bp$stats[3,]==min(bp$stats[3,])]
            Default[[i]]<-PRED

          }
        }
      } else {

        Default<-NULL
        Bondaries<-NULL

      }

    } else {#if(length(FILES)==1) {
      BESTone<-FILES
      print(paste("Selected AMoNet is :", gsub(".Rdata","",BESTone) ))
      load(paste(DIR,BESTone,sep = ""))

      # check it
      NETallpre<-NETallProp$NETallList[[length(NETallProp$NETallList)]]
      NETall<-PlotOptNet(NETall = NETallpre,PDF = F,Optimized = T,
                         NameProj = paste(NameProjbase, "\nmin", ifelse(ValSelect,"Validation","Training"),
                                          "=", round(min(COST),4) ) ,
                         PrintOptNET = T)
    }
  }

  par(mfrow=c(1,1))

  if(!is.null(addOutputs)){
    # addings to NETall
    # outputs
    GenesSelec<-rbind(GenesSelecImmuno,GenesSelecHall)
    MECA<-names_MECA
    NETall<-Outputs(NETall, FamilyGene = GenesSelec[GenesSelec$target_hgnc%in%MECA,],
                    FinalOutput = addOutputs) # FinalOutput=1 for adding a unique output

    # add weigths to new nodes or remove new nodes
    NETallWeights<-addWeights(NETall = NETall[is.na(NETall$Weights),],SdWinit = Default$SdWinit,
                              MeanWinit = Default$MeanWinit ,Scaling = F,Adam = Adam,LSTM = LSTM)
    NETall[is.na(NETall$Weights),]<-NETallWeights

    # how much species are connected to outputs
    CONout<-table(union(NETall$source_hgnc,NETall$target_hgnc)%in%NETall$source_hgnc[NETall$Output])
    print(paste("On",sum(CONout), "species in the net,", as.numeric(CONout["TRUE"]),
                "species are connected to",length(unique(NETall[NETall$Output,3])), "outputs"))
    # unique(NETall[NETall$Output,3])

    # analyse best way to navigate in the network
    # do a sequence of species to update+++
    Analysis<-NetAnalyser(NETall, Propagation = "back",PDF=F) # warnings due to autoregulated nodes (FOXP3): no danger
    NETall$Layer<-Analysis$Layers

    # reload interactions
    NETall$interaction_DIRected_signed<-ifelse(NETall$Weights>0,"ACTIVATE","INHIBIT")

  }


 # PlotOptNet(NETall = NETall,PDF = F,Optimized = F,
 #            NameProj = paste(NameProjbase,ifelse(!is.null(addOutputs), "\n with phenotypes and output", "") ) ,
 #            PrintOptNET = F, LEGEND = F)

  return(list(NETall=NETall, NETallProp=NETallProp, Default=Default, Bondaries=Bondaries))
}
