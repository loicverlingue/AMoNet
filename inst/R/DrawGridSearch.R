#####
# plot the grid search results, estimate best hyperparameters, select best validdation model
# wrapper function to run from command line
# dont run draw if nets of various sizes

#args="--NameProj LUNG_AMoNet_new --Validation T"
# retrieve arguments eithers form command lines or from environments

XV<-try(args)
if("try-error"%in%class(XV)){
  print("Arguments retreived from command line")
  args <- commandArgs(trailingOnly = TRUE)
} else {
  print("Arguments retrieved from R environment")
}

hh <- paste(unlist(args),collapse=' ')
listoptions <- unlist(strsplit(hh,'--'))[-1]
options.args <- sapply(listoptions,function(x){
  unlist(strsplit(x, ' '))[-1]
})
options.names <- sapply(listoptions,function(x){
  option <-  unlist(strsplit(x, ' '))[1]
})
names(options.args) <- unlist(options.names)

#Correct<-options.args[[3]]
options.args<-lapply(options.args,function(Correct){
  if(!is.na(as.numeric(Correct))){
    return(as.numeric(Correct))
  } else if(exists(Correct)){
    return(eval(parse(text = Correct)))
  } else {
    return(Correct)
  }
})

print("options are:")
print(options.args)

if(exists("options.args")){
  list2env(options.args,globalenv())
}
NameProjbase<-NameProj

# set directories - to check
if(!exists("DIR")){
  DIR<-getwd()
}
# check and set directories
if(length(grep("/model",DIR))==0){
    DIR<-file.path(DIR,"model/")
  if(length(list.files(DIR))==0){
    print("Set DIR to root/ or model/")
  }
}

dirData<-DIR
dirPlot<-gsub("/model","/tmp",DIR)

# init matrices
bterms<-list()
Weights<-list()
COSTS<-list()
FILES<-list.files(dirData,pattern = paste("^",NameProjbase, sep = ""))
FILES

### learning curves, cost train, val. No C-index or other metrics used for now

if(TRUE){
  COST<-lapply(FILES,function(f){
    x<-load(file.path(dirData,f))
    net<-get(x)
    COSTsim<-as.numeric(unlist(net$history$Cost))

    MNB<-round(length(unlist(net$TrainSplit$Train))/net$Parameters$Default$MiniBatch)
    #  if(MNB==0){MNB=10}
    if(MNB<1){MNB=1}
    # last cost if batch learning, and median last MNB Costs if mini-batch learning
    COSTmedian<-median(as.numeric(tail(unlist(net$history$Cost),MNB)))

    if("Predict_Val"%in%names(net)){
      #CostVal<-net$Predict_Val$metrics$Cost

      if("Cindex"%in%names(net$Predict_Val$metrics)){
        print("Select on validation 1 - C-index")
        CostVal<-as.numeric(1-net$Predict_Val$metrics$Cindex[1])
      } else {
        print("Select on validation Cost")
        CostVal<-net$Predict_Val$metrics$Cost
      }

    } else {
      CostVal<-0.4
      print(paste("lacks validation cost for", f))
    }

    return(list(COSTsim=COSTsim,COSTmedian=COSTmedian,CostVal=CostVal))
  })

  MedianLast<-unlist(lapply(COST,function(x)x[["COSTmedian"]]))
  names(MedianLast)<-FILES

  CostVal<-unlist(lapply(COST,function(x)x[["CostVal"]]))
  names(CostVal)<-FILES

  COST<-lapply(COST,function(x)x[["COSTsim"]])
  names(COST)<-FILES

}


# if length not the same, harmonize length before do.call:
MAX=max(unlist(lapply(COST,length))) # for weights after
MIN=min(unlist(lapply(COST,length))) # for weights after
if(length(unique(unlist(lapply(COST,length))))>1){
  COST<-lapply(COST,function(x)x[1:MAX])
}

COST<-do.call("rbind",COST)
dim(COST)

# if length not the same:
MINMAX<-apply(COST,1,function(x){max(which(!is.na(x)))})
MIN<-min(MINMAX);MAX<-max(MINMAX)

#Last<-apply(COST,1,function(x){x[max(which(!is.na(x)))]})
#sort(Last,decreasing = T)

### visualization
pdf(paste(dirPlot,"/learningCurves_", NameProjbase, ".pdf",sep = ""))

if(ncol(COST)==1){
  hist(COST,breaks = 100, main="Training cost")
} else {

  # grouping
  Group<-gsub("GridSearch","",gsub(".Rdata","",FILES))
  #Group1<-gsub("_$|_[0-9][0-9]$|_[0-9]$","", gsub("GridSearch","",FILES))
  Group2<-gsub("GridSearch","",gsub(".Rdata","",gsub("\\.[0-9]|[0-9]","", FILES) ))

  # coloring
  RAIN<-rainbow(length(table(Group2)))
  RAIN<-unlist(sapply(seq(length(table(Group2))),function(x){rep(RAIN[x],table(Group2)[x])}))

  # plotting
  par(mfrow=c(1,1))
  par(mar=c(5,4,4,2)+0.1)
  matplot(t(COST),type = 'l',xlim = c(0,MAX), main="Learning curves\n for multiple runs",
          ylab="Training Costs",xlab="iterations*epochs", col=RAIN, lty=1) # , ylim=c(0,0.01)
  legend("topright", legend = names(table(Group2)), lty=1,col=unique(RAIN),cex=0.5)

}

dev.off()

if(FALSE){
# selecting
if(Validation){
  BESTone<-names(CostVal[CostVal==max(CostVal)])
  print(paste("Best training model is:",BESTone))
}else{
  BESTone<-names(MedianLast[MedianLast==min(MedianLast)])
  print(paste("Best training model is:",BESTone))
}
}

######################
# predict the good hyperparameters
### if there are 2 parameters in the grid search
# with extensible list
par(mfrow=c(1,1))

# store good hyperparameters for next steps
#Default<-list(learningrate=0.05,MeanWinit=0.1, SdWinit = 1, MiniBatch=8,SimAnnealMaxSd=0.1,
#              beta1=0.91,beta2=0.91,gradClipping=3,iteration=1,LearningRateDecay="linear")

Default<-AMoNet::Default
Boundaries<-AMoNet::Boundaries

if(TRUE){ # PDF
  pdf(paste(dirPlot,"/PredHyperparamaters_",NameProjbase,".pdf",sep = ""))
  par(mfrow=c(2,2))
  #FILES<-gsub("learning_rate","learningrate",FILES)

  VARname<-unique(gsub(".Rdata","", unlist(strsplit(gsub(NameProjbase,"", FILES),"_"))))
  VARname<-VARname[is.na(as.numeric(VARname))]
  VARname<-VARname[!VARname%in%""]

  RES<-list()
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
  for(i in VARname){
    VAR<-RES[[i]]

    if(!Validation){
      CO<-MedianLast[as.numeric(names(RES[[i]]))]
    }  else {
      CO<-CostVal[as.numeric(names(RES[[i]]))]
    }

    if(is.numeric(RES[[i]])){
      plot(CO,if(all(VAR<1)){log2(VAR)}else{VAR},
           main=i,ylab=if(all(VAR<1)){paste("log",i)}else{i},
           xlab=if(Validation){"Validation"}else{"Training Cost"} )
      #plot(CO, VAR, main=i,ylab=i,xlab="Validation Cost")
      VARmod<-if(all(VAR<1)){log2(VAR)}else{VAR}
      MODEL<-lm(VARmod~CO) # bof
      if(!any(is.na(MODEL$coefficients))){
        abline(MODEL)
        PRED<-MODEL$coefficients[1]+MODEL$coefficients[2]*10^-8
        PVALUE<-summary(MODEL)$coeff[2,4]
        COR<-cor(CO,VAR, method = "spearman")
        legend("topright",legend = c(paste("pred best =", round(PRED,3)),paste("R2 =",round(COR,3) ),
                                     paste("p =",round(PVALUE,3) )), cex=0.7)
      } else {
        PVALUE=1
        COR=0
        PRED=0
      }
      if((PVALUE<0.05|abs(COR)>0.5)&PRED>0){
        Default[[i]]<-as.numeric(PRED)
        Boundaries[[i]]<-as.numeric(confint(MODEL,level = 0.2)[1,])

      } else{
        # select first one to debug when all values are the same, eg beta with minibatch all > than batch size
        Default[[i]]<-as.numeric(VAR[CO==min(CO)])[1]

        MedianBest<-median(VAR[tail(order(CO,decreasing = T), 3)])-VAR[CO==min(CO)]
        MedianBest<-MedianBest[1]
        Boundaries[[i]]<-sort(VAR[CO==min(CO)]+c(MedianBest,-MedianBest))
        Boundaries[[i]]<-Boundaries[[i]][1:2]
      }

    } else {
      CO<-COST[as.numeric(names(RES[[i]])),MIN]
      bp<-boxplot(CO~RES[[i]], ylab="Cost", main = i)
      PRED<-bp$names[bp$stats[3,]==min(bp$stats[3,])]
      Default[[i]]<-PRED
    }
  }

  # by groups / grid

  # discard non numeric Var and associated
  NAS<-unlist(lapply(RES,function(x){all(is.na(as.numeric(x)))}))
  VARname<-VARname[!NAS]
  #VARname<-VARname[!VARname%in%"iteration"]

  # if(length(VARname)>1){

  for(i in VARname){
    # this function retrieve the ith's tween variable to plot it in a grid search mode
    #x=strsplit(gsub(".Rdata","", gsub(NameProjbase,"",FILES)) ,"_")[[32]]
    COMBVAR<-unique(unlist(lapply(strsplit(gsub(".Rdata","", gsub(NameProjbase,"",FILES)) ,"_"),function(x){
      if(any(x%in%i)){
        return(x[x%in%VARname&!x%in%i])
      }
    })))

    if(length(COMBVAR)==1){
      #SIMI<-as.character(intersect(grep(i,FILES),grep(COMBVAR,FILES)))
      SIMI<-as.character(intersect(grep(i,FILES),sapply(COMBVAR,function(z){grep(z,FILES)})))

      plot(RES[[i]][SIMI],RES[[COMBVAR]][SIMI],
           xlab=i, ylab=COMBVAR,
           cex=(1/COST[as.numeric(SIMI),MIN])/10)
    }

  }
  #dev.off()
  write.csv2(t(as.data.frame(unlist(Default))),paste(dirPlot,"/BestHyperparameters_",NameProjbase,".csv",sep = ""), row.names = F)

  Bbond<-t(as.data.frame(unlist(Boundaries)))
  colnames(Bbond)<-reshape2::melt(Boundaries)[,2]
#  colnames(Bbond)<-unlist(lapply(names(Boundaries),function(R)rep(R,2)))
  write.csv2(Bbond,paste(dirPlot,"/BestBoundaries_",NameProjbase,".csv",sep = ""), row.names = F)

}

############
if("Npat"%in%VARname){
  ORD1<-order(RES$Npat,decreasing = F)

  plot(RES$Npat[ORD1],MedianLast[ORD1],ylim=c(0,0.4),
       main = "Learning curve", ylab="Cost", xlab = "N patients")
  points(RES$Npat[ORD1],CostVal[ORD1],col=2)


  Ypredtrain<-predict(loess(MedianLast[ORD1]~RES$Npat[ORD1], span = 1), RES$Npat[ORD1] )
  Xpredtrain<-RES$Npat[ORD1]
  lines(Xpredtrain, Ypredtrain)

  if(length(CostVal[!is.na(CostVal)])>1){
    Ypredval<-predict(loess(CostVal[ORD1][!is.na(CostVal[ORD1])]~RES$Npat[ORD1][!is.na(CostVal[ORD1])]), RES$Npat[ORD1][!is.na(CostVal[ORD1])])
    Xpredval<-RES$Npat[ORD1][!is.na(CostVal[ORD1])]
    lines(Xpredval, Ypredval, col=2)
  }

  legend("topright", legend = c("Train", "Val"), lwd=1, col=c(1,2))

}

dev.off()


####################
# Weights and b terms
# if needed discard worst

#UnWanted<-c("nblayers","MinConnect")
if(any(c("nblayers","MinConnect")%in%VARname)){
  FILES<-FILES[-grep("nblayers|MinConnect", FILES)]
  MedianLast<-MedianLast[-grep("nblayers|MinConnect", names(MedianLast))]

}

if(length(FILES)==0){
#if(FALSE){
#any(c("nblayers","MinConnect")%in%VARname)){
print("Ends before plotting PCA")

} else {

#MAX<-(length(NETallProp$NETallList))
i=1
TwoPoints=F # if you want to follow all the way of learning (F) or just start and end

NETallList<-0
for(f in FILES){
  rm(NETallList)
  x<-load(file.path(dirData,f))
  net<-get(x)
  #select the best net and then select B and W

  NETallList<-net$history$NETallList

  # choose which optimized net you take for Weigths values

  if(TwoPoints){
    for(y in c(1, length(NETallList))){
      bterms[[i]]<-NETallList[[y]]$bterm
      names(bterms[[i]])<-paste(NETallList[[y]]$source_hgnc,NETallList[[y]]$target_hgnc,sep = "_")
      names(bterms)[i]<-paste(gsub(".Rdata","",f),y,sep = "_")

      Weights[[i]]<-NETallList[[y]]$Weights
      names(Weights[[i]])<-paste(NETallList[[y]]$source_hgnc,NETallList[[y]]$target_hgnc,sep = "_")
      names(Weights)[i]<-paste(gsub(".Rdata","",f),y,sep = "_")

      i=i+1
    }
  } else {

    for(y in seq(1,min(length(NETallList),MAX))){ # depending on the learning curve
      bterms[[i]]<-NETallList[[y]]$bterm
      names(bterms[[i]])<-paste(NETallList[[y]]$source_hgnc,NETallList[[y]]$target_hgnc,sep = "_")
      names(bterms)[i]<-paste(gsub(".Rdata","",f),y,sep = "_")

      Weights[[i]]<-NETallList[[y]]$Weights
      names(Weights[[i]])<-paste(NETallList[[y]]$source_hgnc,NETallList[[y]]$target_hgnc,sep = "_")
      names(Weights)[i]<-paste(gsub(".Rdata","",f),y,sep = "_")

      i=i+1
    }
  }
}
#lapply(Weights,length)
WEIGHTS<-t(do.call("rbind",Weights))
bterm<-t(do.call("rbind",bterms))

#############

PCB<-FactoMineR::PCA(t(WEIGHTS),graph = F)

# tag the best (if using the whole FILES)

BEST<-names(tail(sort(MedianLast,decreasing = T),5))
BEST<-gsub(".Rdata","", gsub("GridSearch","",BEST))

# discard non graphic ones

pdf(paste(dirPlot,"/PCAweightsGoodCol_", NameProjbase, ".pdf",sep = ""))

if(TwoPoints){
  COL<-as.vector(sapply(1:(dim(WEIGHTS)[2]/2),function(x)rep(x,2)))
  RAIN<-rainbow(dim(WEIGHTS)[2]/2)

  plot(PCB$ind$coord[,1],PCB$ind$coord[,2],col=RAIN[COL],pch=16,main='PCA Weights',
       xlab=paste("dim1 = ",round(PCB$eig[1,2],2),"%",sep = ""),ylab=paste("dim2 = ",round(PCB$eig[2,2],2),"%",sep = ""))
  arrows(x0 =  PCB$ind$coord[seq(1,length(PCB$ind$coord[,1]),2),1],x1 =  PCB$ind$coord[seq(2,length(PCB$ind$coord[,1]),2),1],
         y0 =  PCB$ind$coord[seq(1,length(PCB$ind$coord[,1]),2),2],y1 =  PCB$ind$coord[seq(2,length(PCB$ind$coord[,2]),2),2],
         col=RAIN,length = 0.1)

  # zoom in
  rownames(PCB$ind$coord)<-gsub("_$|_[0-9][0-9]$|_[0-9]$","", gsub("GridSearch","",rownames(PCB$ind$coord)))
  plot(PCB$ind$coord[,1],PCB$ind$coord[,2],col=RAIN[COL],pch=16,main='PCA Weights',
       xlab=paste("dim1 = ",round(PCB$eig[1,2],2),"%",sep = ""),ylab=paste("dim2 = ",round(PCB$eig[2,2],2),"%",sep = ""),
       xlim=range(PCB$ind$coord[rownames(PCB$ind$coord)%in%BEST,1]), ylim=range(PCB$ind$coord[rownames(PCB$ind$coord)%in%BEST,2]))
  arrows(x0 =  PCB$ind$coord[seq(1,length(PCB$ind$coord[,1]),2),1],x1 =  PCB$ind$coord[seq(2,length(PCB$ind$coord[,1]),2),1],
         y0 =  PCB$ind$coord[seq(1,length(PCB$ind$coord[,1]),2),2],y1 =  PCB$ind$coord[seq(2,length(PCB$ind$coord[,2]),2),2],
         col=RAIN,length = 0.1)

  ORDPCBcoord<-PCB$ind$coord[rownames(PCB$ind$coord)%in%BEST,]
  points(ORDPCBcoord[,1], ORDPCBcoord[,2], cex=2, pch=1)

} else{

  # define 1 group per optimisation
  #Group1<-gsub("_$","", gsub("[0-9]$","", gsub("[0-9]$","", gsub("GridSearch","",rownames(PCB$ind$coord)))))
  Group1<-gsub("_$|_[0-9][0-9]$|_[0-9]$||_[0-9][0-9][0-9]$","", gsub("GridSearch","",rownames(PCB$ind$coord)))
  Group<-c(unlist(sapply(seq(length(table(Group1))),function(REP) rep(REP,length(Group1[Group1%in%names(table(Group1))[REP]]) ))))

  # order by groups
  ORDPCBcoord<-PCB$ind$coord[order(Group),]

  Group<-sort(Group)

  Step<-as.numeric(gsub(".*_","", rownames(ORDPCBcoord)))

  RAIN<-rainbow(max(Group))

  plot(ORDPCBcoord[,1],ORDPCBcoord[,2],col=RAIN[Group],pch=16,main='PCA Weights', type='n',
       xlab=paste("dim1 = ",round(PCB$eig[1,2],2),"%",sep = ""),ylab=paste("dim2 = ",round(PCB$eig[2,2],2),"%",sep = ""))

  sapply((unique(Group)),function(G){
    points( ORDPCBcoord[min(which(Group==G)),1], ORDPCBcoord[min(which(Group==G)),2],col=RAIN[G],pch = 16)
    points( ORDPCBcoord[max(which(Group==G)),1], ORDPCBcoord[max(which(Group==G)),2],col=RAIN[G],pch = 1,
            cex=ifelse(any(BEST%in%Group1[Group==G]),3,2), lwd=ifelse(any(BEST%in%Group1[Group==G]),3,1))
  })

  sapply((unique(Group)),function(G){
    arrows(x0 =  ORDPCBcoord[which(Group==G)[-1]-1,1], x1 = ORDPCBcoord[which(Group==G)[-1],1],
           y0 =  ORDPCBcoord[which(Group==G)[-1]-1,2], y1 = ORDPCBcoord[which(Group==G)[-1],2],
           col=RAIN[G],length = 0.04,
           lwd=ifelse(any(BEST%in%Group1[Group==G]),3,1) )

  })

  #### zoom in

  plot(ORDPCBcoord[,1],ORDPCBcoord[,2],col=RAIN[Group],pch=16,main='PCA Weights', type='n',
       xlim=range(ORDPCBcoord[Group1%in%BEST,1]),ylim=range(ORDPCBcoord[Group1%in%BEST,2]),
       xlab=paste("dim1 = ",round(PCB$eig[1,2],2),"%",sep = ""),ylab=paste("dim2 = ",round(PCB$eig[2,2],2),"%",sep = ""))

  sapply((unique(Group)),function(G){
    points( ORDPCBcoord[min(which(Group==G)),1], ORDPCBcoord[min(which(Group==G)),2],col=RAIN[G],pch = 16)
    points( ORDPCBcoord[max(which(Group==G)),1], ORDPCBcoord[max(which(Group==G)),2],col=RAIN[G],pch = 1,
            cex=ifelse(any(BEST%in%Group1[Group==G]),3,2), lwd=ifelse(any(BEST%in%Group1[Group==G]),3,1))
  })

  sapply((unique(Group)),function(G){
    arrows(x0 =  ORDPCBcoord[which(Group==G)[-1]-1,1], x1 = ORDPCBcoord[which(Group==G)[-1],1],
           y0 =  ORDPCBcoord[which(Group==G)[-1]-1,2], y1 = ORDPCBcoord[which(Group==G)[-1],2],
           col=RAIN[G],length = 0.04,
           lwd=ifelse(any(BEST%in%Group1[Group==G]),3,0.5))
  })

  sapply((unique(Group)),function(G){if(any(BEST%in%Group1[Group==G])){
    COMAX<-MedianLast[gsub(".Rdata","", names(MedianLast)) %in%unique(Group1[Group==G])]
    COVal<-CostVal[gsub(".Rdata","", names(CostVal)) %in%unique(Group1[Group==G])]
    text(ORDPCBcoord[max(which(Group==G)),1],ORDPCBcoord[max(which(Group==G)),2],pos = 1, cex=0.7,
         paste("Training cost =", round(COMAX,3), "\n Validation cost =", round(COVal,3) ))
  }
  })

  #### with density contour plot

  #Group1<-gsub("_$|_[0-9][0-9]$|_[0-9]$||_[0-9][0-9][0-9]$","", gsub("GridSearch","",rownames(ORDPCBcoord)))
  #ORDPCBcoordZOOM<-ORDPCBcoord[Group1%in%BEST,]
  ZOOM<-ORDPCBcoord[,1]>min(ORDPCBcoord[Group1%in%BEST,1])&ORDPCBcoord[,1]<max(ORDPCBcoord[Group1%in%BEST,1])
  ZOOM<-ZOOM|(ORDPCBcoord[,2]>min(ORDPCBcoord[Group1%in%BEST,2])&ORDPCBcoord[,2]<max(ORDPCBcoord[Group1%in%BEST,2]))
  ORDPCBcoordZOOM<-ORDPCBcoord[ZOOM,]
  GroupZOOM<-gsub("_$|_[0-9][0-9]$|_[0-9]$||_[0-9][0-9][0-9]$","", gsub("GridSearch","",rownames(ORDPCBcoordZOOM)))
  ORDPCBcoordZOOM<-ORDPCBcoordZOOM[!duplicated(GroupZOOM,fromLast = T),] # only final

  #dim(ORDPCBcoord);dim(ORDPCBcoordZOOM)
  A<-MASS::kde2d(ORDPCBcoordZOOM[,1],ORDPCBcoordZOOM[,2],n = 500 )
  contour(A,add =F,drawlabels = F,nlevels = 30, col = adjustcolor(4,0.4), axes=F, lwd=1, lty=1,
          xlim=range(ORDPCBcoord[Group1%in%BEST,1]),ylim=range(ORDPCBcoord[Group1%in%BEST,2]))
  #A<-MASS::kde2d(ORDPCBcoordZOOM[,1],ORDPCBcoordZOOM[,2],n = 2 )
  #contour(A,add =T,drawlabels = F,nlevels = 30, col = adjustcolor(3,0.05),axes=F)

  par(new=T)
  plot(ORDPCBcoord[,1],ORDPCBcoord[,2],col=RAIN[Group],pch=16,main='PCA Weights', type='n',
       xlim=range(ORDPCBcoord[Group1%in%BEST,1]),ylim=range(ORDPCBcoord[Group1%in%BEST,2]),
       xlab=paste("dim1 = ",round(PCB$eig[1,2],2),"%",sep = ""),ylab=paste("dim2 = ",round(PCB$eig[2,2],2),"%",sep = ""))

  sapply((unique(Group)),function(G){
    points( ORDPCBcoord[min(which(Group==G)),1], ORDPCBcoord[min(which(Group==G)),2],col=RAIN[G],pch = 16)
    points( ORDPCBcoord[max(which(Group==G)),1], ORDPCBcoord[max(which(Group==G)),2],col=RAIN[G],pch = 1,
            cex=ifelse(any(BEST%in%Group1[Group==G]),3,2), lwd=ifelse(any(BEST%in%Group1[Group==G]),3,1))
  })

  sapply((unique(Group)),function(G){
    arrows(x0 =  ORDPCBcoord[which(Group==G)[-1]-1,1], x1 = ORDPCBcoord[which(Group==G)[-1],1],
           y0 =  ORDPCBcoord[which(Group==G)[-1]-1,2], y1 = ORDPCBcoord[which(Group==G)[-1],2],
           col=RAIN[G],length = 0.04,
           lwd=ifelse(any(BEST%in%Group1[Group==G]),3,0.5))
  })
  mtext("with end points densities",side = 3,line = 0.3)
  sapply((unique(Group)),function(G){if(any(BEST%in%Group1[Group==G])){
    COMAX<-MedianLast[gsub(".Rdata","", names(MedianLast)) %in%unique(Group1[Group==G])]
    COVal<-CostVal[gsub(".Rdata","", names(CostVal)) %in%unique(Group1[Group==G])]
    text(ORDPCBcoord[max(which(Group==G)),1],ORDPCBcoord[max(which(Group==G)),2],pos = 1, cex=0.7,
         paste("Training cost =", round(COMAX,3), "\n Validation cost =", round(COVal,3) ))

    #text(ORDPCBcoord[max(which(Group==G)),1],ORDPCBcoord[max(which(Group==G)),2],pos = 1, cex=0.7,
    #    paste("median cost =", round(COMAX,3) ))
  }
  })
}

#legend("topleft",legend = unique(Group1[order(Group)]), col = RAIN, lty = 1,cex=0.4 )
dev.off()

#############
## if pooled representation, to check how is who
#sapply((unique(Group)),function(G){
#    text(ORDPCBcoord[max(which(Group==G)),1],ORDPCBcoord[max(which(Group==G)),2],pos = 4, cex=0.6,
#         gsub("_.*","", unique(Group1[which(Group==G)])), col=RAIN[G],  srt=90)
#})

##############
### select ~10 (or <x) best ones and find the most clustered ones
if(!Validation){
  SELECT=min(4,length(MedianLast))
  BEST<-names(tail(sort(MedianLast,decreasing = T),SELECT))
  VeryBEST<-names(tail(sort(MedianLast,decreasing = T),1))
}else{
  SELECT=min(4,length(CostVal))
  BEST<-names(tail(sort(CostVal,decreasing = T),SELECT))
  VeryBEST<-names(tail(sort(CostVal,decreasing = T),1))
}

BEST<-gsub(".Rdata","", gsub("GridSearch","",BEST))
VeryBEST<-gsub(".Rdata","", gsub("GridSearch","",VeryBEST))
FILES<-FILES[FILES%in%paste(BEST,".Rdata",sep = "")]

bterms<-list()
Weights<-list() # and redo file with weigths


i=1
TwoPoints=F # if you want to follow all the way of learning (F) or just start and end
NETallList<-0
for(f in FILES){
  rm(NETallList)
  x<-load(file.path(dirData,f))
  net<-get(x)

  NETallList<-net$history$NETallList

  # choose which optimized net you take for Weigths values

  if(TwoPoints){
    for(y in c(1, length(NETallList))){ # c(1, length(NETallList))
      bterms[[i]]<-NETallList[[y]]$bterm
      names(bterms[[i]])<-paste(NETallList[[y]]$source_hgnc,NETallList[[y]]$target_hgnc,sep = "_")
      names(bterms)[i]<-paste(gsub(".Rdata","",f),y,sep = "_")

      Weights[[i]]<-NETallList[[y]]$Weights
      names(Weights[[i]])<-paste(NETallList[[y]]$source_hgnc,NETallList[[y]]$target_hgnc,sep = "_")
      names(Weights)[i]<-paste(gsub(".Rdata","",f),y,sep = "_")

      #    COSTS[[i]]<-as.numeric(unlist(NETallList[[y]]$COST))
      i=i+1
    }
  } else {

    for(y in seq(1,min(length(NETallList),MAX))){ # depending on the learning curve
      bterms[[i]]<-NETallList[[y]]$bterm
      names(bterms[[i]])<-paste(NETallList[[y]]$source_hgnc,NETallList[[y]]$target_hgnc,sep = "_")
      names(bterms)[i]<-paste(gsub(".Rdata","",f),y,sep = "_")

      Weights[[i]]<-NETallList[[y]]$Weights
      names(Weights[[i]])<-paste(NETallList[[y]]$source_hgnc,NETallList[[y]]$target_hgnc,sep = "_")
      names(Weights)[i]<-paste(gsub(".Rdata","",f),y,sep = "_")

      #    COSTS[[i]]<-as.numeric(unlist(NETallList[[y]]$COST))
      i=i+1
    }
  }
}
WEIGHTS<-t(do.call("rbind",Weights))
bterm<-t(do.call("rbind",bterms))
dim(WEIGHTS)

# a function to select the larger cluster in term of density
CLUSTCUT<-function(DEN=DENX, Cutoffdens=0.8){

  #  abline(h=quantile(DEN$y,Cutoffdens))
  POS<-which(DEN$y>quantile(DEN$y,Cutoffdens))

  POScutup<-POS[seq(2,length(POS)+1,length.out = length(POS))]-POS
  POScutdn<-POS-c(NA,POS[-length(POS)])

  GPPOS<-sort(POS[c(which(POScutup!=1),which(POScutdn!=1))])
  GPPOS<-c(POS[1],GPPOS,POS[length(POS)])
  plot(DEN)
  points(DEN$x[GPPOS],rep(quantile(DEN$y,Cutoffdens),length(GPPOS)) )

  if(length(GPPOS)>2){
    CLUST<-sapply(GPPOS[seq(1,length(GPPOS),2)],function(x){
      DEN$y[x:GPPOS[which(GPPOS==x)+1]]
    })

    # select the highest density
    POScutdn<-GPPOS[seq(1,length(GPPOS),2)][lapply(CLUST,max)==max(as.numeric(lapply(CLUST,max)))]
    POScutup<-GPPOS[which(GPPOS%in%POScutdn)+1]
    # or select the one with lowest cost?


    abline(v=DEN$x[POScutdn],col=2)
    abline(v=DEN$x[POScutup],col=2)
    return(c(DEN$x[POScutdn],DEN$x[POScutup]))
  } else {
    return(DEN$x[range(POS)])
  }
}

# do the PCA and grouping
PCB<-FactoMineR::PCA(t(WEIGHTS),graph = F)
Group1<-gsub("_$|_[0-9][0-9]$|_[0-9]$","", gsub("GridSearch","",rownames(PCB$ind$coord)))
Group<-c(unlist(sapply(seq(length(table(Group1))),function(REP) rep(REP,length(Group1[Group1%in%names(table(Group1))[REP]]) ))))
ORDPCBcoord<-PCB$ind$coord[order(Group),]
Group<-sort(Group)
Step<-as.numeric(gsub(".*_","", rownames(ORDPCBcoord)))
RAIN<-rainbow(max(Group))

# X
DENX<-density(ORDPCBcoord[,1], adjust=1)
Xselect<-CLUSTCUT(DENX,Cutoffdens = 0.9)

DENY<-density(ORDPCBcoord[,2], adjust=1)
Yselect<-CLUSTCUT(DENY,Cutoffdens = 0.8)

#unique(c(unique(Group1[Xselect[1]:Xselect[2]]),unique(Group1[Yselect[1]:Yselect[2]])))
if(TRUE){
  pdf(paste(dirPlot,"/PCAweights_BestTrain_",NameProjbase,".pdf",sep = ""))

  plot(ORDPCBcoord[,1],ORDPCBcoord[,2],col=RAIN[Group],pch=16,
       main='PCA Weights to select best & clustered nets', type='n',
       xlab=paste("dim1 = ",round(PCB$eig[1,2],2),"%",sep = ""),ylab=paste("dim2 = ",round(PCB$eig[2,2],2),"%",sep = ""))

  sapply((unique(Group)),function(G){
    points( ORDPCBcoord[min(which(Group==G)),1], ORDPCBcoord[min(which(Group==G)),2],col=RAIN[G],pch = 16)
    points( ORDPCBcoord[max(which(Group==G)),1], ORDPCBcoord[max(which(Group==G)),2],col=RAIN[G],pch = 1,
            cex=ifelse(any(BEST%in%Group1[Group==G]),3,2), lwd=ifelse(any(BEST%in%Group1[Group==G]),3,1))
  })

  sapply((unique(Group)),function(G){
    arrows(x0 =  ORDPCBcoord[which(Group==G)[-1]-1,1], x1 = ORDPCBcoord[which(Group==G)[-1],1],
           y0 =  ORDPCBcoord[which(Group==G)[-1]-1,2], y1 = ORDPCBcoord[which(Group==G)[-1],2],
           col=RAIN[G],length = 0.04,
           lwd=ifelse(any(BEST%in%Group1[Group==G]),3,1) )

  })

  SELEC<-intersect(unique(Group1[ORDPCBcoord[,1]>Xselect[1]&ORDPCBcoord[,1]<Xselect[2]]),
                   unique(Group1[ORDPCBcoord[,2]>Yselect[1]&ORDPCBcoord[,2]<Yselect[2]]))

  if(!VeryBEST%in%SELEC){
    SELEC<-VeryBEST
  } else{
    rect(xleft = Xselect[1],ybottom = Yselect[1],xright = Xselect[2],ytop = Yselect[2],border = 2)
  }

  sapply((unique(Group)),function(G){
    if(any(SELEC%in%Group1[Group==G])){
      text(ORDPCBcoord[max(which(Group==G)),1],ORDPCBcoord[max(which(Group==G)),2],pos = 2, cex=0.7,
           unique(Group1[which(Group==G)]), col=RAIN[G])

      COMAX<-MedianLast[ gsub(".Rdata","",names(MedianLast)) %in%unique(Group1[Group==G])] #gsub(NameProjbase,"", rownames(COST)))
      COVal<-CostVal[gsub(".Rdata","", names(CostVal)) %in%unique(Group1[Group==G])]
      text(ORDPCBcoord[max(which(Group==G)),1],ORDPCBcoord[max(which(Group==G)),2],pos = 1, cex=0.7,
           paste("Training cost =", round(COMAX,3), "\n Validation cost =", round(COVal,3) ))

      #text(ORDPCBcoord[max(which(Group==G)),1],ORDPCBcoord[max(which(Group==G)),2],pos = 4, cex=0.7,
      #     paste("Cost =", round(COMAX,3) ))

    }
  })

  # plot the closest nets
  par(mfrow=c(2,2))
  par(mar=c(rep(1,4)))
  FILES[gsub(".Rdata","",FILES)%in%SELEC]
  #f<-FILES[gsub(".Rdata","",FILES)%in%SELEC][1]
  for(f in FILES[gsub(".Rdata","",FILES)%in%SELEC]){
    rm(NETallList)
    load(file.path(dirData,f))
    NETallList<-net$history$NETallList

    NETall<-NETallList[[length(NETallList)]]
    NETallopt<-PlotOptNet(NETall,PDF = F,Optimized = T,NameProj = gsub(".Rdata" ,"",f),PrintOptNET = T,LEGEND = F)

    # check edges that are reclassified
    # final format edges:
    INT<-paste(NETallopt$source_hgnc,NETallopt$interaction_directed_signed, NETallopt$target_hgnc,sep = "_")

    INT[NETall$interaction_directed_signed%in%"ACTIVATE"&NETall$Weights<0]
    INT[NETall$interaction_directed_signed%in%"INHIBIT"&NETall$Weights>0]

  }


  dev.off()
}
}

print("end")
