#TotAttractors<-TotAttractorsMult
DrawPath<-function(TotAttractors, Reduction=50, NameProj="transition_graph", PDF=T){

  # if state transition too big, reduce it
  # but can impare the detection of last species before attraction
#    if(min(do.call("rbind",lapply(TotAttractors$LAttractors,dim))[,1])>Reduction){
#      TotAttractors$LAttractors<-lapply(TotAttractors$LAttractors,function(x){
#        x<-rbind(x[seq(1,nrow(x)-20,length.out = 20),],x[(nrow(x)-20):nrow(x),])
#        return(x)
#      })
#    }

  # true attractors
  # all unique possible states
  #LAttractors<-lapply(TotAttractors$LAttractors,function(x){
  #  apply(x,1,function(y){paste(y,collapse = "")})
  #})

  LAttractors<-lapply(seq(dim(TotAttractors)[4]), function(x){
      apply(TotAttractors[,"A",,x],1,function(y){paste(y,collapse = "")})
    })
    length(LAttractors)

  # give a number to each state
  TotiStateUnique<-unique(as.character(unlist(LAttractors)))
  TotiStateUnique<-data.frame(TotiStateUnique,NUM=seq(length(TotiStateUnique)))

  #nrow(TotiStateUnique)
  #length(unique(TotiStateUnique$TotiStateUnique))

  # number of the best attractors
  #  NumBestAttractors<-TotiStateUnique$NUM[TotiStateUnique$TotiStateUnique%in%rownames(BestAttract)]

  #x<-TotAttractors1$LAttractors$`1`

  # compute all the path
  # put all path the same number of steps by repeting the last state for those who finished early
  MAX<-max(unlist(lapply(seq(dim(TotAttractors)[4]),function(x){
    nrow(TotAttractors[,"A",,x])
    })))

  # parallelise!
 # x<-TotAttractors$LAttractors[[1]]
#  Path<-lapply(TotAttractors$LAttractors,function(x){
#    PathAtname<-as.character(apply(x,1,function(y){paste(y,collapse = "")}))
#    PathAt<-as.numeric(sapply(PathAtname,function(h){which(TotiStateUnique$TotiStateUnique%in%h)}))
#    Attract<-names(table(PathAt))[table(PathAt)%in%max(table(PathAt))]
#    Attract<-names(table(PathAt))[table(PathAt)>1]
#    Path1<-c(PathAt,rep(PathAt[length(PathAt)], MAX-length(PathAt))) # repeat the last state
#    GeneTransition=c(rownames(x),rep(tail(rownames(x),1),MAX-length(rownames(x)) ))
#    list(Path1=Path1, GeneTransition=GeneTransition, Attract=Attract)
#  })

  Path<-lapply(seq(dim(TotAttractors)[4]),function(x){
    PathAtname<-as.character(apply(TotAttractors[,"A",,x],1,function(y){paste(y,collapse = "")}))
    PathAt<-as.numeric(sapply(PathAtname,function(h){which(TotiStateUnique$TotiStateUnique%in%h)}))
    Attract<-names(table(PathAt))[table(PathAt)%in%max(table(PathAt))]
    Attract<-names(table(PathAt))[table(PathAt)>1]

    Path1<-c(PathAt,rep(PathAt[length(PathAt)], MAX-length(PathAt))) # repeat the last state
    GeneTransition=c(rownames(TotAttractors[,"A",,x]),rep(tail(rownames(TotAttractors[,"A",,x]),1),MAX-length(rownames(TotAttractors[,"A",,x])) ))
    list(Path1=Path1, GeneTransition=GeneTransition, Attract=Attract)
  })

  # check species activities for each attractors == already done in PlotAttractors
#  ATT<-which( PathAtname%in%TotiStateUnique$TotiStateUnique[TotiStateUnique$NUM%in%as.numeric(Attract)])
#  SpeciesActivityATT<-TotAttractors$LAttractors[[length(TotAttractors$LAttractors)]][ATT,]
#  apply(SpeciesActivityATT,2,table)

  DrawPath1<-function(Path=Path, NameProj=NameProj){
   # setwd(gsub("/tmp","",getwd()))
  #  setwd(paste(getwd(),"/tmp/",sep = ""))
  #  getwd()


    #  GeneTransition<-lapply(Path,function(x)x[["GeneTransition"]])
    #  Attract<-lapply(Path,function(x)as.numeric(x[["Attract"]]))
    Path<-lapply(Path,function(x)x[["Path1"]])

    # store all path in data frame
    Path<-as.data.frame(Path)
    Path[,ncol(Path)+1]<-seq(nrow(Path))
    Path<-as.matrix(Path)

    ##### do the representation of the trajectories
    if(PDF){
    Latt<-length(list.files(paste(getwd(),"/tmp/",sep = ""),
                            pattern =  paste( NameProj, "_transition_graph", sep = "") ))
    NAME<-paste(NameProj, "_transition_graph",Latt+1,sep = "")

    pdf(paste(getwd(),"/tmp/",NAME,".pdf",sep = ""))
    }

    par(mar=c(5, 4, 4, 2) + 0.1)
    matplot(Path[,-ncol(Path)],type="l",
            main="Trajectories",xlab="Time steps", ylab="State trasitions",axes=F)
    arrows(x0 = 0,y0 = -max(Path)/10,x1 = nrow(Path),y1 = -max(Path)/10,length = 0.2,xpd=T)
    axis(2,las=2)

    matplot(Path,type="l",xlim = c(13,18),
            main="Trajectories",xlab="Time steps", ylab="State trasitions",axes=F)

    if(PDF){
    dev.off()
    }
    # can be used to locate the most important transitions / the contour of the bassins

  }
 # if(PDF){
    DrawPath1(Path)
#  }
  return(Path)
}

# dissect path and do 3D plot
DissectPath<-function(Path, NTrans=3, Smooth=5, Cordplot=T, plot3D=T, PDF=F){
#  setwd(gsub("/tmp","",getwd()))
#  setwd(paste(getwd(),"/tmp/",sep = ""))
#  getwd()

  # load function for 3D
  Plot3DAttract<-function(TopTrans=TopTrans,Smooth=Smooth){

    TopT<-TopTrans[TopTrans$OverallFreq>0,]
    TopT<-TopT[order(TopT$To),] # try to give meaning to the coordinates....

    #TopT<-TopT[order(TopT$To),]

    kd<-MASS::kde2d(x = TopT$OverallFreq,y=TopT$To,h = Smooth)
    #kd<-kde2d(x = TopT$OverallFreq,y=seq(nrow(TopT)),h = 5)

    # allTrans1<-t(matrix(as.numeric(unlist(strsplit( as.character(Tab$allTrans)," " ) )),nrow = 2))
  #  kd<-kde2d(x = allTrans1[,1],y=allTrans1[,2],h = 500)

    # kd<-with(TopT,MASS::kde2d(seq(nrow(TopT)),OverallFreq))

    #kd <- with(MASS::geyser, MASS::kde2d(duration, waiting, n = 50))

    zlim <- range(kd$y)
    zlen <- zlim[2] - zlim[1] + 1
    colorlut <- terrain.colors(10) # height color lookup table
    col <- cm.colors(20)[1 + round(19*(kd$z - min(kd$z))/diff(range(kd$z)))]

    rgl::open3d()
    rgl::surface3d(x = kd$x,y = kd$y,z = -kd$z, col=col)
    rgl::aspect3d(1,1,1)
    rgl::decorate3d(main = "Best attractors", xlab = "Frequency", ylab = "Proximity", zlab ="Kernel Density")

    Latt<-length(list.files(getwd(),pattern =  "BestAttractLog"))
    NAME<-paste("BestAttractLog",Latt+1,sep = "")
    rgl::movie3d(rgl::spin3d(axis = c(0.5,0,1), rpm = 10), duration=6, movie = NAME,dir = getwd())
  }

  #

  GeneTransition<-lapply(Path,function(x)x[["GeneTransition"]])
  Attract<-lapply(Path,function(x)as.numeric(x[["Attract"]]))
  Path<-lapply(Path,function(x)x[["Path1"]])

  #
  #dist(do.call("rbind",GeneTransition))

  #hm<-hclust(dist(as.matrix(Path)))

  # store all path in data frame
  Path<-as.data.frame(Path)
  Path[,ncol(Path)+1]<-seq(nrow(Path))
  Path<-as.matrix(Path)

  #D<-dist(as.matrix(Path))
  #D<-as.matrix(D)
  #col <- cm.colors(20)[1 + round(19*(D - min(D))/diff(range(D)))]
  #surface3d(x = seq(ncol(D)),y = seq(nrow(D)),z = D,col=col)
  #aspect3d(1,1,1)

  # compute list of vectors of all possible transitions from A to B
  ROWPath<-apply(Path,1,function(x)paste(x,collapse = ""))
  AllPath<-apply(Path,1,function(x){
    rbind(Path[grep(paste(x,collapse = ""),ROWPath)-1,-ncol(Path)],x[-ncol(Path)])
  })
  allTrans<-lapply(AllPath,function(x){
    apply(x,2,function(y){paste(y,collapse = " ")})
  })
  allTrans<-unlist(allTrans)
  allTrans<-allTrans[grep(" ",allTrans)]

  # do a table, estimate frequencies by transition
  allTrans<-data.frame(allTrans,WhichTotPath = names(allTrans))
  Tab<-table(allTrans)
  Tab<-data.frame(Tab)
  TabTrans<-table(allTrans$allTrans)
  Tab$Freq<-TabTrans[Tab$allTrans]

  Tab<-cbind(t(matrix(as.numeric(unlist(strsplit( as.character(Tab$allTrans)," " ) )),nrow = 2)),Tab$Freq)
  colnames(Tab)<-c("From","To","OverallFreq")

  # the frequent transitions and the attractors:

  # attractors sorted by frequency
  #AttractNames<-unique(unlist(Attract))
  AttractNames<-as.numeric(names(sort(table(unlist(Attract)),decreasing = T)))

  TopTrans<-Tab[Tab[,2]%in%AttractNames,]

  for(i in seq(NTrans)){
    TopTrans<-rbind(TopTrans, Tab[Tab[,2]%in%TopTrans[,1],]) # identify and select thoses just before
  }

  TopTrans<-TopTrans[!duplicated(paste(TopTrans[,1],TopTrans[,2],TopTrans[,3])),] #remove duplicates
  TopTrans<-TopTrans[!TopTrans[,1]==TopTrans[,2],] # remove the transition within a single path

  TopTrans<-as.data.frame(TopTrans)
  TopTrans$OverallFreq<-log(TopTrans$OverallFreq) # only none zero relations
  #TopTrans$Attract<-0
  #TopTrans$Attract[TopTrans$To%in%AttractNames]<-1

  # cord plot
  if(Cordplot){
    library(circlize)
    if(PDF){
    Latt<-length(list.files(getwd(),pattern =  "CordPlotTransitionStates"))
    NAME<-paste("CordPlotTransitionStates",Latt+1,sep = "")
    pdf(paste(getwd(),"/",NAME,".pdf",sep = ""))
    }
    par(mar=c(5, 4, 4, 2) + 0.1)
    circlize::chordDiagram(TopTrans,directional = 1,direction.type = "arrows",col=ifelse(TopTrans$To%in%AttractNames,3,4))
    legend("bottomright",legend = c("Attractors","Transitions"),fill = c(3,4),cex = 0.8)
    if(PDF){
    dev.off()
    }
  }

  if(plot3D){
    # do 3D plot if there are various frequencies
    if(length(unique(TopTrans[TopTrans$OverallFreq>0,"OverallFreq"]))!=1){
      library(MASS);library(rgl)
      Plot3DAttract(TopTrans)
    } else{
      print("Frequencies of transition are all the same: no 3D plot feasible")
    }
   }

  #    TopTrans[TopTrans$OverallFreq>0,"From"]
  #    TopTrans[TopTrans$OverallFreq>,"To"]
  #identify which species are involved in last paths
  GeneTransition<-as.data.frame(GeneTransition)
  #x<-Path[,98]
  TopGeneTransition<-apply(Path[,-ncol(Path)],2,function(x){
    NCOL<-grep(paste(x,collapse = ""),apply(Path,2,function(h)paste(h,collapse = "")))
    #AttractNames<-unique(unlist(Attract))
    Attractor<-AttractNames[AttractNames%in%x][1]

    # positions of unique states and attractors, but some >1
    States<-which(!duplicated(x))
    #PosAttract<-grep( Attractor,x)[1]
    PosAttract<-which(x%in%Attractor)[1]

    if(grep(PosAttract,States)>NTrans){
      Last<-gsub("_[0-9].*","",GeneTransition[States[seq(grep(PosAttract,States)-NTrans,grep(PosAttract,States))],NCOL])
      Attractor<-x[States[seq(grep(PosAttract,States)-NTrans,grep(PosAttract,States))]]
    } else {
      Last<-gsub("_[0-9].*","",GeneTransition[States[seq(0,grep(PosAttract,States))],NCOL])
      Attractor<-x[States[seq(0,grep(PosAttract,States))]]
    }
    #  Attractor<-Path[States[seq(grep(PosAttract,States)-3,grep(PosAttract,States))],NCOL]
    #  Last<-gsub("_[0-9].*","",GeneTransition[seq(grep(Attractor,x)[1]-3,grep(Attractor,x)[1]),NCOL])
    #  Attractor<-Path[seq(grep(Attractor,x)[1]-3,grep(Attractor,x)[1]),NCOL]
    list(Last=Last,Attractor=Attractor)
  })

  TopGeneTransition<-as.data.frame(do.call("rbind",TopGeneTransition))
  Transition<-as.data.frame(do.call("rbind",TopGeneTransition$Last))
  Transition<-cbind(Transition,as.data.frame(do.call("rbind",TopGeneTransition$Attractor)))
  colnames(Transition)<-paste("V",seq(ncol(Transition)),sep = "")
  #Transition[1:2,]

  SpeciesImportance<-head(sort(table(as.character(unlist(
    Transition[,which(!sapply(Transition[1,],is.numeric))]))),
    decreasing = T),4) # or NTrans?

  # careful, 2 best attractors maybe too low in some cases
  BestTransition<-Transition[Transition[,ncol(Transition)]%in%AttractNames[1:2],]
  SpeciesImportanceAttract<-head(sort(table(as.character(unlist(
    BestTransition[,which(!sapply(BestTransition[1,],is.numeric))]))),
    decreasing = T),4)

  return(list(Transition=Transition,SpeciesImportance=SpeciesImportance,SpeciesImportanceAttract=SpeciesImportanceAttract))
}


