
# compute descrete survival variable for learning

#Ytrain<-YTCGA
#colnames(Ytrain)<-c("time","status")
#XtrainN<-XTCGA
#Interval=5

#time=Y_base_raw[1,];status=Y_base_raw[2,]
#time=Ytrain$time;status=Ytrain$status
#plot(survfit(Surv(time,status)~NULL), mark.time = T)
#Interval=nrow(Ytrain)
SurvEstimator<-function(Interval=30, time=time, status=status, PDF=F){
  
  Interval=Interval+1
  L=seq(min(time),max(time),length.out = (Interval)-1)
#  L=c(1,L) # add a first interval
  L=c(0,L) # add a first interval
  ProbTable<-data.frame(TimeInt=L)
  
  # cumulative number of events to calculate the nb of subject at risk remaining
  Nevent<-sapply(seq(Interval),function(l){
    subjectRisk<-time>=L[l]
    Nevent<-as.numeric(table(subjectRisk)["TRUE"])
    return(Nevent)  
  })
  
  ProbTable[,"NatRisk"]<-Nevent
  
  # nb of subject censored
  Ncensored<-sapply(seq(Interval),function(l){
    MAX<-ifelse(is.na(L[l+1]),L[l],L[l+1])
    
    Int<-time<MAX &  time>=L[l]
    Ncensored<-ifelse(any(status[Int]==0), as.numeric(table(status[Int]==0)["TRUE"]), 0)
    ProbTable[l,"Ncensored"]<<-Ncensored
    Ndeath<-ifelse(any(status[Int]==1), as.numeric(table(status[Int]==1)["TRUE"]), 0)
    ProbTable[l,"Ndeath"]<<-Ndeath
    
    return(Ncensored)  
  })
  ProbTable[,"Ncensored"]<-Ncensored
#  head(ProbTable)
  
  # probability of survival relative to intervals
  ProbTable[1,"Prob"]<-1
  for(l in seq(Interval)[-1]){
    
    Prob<-ProbTable$Prob[l-1]*( (ProbTable$NatRisk[l]-ProbTable$Ndeath[l]) / ProbTable$NatRisk[l] )
    ProbTable[l,"Prob"]<-Prob
  }
  
  # if more intervals than events, no because not real time scale then with matrix
 # ProbTable<-ProbTable[!duplicated(ProbTable$NatRisk),]
  
  if(PDF){
    plot(ProbTable$TimeInt, ProbTable$Prob, type = 'l', main="Kaplan Meier curve",
         xlab="time", ylab="Probability of survival", ylim=c(0,1))
    plot(ProbTable$Prob, type = 'l', main="Kaplan Meier curve", 
         xlab="Descrete time intervals", ylab="Probability of survival", ylim=c(0,1))
    #lines(ProbTable$Prob,col=2)
  }
  
  # survival function
  St<-ProbTable$Prob
  
  #discrete probability function
  Ft<-St[-length(St)]-St[-1]
  Ft=c(0,Ft)
  
  if(PDF){
    plot(Ft, main="Discrete probability function", xlab="Descrete time intervals")
    #plot(Ft,St[-1])
  }
  
  #discrete hazard rate = conditional failure probability
  Ht=Ft/St
  ProbTable$Ht<-Ht
  
  # St should be = prod 1-ht at interval l
 # Ht<-c(0,Ht)
  if(PDF){
    plot(Ht,main="Ht")
  }
  #plot(St,cumprod(1-Ht)) ##############
  #plot(Ft,Ht*cumprod(1-Ht)) #############
  # ok!
  
  ProbTable<-ProbTable[-1,]
  
  ############ patients' based S(t)
  PatProb<-sapply(time , function(Time){
    1-ProbTable$Prob[max(which(ProbTable$TimeInt<=Time))]
  })
  
  ############# survival estimator
  i=which(time==1)[1]
  
  IntPat<-sapply(seq(length(time)),function(i){  
    Proba<-vector()
    time_interval<-max(which(ProbTable$TimeInt<=time[i]))
    
    #for( l in seq(Interval)){
     for( l in seq(nrow(ProbTable))){
        # uncensored
      Pu<-ifelse(time_interval>=l,1,0)
      
      # censored
      Pc<-ifelse(time_interval>=l,1,St[l])
      Proba<-c(Proba,ifelse(status[i]==0,Pc,Pu))
    }
    return(Proba)  
  })
  #dim(IntPat)
  #image(t(IntPat))
  
  
  ######### in-house approx
  normalized<-function(x) {(x-min(x))/(max(x)-min(x))}
  normTime<-normalized(time)
  
  IntPatC3<-sapply(seq(length(time)),function(x){
    if(status[x]==1){
      Pc = normTime[x]
    } else {
      time_interval=max(which(ProbTable$TimeInt<=time[x]))
#      Pc = normTime[x] + prod(1-Ht[time_interval:length(Ht)]) ###?
      Pc = normTime[x] + Ft[time_interval]
    }
    return(Pc)
  })
  #plot(IntPatC3,time,pch=ifelse(status==0,3,1))

  ############### mean of intervals status
  MeanProbMat<-colMeans(IntPat)
  
  
  
  return(list(ProbTable=ProbTable,ProbMat=IntPat,PatProb=PatProb,ProbVec=IntPatC3,MeanProbMat=MeanProbMat))
}

# compute descrete survival variable for learning

#Ytrain<-YTCGA
#colnames(Ytrain)<-c("time","status")
#XtrainN<-XTCGA
#Interval=5

#time=Y_base_raw[1,];status=Y_base_raw[2,]
#time=Ytrain$time;status=Ytrain$status
#plot(survfit(Surv(time,status)~NULL), mark.time = T)
#Interval=nrow(Ytrain)
SurvEstimator<-function(Interval=30, time=time, status=status, PDF=F){
  
  Interval=Interval+1
  L=seq(min(time),max(time),length.out = (Interval)-1)
#  L=c(1,L) # add a first interval
  L=c(0,L) # add a first interval
  ProbTable<-data.frame(TimeInt=L)
  
  # cumulative number of events to calculate the nb of subject at risk remaining
  Nevent<-sapply(seq(Interval),function(l){
    subjectRisk<-time>=L[l]
    Nevent<-as.numeric(table(subjectRisk)["TRUE"])
    return(Nevent)  
  })
  
  ProbTable[,"NatRisk"]<-Nevent
  
  # nb of subject censored
  Ncensored<-sapply(seq(Interval),function(l){
    MAX<-ifelse(is.na(L[l+1]),L[l],L[l+1])
    
    Int<-time<MAX &  time>=L[l]
    Ncensored<-ifelse(any(status[Int]==0), as.numeric(table(status[Int]==0)["TRUE"]), 0)
    ProbTable[l,"Ncensored"]<<-Ncensored
    Ndeath<-ifelse(any(status[Int]==1), as.numeric(table(status[Int]==1)["TRUE"]), 0)
    ProbTable[l,"Ndeath"]<<-Ndeath
    
    return(Ncensored)  
  })
  ProbTable[,"Ncensored"]<-Ncensored
#  head(ProbTable)
  
  # probability of survival relative to intervals
  ProbTable[1,"Prob"]<-1
  for(l in seq(Interval)[-1]){
    
    Prob<-ProbTable$Prob[l-1]*( (ProbTable$NatRisk[l]-ProbTable$Ndeath[l]) / ProbTable$NatRisk[l] )
    ProbTable[l,"Prob"]<-Prob
  }
  
  # if more intervals than events, no because not real time scale then with matrix
 # ProbTable<-ProbTable[!duplicated(ProbTable$NatRisk),]
  
  if(PDF){
    plot(ProbTable$TimeInt, ProbTable$Prob, type = 'l', main="Kaplan Meier curve",
         xlab="time", ylab="Probability of survival", ylim=c(0,1))
    plot(ProbTable$Prob, type = 'l', main="Kaplan Meier curve", 
         xlab="Descrete time intervals", ylab="Probability of survival", ylim=c(0,1))
    #lines(ProbTable$Prob,col=2)
  }
  
  # survival function
  St<-ProbTable$Prob
  
  #discrete probability function
  Ft<-St[-length(St)]-St[-1]
  Ft=c(0,Ft)
  
  if(PDF){
    plot(Ft, main="Discrete probability function", xlab="Descrete time intervals")
    #plot(Ft,St[-1])
  }
  
  #discrete hazard rate = conditional failure probability
  Ht=Ft/St
  ProbTable$Ht<-Ht
  
  # St should be = prod 1-ht at interval l
 # Ht<-c(0,Ht)
  if(PDF){
    plot(Ht,main="Ht")
  }
  #plot(St,cumprod(1-Ht)) ##############
  #plot(Ft,Ht*cumprod(1-Ht)) #############
  # ok!
  
  ProbTable<-ProbTable[-1,]
  
  ############ patients' based S(t)
  PatProb<-sapply(time , function(Time){
    1-ProbTable$Prob[max(which(ProbTable$TimeInt<=Time))]
  })
  
  ############# survival estimator
  i=which(time==1)[1]
  
  IntPat<-sapply(seq(length(time)),function(i){  
    Proba<-vector()
    time_interval<-max(which(ProbTable$TimeInt<=time[i]))
    
    #for( l in seq(Interval)){
     for( l in seq(nrow(ProbTable))){
        # uncensored
      Pu<-ifelse(time_interval>=l,1,0)
      
      # censored
      Pc<-ifelse(time_interval>=l,1,St[l])
      Proba<-c(Proba,ifelse(status[i]==0,Pc,Pu))
    }
    return(Proba)  
  })
  #dim(IntPat)
  #image(t(IntPat))
  
  
  ######### in-house approx
  normalized<-function(x) {(x-min(x))/(max(x)-min(x))}
  normTime<-normalized(time)
  
  IntPatC3<-sapply(seq(length(time)),function(x){
    if(status[x]==1){
      Pc = normTime[x]
    } else {
      time_interval=max(which(ProbTable$TimeInt<=time[x]))
#      Pc = normTime[x] + prod(1-Ht[time_interval:length(Ht)]) ###?
      Pc = normTime[x] + Ft[time_interval]
    }
    return(Pc)
  })
  #plot(IntPatC3,time,pch=ifelse(status==0,3,1))

  ############### mean of intervals status
  MeanProbMat<-colMeans(IntPat)
  
  
  
  return(list(ProbTable=ProbTable,ProbMat=IntPat,PatProb=PatProb,ProbVec=IntPatC3,MeanProbMat=MeanProbMat))
}

# compute descrete survival variable for learning

#Ytrain<-YTCGA
#colnames(Ytrain)<-c("time","status")
#XtrainN<-XTCGA
#Interval=5

#time=Y_base_raw[1,];status=Y_base_raw[2,]
#time=Ytrain$time;status=Ytrain$status
#plot(survfit(Surv(time,status)~NULL), mark.time = T)
#Interval=nrow(Ytrain)
SurvEstimator<-function(Interval=30, time=time, status=status, PDF=F){
  
  Interval=Interval+1
  L=seq(min(time),max(time),length.out = (Interval)-1)
#  L=c(1,L) # add a first interval
  L=c(0,L) # add a first interval
  ProbTable<-data.frame(TimeInt=L)
  
  # cumulative number of events to calculate the nb of subject at risk remaining
  Nevent<-sapply(seq(Interval),function(l){
    subjectRisk<-time>=L[l]
    Nevent<-as.numeric(table(subjectRisk)["TRUE"])
    return(Nevent)  
  })
  
  ProbTable[,"NatRisk"]<-Nevent
  
  # nb of subject censored
  Ncensored<-sapply(seq(Interval),function(l){
    MAX<-ifelse(is.na(L[l+1]),L[l],L[l+1])
    
    Int<-time<MAX &  time>=L[l]
    Ncensored<-ifelse(any(status[Int]==0), as.numeric(table(status[Int]==0)["TRUE"]), 0)
    ProbTable[l,"Ncensored"]<<-Ncensored
    Ndeath<-ifelse(any(status[Int]==1), as.numeric(table(status[Int]==1)["TRUE"]), 0)
    ProbTable[l,"Ndeath"]<<-Ndeath
    
    return(Ncensored)  
  })
  ProbTable[,"Ncensored"]<-Ncensored
#  head(ProbTable)
  
  # probability of survival relative to intervals
  ProbTable[1,"Prob"]<-1
  for(l in seq(Interval)[-1]){
    
    Prob<-ProbTable$Prob[l-1]*( (ProbTable$NatRisk[l]-ProbTable$Ndeath[l]) / ProbTable$NatRisk[l] )
    ProbTable[l,"Prob"]<-Prob
  }
  
  # if more intervals than events, no because not real time scale then with matrix
 # ProbTable<-ProbTable[!duplicated(ProbTable$NatRisk),]
  
  if(PDF){
    plot(ProbTable$TimeInt, ProbTable$Prob, type = 'l', main="Kaplan Meier curve",
         xlab="time", ylab="Probability of survival", ylim=c(0,1))
    plot(ProbTable$Prob, type = 'l', main="Kaplan Meier curve", 
         xlab="Descrete time intervals", ylab="Probability of survival", ylim=c(0,1))
    #lines(ProbTable$Prob,col=2)
  }
  
  # survival function
  St<-ProbTable$Prob
  
  #discrete probability function
  Ft<-St[-length(St)]-St[-1]
  Ft=c(0,Ft)
  
  if(PDF){
    plot(Ft, main="Discrete probability function", xlab="Descrete time intervals")
    #plot(Ft,St[-1])
  }
  
  #discrete hazard rate = conditional failure probability
  Ht=Ft/St
  ProbTable$Ht<-Ht
  
  # St should be = prod 1-ht at interval l
 # Ht<-c(0,Ht)
  if(PDF){
    plot(Ht,main="Ht")
  }
  #plot(St,cumprod(1-Ht)) ##############
  #plot(Ft,Ht*cumprod(1-Ht)) #############
  # ok!
  
  ProbTable<-ProbTable[-1,]
  
  ############ patients' based S(t)
  PatProb<-sapply(time , function(Time){
    1-ProbTable$Prob[max(which(ProbTable$TimeInt<=Time))]
  })
  
  ############# survival estimator
  i=which(time==1)[1]
  
  IntPat<-sapply(seq(length(time)),function(i){  
    Proba<-vector()
    time_interval<-max(which(ProbTable$TimeInt<=time[i]))
    
    #for( l in seq(Interval)){
     for( l in seq(nrow(ProbTable))){
        # uncensored
      Pu<-ifelse(time_interval>=l,1,0)
      
      # censored
      Pc<-ifelse(time_interval>=l,1,St[l])
      Proba<-c(Proba,ifelse(status[i]==0,Pc,Pu))
    }
    return(Proba)  
  })
  #dim(IntPat)
  #image(t(IntPat))
  
  
  ######### in-house approx
  normalized<-function(x) {(x-min(x))/(max(x)-min(x))}
  normTime<-normalized(time)
  
  IntPatC3<-sapply(seq(length(time)),function(x){
    if(status[x]==1){
      Pc = normTime[x]
    } else {
      time_interval=max(which(ProbTable$TimeInt<=time[x]))
#      Pc = normTime[x] + prod(1-Ht[time_interval:length(Ht)]) ###?
      Pc = normTime[x] + Ft[time_interval]
    }
    return(Pc)
  })
  #plot(IntPatC3,time,pch=ifelse(status==0,3,1))

  ############### mean of intervals status
  MeanProbMat<-colMeans(IntPat)
  
  
  
  return(list(ProbTable=ProbTable,ProbMat=IntPat,PatProb=PatProb,ProbVec=IntPatC3,MeanProbMat=MeanProbMat))
}

# compute descrete survival variable for learning

#Ytrain<-YTCGA
#colnames(Ytrain)<-c("time","status")
#XtrainN<-XTCGA
#Interval=5

#time=Y_base_raw[1,];status=Y_base_raw[2,]
#time=Ytrain$time;status=Ytrain$status
#plot(survfit(Surv(time,status)~NULL), mark.time = T)
#Interval=nrow(Ytrain)
SurvEstimator<-function(Interval=30, time=time, status=status, PDF=F){
  
  Interval=Interval+1
  L=seq(min(time),max(time),length.out = (Interval)-1)
#  L=c(1,L) # add a first interval
  L=c(0,L) # add a first interval
  ProbTable<-data.frame(TimeInt=L)
  
  # cumulative number of events to calculate the nb of subject at risk remaining
  Nevent<-sapply(seq(Interval),function(l){
    subjectRisk<-time>=L[l]
    Nevent<-as.numeric(table(subjectRisk)["TRUE"])
    return(Nevent)  
  })
  
  ProbTable[,"NatRisk"]<-Nevent
  
  # nb of subject censored
  Ncensored<-sapply(seq(Interval),function(l){
    MAX<-ifelse(is.na(L[l+1]),L[l],L[l+1])
    
    Int<-time<MAX &  time>=L[l]
    Ncensored<-ifelse(any(status[Int]==0), as.numeric(table(status[Int]==0)["TRUE"]), 0)
    ProbTable[l,"Ncensored"]<<-Ncensored
    Ndeath<-ifelse(any(status[Int]==1), as.numeric(table(status[Int]==1)["TRUE"]), 0)
    ProbTable[l,"Ndeath"]<<-Ndeath
    
    return(Ncensored)  
  })
  ProbTable[,"Ncensored"]<-Ncensored
#  head(ProbTable)
  
  # probability of survival relative to intervals
  ProbTable[1,"Prob"]<-1
  for(l in seq(Interval)[-1]){
    
    Prob<-ProbTable$Prob[l-1]*( (ProbTable$NatRisk[l]-ProbTable$Ndeath[l]) / ProbTable$NatRisk[l] )
    ProbTable[l,"Prob"]<-Prob
  }
  
  # if more intervals than events, no because not real time scale then with matrix
 # ProbTable<-ProbTable[!duplicated(ProbTable$NatRisk),]
  
  if(PDF){
    plot(ProbTable$TimeInt, ProbTable$Prob, type = 'l', main="Kaplan Meier curve",
         xlab="time", ylab="Probability of survival", ylim=c(0,1))
    plot(ProbTable$Prob, type = 'l', main="Kaplan Meier curve", 
         xlab="Descrete time intervals", ylab="Probability of survival", ylim=c(0,1))
    #lines(ProbTable$Prob,col=2)
  }
  
  # survival function
  St<-ProbTable$Prob
  
  #discrete probability function
  Ft<-St[-length(St)]-St[-1]
  Ft=c(0,Ft)
  
  if(PDF){
    plot(Ft, main="Discrete probability function", xlab="Descrete time intervals")
    #plot(Ft,St[-1])
  }
  
  #discrete hazard rate = conditional failure probability
  Ht=Ft/St
  ProbTable$Ht<-Ht
  
  # St should be = prod 1-ht at interval l
 # Ht<-c(0,Ht)
  if(PDF){
    plot(Ht,main="Ht")
  }
  #plot(St,cumprod(1-Ht)) ##############
  #plot(Ft,Ht*cumprod(1-Ht)) #############
  # ok!
  
  ProbTable<-ProbTable[-1,]
  
  ############ patients' based S(t)
  PatProb<-sapply(time , function(Time){
    1-ProbTable$Prob[max(which(ProbTable$TimeInt<=Time))]
  })
  
  ############# survival estimator
  i=which(time==1)[1]
  
  IntPat<-sapply(seq(length(time)),function(i){  
    Proba<-vector()
    time_interval<-max(which(ProbTable$TimeInt<=time[i]))
    
    #for( l in seq(Interval)){
     for( l in seq(nrow(ProbTable))){
        # uncensored
      Pu<-ifelse(time_interval>=l,1,0)
      
      # censored
      Pc<-ifelse(time_interval>=l,1,St[l])
      Proba<-c(Proba,ifelse(status[i]==0,Pc,Pu))
    }
    return(Proba)  
  })
  #dim(IntPat)
  #image(t(IntPat))
  
  
  ######### in-house approx
  normalized<-function(x) {(x-min(x))/(max(x)-min(x))}
  normTime<-normalized(time)
  
  IntPatC3<-sapply(seq(length(time)),function(x){
    if(status[x]==1){
      Pc = normTime[x]
    } else {
      time_interval=max(which(ProbTable$TimeInt<=time[x]))
#      Pc = normTime[x] + prod(1-Ht[time_interval:length(Ht)]) ###?
      Pc = normTime[x] + Ft[time_interval]
    }
    return(Pc)
  })
  #plot(IntPatC3,time,pch=ifelse(status==0,3,1))

  ############### mean of intervals status
  MeanProbMat<-colMeans(IntPat)
  
  
  
  return(list(ProbTable=ProbTable,ProbMat=IntPat,PatProb=PatProb,ProbVec=IntPatC3,MeanProbMat=MeanProbMat))
}
