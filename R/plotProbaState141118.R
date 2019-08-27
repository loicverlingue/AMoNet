#
plotProbaState<-function(Species=Species,TotAttractors=TotAttractors, Level=c("A","Z"), 
                         Plot=T, PDF=F, Call=NULL, NameProj=NameProj, Legend=F){
#  setwd(gsub("/tmp","",getwd()))
#  setwd(paste(getwd(),"/tmp/",sep = ""))
#  getwd()
  
 # dim(TotAttractors ) 
  
  # put the call of the simulation to put it into the title of the graph
  #Details=paste(names(unlist(TotAttractors$Call)),"=",unlist(TotAttractors$Call),collapse = "; ")
  #Logic=TotAttractors$Call$Logic
  if(!is.null(Call)){
    Details=paste(names(unlist(Call)),"=",unlist(Call),collapse = "; ")
    Logic=Call$Logic
  } else {
    Details=NULL
    Logic=NULL
  }
  
#  if(!all(Species%in%colnames(TotAttractors$TotAttractors))){
#    stop("Species selected not in the network")
#  }
  
  #lapply(TotAttractors$LAttractors,dim)
  # function to call
  
 # dim(TotAttractors)
  SpeciesActivity<-sapply(Species,function(SP){
    if(dim(TotAttractors)[1]>1){
    colMeans(cbind(TotAttractors[,"iStates",SP,1],TotAttractors[,Level,SP,]))
    } else {
      t(c(TotAttractors[,"iStates",SP,1],TotAttractors[,Level,SP,]))
    }
    
  })
  #    colMeans(TotAttractors[,Level,SP,])    
  
  if(FALSE){
  #LA=TotAttractorsMult$LAttractors
    ProbaState<-function(LA=TotAttractors$LAttractors,Ospecies=Ospecies){
      MAX<-max(unlist(lapply(LA,nrow)))
    #Ospecies="FOXP3"
      
      LAttractors<-lapply(LA,function(x){
        x<-x[,Ospecies]
        x<-c(x,rep(x[length(x)], MAX-length(x))) # repeat the last state
      })
      LAttractors<-do.call("rbind",LAttractors)
     # table(LAttractors)
      MEAN<-colMeans(LAttractors)
      return(MEAN)
  }
    if("A"%in%Level){
      SpeciesActivity<-sapply(Species,function(x){
        ProbaState(LA = TotAttractors$LAttractors,Ospecies = x)
      })
    } else if("Z"%in%Level){
      SpeciesActivity<-sapply(Species,function(x){
        ProbaState(LA = TotAttractors$LAttractorsZ,Ospecies = x)
      })
    }
  }
  
  
    if(Plot){
      
      if(nrow(SpeciesActivity)>4){
        z2<-apply(SpeciesActivity,2,function(h)smooth.spline(h,df = length(h))$y)
      } else{
        z2<-SpeciesActivity
      }
      
      #matplot(z2,type='l')
      
      if(PDF){
        NameProj=paste(NameProj,"Species_activities",sep = "_")
        
        Latt<-length(list.files(paste(getwd(),"/tmp/",sep = ""),pattern =  NameProj))
        if(Latt>0){
          Latt<-max(as.numeric(gsub(".pdf","", gsub(NameProj,"", list.files(paste(getwd(),"/tmp/",sep = ""),pattern =  NameProj)))))
        } 
        
        NAME<-paste("/tmp/",NameProj, Latt+1,sep = "") #ifelse(Optimized,"Opt","Init"),
        pdf(paste(getwd(),"/",NAME,".pdf",sep = ""))
      }
      
      if(FALSE){
        Latt<-length(list.files(paste(getwd(),"/tmp/",sep = ""),pattern =  "Species_activities"))
        if(Latt>0){
          Latt<-max(as.numeric(gsub(".pdf","", gsub("Species_activities","", list.files(paste(getwd(),"/tmp/",sep = ""),pattern =  "Species_activities")))))
          NAME<-paste("Species_activities",Latt+1,sep = "")
          pdf(paste(getwd(),"/tmp/",NAME,".pdf",sep = ""))
        }
      }
  
      if(Legend){
        par(mfrow=c(1,1))
        par(mar=c(5,5,4,10)+0.1) #No margin on the right side
      } else {
       # plot.new()
        par(mfrow=c(1,1))
        par(mar=c(5,5,4,2)+0.1)
      }
      matplot(z2,type='l',axes=F,ylab = paste("Mean species activities \n on", dim(TotAttractors)[1],"simulations" ) ,xlab = "Steps",
              main=paste("Simulations with",Details,"\n for level:",Level) )
      #  arrows(x0 = 0,y0 = -0.05,x1 = nrow(z2),y1 = -0.05,length = 0.2,xpd=T)
      axis(2,las=2)
      axis(1)
      
      if(Legend){
        par(new=T)
        par(mar=c(5,8,4,2)) #No margin on the left side
        plot(c(0,1),type="n", axes=F, xlab="", ylab="")
        legend("right", colnames(z2),col=seq_len(ncol(z2)),lty=seq_len(ncol(z2)),cex=0.4)
        
      }
      
      if(PDF){ 
        dev.off()
      }
    }
  return(SpeciesActivity)
  
  if(FALSE){
    ###### problem of zeros
    LA=TotAttractors$LAttractors
      MAX<-max(unlist(lapply(LA,nrow)))
     # Ospecies="CDKN1A_Input"
      LAttractors<-lapply(LA,function(x){
        x<-x[,Ospecies]
        x<-c(x,rep(x[length(x)], MAX-length(x))) # repeat the last state
      })
      LAttractors<-do.call("rbind",LAttractors)
      MEAN<-colMeans(LAttractors)
      plot(MEAN, type='l')
      dim(LAttractors)
      matplot(t(LAttractors),type='l', main=paste(Ospecies,"activity"))
      
      table(LAttractors[,ncol(LAttractors)])/nrow(LAttractors)
      (sort(table(LAttractors),decreasing = T)/length(as.numeric(LAttractors)))
      
      KM<-kmeans(LAttractors[,ncol(LAttractors)],centers = 2)
      KM$centers
      
      DENS<-density(t(LAttractors))
      plot(DENS)
      DENS$x[DENS$y==max(DENS$y)] 
      DENS$x[DENS$y==min(DENS$y)] 
      hm<-hclust(dist(LAttractors[,ncol(LAttractors)]))
      
    }
    
}
