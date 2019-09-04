#' This function builds a network and clean it from genes' query and interaction database
#' @param GENESman character vector. A vector of gene(s) selected to build the *AMoNet* object.
#'
#' @param treatmt character vector. A vector of gene(s) targeted by treatment(s). Used in to build the *AMoNet* object.
#' @param OMNI data frame. Interaction database with columns ordered as: source_hgnc, interaction_directed_signed, target_hgnc corresponding to source species, interation type (ACTIVATE or INHIBIT), target species, respectively.
#' @param VeryRestricted boolean. A restriction during the search in the PPI database to keep only the new nodes that have already a partner in the network at the previous time step.
#' @param nblayers integer. Number of layers from the gene queries to explore within the InteractionBase (OMNI by default), i.e number of search steps. Should be >=1.
#' @param MinConnect integer. During the filtering phase, what is the minimal value of connections to keep a node. If MinConnect=1, you'll discard all nodes that are only "passengers".
#' @param WRITE boolean. Should the network be saved as csv?
#' @param no_cores numeric. Set the number of cores to parallelize on. Package parallel required.
#' @param FamilyGene data frame. Restrict the building of the nets to genes related to the phenotypes selected with MECA. Default are \code{GenesSelectImmuno} and/or \code{GenesSelecHall}
NetBuilding <- function(GENESman=GENESman,treatmt=treatmt, OMNI=OMNI, nblayers=2,
                        FamilyGene=NULL, WRITE=F, MinConnect=1, VeryRestricted=F, no_cores=4){

   ################
  # keep all genes at each layers and then filter those with feedbacks if VeryRestricted = T

  i=0
  Nspecies<-vector()

  if(!is.null(OMNI)){

    for(y in seq(nblayers) ){
      i=i+1
      if(y==1){
        NET<-OMNI[OMNI$source_hgnc%in%union(GENESman,treatmt)|OMNI$target_hgnc%in%union(GENESman,gsub("_TTT","",treatmt)),]
        NETall<-NET
      } else {

        Species<-union( NETall$source_hgnc, NETall$target_hgnc)

      NETall<-OMNI[OMNI$source_hgnc%in%Species|OMNI$target_hgnc%in%Species,]

        if(!is.null(FamilyGene)){
          NETall<-NETall[NETall$source_hgnc%in%FamilyGene|NETall$target_hgnc%in%FamilyGene,] # keep FamilyGene only
        }

      if(VeryRestricted){
        NETall<-NETall[NETall$target_hgnc%in%Species,]
        NETall<-rbind(NETall,NET) # restaure initial genes
      }

      }

      NETall<-NETall[!duplicated(paste(NETall$source_hgnc,NETall$target_hgnc)),] # no dups

      Nspecies[i]<-length(unique(NETall$source_hgnc))[1]
      print(paste( "Round",i,"=",length(unique(NETall$source_hgnc))[1], "species"))
      if(length(unique(NETall$source_hgnc))==0){
        print("not enought species")
        return(NETall=NULL)
      }
    }

  } else {
    print("Necessary to provide a base of directed protein interactions")
  }

  #write.csv2(NETall,"C:/Users/L_VERLINGUE/Desktop/ModelK/Roma/networkexpansion3_PDL1_PD1_CTLA4_flash.csv")
  #write.csv2(NETall,"C:/Users/L_VERLINGUE/Desktop/ModelK/Roma/networkexpansion3_PDL1_PD1_CTLA4_flash_restricted.csv")

  #no_core<-detectCores()
  cl<-parallel::makeCluster(no_cores)

  Nspecies<-c(1,2)
  while( length(unique(tail(round(Nspecies),n = 2)))!=1 )  {
    i=i+1

    # then do a restriciton to fb edges
    TS<-NETall$target_hgnc%in%NETall$source_hgnc # node without output
    ST<-NETall$source_hgnc%in%NETall$target_hgnc # node without input
    table(TS&ST)


   # while(length(unique(TS&ST))!=1){
    while(!all(TS&ST)){
      TS<-NETall$target_hgnc%in%NETall$source_hgnc # node without output
      ST<-NETall$source_hgnc%in%NETall$target_hgnc # node without input
      #data.frame(NETall[,c(1:3)],TS,ST)
      NETall<-NETall[TS&ST,]

      # dim(NETall)
      #if(dim(NETall)[1]<10) break
    }

    #table(NETall$source_hgnc%in%NETall$target_hgnc)
    #dim(NETall)

    #GENESman%in%union(NETall$source_hgnc,NETall$target_hgnc)

    print(paste("post restriction to FB N=",length(unique(NETall$source_hgnc)) ))

    # are nodes only "passengers": discard them
    #TAB<-sort(table(c(as.character(NETall$source_hgnc),as.character(NETall$target_hgnc))))
    TAB<-table(NETall$source_hgnc)
    #tail(sort(TAB))
    #NETall[NETall$source_hgnc%in%"FASL"|NETall$target_hgnc%in%"FASL",1:3]
    #TAB[names(TAB)%in%unique(NET$source_hgnc)]
    PASSENGERS<-names(TAB[TAB<=MinConnect])
    # dont take the queried genes
    PASSENGERS<-setdiff(PASSENGERS,GENESman)

    #length(PASSENGERS)
    #MinConnect=4
    # function to store the passenger inputs and outputs
    #x<-PASSENGERS[3]
    #PAS<-PASSENGERS[3]

    # check the type of relation
    #PAS<-PASSENGERS[1]
    #PAS="WWOX"
    print("Imputing connections for reduction")

  #  if(length(PASSENGERS)>10){

      EXPORT<-c("PASSENGERS","NETall")
      parallel::clusterExport(cl,EXPORT, envir = environment())

      PASSV<-parallel::parLapply(cl,PASSENGERS,function(PAS){
        SOU<-NETall[NETall$target_hgnc%in%PAS,"source_hgnc"]
        ISOU<-NETall[NETall$target_hgnc%in%PAS,"interaction_directed_signed"]
        TAR<-NETall[NETall$source_hgnc%in%PAS,"target_hgnc"]
        ITAR<-NETall[NETall$source_hgnc%in%PAS,"interaction_directed_signed"]

        if(length(TAR)==1&length(SOU)==1){
          V<-c(ISOU,ITAR)
          names(V)<-c(SOU,TAR)
          Xlist<-list(V)
          names(Xlist)<-PAS
          #          PASSV[[length(PASSV)+1]]<-V
          #          names(PASSV)[length(PASSV)]<-PAS
        } else {
          names(ISOU)<-SOU
          names(ITAR)<-TAR

          if(length(SOU)>1){
            TAR<-rep(TAR,length(SOU))
            ITAR<-rep(ITAR,length(ISOU))
            #TAR;ITAR
          }

          V<-rbind(SOU,TAR)
          #  print( c(PAS,length(SOU),length(TAR),dim(V)))
          # Xlist<-t(V)
          # c(ISOU[V["SOU",NN]],ITAR[V["TAR",NN]])

          Xlist<-lapply(seq(ncol(V)),function(NN){
            c(ISOU[V["SOU",NN]],ITAR[V["TAR",NN]])
          })
          names(Xlist)<-rep(PAS,length(Xlist))
        }


        #x<-Xlist[[3]]
        RESlist<- lapply(Xlist,function(x){

          # auto loop autorized!
#          if( length(unique(paste(names(x),x,sep = "_")))==1 ){
#            RES<-NULL
#          }else
          if(length(unique(x))==1){ # if only activ or 2 inhib
            RES<-c(names(x)[1],"ACTIVATE",names(x)[2])
            #          NETall[,"source_hgnc"]<<-names(x)[1]
            #          NETall[nrow(NETall),"target_hgnc"]<<-names(x)[2]
            #          NETall[nrow(NETall),"interaction_directed_signed"]<<-"ACTIVATE"
          } else {
            RES<-c(names(x)[1],"INHIBIT",names(x)[2])
            #          NETall[nrow(NETall)+1,"source_hgnc"]<<-names(x)[1]
            #          NETall[nrow(NETall),"target_hgnc"]<<-names(x)[2]
            #          NETall[nrow(NETall),"interaction_directed_signed"]<<-"INHIBIT"
          }

        })

        RES<-do.call('rbind',RESlist)
        return(RES)
      })

      PASSV<-do.call('rbind',PASSV)

      if(!is.null(PASSV)){
       # PASSV<-as.data.frame(PASSV[,1:3],stringsAsFactors = F,row.names = NULL)
        colnames(PASSV)<-colnames(NETall)[1:3]
        #class(PASSV[,1])
        #table(PASSV[,1]%in%PASSENGERS) # possible
        #length(union(PASSV[,1],PASSV[,3]))
        # remove passengers from the network table
        # and add their parent and children nodes as a direct interaction
        # taking into account the type of interaction
        #length((PASSENGERS))
        #table(!(NETall$source_hgnc%in%PASSENGERS|NETall$target_hgnc%in%PASSV))
        NETall<-NETall[!(NETall$source_hgnc%in%PASSENGERS|NETall$target_hgnc%in%PASSENGERS),1:3]
        #length(union(NETall$source_hgnc,NETall$target_hgnc))
        NETall<-rbind(NETall,PASSV)
      }

    #   } else {

      if(FALSE){
      PASSV<-list()

      for(PAS in PASSENGERS){ # slow loop
        SOU<-NETall[NETall$target_hgnc%in%PAS,"source_hgnc"]
        ISOU<-NETall[NETall$target_hgnc%in%PAS,"interaction_directed_signed"]
        TAR<-NETall[NETall$source_hgnc%in%PAS,"target_hgnc"]
        ITAR<-NETall[NETall$source_hgnc%in%PAS,"interaction_directed_signed"]

        if(length(TAR)==1&length(SOU)==1){
          V<-c(ISOU,ITAR)
          names(V)<-c(SOU,TAR)
          PASSV[[length(PASSV)+1]]<-V
          names(PASSV)[length(PASSV)]<-PAS
        } else {
          names(ISOU)<-SOU
          names(ITAR)<-TAR

          if(length(SOU)>1){
            TAR<-rep(TAR,length(SOU))
            ITAR<-rep(ITAR,length(ISOU))
            #TAR;ITAR
          }

          V<-rbind(SOU,TAR)
          #  print( c(PAS,length(SOU),length(TAR),dim(V)))


          for(NN in seq(ncol(V))){
            PASSV[[length(PASSV)+1]]<-c(ISOU[V["SOU",NN]],ITAR[V["TAR",NN]])
            names(PASSV)[length(PASSV)]<-PAS
          }

          if(FALSE){
            SOU<-rep(SOU,length(TAR))
            ISOU<-rep(ISOU,length(TAR))

            ITAR<-unlist(lapply(ITAR,function(x)rep(x,length(TAR))))
            TAR<-unlist(lapply(TAR,function(x)rep(x,length(TAR))))

            for(NN in seq(length(SOU))){
              V<-c(ISOU[NN],ITAR[NN])
              names(V)<-c(SOU[NN],TAR[NN])
              PASSV[[length(PASSV)+1]]<-V
              names(PASSV)[length(PASSV)]<-PAS
            }
          }
        }
      } # en for


      #table(unlist(lapply(PASSV, length))==2)
      #KEEP<-as.numeric(which((lapply(PASSV,length)==MinConnect))) # should be everyone, just a check
      #PASSV<-PASSV[KEEP]
      #rm(KEEP)

      # remove passengers from the network table
      # and add their parent and children nodes as a direct interaction
      # taking into account the type of interaction
      #table((NETall$source_hgnc%in%names(PASSV)|NETall$target_hgnc%in%names(PASSV)))
      NETall<-NETall[!(NETall$source_hgnc%in%names(PASSV)|NETall$target_hgnc%in%names(PASSV)),]
      dim(NETall)


      #x<-PASSV[[1]]
      if(length(PASSV)>0){
        for(y in seq(length(PASSV))){
          x<-PASSV[[y]]
          if( length(unique(paste(names(x),x,sep = "_")))==1 ){
            next
          }else if(length(unique(x))==1){ # if only activ or 2 inhib
            NETall[nrow(NETall)+1,"source_hgnc"]<-names(x)[1]
            NETall[nrow(NETall),"target_hgnc"]<-names(x)[2]
            NETall[nrow(NETall),"interaction_directed_signed"]<-"ACTIVATE"
          } else {
            NETall[nrow(NETall)+1,"source_hgnc"]<-names(x)[1]
            NETall[nrow(NETall),"target_hgnc"]<-names(x)[2]
            NETall[nrow(NETall),"interaction_directed_signed"]<-"INHIBIT"
          }
        }
      }

      }

      #  print(paste("post passenger filtering N sources =",length(unique(NETall$source_hgnc)) ))
  #  } # end if

    # discard duplicates
    NETall<-NETall[!duplicated(paste(NETall$source_hgnc,NETall$target_hgnc)),] # no dups

    # discard auto regulated node only
    AUTO<-table(NETall$source_hgnc[NETall$source_hgnc%in%NETall$source_hgnc[NETall$source_hgnc==NETall$target_hgnc]])
    if(any(AUTO==1)){
      NETall<-NETall[!NETall$source_hgnc%in%names(AUTO[AUTO==1]),]
    }

    Nspecies[i]<-length(unique(NETall$source_hgnc))[1]
    print(paste( "Round",i,"=",length(unique(NETall$source_hgnc))[1], "species"))

    if(length(unique(NETall$source_hgnc))==0){
      print("not enought species")
      NETall<-NULL
      break
      }
  }

  parallel::stopCluster(cl)

  if(WRITE){
    colnames(NETall)
    write.csv2(NETall[,c(1:3,9,10)],paste(getwd(),"/model/NETall.csv",sep = ""))
  }
  return(NETall)
}
