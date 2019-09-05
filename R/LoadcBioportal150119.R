

#### load TCGA data from cbioportal
# needs several docs in data: clinicTCGA.csv (TCGA studies with OS) ; GTExTCGAcode.txt
# check that pooling studies desnt duplicate the data!

# parameters
#Organ="luad"
#setwd("C:/Users/L_VERLINGUE/Desktop/ModelK/Rpack/ArtMolNet/")

LoadcBioportal<-function(Genes=c("TP53","KRAS"), ClinicNeeded=T,
                       RNANeeded=F, MutNeeded=T, FunctionalAnnot=T,
                       Organ="lung", NormalizeRNA=T,
                       PDF=F,Tests=F){

  ### define function to normalise expression data
  normalized = function(x) {(x-min(x))/(max(x)-min(x))}

  #### load cbioportal
  #library(cgdsr)
  mycgds = cgdsr::CGDS("http://www.cbioportal.org/")
  if(Tests){
    cgdsr::test(mycgds)
  }

  # these are the tumor types in TCGA nomenclature, taken from NCI
  OfficialNamesTCGA=tolower(unique(scan(text="GBM	GBM	GBM	GBM	GBM	GBM	GBM	GBM	OV	GBM	OV	OV	OV	OV	OV	GBM	OV	OV	OV	GBM	OV	OV	LUSC	OV	LAML	GBM	OV	COAD	COAD	COAD	LUSC	KIRC	COAD	LUAD	COAD	COAD	LUAD	GBM	LUSC	OV	COAD	READ	READ	OV	COAD	READ	BRCA	STAD	UCEC	KIRC	KIRP	LUAD	LUSC	HNSC	LIHC	BRCA	STAD	LUAD	UCEC	LUSC	BRCA	GBM	KIRC	KIRC	KIRC	COAD	READ	KIRC	KIRC	KIRC	KIRP	BRCA	UCEC	BRCA	UCEC	COAD	LUSC	LGG	GBM	BRCA	UCEC	KIRC	HNSC	LUAD	BRCA	BLCA	THCA	CESC	COAD	KIRC	PRAD	UCEC	BRCA	UCEC	STAD	BRCA	UCEC	COAD	READ	LIHC	LUSC	READ	BRCA	UCEC	KIRC	PAAD	HNSC	PRAD	BRCA	UCEC	GBM	LGG	BLCA	CESC	THCA	COAD	BRCA	UCEC	LUAD	BRCA	UCEC	READ	COAD	BRCA	UCEC	THCA	CESC	BLCA	STAD	GBM	LIHC	COAD	READ	PAAD	DLBC	BRCA	UCEC	COAD	READ	LUSC	BRCA	UCEC	LUAD	HNSC	LGG	BRCA	CESC	THCA	BLCA	HNSC	STAD	LIHC	COAD	BRCA	UCEC	COAD	READ	LUSC	LUAD	PRAD	KIRP	LGG	HNSC	STAD	LUAD	BRCA	UCEC	CESC	BLCA	THCA	COAD	LIHC	GBM	BLCA	THCA	BRCA	UCEC	CESC	SKCM	LUSC	THCA	LUAD	PRAD	BRCA	UCEC	CESC	HNSC	LGG	HNSC	THCA	BLCA	LUSC	KIRP	PAAD	LUAD	THCA	SKCM	BLCA	CESC	UCEC	BRCA	LIHC	LUAD	THCA	SKCM	BLCA	LUSC	KIRP	PAAD	STAD	DLBC	LUAD	LUSC	HNSC	BRCA	CESC	SARC	LGG	STAD	PRAD	LUAD	BLCA	THCA	LUSC	KICH	BRCA	UCEC	SARC	THCA	LIHC	LUAD	LUSC	BRCA	BLCA	CESC	THCA	LUAD	BRCA	SKCM	HNSC	STAD	LUSC	PRAD	LGG	KIRP	PAAD	DLBC	BLCA	THCA	KICH	BLCA	THCA	ESCA	BRCA	CESC	STAD	LUAD	LUSC	HNSC	KIRP	SKCM	SARC	LUAD	HNSC	KIRP	LUSC	PRAD	STAD	THCA	BRCA	ESCA	HNSC	KIRC	LIHC	LUSC	SKCM	LUAD	BLCA	CESC	KIRP	LGG	LUSC	PAAD	PRAD	SARC	LIHC	BRCA	UCEC	SARC	SKCM	LGG	LUSC	BLCA	LGG	BRCA	CESC	ESCA	KIRP	COAD	HNSC	UCS	LIHC	ACC	BRCA	LGG	SARC	PRAD	BLCA	SARC	HNSC	PRAD	ACC	LIHC	PRAD	SKCM	BLCA	HNSC	LGG	PRAD	STAD	BRCA	ESCA	UCEC	KIRP	DLBC	LIHC	BLCA	LGG	PAAD	PRAD	SKCM	STAD	BRCA	CESC	SARC	BLCA	BRCA	ESCA	KIRC	LIHC	SARC	STAD	KIRP	LIHC	GBM	LGG	PRAD	PAAD	CESC	KICH	LGG	MESO	UCEC	HNSC	BLCA	PRAD	SKCM	STAD	BRCA	CESC	ESCA	SARC	KIRP	LIHC	PCPG	BLCA	LGG	PAAD	PRAD	STAD	BRCA	CESC	ESCA	SARC	KIRP	LIHC	GBM	BRCA	CESC	UCEC	ESCA	SARC	LIHC	COAD	KIRP	KIRC	SKCM	PRAD	PAAD	BLCA	PAAD	SKCM	CESC	ESCA	SARC	DLBC	KIRP	LIHC	BLCA	GBM	HNSC	KIRC	LIHC	LUAD	DLBC	SKCM	OV	PAAD	READ	SARC	STAD	THCA	LUSC	SKCM	UVM	MESO	ESCA	CESC	UCEC	COAD	LUAD	PAAD	LIHC	SARC	STAD	CHOL	PRAD	TGCT	THYM	KIRP	BLCA	PAAD	CHOL",what = "")))

  if(!is.null(Organ)){
    # table to translate TCGA tumor types to organ names (inspired from GTEx)
    # either you have the file or directly build the table

    TCGAcode=scan(text= "sarc	NA	skcm	sarc	NA	thca	lu	NA	esca	NA	NA	esca	NA	NA	NA	brca	coad	NA	stad	tgct	paad	esca	coad	acc	NA	gbm	lihc	dlbc	gbm	gbm	gbm	gbm	prad	gbm	NA	NA	gbm	ov	NA	NA	NA	NA	gbm	ucec	NA	NA	NA	NA	kirc	blca	cesc	NA	cesc",what = "")

    Tissue=scan(text = "Muscle	Blood	Skin	Adipose	Artery	Thyroid	Lung	Nerve	Esophagus	fibroblasts	Skin	Esophagus	Adipose	Artery	Heart	Breast	Colon	Heart	Stomach	Testis	Pancreas	Esophagus	Colon	Adrenal	Artery	Brain	Liver	lymphocytes	Brain	Brain	Brain	Brain	Prostate	Brain	Spleen	Pituitary	Brain	Ovary	Hypothalamus	Vagina	Hippocampus	Ileum	Brain	Uterus	Brain	Spinal_cord	Brain	Salivary_Gland	Kidney	Bladder	Ectocervix	Fallopian_Tube	Endocervix", what = "")

    GTExTCGAcode<-t(data.frame(Tissue=Tissue,TCGAcode=TCGAcode))
    colnames(GTExTCGAcode)<-Tissue
  } else {
    print("Organ = NULL --> load all data")
    ANS <- readline("do you want to continue - it maybe long and heavy [Y/N]?")
    if("y"%in%ANS){
      NUMstudy<-unique(cgdsr::getCancerStudies(mycgds)$cancer_study_id)
    } else {
      print(paste("You should choose organ within:",paste(OfficialNamesTCGA,collapse = ", ") ))
      stop()
    }
  }
  #grep("impact",tolower(NUMstudy),value = T)

  # OS_TCGA = studies with RNAseq and OS
  OS_TCGA<-scan(text = "laml_tcga_pub	laml_tcga	acyc_mskcc_2013	acbc_mskcc_2015	acc_tcga	ampca_bcm_2016	blca_mskcc_solit_2014	blca_mskcc_solit_2012	blca_tcga_pub_2017	blca_plasmacytoid_mskcc_2016	blca_tcga_pub	blca_tcga	lgg_tcga	brca_metabric	brca_tcga_pub2015	brca_tcga_pub	brca_tcga	cesc_tcga	chol_nus_2012	chol_tcga	coadread_tcga_pub	coadread_tcga	coadread_mskcc	cscc_hgsc_bcm_2014	cscc_dfarber_2015	esca_tcga	escc_icgc	es_iocurie_2014	gbc_shanghai_2014	egc_tmucih_2015	prad_cpcg_2017	gct_msk_2016	mixed_allen_2018	gbm_tcga_pub2013	gbm_tcga_pub	gbm_tcga	hnsc_tcga_pub	hnsc_tcga	kich_tcga_pub	kich_tcga	kirc_tcga_pub	kirc_tcga	kirp_tcga	lihc_amc_prv	lihc_tcga	lgg_ucsf_2014	luad_tcga_pub	luad_tcga	lusc_tcga	dlbc_tcga	plmeso_nyu_2015	mbl_icgc	mbl_pcgp	mbl_sickkids_2016	skcm_broad_dfarber	lgggbm_tcga_pub	meso_tcga	egc_msk_2017	prad_mich	mixed_pipseq_2017	odg_msk_2017	npc_nusingapore	nbl_amc_2012	nbl_broad_2013	skcm_vanderbilt_mskcc_2015	nhl_bcgsc_2011	hnsc_mdanderson_2013	ov_tcga_pub	ov_tcga	mel_tsam_liang_2017	paad_tcga	thca_tcga_pub	all_phase2_target_2018_pub	es_dfarber_broad_2014	nbl_target_2018_pub	rt_target_2018_pub	wt_target_2018_pub	pcpg_tcga	thyroid_mskcc_2016	pcnsl_mayo_2015	prad_tcga	prad_mskcc_2014	hnc_mskcc_2016	sarc_tcga	skcm_tcga	sclc_ucologne_2015	stad_tcga_pub	stad_tcga	stad_uhongkong	stes_tcga_pub	urcc_mskcc_2016	crc_msk_2018	utuc_mskcc_2013	tgct_tcga	thym_tcga	thca_tcga	ucs_tcga	ucec_tcga_pub	ucec_tcga	uvm_tcga	panet_arcnet_2017	skcm_ucla_2016	past_dkfz_heidelberg_2013", what = "")
  NUMstudybase<-unique(cgdsr::getCancerStudies(mycgds)$cancer_study_id)

  if(ClinicNeeded){


    NUMstudy<-NULL
    if(!is.null(Organ)){
      CODE<-unique((GTExTCGAcode[2,grep(tolower(Organ),tolower(colnames(GTExTCGAcode)))]))

      if(length(CODE)>0){
        NUMstudy<-sapply(CODE[!is.na(CODE)], function(x){grep(x, OS_TCGA, value = T)})
        NUMstudy<-unique(unlist(NUMstudy))
      } else {
        NUMstudy<-grep(tolower(Organ), tolower(OS_TCGA),value = T)
      }
    } else if(is.null(Organ)){
      NUMstudy<-unique(OS_TCGA)
    }
  } else if(!ClinicNeeded){
      if(!is.null(Organ)){
       CODE<-(GTExTCGAcode[2,grep(tolower(Organ),tolower(names(GTExTCGAcode)))])
       if(length(CODE)>0){
        NUMstudy<-sapply(CODE[!is.na(CODE)], function(x){grep(x, NUMstudybase, value = T)})
        NUMstudy<-unique(unlist(NUMstudy))
      } else {
        NUMstudy<-grep(tolower(Organ), NUMstudybase,value = T)
      }

     } else if(is.null(Organ)){
      NUMstudy<-unique(NUMstudybase)
    }
  }

  if(is.null(NUMstudy)) {
    print("Check 'Organ' name is in TCGA diseases' format, eg: lung adenocarcinoma = luad")
    NUMstudy<-grep(tolower(Organ), tolower(NUMstudybase),value = T)
  }
  if(is.null(NUMstudy)) {
    print("Didn't find your TCGA study, revise the 'organ' parameter")
    stop()
  } else {
    # remove duplicated studies
    NUMstudy<-NUMstudy[!duplicated(gsub("_pub","",NUMstudy))]
  }

  print(paste("Pooling :", paste(NUMstudy, collapse = " & ")))

  MUT<-data.frame()
  EXP<-data.frame()
  CLINIC<-data.frame()
  GenesNAS<-vector()
  Study<-data.frame(row.names = NUMstudy)
  #mycancerstudy<-NUMstudy[1]

  for(mycancerstudy in NUMstudy){

    #############
    # get study
    # mycancerstudy = getCancerStudies(mycgds)[NUMstudy,1]
    mycaselist = grep("complete|all", tolower(cgdsr::getCaseLists(mycgds,mycancerstudy)[,1]),value = T)
    if(length(mycaselist)>1){
      mycaselist = grep("complete", tolower(cgdsr::getCaseLists(mycgds,mycancerstudy)[,1]),value = T)
    }
    #mycaselist="luad_tcga_pub_all"

    ################
    # clinic
    if(ClinicNeeded){
      myclinicaldata = cgdsr::getClinicalData(mycgds,mycaselist)

      #colnames(myclinicaldata)
      #ClinVar<-c("OS_MONTHS", "OS_STATUS", "PFS_MONTHS","PFS_STATUS","DFS_MONTHS","DFS_STATUS","CANCER_TYPE","AJCC_PATHOLOGIC_TUMOR_STAGE")
      ClinVar<-c("OS_MONTHS", "OS_STATUS")
      #myclinicaldata[,"CANCER_TYPE"]

      ClinVar<-ClinVar[ClinVar%in%colnames(myclinicaldata)]

    #  ClinVar<-c(ClinVar[which(ClinVar%in%colnames(myclinicaldata))],grep("TUMOR_STAGE",colnames(myclinicaldata),value = T)[1])

      if(!length(ClinVar)%in%2){
        print(paste(mycancerstudy, "does't have OS and Status recorded"))
        next
      }

      Study[mycancerstudy,1]<-nrow(myclinicaldata)

      CLINIC1<-myclinicaldata[,ClinVar]
      CLINIC1[,"study"]<-gsub("_.*","", mycancerstudy)
      CLINIC<-rbind(CLINIC,CLINIC1)
    } else {
      CLINIC=NULL
    }

    ###############
    # for expression
    if(RNANeeded){
      mygeneticprofile = grep("mrna", cgdsr::getGeneticProfiles(mycgds,mycancerstudy)[,1],value = T)
      if(length(setdiff(grep("v2",mygeneticprofile),grep("Zscores",mygeneticprofile)))==1){
        mygeneticprofile =mygeneticprofile[setdiff(grep("v2",mygeneticprofile),grep("Zscores",mygeneticprofile))]

        EXP1 <- cgdsr::getProfileData(mycgds, genes = Genes ,mygeneticprofile,mycaselist)

        # tag NAs
        NAS<-apply(EXP1,2,function(Exp)any(is.na(Exp)))
        GenesNAS<-union(GenesNAS,Genes[NAS])

        if(NormalizeRNA){
          EXP1<-normalized(log10(as.matrix(EXP1[,!NAS])+1))
        } else {
          EXP1<-EXP1[,!NAS]
        }

        Study[mycancerstudy,1]<-nrow(EXP1)

      } else {
        print(paste("No RNAseq_V2 for", mycancerstudy))
        next()
      }

      if(nrow(EXP)>1){
        GeneID<-intersect(colnames(EXP),colnames(EXP1))
        EXP<-rbind(EXP[,GeneID],EXP1[,GeneID])
      } else{
        EXP<-rbind(EXP,EXP1)
      }
      #dim(EXP)
    } else {
      EXP<- NULL
    }

    #### for methylation
    if(FALSE){
      mygeneticprofile = grep("methylation", getGeneticProfiles(mycgds,mycancerstudy)[,1],value = T)
      mygeneticprofile="luad_tcga_pub_methylation_hm450"


        mETHYL <- cgdsr::getProfileData(mycgds, genes = Genes ,mygeneticprofile,mycaselist)

    }



    ###############
    # for mutations
    if(MutNeeded){
      mygeneticprofile = grep("mut", getGeneticProfiles(mycgds,mycancerstudy)[,1],value = T)
      MUT1 <- cgdsr::getProfileData(mycgds, genes = Genes ,mygeneticprofile,mycaselist)

      if(!FunctionalAnnot){
        if(nrow(MUT)>1){
          GeneID<-intersect(colnames(MUT),colnames(MUT1))
          MUT<-rbind(MUT[,GeneID],MUT1[,GeneID])
        } else{
          MUT<-rbind(MUT,MUT1)
        }
      } else{

       # MUT1[MUT1=="NaN"]<-NA
        MUT1[is.na(MUT1)]<-NaN

        #table(is.na(MUT1))
        # functional impact: from cgdsr, very strage...
       MUTa<-cgdsr::getMutationData(mycgds,mycaselist,mygeneticprofile, genes = Genes )
       MUTa$case_id<-gsub("-",".",MUTa$case_id)
       #table(MUTa$mutation_type);table(MUTa$functional_impact_score)
       MUTa[grep("Not Available",MUTa$functional_impact_score),"functional_impact_score"]<-"N"
       MUTa[is.na(MUTa$functional_impact_score),"functional_impact_score"]<-"N"
       MUTa[grep("Frame_Shift",MUTa$mutation_type),"functional_impact_score"]<-"H"
       MUTa[grep("fs",MUTa$amino_acid_change),"functional_impact_score"]<-"H"
       MUTa[grep("Nonsense",MUTa$mutation_type),"functional_impact_score"]<-"H"
       #GeNe<-MUTa$gene_symbol[1]
       ID<-rownames(MUT1)
       #MUTa[,c("case_id","gene_symbol","functional_impact_score")]
       # to change once....
       print("warning message ok")
       #GeNe<-colnames(MUT1)[ colnames(MUT1)%in%MUTa$gene_symbol][2]

       ImpBase<-c("N","L","M","H")
       ImpVal<-c(NaN,0.3,0.6,1)
      #MUT1[BIZ[,1],BIZ[,2]]

      for(GeNe in colnames(MUT1)[ colnames(MUT1)%in%MUTa$gene_symbol]){
         MUTf<-reshape2::dcast(data =  MUTa[MUTa$gene_symbol%in%GeNe,c("case_id","gene_symbol","functional_impact_score")],
                     formula = case_id~functional_impact_score,
                     value.var = "functional_impact_score" )

         levels(MUT1[,GeNe])<-c(levels(MUT1[,GeNe]),ImpVal)
         #MA<-c(1,match( ImpBase,colnames(MUTf) ))

         Imp<-colnames(MUTf[,-1,drop=F])[as.numeric(unlist(apply(MUTf[,-1,drop=F],1,function(m){
           L<-which(m!=0)
           if(length(L)>1){
             L<-max(L)
           } else if (length(L)==0){
             L<-"N"
           }
           return(L)
         }) ))]

         MUT1[MUTf$case_id,GeNe]<-ImpVal[match(Imp,ImpBase)]
         MUT1[ MUT1[MUTf$case_id,GeNe]==0&!is.na(MUT1[MUTf$case_id,GeNe]),GeNe]<-NaN
       #  MUT1[,GeNe]
       } # end loop

       #BIZ<-which(is.na(MUT1),arr.ind = T)
       MUT1[is.na(MUT1)]<-NaN

       #table(as.numeric(as.character(MUT1[,3])))
       # pool
       MUT1<-apply(MUT1,2,function(x) as.numeric(as.character(x)))
       rownames(MUT1)<-ID

       if(nrow(MUT)>1){
         GeneID<-intersect(colnames(MUT),colnames(MUT1))
         MUT<-rbind(MUT[,GeneID],MUT1[,GeneID])
       } else{
         MUT<-rbind(MUT,MUT1)
       }
      }

      Study[mycancerstudy,1]<-nrow(MUT1)

    } else {
      MUT<-NULL
    }

  }

  ####################
  # remove NAs
  # from RNaseq
  if(RNANeeded){
    GenesNAS<-union(GenesNAS, setdiff(colnames(EXP),Genes) )

  #if(RNANeeded){
    if(ClinicNeeded){
      CLINIC<-CLINIC[rownames(EXP),] # for the patients without RNAseq_V2
  #  }
  }

  if(length(GenesNAS)>0){
    print(paste("Discard", paste(GenesNAS,collapse = ", "),"because NAs in RNAseq"))
   # if(RNANeeded){
     EXP<-EXP[,!colnames(EXP)%in%GenesNAS]
   # }
    if(MutNeeded){
     # MUT<-MUT[,!colnames(MUT)%in%GenesNAS]
      MUT<-MUT[,colnames(EXP)]
    }
  }
}
  if(ClinicNeeded){
    # reorder the patients
    MUT=MUT[rownames(CLINIC),]
    EXP=EXP[rownames(CLINIC),]
  }

  # vizualization
  if(PDF){
    pdf("Plots.pdf")
    #heatmap(as.matrix(EXP))
    #VAR<-apply(EXP,2,var)
    #barplot(sort(VAR),las=2,cex.names=0.5)

    MED<-apply(EXP,2,median)
    boxplot(EXP[,order(MED)],las=2, cex.axis=0.5,
            main=paste("Models' genes expression \n across",Organ, "cBioportal TCGA samples"),
            ylab="normalized log(x+1) expression")
    MUT1<-MUT
    MUT1[which(MUT=="NaN",arr.ind = T)]<-NA
    MUT1<-apply(MUT1,1,is.na)
    GN<-apply(MUT1,1,function(x){table(x)["FALSE"]})
    MUT1<-MUT1[order(GN, na.last = F),]
    PN<-apply(MUT1,2,function(x){table(x)["FALSE"]})
    MUT1<-MUT1[,order(PN)]
    image(t(MUT1),col=c(1,0), main ="Pattern of mutations", axes=F)
    axis(2,at = seq(0,1,length.out = nrow(MUT1)),labels = rownames(MUT1),cex.axis=0.3,las=2,tick = F)

    dev.off()
  }

  return(list(MUT=MUT,EXP=EXP,CLINIC=CLINIC, STUDY=Study))
}

