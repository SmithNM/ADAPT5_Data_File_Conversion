#~~~~~~~~~~~~~~~~~~~~~~USER VARIABLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Data file 
DATA_FILE = "RAW_DATA.csv"
#Dosing File
DOSE_FILE = "Dosing_information_EDITED.csv"
Number_Doses_Per_Subject = 2                     #This input is only used if no data file is availble to create a data sheet
#Output File
Output_File = "ADAPT5_DATA.dat"
#Data file format & Components
COVARIATE_LIST = c("PB_0h",
                   "PB_Mtn",
                   "TK_Conc",
                   "TK",
                   "TAU_PB")
#THIS SHOULD ALWAYS HAVE DVID AT LEAST
DV_Labels      = c("DVID",
                   "DVIDTXT")

#Bacterial Experiments?
Bacterial_Experiment_1 = 1
#Plating Volume in microliters (uL)
Plating_Volume_uL = 50



#~~~~~~~~~~~~~~~~~~~~DON'T TOUCH BELOW THIS LINE~~~~~~~~~~~~~~~~~~~~~~~ 
{required.packages=c("dplyr",
                     "reshape2")
  
  is.installed = function(package){
    is.element(package,installed.packages()[,1])
  }
  
  for (i in 1:length(required.packages)){
    if (!is.installed(required.packages[i])){
      install.packages(required.packages[i],repos='https://cran.cnr.berkeley.edu/')
    }
  }
  
  require("dplyr")
  require("reshape2")
  rm(required.packages,is.installed)
}
#~~~~~~~~~~~~~~~~~~~~~~Data Reading~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
n.covariates = length(COVARIATE_LIST)
n.DV.labels = length(DV_Labels)
#COVARIATES
if(n.covariates>0){
  Covariates = read.csv(DATA_FILE,na.strings = ".", header = FALSE, nrows = (n.covariates))
  Covariates = Covariates[-c(seq(2,(1+n.DV.labels)))]
  n.subjects.covariates = ncol(Covariates) - 1
  colnames(Covariates) = as.character(c("CovariateID",seq(1,n.subjects.covariates)))
  Covariates = melt(Covariates,variable.name = "ID",value.name = "COV")
}else{
  Covariates=NULL
}
#DVs
DV = read.csv(DATA_FILE,na.strings = ".",header = FALSE, skip = n.covariates+1)
n.subjects.DV = ncol(DV)-length(DV_Labels) - 1
colnames(DV) = as.character(c("Time",DV_Labels,seq(1,n.subjects.DV)))
DV = melt(DV,variable.name = "ID",id.vars = c("Time",DV_Labels),value.name = "DV")
DV[is.na(DV)]=-1 #Want to keep NAs in the Covariate file, but can remove for the DVs
names(DV) = toupper(names(DV))
DV$DVID = as.factor(DV$DVID)
}
#~~~~~~~~~~~~~~~~~~~~~~Dose Reading~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
{
  
Read.Universal.Dose = function(DOSE_FILE){
  if (file.exists(DOSE_FILE)){
    Headers = try(read.table(DOSE_FILE,header=FALSE,na.strings = ".",nrows = 1,sep=",",comment.char = ""),silent=FALSE)
    Headers[,]=sub("# ","",as.matrix(Headers[1,]))
    Headers = (as.vector(Headers))
    Dose = try(read.table(DOSE_FILE,header = FALSE, comment.char="#", na.strings = ".",sep=",",skip=1),silent=TRUE)
    if ("try-error" %in% Dose){
      Dose=NA
      return(Dose)
    }else{
      colnames(Dose) = Headers
      return(Dose)
    }
  }else{
    return(FALSE)
  }
}
Create.Universal.Dose = function(SUBJECT_NUMBER,DOSING_NUMBER){
  #SUBJECT_NUMBER = 2
  #DOSING_NUMBER = 5
  if(missing(SUBJECT_NUMBER)){stop("NO SUBJECTS WERE DETECTED IN THE DATA FILE. UNABLE TO CREATE A DOSING TEMPLATE")}
  if(missing(DOSING_NUMBER)){DOSING_NUMBER=1}
  DOSE_FORM = data.frame("# ID"=rep(seq(1:SUBJECT_NUMBER),each=DOSING_NUMBER),
                         CMT=NA,
                         TIME=NA,
                         AMT=NA,
                         DURATION=NA,
                         RATE=NA,
                         TARGET_PEAK=NA,
                         VOLUME=NA,
                         TARGET_CONSTANT_CONC=NA,
                         ELIMINATION_HALF_LIFE=NA,
                         check.names = FALSE)
  write.csv(DOSE_FORM,file="Dosing_information.csv",na=".",row.names = FALSE)
  rm(list=ls())
  stop("DATA FORM CREATED")
}

Dose= try(Read.Universal.Dose(DOSE_FILE),silent=TRUE)
if(Dose==FALSE){Create.Universal.Dose(n.subjects.DV,Number_Doses_Per_Subject)}

}
#~~~~~~~~~~~~~~~~~~Writing the output table~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Make_Subject_Data = function(Dose, Covariates, DV,n){
  subj.list = list()
  #n=5
  
  Dose.Test = if(length(Dose)>0){filter(Dose,ID==n)}else{NULL}
  Cov.Test = if(length(Covariates)>0){(filter(Covariates,ID==n))}else{NULL}
  Cov.Test[is.na(Cov.Test)] = -1
  DV.Test = filter(DV,ID==n)
  DV.Test[is.na(DV.Test)] = -1
  DV_Labels = levels(DV$DVID)
  COVARIATE_LIST=if(length(Covariates)>0){try(levels(Cov.Test$CovariateID),silent=TRUE)}else{NULL}
  if (all(!is.na(Dose.Test$RATE)) & all(!is.na(Dose.Test$DURATION))){
    if(!all(Dose.Test$AMT/Dose.Test$DURATION==Dose.Test$RATE)){
      print("There is a mismatch between the dosing amount, rate, and duration.")
    }
  } else if(all(!is.na(Dose.Test$RATE))){
    Dose.Test$DURATION=Dose.Test$AMT/Dose.Test$RATE
  } else if(all(!is.na(Dose.Test$DURATION))){
    if(all(Dose.Test$DURATION==0)){}
    else {
      Dose.Test$RATE[Dose.Test$DURATION!=0] = Dose.Test$AMT[Dose.Test$DURATION!=0]/Dose.Test$DURATION[Dose.Test$DURATION!=0]    }
  } else{
    #print("YOU NEED To havE either A RATE OR DURATION; USE DURATION OF 0 TO INDICATE BOLUS")
  }
  
  n.bolus.input = nlevels(as.factor(Dose.Test$CMT[Dose.Test$DURATION==0]))
  n.infus.input = nlevels(as.factor(Dose.Test$CMT[Dose.Test$DURATION!=0]))
  n.bolus.events= length(Dose.Test$AMT[Dose.Test$DURATION==0])
  n.infus.events= length(Dose.Test$AMT[Dose.Test$DURATION!=0])
  #PArt1: SUBJECT NAME
  {
    sub.name=paste("SUBJECT",n,sep="_")
  }
  #PART2: INPUT NUMBERS
  {
    #Number Infusion/Covariate Inputs
    n.model.input=n.infus.input + length(COVARIATE_LIST)
    #Number Bolus Inputs
    n.bolus.input=n.bolus.input
    #Number Input Events
    Infusion_Start_Times=unique(Dose.Test$TIME[Dose.Test$DURATION!=0])
    Bolus_Times=unique(Dose.Test$TIME[Dose.Test$DURATION==0])
    Infusion_End_Times=sort(Infusion_Start_Times+Dose.Test$DURATION[Dose.Test$DURATION!=0])
    Input_Events=sort(unique(c(Bolus_Times,Infusion_Start_Times,Infusion_End_Times)))
    n.Input.Events=0
    if(length(Input_Events)>0){
      n.Input.Events=length(Input_Events)
    }else if (length(COVARIATE_LIST>0)){
      n.Input.Events=1
    }
  }
  #PART3: Vector of Bolus and Infusion Names
  {
    Input.Names = c("#IVAR")
    if(length(COVARIATE_LIST)>0){for (i in 1:length(COVARIATE_LIST)){
      Input.Names = c(Input.Names, paste("#",COVARIATE_LIST[i],sep = ""))
    }}
    if(n.infus.input>0){for (i in 1:(n.infus.input)){
      Input.Names = c(Input.Names,paste("#Infusion",i,sep=""))
    }}
    if(n.bolus.input>0){for (i in 1:(n.bolus.input)){
      Input.Names = c(Input.Names,paste("#Bolus",i,sep=""))
    }}
    Input.Names = t(data.frame(Input.Names))
  }
  #PART4: Matrix of Infusions, Covariates, and Infusions
  {
    #Part4a Infusions
    INPUT=data.frame()
    if (n.infus.input>0){
      Dose.Inf.St=filter(Dose.Test, DURATION>0)
      Dose.Inf.St=select(Dose.Inf.St,ID,TIME,CMT,AMT,DURATION,RATE)
      Dose.Inf.Ed=Dose.Inf.St
      Dose.Inf.Ed$TIME=unique(c(Infusion_End_Times))
      Dose.Inf.Ed$AMT=0
      Dose.Inf.Ed$DURATION=NA
      Dose.Inf.Ed$RATE=0
      Dose.Inf = (bind_rows(Dose.Inf.St,Dose.Inf.Ed))
      Dose.Inf = arrange(Dose.Inf,TIME)
      Dose.Inf=dcast(Dose.Inf,TIME~CMT,value.var="RATE")
      rm(Dose.Inf.St,Dose.Inf.Ed)
      INPUT=(Dose.Inf)
    }
    
    #Part4b Boluses
    if (n.bolus.input>0){
      Dose.Bol=filter(Dose.Test,DURATION==0)
      Dose.Bol=select(Dose.Bol,ID,TIME,CMT,AMT,DURATION,RATE)
      Dose.Bol=dcast(Dose.Bol,TIME~CMT,value.var = 'AMT')
      NonInf_Bolus_Times = Bolus_Times[!(Bolus_Times %in% c(Infusion_Start_Times,Infusion_End_Times))]
      INPUT=full_join(INPUT,Dose.Bol,by="TIME")
      INPUT=arrange(INPUT,TIME)
      
    }
    
    #Part4c Covariates
    if (n.covariates>0 & ((n.bolus.input>0)|(n.infus.input>0))){
      Cov.Test=dcast(Cov.Test,ID~CovariateID,value.var = "COV")
      Cov.Test=Cov.Test[rep(seq_len(nrow(Cov.Test)), each=nrow(INPUT)),]
      INPUT=bind_cols(INPUT[1],Cov.Test[2:ncol(Cov.Test)],INPUT[(2):ncol(INPUT)])
    }else if (n.covariates>0 & n.bolus.input==0 & n.infus.input==0){
      Cov.Test = dcast(Cov.Test,ID~CovariateID,value.var = "COV")
      Cov.Test$ID=0
      INPUT = Cov.Test
    }
  }
  #PART5: DV NUMBERS
  {
    n.DV.Type = length(DV_Labels)
    n.DV.Measure = length(DV.Test$TIME)
  }
  #PART6: Vector of DV Names
  {
    DV.Name=c("#IVAR")
    if (length(DV_Labels)==0){
      DV.Name=c(DV.Name,"#DV_1")
    }
    else{
      for (i in 1:length(DV_Labels)){
        DV.Name = c(DV.Name, paste("#DV_",DV_Labels[i],sep=""))
      }
    }
    DV.Name = t(data.frame(DV.Name))
  }
  #PART7: Matrix of DVs
  {
    if(length(DV_Labels)>0){
      DV.Test = dcast(DV.Test,TIME~DVID,value.var = "DV")
      #DV.Test[is.na(DV.Test)] = -1
    }
    else {
      DV.Test = select(DV.Test,TIME,DV)
    }
  }
  # STORING ALL VALUES FOR RETURN
  if (all(is.na(INPUT))){
    subj.list = list(sub.name,
                     n.model.input,
                     n.bolus.input,
                     n.Input.Events,
                     n.DV.Type,
                     n.DV.Measure,
                     DV.Name,
                     as.matrix(DV.Test)
                    )
  } else{
      subj.list = list(sub.name,
                       n.model.input,
                       n.bolus.input,
                       n.Input.Events,
                       Input.Names,
                       as.matrix(INPUT),
                       n.DV.Type,
                       n.DV.Measure,
                       DV.Name,
                       as.matrix(DV.Test)
      )      
    }
  return(subj.list)
}

write.table(c(),file=Output_File,sep=",",append=FALSE,col.names = FALSE,row.names = FALSE)

if(n.subjects.DV>1){
  for (i in 1:n.subjects.DV){
    Single_Subj = Make_Subject_Data(Dose,Covariates,DV,i)
    for (j in 1:length(Single_Subj)){
    write.table(Single_Subj[j],file=Output_File,sep=",",append = TRUE,col.names =FALSE,row.names = FALSE,na="")
    }
    print(paste("Subject ",i," of ",n.subjects.DV," has been written to the data file",sep=""))
  }
}else if(n.subjects.DV==1){
  Single_Subj = Make_Subject_Data(Dose,Covariates,DV,1)
  for (i in 2:length(Single_Subj)){
    write.table(Single_Subj[i],file=Output_File,sep=",",append = TRUE,col.names =FALSE,row.names = FALSE,na="")
  }
}else{}

print("Operation Complete")
#rm(list=ls())
