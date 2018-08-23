#================================================================================USER VARIABLES===================================================================================
#Script input information=========================================================================================================================================================

#Data file name (as a .csv file)
Data.File                 = "RAW_DATA.csv"

#Dosing File (as a .csv file); if no dosing file is provided one will be created with the number of dosing events as indicated below
#For the dosing file, use should include two of the three: AMT, RATE, OR DURATION. The script will automatically calculate the third. 
#You should also use a DURATION of '0' to indicate a bolus administration
Dose.File                 = "Dosing_information_EDITED.csv"
Number_Doses_Per_Subject  = 10

#THIS SHOULD ALWAYS HAVE DVID and DVIDTXT, AT LEAST
DV.Labels                 = c("DVID",
                              "DVIDTXT")

#Script output information==========================================================================================================================================================
#Output File
Output.File.Name          = "ADAPT5_DATA.dat"
#Folder name to put the Output Files
Output.Folder             = "Results"
#Folder name to put the Diagnostic Plots/Graphs
Plot.Folder               = "Figures"
             


#Bacterial Experiments?=============================================================================================================================================================
Bacterial.Experiment              = 0                     #Set this to 1 if the experiment is using bacterial growth data and you want to create a special data set
Plating.Volume.uL                 = 50                    #The volume of sample used to to determine the cfu/mL
Log.DV.Cut.Off                    = 1.99                  #cfu/mL values at or below this will automatically be converted into colony number on a plate (i.e., 1.60206 is 40 cfu/mL
                                                          # which is converted to 2 colonies in a 50 uL sample)
Automatically.Round.Colony.Count  = 1                     #Set this to 1 if you want to automactially round colony counts (i.e. number of colonys on a plate) to an integer value;
                                                          #the default is to round the data down
DV.Value.For.No.Colonies.Reported = 1                     #This is the value that indicates no counts were observed (Log(cfu/mL) = 1 is preferred)


#=====================================================================General Notes on using this script====================================================================================
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
#
# This script will take a user supplied data file (template will be generated if none is found in the same directory of this script) and a user supplied dosing
# file (template will be generated if none is found in the same directory of this script). 
#
# The script will then convert the data from a wide format (which is easier to input data into for most experimetnal designs) and convert it into a long, time-based
# format--for ADAPT 5. 
#
#
# IMPORTANT: This script comes with absolutely no warranty (neither implicit or explicit) for any task. 
# Every use or applicaiton of this script is the sole responsibility of the user.
#
# THis script was written by Dr. Nicholas Smith. 
# The first version of this script dates from June 1st, 2017
#
# The current version of the script - v 0.6 - was last edited August 1st, 2018
#
# This script can be FREELY DISTRBUTED to other users.  
# This script may not be modified or built into any compercial application without the prior 
# written permission of the author. 
#
#
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#=======================================================================REQUIRED FUNCTIONS FOR SCRIPT=======================================================================================
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
Create.Data.File = function(DATA.FILE="RAW_DATA.csv",...){
  #this function takes a single argument for the output file name and creates a file that is Tsuji Lab standard for data acquisition. This is a wide format, which indicates the independent
  #variable (time, by default), various DV labels (DVID and DVIDTXT, by default), and time-independent covariates. The time-independent covariates can be used to indicate subjects dose group, 
  #Dosing interval, gender, etc.
  
  cat("Creating template data file","\n")
  cat(paste("The file will be named: ",DATA.FILE),"\n")
  
  Col1 = c(paste("COV",seq(1:3)),"#TIME",rep(seq(from=0,to=48,by=2),2))
  Col2 = c(rep(NA,3),"DVID",rep(1,25),rep(2,25))
  Col3 = c(rep(NA,3),"DVIDTXT",rep("Cocnentration",25),rep("Effect",25))
  Col4 = c(rep(sample(2,1),3),rep(NA,51))
  
  Template.File = data.frame(Col1,Col2,Col3,Col4,Col4,Col4)
  if (file.exists(DATA.FILE)){
    write.table(Template.File,file = paste("NEW_",DATA.FILE),quote=FALSE,sep=",",na = ".",row.names = FALSE,col.names = FALSE)
    stop("This message was generated due to an error with the supplied data file. A new one has been generated for you. Because you already had a file present, the new one has been named with a 'NEW_' prefix.")
  }else{
    write.table(Template.File,file = DATA.FILE,quote=FALSE,sep=",",na = ".",row.names = FALSE,col.names = FALSE)
    stop("This message was generated due to an error with the supplied data file. A new one has been generated for you.")
  }
}
Read.Data.File = function(DATA.FILE="RAW_DATA.csv",DVLABS=c("DVID","DVIDTXT"),Bacterial.Experiment=0,DV.Value.For.No.Colonies.Reported=1,Return.Covariates.Only=0,
                          Return.Data.Only = 0,DV.NA.Char=-1, Plating.Volume.uL=50,Log.DV.Cut.Off=1.99,Automatically.Round.Colony.Count=1,...){
  #   This function takes a Microbiology lab-standardized data file (a template Data file will be created if none is detected), and converts it to a long format. The data
  #   file will be organized by time and can be made available for use in other programs. 
  
  
  #DATA.FILE = "RAW_DATA.csv"
  #DVLABS = c("DVID","DVIDTXT")
  
  if(!file.exists(DATA.FILE)){Create.Data.File(DATA.FILE)}
  
  RAW.DATA = read.csv(DATA.FILE,sep=",",header=FALSE,check.names = FALSE,na.strings = c(".","")) #Read all data in order to locate where Covariates and DVs are separated
  
  Data.Label.Row.Number = grep("TIME",RAW.DATA[[1]])
  
  RAW.DATA = read.csv(DATA.FILE,sep=",",header=FALSE,check.names = FALSE, comment.char = "#",na.strings = c(".","")) #Re-read the data, but remove the comment line
  
  #This segment of the script reads the Time-invariate covariates and saves as a separate data frame
  if (Data.Label.Row.Number>1){
      Covariate.Names = try(as.character(RAW.DATA[[1]][1:(Data.Label.Row.Number-1)]))
  } else {
      Covariate.Names = c()
  }
  
  
  if (length(Covariate.Names) > 0) {
    cat("Covariates Found!\n")
    cat(Covariate.Names,"\n")
  } else{
    cat("NO COVARIATES FOUND!\n") 
  }
  n.ID.COV = 0 # This variable will count the number of IDs/subjects based on the covariates provided
  if (length(Covariate.Names)>0){
    #Covariates = filter(RAW.DATA, V1==Covariate.Names)
    Covariates = read.csv(DATA.FILE,sep=",",header=FALSE,check.names = FALSE,colClasses = "character", comment.char = "#",na.strings = c(".",""),nrows = length(Covariate.Names))
    Covariates = t(Covariates)
    Covariates = as.data.frame(Covariates[-seq(1:(1+length(DVLABS))),])
    colnames(Covariates) = Covariate.Names
    rownames(Covariates) = seq(1:nrow(Covariates))
    Covariates = as.data.frame(Covariates)
    Covariates[["ID"]] = seq(1:nrow(Covariates))
    n.ID.COV = nrow(Covariates)
  } else {Covariates = c()}
  
  #This segment of the script reads the dependent variables and associated labels and converts into a long data format
  DV = RAW.DATA[seq(Data.Label.Row.Number,nrow(RAW.DATA)),]
  names(DV) = c("TIME",DVLABS,seq(1,(ncol(DV)-length(DVLABS)-1)))
  n.ID.DV = 0 # This variable will count the number of IDs/subjects based on the DVs provided
  
  DV = melt(DV,id.vars = c("TIME",DVLABS),variable.name = "ID",value.name = "DV")
  
  
  if (Bacterial.Experiment==1){
    
    DV[["DV"]][DV[["DV"]]==DV.Value.For.No.Colonies.Reported] = 0
    
    Temp.DV = DV
    
    Temp.DV = filter(Temp.DV, !is.na(DV) & DV<=Log.DV.Cut.Off)
    Temp.DV[["ColCount"]] = (10^Temp.DV[["DV"]])*(Plating.Volume.uL/1000)
    
    if (Automatically.Round.Colony.Count==1) {
      Temp.DV[["ColCount"]] = round(Temp.DV[["ColCount"]],digits=0)
      warning("You have chosen to round the colony counts automatically. You may wish to re-evaluate how much you trust your data before continuing with this course of action...")
    }
    
    Temp.DV[["DVID"]] = Temp.DV[["DVID"]]+1
    Temp.DV[["DVIDTXT"]] = paste(Temp.DV[["DVIDTXT"]],"_adj",sep="")
    
    Temp.DV[["DV"]] = Temp.DV[["ColCount"]]
    Temp.DV[["ColCount"]]=NULL
    
    DV = rbind(DV,Temp.DV)
    DV = arrange(DV,ID,DVID,TIME)
  }
  
  if (nrow(DV)>0) {
    n.ID.DV=max(as.numeric(DV[["ID"]]))
  }else{n.ID.DV=0}
  
  

  #DV[["DV"]][is.na(DV[["DV"]])]=DV.NA.Char
  if (length(Covariates)>0){
    Final.Data = merge(Covariates,DV, by = "ID")  
  } else {
    Final.Data = DV
  }
  
  Final.Data = Final.Data[!is.na(Final.Data[["DV"]]),]
  #if (Bacterial.Experiment.1 == 1) {
  #  Final.Data[["DV"]][Final.Data[["DV"]] == No.Colonies.Reported.Value] = 0
  #}
  
  if (n.ID.COV!=n.ID.DV) {
    cat("THE NUMBER OF INDIVIDUALS PRESENT IN THE COVARIATE LIST DOES NOT MATCH THOSE IN THE DATA LIST. PLEASE BE CAREFUL AS THE RESULTS MAY BE INCORRECT! \n")
    cat("If the file contained no covariates, you may disregard this message...I guess \n")
  }
  
  write.csv(Final.Data,file=paste(gsub(".csv","",DATA.FILE),"_Long_Format.csv",sep=""),row.names = FALSE,na = ".")
  
  if (Return.Covariates.Only==1) {
    return(Covariates)
  } else if (Return.Data.Only==1) {
    return(DV)  
  } else {
    return(Final.Data)
  }
  
}
#~~~~~~~~~~~~~~~~~~~~~~Dose Reading~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Read.Universal.Dose = function(Dose.File){
  if (file.exists(Dose.File)){
    Headers = try(read.table(Dose.File,header=FALSE,na.strings = ".",nrows = 1,sep=",",comment.char = ""),silent=FALSE)
    Headers[,]=sub("# ","",as.matrix(Headers[1,]))
    Headers = (as.vector(Headers))
    Dose = try(read.table(Dose.File,header = FALSE, comment.char="#", na.strings = ".",sep=",",skip=1),silent=TRUE)
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
  write.csv(DOSE_FORM,file="Dosing_information.csv",na=".",row.names = FALSE,quote=FALSE)
  rm(list=ls())
  stop("DATA FORM CREATED")
}
#~~~~~~~~~~~~~~~~~~~~~~Subject Data Creation~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make_Subject_Data = function(Dose, Covariates, DV,n){
  #This function takes four arguments, all mandatory. It takes a population dosing data frame (long format), a covariate data frame (long format), a dependent
  #variable data frame (contaiing only DV related information), and a subject ID number. This function requires that the subject ID=n must be present in all of
  #the supplied arguments
  
  
  subj.list = list()
  #n=5
  
  Covariates[is.na(Covariates)] = -1
  
  Dose.Test = if(length(Dose)>0){filter(Dose,ID==n)}else{NULL}
  
  Cov.Test = (if(length(Covariates)>0){(filter(Covariates,ID==n))}else{NULL})
  Cov.Test[is.na(Cov.Test)] = -1
  
  DV.Test = filter(DV,ID==n)
  DV.Test[is.na(DV.Test)] = -1
  DV.Labels = levels(as.factor(DV$DVID))
  Covariate.List=if(length(Covariates)>0){try(levels(Cov.Test$CovariateID),silent=TRUE)}else{NULL}
  
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
  
  n.covariates = if (length(Covariates)>0) {nrow(Cov.Test)} else {0}
  
  n.bolus.input.all = nlevels(as.factor(Dose$CMT[ceiling(Dose$DURATION)==0])) 
  n.infus.input.all = nlevels(as.factor(Dose$CMT[Dose$DURATION!=0]))
  n.bolus.input = nlevels(as.factor(Dose.Test$CMT[ceiling(Dose.Test$DURATION)==0])) #number of bolus inputs for this specific (ADAPT 5 mandates that each subject have same number of bolus inputs, even if dose is 0 for that input)
  n.infus.input = nlevels(as.factor(Dose.Test$CMT[Dose.Test$DURATION!=0])) #number of infusion inputs for this specific (see above for rationale)
  n.bolus.events= length(Dose.Test$AMT[Dose.Test$DURATION==0])        #number of bolus events for this specific subject
  n.infus.events= length(Dose.Test$AMT[Dose.Test$DURATION!=0])        #number of infusion START events for this subject
  #PArt1: SUBJECT NAME
  {
    sub.name=paste("SUBJECT",n,sep="_")
  }
  #PART2: INPUT NUMBERS
  {
    #Number Infusion/Covariate Inputs
    n.model.input=n.infus.input.all + length(Covariate.List)
    #Number Bolus Inputs
    #n.bolus.input=n.bolus.input.all #pretty sure this line is not necessary and is causing an error...but who knows...
    #Number Input Events
    Infusion_Start_Times=unique(Dose.Test$TIME[Dose.Test$DURATION!=0])
    Bolus_Times=unique(Dose.Test$TIME[Dose.Test$DURATION==0])
    Infusion_End_Times=sort(Infusion_Start_Times+Dose.Test$DURATION[Dose.Test$DURATION!=0])
    Input_Events=sort(unique(c(Bolus_Times,Infusion_Start_Times,Infusion_End_Times)))
    n.Input.Events=0
    if(length(Input_Events)>0){
      n.Input.Events=length(Input_Events)
    }else if (length(Covariate.List>0)){
      n.Input.Events=1
    }
  }
  #PART3: Vector of Bolus and Infusion Names
  {
    Input.Names = c("#IVAR")
    if(length(Covariate.List)>0){for (i in 1:length(Covariate.List)){
      Input.Names = c(Input.Names, paste("#",Covariate.List[i],sep = ""))
    }}
    if(n.infus.input.all>0){for (i in 1:(n.infus.input.all)){
      Input.Names = c(Input.Names,paste("#Infusion",i,sep=""))
    }}
    if(n.bolus.input.all>0){for (i in 1:(n.bolus.input.all)){
      Input.Names = c(Input.Names,paste("#Bolus",i,sep=""))
    }}
    Input.Names = t(data.frame(Input.Names))
  }
  #PART4: Matrix of Infusions, Covariates, and Infusions
  {
    #Part4a Infusions
    INPUT=data.frame()
    if (n.infus.input.all>0 & n.infus.input==n.infus.input.all){
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
    } else if (n.infus.input>0){
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
      CMT.for.this.subject = names(Dose.Inf)[-1]
      CMT.to.add.to.subject = unique(Dose$CMT)[(!(unique(Dose$CMT) %in% CMT.for.this.subject))]
      for (i in 1: length(CMT.to.add.to.subject)){
        Dose.Inf[[paste(CMT.to.add.to.subject[i])]] = 0
      }
      
      new.order.inf = c(1)
      for (i in 1:n.infus.input.all){
        new.order.inf = c(new.order.inf,which(names(Dose.Inf)==i))
      }
      Dose.Inf = Dose.Inf[new.order.inf]
      INPUT=(Dose.Inf)
    } else if (n.infus.input.all>0){
      Dose.Inf = data.frame("TIME"=0)
      for (i in 1:n.infus.input.all){
        Dose.Inf[[paste(i)]]=0
      }
      INPUT=(Dose.Inf)
    }
    
    #Part4b Boluses
    if (n.bolus.input.all>0 & n.bolus.input==n.bolus.input.all){
      Dose.Bol=filter(Dose.Test,DURATION==0)
      Dose.Bol=select(Dose.Bol,ID,TIME,CMT,AMT,DURATION,RATE)
      Dose.Bol=dcast(Dose.Bol,TIME~CMT,value.var = 'AMT')
      NonInf_Bolus_Times = Bolus_Times[!(Bolus_Times %in% c(Infusion_Start_Times,Infusion_End_Times))]
      INPUT=full_join(INPUT,Dose.Bol,by="TIME")
      INPUT=arrange(INPUT,TIME)
      
    } else if (n.bolus.input>0) {
      Dose.Bol=filter(Dose.Test,DURATION==0)
      Dose.Bol=select(Dose.Bol,ID,TIME,CMT,AMT,DURATION,RATE)
      Dose.Bol=dcast(Dose.Bol,TIME~CMT,value.var = 'AMT')
      NonInf_Bolus_Times = Bolus_Times[!(Bolus_Times %in% c(Infusion_Start_Times,Infusion_End_Times))]
      
      CMT.for.this.subject = names(Dose.Bol)[-1]
      CMT.to.add.to.subject = unique(Dose$CMT)[(!(unique(Dose$CMT) %in% CMT.for.this.subject))]
      for (i in 1: length(CMT.to.add.to.subject)){
        Dose.Bol[[paste(CMT.to.add.to.subject[i])]] = 0
      }
      
      new.order.bol = c(1)
      for (i in 1:n.bolus.input.all){
        new.order.bol = c(new.order.bol,which(names(Dose.Bol)==i))
      }
      Dose.Bol = Dose.Bol[new.order.bol]
      
      
      INPUT=full_join(INPUT,Dose.Bol,by="TIME")
      INPUT=arrange(INPUT,TIME)
      
      
    } else if (n.bolus.input.all>0){
      Dose.Bol = data.frame("TIME"=0)
      for (i in 1:n.bolus.input.all){
        Dose.Bol[[paste(i)]]=0
      }
      INPUT=full_join(INPUT,Dose.Bol,by="TIME")
      INPUT=arrange(INPUT,TIME)
    }
    
    #Part4c Covariates
    if (n.covariates>0 & ((n.bolus.input.all>0)|(n.infus.input.all>0))){
      Cov.Test=dcast(Cov.Test,ID~CovariateID,value.var = "COV")
      Cov.Test=Cov.Test[rep(seq_len(nrow(Cov.Test)), each=nrow(INPUT)),]
      INPUT=bind_cols(INPUT[1],Cov.Test[2:ncol(Cov.Test)],INPUT[(2):ncol(INPUT)])
    }else if (n.covariates>0 & n.bolus.input.all==0 & n.infus.input.all==0){
      Cov.Test = dcast(Cov.Test,ID~CovariateID,value.var = "COV")
      Cov.Test$ID=0
      INPUT = Cov.Test
    }
    
    INPUT[is.na(INPUT)] = 0
  }
  #PART5: DV NUMBERS
  {
    n.DV.Type = length(DV.Labels)
    n.DV.Measure = length(DV.Test$TIME)
  }
  #PART6: Vector of DV Names
  {
    DV.Name=c("#IVAR")
    if (length(DV.Labels)==0){
      DV.Name=c(DV.Name,"#DV_1")
    }
    else{
      for (i in 1:length(DV.Labels)){
        DV.Name = c(DV.Name, paste("#DV_",DV.Labels[i],sep=""))
      }
    }
    DV.Name = t(data.frame(DV.Name))
  }
  #PART7: Matrix of DVs
  {
    if(length(DV.Labels)>0){
      DV.Test = dcast(DV.Test,TIME~DVID,value.var = "DV")
      DV.Test$TIME = as.numeric(as.character(DV.Test$TIME))
      DV.Test = arrange(DV.Test,TIME)
      #DV.Test[is.na(DV.Test)] = -1
      DV.Test[is.na(DV.Test)]=-1
      DVID.for.this.subject = names(DV.Test)[-1]
      DVID.to.add.to.subject = unique(DV$DVID)[(!(unique(DV$DVID) %in% DVID.for.this.subject))]
      
      if (length(DVID.to.add.to.subject)>0){
        for (i in 1: length(DVID.to.add.to.subject)){
          DV.Test[[paste(DVID.to.add.to.subject[i])]] = -1
        }
        
        new.order.DV = c(1)
        for (i in 1:n.DV.Type){
          new.order.DV = c(new.order.DV,which(names(DV.Test)==i))
        }
        DV.Test = DV.Test[new.order.DV]
      }
    }
    else {
      DV.Test = select(DV.Test,TIME,DV)
    }
  }
  # STORING ALL VALUES FOR RETURN
  if (all(is.na(INPUT))){
    subj.list = list(sub.name,
                     n.model.input,
                     n.bolus.input.all,
                     n.Input.Events,
                     n.DV.Type,
                     n.DV.Measure,
                     DV.Name,
                     as.matrix(DV.Test)
    )
  } else{
    subj.list = list(sub.name,
                     n.model.input,
                     n.bolus.input.all,
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



#=============================================================================MAIN========================================================================================

Covariates = Read.Data.File(Data.File,DVLABS= DV.Labels,Bacterial.Experiment,Return.Covariates.Only = 1)
if (length(Covariates)>0){Covariates = melt(Covariates, id.vars =  "ID",variable.name="CovariateID",value.name = "COV")}
DV = Read.Data.File(Data.File,DVLABS = DV.Labels,Bacterial.Experiment,Return.Data.Only = 1)
n.subjects.DV = nlevels(DV[["ID"]])
Dose= try(Read.Universal.Dose(Dose.File),silent=TRUE)
if(Dose==FALSE){Create.Universal.Dose(n.subjects.DV,Number_Doses_Per_Subject)}
#~~~~~~~~~~~~~~~~~~Writing the output table~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(file.exists(Output.Folder)){unlink(Output.Folder,recursive = TRUE)}
Main.WD = getwd()
dir.create(Output.Folder)
Results.WD = file.path(Main.WD,Output.Folder)
setwd(Results.WD)

file.rename(file.path(Main.WD,paste(sub(".csv","",Data.File),"_Long_Format.csv",sep="")),
            file.path(Results.WD,paste(sub(".csv","",Data.File),"_Long_Format.csv",sep=""))
)

write.table(c(),file=Output.File.Name,sep=",",append=FALSE,col.names = FALSE,row.names = FALSE,quote=FALSE)
 
if(nlevels(DV[["ID"]])>1){
  for (i in 1:nlevels(DV[["ID"]])){
    Single_Subj = Make_Subject_Data(Dose,Covariates,DV,i)
    for (j in 1:length(Single_Subj)){
    write.table(Single_Subj[j],file=Output.File.Name,sep=",",append = TRUE,col.names =FALSE,row.names = FALSE,na="",quote=FALSE)
    }
    cat(paste("Subject ",i," of ",n.subjects.DV," has been written to the data file \n",sep=""))
  }
}else if(nlevels(DV[["ID"]])==1){
  Single_Subj = Make_Subject_Data(Dose,Covariates,DV,1)
  for (i in 2:length(Single_Subj)){
    write.table(Single_Subj[i],file=Output.File.Name,sep=",",append = TRUE,col.names =FALSE,row.names = FALSE,na="",quote=FALSE)
  }
}else{stop("The number of subjects was found to be 0. This is because of incorrect DV Labels or improper data file creation. Please confirm data file format.")}

cat("Operation Complete  \n")
setwd(Main.WD)
#rm(list=ls())
