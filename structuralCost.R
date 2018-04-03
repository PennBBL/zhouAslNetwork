##############################################################################################
###   This script will construct a subject cohort for CBF-DTI analyses, retaining   	     ###
## 			subjects who passed quality assurance protocols for T1, DTI, and PCASL      	      ##
##############################################################################################
# T1
t1_qa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/t1struct/n1601_t1QaData_20170306.csv")
# DTI 
dti_qa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/n1601_dti_qa_20170301.csv")
# Resting-state QA
rest_qa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/rest/n1601_RestQAData_20170714.csv")
# ASL QA
asl_qa <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/asl/n1601_PcaslQaData_20170403.csv")
# B0 Acquisition
protVal <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/n1601_pnc_protocol_validation_params_status_20161220.csv")
# LTN Status
health <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/health/n1601_health_20161214.csv")
# Demographics
demo <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv")
# Cognitive Scores
cognitive <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/cnb/n1601_cnb_factorscores.csv")

# Subject identifier
tracker <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/demographics/n1601_tracker_go1_20161212.csv")
tracker <- tracker[c("bblid","scanid")]

############################
## Merge demographic data ##
############################
demo <- merge(demo, health, by=c("bblid","scanid"))
demo <- merge(demo, t1_qa, by=c("bblid","scanid"))
demo <- merge(demo, dti_qa, by=c("bblid","scanid"))
demo <- merge(demo, rest_qa, by=c("bblid","scanid"))
demo <- merge(demo, asl_qa, by=c("bblid","scanid"))
demo <- merge(demo, protVal, by=c("bblid","scanid"))
demo <- merge(demo, cognitive, by=c("bblid","scanid"))

############################
## Sample Exclusions ##
############################
QA_df <- demo

# DTI exclusions
QA_df <- QA_df[which(QA_df$dti64ProtocolValidationStatus == 1), ]
QA_df <- QA_df[which(QA_df$b0ProtocolValidationStatus==1), ]
QA_df <- QA_df[which(QA_df$dti64Exclude == 0), ]

# Health exclusion
QA_df <- QA_df[which(QA_df$healthExcludev2==0), ]

# PCASL Exclusion
QA_df <- QA_df[which(QA_df$pcaslExclude==0), ]
QA_df <- QA_df[which(QA_df$pcaslNoDataExclude==0), ]
QA_df <- QA_df[which(QA_df$pcaslRelMeanRMSMotionExclude==0), ]
QA_df <- QA_df[which(QA_df$pcaslTSNRExclude==0), ]
QA_df <- QA_df[which(QA_df$pcaslNVolumesAcquiredExclude==0), ]
QA_df <- QA_df[which(QA_df$pcaslMeanGMValueExclude==0), ]
QA_df <- QA_df[which(QA_df$pcaslCoverageExclude==0), ]

# T1 Exclusion
QA_df <- QA_df[which(QA_df$t1Exclude==0), ]

###########################################
## Global, module, and regional analysis ##
###########################################
library(ppcor)
library(mgcv)
library(visreg)
library(akima)
library(scatterplot3d)

## To do: control for GM density for CBF? Subset by euclidian distance bins?
## Perfusion-weighted structural networks?

# Calculate regional mean CBF

#########
#Glasser#
#########

# Merge final selection 
pcasl_df <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/asl/n1601_glasserPcaslValues_20161014.csv")

cost_df_proc <- merge(QA_df,pcasl_df,by=c("bblid","scanid"))
vars <- ls(cost_df_proc)
vars <- vars[215:574]
cost_df <- cbind(cost_df_proc[vars])
scanid <- cost_df_proc$scanid
cost_df <- cbind(scanid, cost_df)

# Count number of NAs-- 118 total, see na_count for table of NAs regionally
na_count <-sapply(cost_df, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)

# Subject scan IDs with NAs
na_index <- which(is.na(cost_df), arr.ind=TRUE)
na_index <- na_index[,1]

#Modularity metrics
modularity <- read.csv("/data/joy/BBL/projects/zhouCbfNetworks/results/modularity.txt",header=F,sep=" ")
mod_df_colnames <- c("scanid","Q",sprintf("%01d", seq(1,360)))
colnames(modularity) <- mod_df_colnames
mod_df <- merge(QA_df,modularity,by=c("scanid"))[names(modularity)]
cost_df<-merge(cost_df,mod_df,by=c("scanid"))[names(cost_df)]
QA_df <- merge(QA_df,cost_df,by=c("scanid"))

# Load in Glasser index
glasser_index <- read.csv('/home/rciric/xcpAccelerator/xcpEngine/atlas/glasser360/glasser360NodeIndex.1D')
glasser_names <- read.csv('/home/rciric/xcpAccelerator/xcpEngine/atlas/glasser360/glasser360NodeNames.txt')

# Calculate subject global mean CBF, removing NAs
globalCBF <- cbind(cost_df$scanid, rowMeans(cost_df, na.rm=T))

colnames(globalCBF) <- c("scanid","globalCBF")
cor.test(globalCBF[,2],mod_df$Q,method = "spearman",na.rm=T)

merge(globalCBF, mod_df, by=c("scanid"))[colnames(globalCBF)]
covs <- cbind(QA_df$ageAtScan1,QA_df$sex, QA_df$dti64MeanRelRMS, QA_df$pcaslRelMeanRMSMotion)
covs <- na.omit(covs)
pcor.test(na.omit(globalCBF), na.omit(mod_df$Q),covs)

QA_df<- merge(mod_df,QA_df,by=c("scanid"))
QA_df<- merge(globalCBF,QA_df,by=c("scanid"))

modularityByAge<- gam(Q ~ s(ageAtScan1,k=4) + dti64MeanRelRMS + sex + pcaslRelMeanRMSMotion,fx=TRUE,data=QA_df)
visreg(modularityByAge,"ageAtScan1", xlab="Age (months)", ylab="Modularity")

cbfByAge<- gam(globalCBF ~ s(ageAtScan1,k=4) + dti64MeanRelRMS + sex + pcaslRelMeanRMSMotion,fx=TRUE,data=QA_df)
visreg(cbfByAge,"ageAtScan1", xlab="Age (months)", ylab="Global CBF")

globalModularityCbf<- gam(globalCBF ~ Q + s(ageAtScan1,k=4) + dti64MeanRelRMS + sex + pcaslRelMeanRMSMotion,fx=TRUE,data=QA_df)
visreg(globalModularityCbf,"Q", xlab="Modularity", ylab="Global CBF")

globalModularityCbf_noAge<- gam(globalCBF ~ Q + dti64MeanRelRMS + sex + pcaslRelMeanRMSMotion,fx=TRUE,data=QA_df)
visreg(globalModularityCbf,"Q", xlab="Modularity", ylab="Global CBF")

scatterplot3d(QA_df$sex,mod_df$Q,QA_df$globalCBF,xlab="Sex (1=Male; 2=Female)",ylab="Glasser Data-driven Modularity (Q)",zlab="global mean CBF (ml/100g/min)",main="Metabolic cost of Q by sex",pch=19,color="blue")

vis.gam(globalModularityCbf,view=c("sex","Q"),type="response",theta=115,phi=0)
title("Sex differences in metabolic cost of modularity (data-driven Glasser)")

vis.gam(globalModularityCbf,view=c("ageAtScan1","Q"),type="response",theta=10,phi=0)
title("Renegotiation over ages 8 to 23 for \n metabolic cost of modularity by age (data-driven Glasser)")

##########
#Schaefer#
##########

# Work in progress

# Merge final selection 
pcasl_df <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/asl/n1601_glasserPcaslValues_20161014.csv")

cost_df_proc <- merge(QA_df,pcasl_df,by=c("bblid","scanid"))
vars <- ls(cost_df_proc)
vars <- vars[215:574]
cost_df <- cbind(cost_df_proc[vars])
scanid <- cost_df_proc$scanid
cost_df <- cbind(scanid, cost_df)

# Count number of NAs-- 118 total, see na_count for table of NAs regionally
na_count <-sapply(cost_df, function(y) sum(length(which(is.na(y)))))
na_count <- data.frame(na_count)

# Subject scan IDs with NAs
na_index <- which(is.na(cost_df), arr.ind=TRUE)
na_index <- na_index[,1]

#Modularity metrics
modularity <- read.csv("/data/joy/BBL/projects/zhouCbfNetworks/results/modularitySchaefer.txt",header=F,sep=" ")
mod_df_colnames <- c("scanid","Q")
colnames(modularity) <- mod_df_colnames
mod_df <- merge(QA_df,modularity,by=c("scanid"))[names(modularity)]
cost_df<-merge(cost_df,mod_df,by=c("scanid"))[names(cost_df)]
QA_df <- merge(QA_df,cost_df,by=c("scanid"))