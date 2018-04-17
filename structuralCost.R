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

QA_df$Q <- (mod_df$Q- mean(mod_df$Q))/sd(mod_df$Q)
QA_df$globalCBF <- globalCBF[,2]

# Residualized Q linearly regressing out age
lm1 <- lm(QA_df$Q ~ QA_df$ageAtScan1)
QA_df$Q_resid <- residuals(lm1)

modularityByAge<- gam(Q ~ s(ageAtScan1,k=4) + dti64MeanRelRMS + sex + pcaslRelMeanRMSMotion,fx=TRUE,data=QA_df)
visreg(modularityByAge,"ageAtScan1", xlab="Age (months)", ylab="Modularity")

cbfByAge<- gam(globalCBF ~ s(ageAtScan1,k=4) + dti64MeanRelRMS + sex + pcaslRelMeanRMSMotion,fx=TRUE,data=QA_df)
visreg(cbfByAge,"ageAtScan1", xlab="Age (months)", ylab="Global CBF")
                  
globalModularityCbf<- gam(globalCBF ~ Q_resid + s(ageAtScan1,k=4) + dti64MeanRelRMS + sex + pcaslRelMeanRMSMotion,fx=TRUE,data=QA_df)
visreg(globalModularityCbf,"Q_resid", xlab="Modularity", ylab="Global CBF")

globalModularityCbf_noAge<- gam(globalCBF ~ Q + dti64MeanRelRMS + sex + pcaslRelMeanRMSMotion,fx=TRUE,data=QA_df)
visreg(globalModularityCbf_noAge,"Q", xlab="Modularity", ylab="Global CBF")

scatterplot3d(QA_df$sex,mod_df$Q,QA_df$globalCBF,xlab="Sex (1=Male; 2=Female)",ylab="Glasser Data-driven Modularity (Q)",zlab="global mean CBF (ml/100g/min)",main="Metabolic cost of Q by sex",pch=19,color="blue")

vis.gam(globalModularityCbf,view=c("sex","Q_resid"),type="response",theta=135,phi=0)
title("Sex differences in metabolic cost of modularity (data-driven Glasser)")

vis.gam(globalModularityCbf,view=c("ageAtScan1","Q_resid"),type="response",theta=150,phi=0)
title("Renegotiation over ages 8 to 23 for \n metabolic cost of modularity by age (data-driven Glasser)")

##########
#Schaefer#
##########

# load pcasl data

#Schaefer 200
schaefer_index <- read.csv('/home/rciric/xcpAccelerator/xcpEngine/atlas/schaefer200/schaefer200x7CommunityAffiliation.1D',header=FALSE)
schaefer_names <- read.csv('/home/rciric/xcpAccelerator/xcpEngine/atlas/schaefer200/schaefer200x7CommunityNames.txt',header=F)
pcasl_df <- read.csv("/data/jux/BBL/projects/ASLnetwork/results/all_roiquant.txt",header = TRUE,sep="\t")

colnames(pcasl_df)[3:204]=colnames(pcasl_df)[5:206]
colnames(pcasl_df)[1]="bblid"
pcasl_df<- as.data.frame(pcasl_df[1:202])
pcasl_df<-pcasl_df[-2]

# remove time 2 data
sum(duplicated(pcasl_df[,1]))
pcasl_df<- pcasl_df[!duplicated(pcasl_df[,1]),]

# merge final selection
cost_df <- merge(QA_df,pcasl_df,by=c("bblid"))
cost_df<-merge(dti_qa,cost_df,by=c("bblid"))

#negatives <- cost_df_proc[298:497]<0
#negatives_sum<- as.data.frame(apply(negatives,1,sum))

#Schaefer 400
pcasl_df <- read.csv("/data/jux/BBL/projects/ASLnetwork/data/cbfRoiquant400/all_roiquant.txt",header = F,sep="\t")
schaefer_index <- read.csv('/home/rciric/xcpAccelerator/xcpEngine/atlas/schaefer400/schaefer400x7CommunityAffiliation.1D',header=FALSE)
schaefer_names <- read.csv('/home/rciric/xcpAccelerator/xcpEngine/atlas/schaefer200/schaefer200x7CommunityNames.txt',header=F)

colnames(pcasl_df)[3:402]=colnames(pcasl_df)[1:400]
colnames(pcasl_df)[1]="bblid"
pcasl_df<-pcasl_df[-2]
pcasl_df<- as.data.frame(pcasl_df[1:401])

# remove time 2 data
sum(duplicated(pcasl_df[,1]))
pcasl_df<- pcasl_df[!duplicated(pcasl_df[,1]),]

# merge final selection
cost_df <- merge(QA_df,pcasl_df,by=c("bblid"))
cost_df<-merge(dti_qa,cost_df,by=c("bblid"))

#negatives <- cost_df_proc[298:497]<0
#negatives_sum<- as.data.frame(apply(negatives,1,sum))

#Modularity metrics
schaefer_index<- as.data.frame(schaefer_index)
modularity <- read.csv("/data/jux/BBL/projects/ASLnetwork/results/modularitySchaefer.txt",header=F,sep=" ")
modularity <-modularity[,1:2]
mod_df_colnames <- c("bblid","Q")
colnames(modularity) <- mod_df_colnames
mod_df <- merge(QA_df,modularity,by=c("bblid"))[names(modularity)]
cost_df<-merge(cost_df,mod_df,by=c("bblid"))[names(cost_df)]
QA_df <- merge(QA_df,cost_df,by=c("bblid"))


# Calculate subject global mean CBF, removing NAs
globalCBF <- cbind(cost_df$bblid, rowMeans(cost_df[346:745], na.rm=T))

colnames(globalCBF) <- c("bblid","globalCBF")
cor.test(globalCBF[,2],mod_df$Q,method = "spearman",na.rm=T)

QA_df$Q <- mod_df$Q
QA_df$globalCBF <- globalCBF[,2]

# Residualized Q linearly regressing out age
#lm1 <- lm(QA_df$Q ~ QA_df$ageAtScan1)
#QA_df$Q_resid <- residuals(lm1)

modularityByAge<- gam(Q ~ s(ageAtScan1.x,k=4) + dti64MeanRelRMS.x + sex.x + pcaslRelMeanRMSMotion.x,fx=TRUE,data=QA_df)
visreg(modularityByAge,"ageAtScan1.x", xlab="Age (months)", ylab="Modularity")

cbfByAge<- gam(globalCBF ~ s(ageAtScan1.x,k=4) + dti64MeanRelRMS.x + sex.x + pcaslRelMeanRMSMotion.x,fx=TRUE,data=QA_df)
visreg(cbfByAge,"ageAtScan1.x", xlab="Age (months)", ylab="Global CBF")

globalModularityCbf<- gam(globalCBF ~ Q + s(ageAtScan1.x,k=4) + dti64MeanRelRMS.x + sex.x + pcaslRelMeanRMSMotion.x,fx=TRUE,data=QA_df)
visreg(globalModularityCbf,"Q", xlab="Modularity", ylab="Global CBF")

vis.gam(globalModularityCbf,view=c("sex.x","Q"),type="response",theta=135,phi=0)
title("Sex differences in metabolic cost of modularity (Schaefer 200x7")

vis.gam(globalModularityCbf,view=c("ageAtScan1.x","Q"),type="response",theta=150,phi=0)
title("Renegotiation over ages 8 to 23 for \n metabolic cost of modularity by age (Schaefer 200x7)")


# Module level
#Visual
visual_CBF <- cost_df[346:745]
mindex<- t(schaefer_index==1)
visual_CBF<- as.data.frame(visual_CBF[,mindex])
globalCBF <- cbind(cost_df$bblid, rowMeans(visual_CBF, na.rm=T))
QA_df$globalCBF <- globalCBF[,2]
globalModularityCbf<- gam(globalCBF ~ Q + s(ageAtScan1.x,k=4) + dti64MeanRelRMS.x + sex.x + pcaslRelMeanRMSMotion.x,fx=TRUE,data=QA_df)
visreg(globalModularityCbf,"Q", xlab="Modularity", ylab="Global CBF")
summary(globalModularityCbf)

#Somatomotor
SM_CBF <- cost_df[346:745]
mindex<- t(schaefer_index==2)
SM_CBF<- as.data.frame(SM_CBF[,mindex])
globalCBF <- cbind(cost_df$bblid, rowMeans(SM_CBF, na.rm=T))
QA_df$globalCBF <- globalCBF[,2]
globalModularityCbf<- gam(globalCBF ~ Q + s(ageAtScan1.x,k=4) + dti64MeanRelRMS.x + sex.x + pcaslRelMeanRMSMotion.x,fx=TRUE,data=QA_df)
visreg(globalModularityCbf,"Q", xlab="Modularity", ylab="Global CBF")
summary(globalModularityCbf)

#Dorsal attention
dAtt_CBF <- cost_df[346:745]
mindex<- t(schaefer_index==3)
dAtt_CBF<- as.data.frame(dAtt_CBF[,mindex])
globalCBF <- cbind(cost_df$bblid, rowMeans(dAtt_CBF, na.rm=T))
QA_df$globalCBF <- globalCBF[,2]
globalModularityCbf<- gam(globalCBF ~ Q + s(ageAtScan1.x,k=4) + dti64MeanRelRMS.x + sex.x + pcaslRelMeanRMSMotion.x,fx=TRUE,data=QA_df)
visreg(globalModularityCbf,"Q", xlab="Modularity", ylab="Global CBF")
summary(globalModularityCbf)

#salience/Ventral attention
salience_CBF <- cost_df[346:745]
mindex<- t(schaefer_index==4)
salience_CBF<- as.data.frame(salience_CBF[,mindex])
globalCBF <- cbind(cost_df$bblid, rowMeans(salience_CBF, na.rm=T))
QA_df$globalCBF <- globalCBF[,2]
globalModularityCbf<- gam(globalCBF ~ Q + s(ageAtScan1.x,k=4) + dti64MeanRelRMS.x + sex.x + pcaslRelMeanRMSMotion.x,fx=TRUE,data=QA_df)
visreg(globalModularityCbf,"Q", xlab="Modularity", ylab="Global CBF")
summary(globalModularityCbf)

#limbic
limbic_CBF <- cost_df[346:745]
mindex<- t(schaefer_index==5)
limbic_CBF<- as.data.frame(limbic_CBF[,mindex])
globalCBF <- cbind(cost_df$bblid, rowMeans(limbic_CBF, na.rm=T))
QA_df$globalCBF <- globalCBF[,2]
globalModularityCbf<- gam(globalCBF ~ Q + s(ageAtScan1.x,k=4) + dti64MeanRelRMS.x + sex.x + pcaslRelMeanRMSMotion.x,fx=TRUE,data=QA_df)
visreg(globalModularityCbf,"Q", xlab="Modularity", ylab="Global CBF")
summary(globalModularityCbf)

#FPC
FPC_CBF <- cost_df[346:745]
mindex<- t(schaefer_index==6)
FPC_CBF<- as.data.frame(FPC_CBF[,mindex])
globalCBF <- cbind(cost_df$bblid, rowMeans(FPC_CBF, na.rm=T))
QA_df$globalCBF <- globalCBF[,2]
globalModularityCbf<- gam(globalCBF ~ Q + s(ageAtScan1.x,k=4) + dti64MeanRelRMS.x + sex.x + pcaslRelMeanRMSMotion.x,fx=TRUE,data=QA_df)
visreg(globalModularityCbf,"Q", xlab="Modularity", ylab="Global CBF")
summary(globalModularityCbf)

#DMN
DMN_CBF <- cost_df[346:745]
mindex<- t(schaefer_index==7)
DMN_CBF<- as.data.frame(DMN_CBF[,mindex])
globalCBF <- cbind(cost_df$bblid, rowMeans(DMN_CBF, na.rm=T))
QA_df$globalCBF <- globalCBF[,2]
globalModularityCbf<- gam(globalCBF ~ Q + s(ageAtScan1.x,k=4) + dti64MeanRelRMS.x + sex.x + pcaslRelMeanRMSMotion.x,fx=TRUE,data=QA_df)
visreg(globalModularityCbf,"Q", xlab="Modularity", ylab="Global CBF")
summary(globalModularityCbf)