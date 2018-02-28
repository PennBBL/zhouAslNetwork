ASL_sample <- read.csv("/Users/student/Desktop/dale/outliers/quality_group.csv")
data <- data.frame(c(ASL_sample["relMeanRMSMotion"],ASL_sample["negativeVoxels"],
                     ASL_sample["normCoverage"],ASL_sample["coregCrossCorr"], ASL_sample["coregJaccard"], 
                     ASL_sample["coregDice"], ASL_sample["coregCoverage"]))

# findOutlier <- function(data, cutoff =3) {
#   sds <- apply(data, 2, sd)
#   result <- mapply(function(d, s) {
#     which(d > cutoff * s)
#   }, data, sds)
#   result
# }
# 
# outliers <- findOutlier(data)

#outliers for relMeanRMSMotion
motion_sd <- sd(unlist(data[1]))*2
motion_mean <- mean(unlist(data[1]))
motion_cutoff <- motion_mean + motion_sd
index1<-which(data[1]>motion_cutoff)
ASL_sample[index1,1]

#outliers for negative voxels
NV_sd <- sd(unlist(data["negativeVoxels"]))*2
NV_mean <- mean(unlist(data["negativeVoxels"]))
NV_cutoff <- NV_mean + NV_sd
index2<-which(data[2]>NV_cutoff)
ASL_sample[index2,1]

#outliers for normCoverage
coverage_sd <- sd(unlist(data[3]))*2
coverage_mean <- mean(unlist(data[3]))
coverage_cutoff <- coverage_mean - coverage_sd
index3<-which(data[3]<coverage_cutoff)
ASL_sample[index3,1]

#outliers for coreg cross correlation
coregcc_sd <- sd(unlist(data[4]))*2
coregcc_mean <- mean(unlist(data[4]))
coregcc_cutoff <- coregcc_mean - coregcc_sd
index4<-which(data[4]<coregcc_cutoff)
ASL_sample[index4,1]

#outliers for coreg Jaccard
coregjacc_sd <- sd(unlist(data[5]))*2
coregjacc_mean <- mean(unlist(data[5]))
coregjacc_cutoff <- coregjacc_mean - coregjacc_sd
index5<-which(data[5]<coregjacc_cutoff)
ASL_sample[index5,1]

#outliers for coreg Jaccard
coregdice_sd <- sd(unlist(data[6]))*2
coregdice_mean <- mean(unlist(data[6]))
coregdice_cutoff <- coregdice_mean - coregdice_sd
index6<-which(data[6]<coregdice_cutoff)
ASL_sample[index6,1]

#outliers for coreg cov
coregcov_sd <- sd(unlist(data[7]))*2
coregcov_mean <- mean(unlist(data[7]))
coregcov_cutoff <- coregcov_mean - coregcov_sd
index7<-which(data[7]<coregcov_cutoff)
ASL_sample[index7,1]