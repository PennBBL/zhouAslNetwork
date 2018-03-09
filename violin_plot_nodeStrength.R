############################################################
## Code to create violin plots for within-module coupling ##
############################################################
require(ggplot2)

##############
## CBF - FA ##
##############
## Read in module-level coupling measures (Number of Subjects x Number of Modules)
withinMod_coupling <- read.csv("/Users/student/Desktop/dale/faCbf_commComm_nodeStrength.txt", header=FALSE, sep = ' ')
withinMod_coupling <- as.data.frame(withinMod_coupling[-1,])
withinMod_coupling <- as.data.frame(withinMod_coupling[-c(8:49)])
colnames(withinMod_coupling) <- c("within_VIS_coupling", "within_SOM_coupling", "within_DORS_coupling", "within_VENT_coupling", "within_LIM_coupling", "within_FPC_coupling", "within_DMN_coupling")

## Transpose dataframe
tdata <- as.data.frame(t(withinMod_coupling))
tdata$module <- as.factor(c(1:7))

## Melt dataframe
require(reshape)
mdata <- melt(tdata, id="module")

## Create violin plots

withinMod_coupling_plot <- ggplot(data=mdata, aes(x=factor(module), y=value)) + ggtitle("CBF-FA node strength coupling")+ theme(plot.title = element_text(size = 25,face="bold", hjust=0.5))
withinMod_coupling_plot + geom_violin(aes(y=value, 
fill=factor(module))) + stat_summary(fun.y=mean, geom="point", size=4, color="white") + theme(panel.grid.minor = element_blank(),
panel.background = element_blank()) + scale_fill_manual(values = c("purple3", "royalblue1", "green3", "violet", 
"yellow", "darkorange1", "firebrick3")) + theme(axis.line = element_line(colour = 'black', size = 1.5), 
axis.ticks.length = unit(.25, "cm")) + theme(legend.position="none") + scale_x_discrete(breaks=1:7, 
labels=c("VIS", "SOM", "DORS", "VENT", "LIM", "FPC", "DMN")) + xlab(NULL) + ylab("Spearman's rho") + ylim(-.2,0.5)

###############
## CBF - ODI ##
###############
## Read in module-level coupling measures (Number of Subjects x Number of Modules)
withinMod_coupling <- read.csv("/Users/student/Desktop/dale/odiCbf_commComm_nodeStrength.txt", header=FALSE, sep = ' ')
withinMod_coupling <- as.data.frame(withinMod_coupling)
withinMod_coupling <- as.data.frame(withinMod_coupling[-1,])
withinMod_coupling <- as.data.frame(withinMod_coupling[-c(8:49)])
colnames(withinMod_coupling) <- c("within_VIS_coupling", "within_SOM_coupling", "within_DORS_coupling", "within_VENT_coupling", "within_LIM_coupling", "within_FPC_coupling", "within_DMN_coupling")

## Transpose dataframe
tdata <- as.data.frame(t(withinMod_coupling))
tdata$module <- as.factor(c(1:7))

## Melt dataframe
require(reshape)
mdata <- melt(tdata, id="module")

## Create violin plots

withinMod_coupling_plot <- ggplot(data=mdata, aes(x=factor(module), y=value)) + ggtitle("CBF-ODI node strength coupling")+ theme(plot.title = element_text(size = 25,face="bold", hjust=0.5))
withinMod_coupling_plot + geom_violin(aes(y=value, 
fill=factor(module))) + stat_summary(fun.y=mean, geom="point", size=4, color="white") + theme(panel.grid.minor = element_blank(),
panel.background = element_blank()) + scale_fill_manual(values = c("purple3", "royalblue1", "green3", "violet", 
"yellow", "darkorange1", "firebrick3")) + theme(axis.line = element_line(colour = 'black', size = 1.5), 
axis.ticks.length = unit(.25, "cm")) + theme(legend.position="none") + scale_x_discrete(breaks=1:7, 
labels=c("VIS", "SOM", "DORS", "VENT", "LIM", "FPC", "DMN")) + xlab(NULL) + ylab("Spearman's rho") + ylim(-.2,0.5)

################
## CBF - ICVF ##
################
## Read in module-level coupling measures (Number of Subjects x Number of Modules)
withinMod_coupling <- read.csv("/Users/student/Desktop/dale/icvfCbf_commComm_nodeStrength.txt", header=FALSE, sep = ' ')
withinMod_coupling <- as.data.frame(withinMod_coupling)
withinMod_coupling <- as.data.frame(withinMod_coupling[-1,])
withinMod_coupling <- as.data.frame(withinMod_coupling[-c(8:49)])
colnames(withinMod_coupling) <- c("within_VIS_coupling", "within_SOM_coupling", "within_DORS_coupling", "within_VENT_coupling", "within_LIM_coupling", "within_FPC_coupling", "within_DMN_coupling")

## Transpose dataframe
tdata <- as.data.frame(t(withinMod_coupling))
tdata$module <- as.factor(c(1:7))

## Melt dataframe
require(reshape)
mdata <- melt(tdata, id="module")

## Create violin plots

withinMod_coupling_plot <- ggplot(data=mdata, aes(x=factor(module), y=value)) + ggtitle("CBF-ICVF node strength coupling")+ theme(plot.title = element_text(size = 25,face="bold", hjust=0.5))
withinMod_coupling_plot + geom_violin(aes(y=value, 
fill=factor(module))) + stat_summary(fun.y=mean, geom="point", size=4, color="white") + theme(panel.grid.minor = element_blank(),
panel.background = element_blank()) + scale_fill_manual(values = c("purple3", "royalblue1", "green3", "violet", 
"yellow", "darkorange1", "firebrick3")) + theme(axis.line = element_line(colour = 'black', size = 1.5), 
axis.ticks.length = unit(.25, "cm")) + theme(legend.position="none") + scale_x_discrete(breaks=1:7, 
labels=c("VIS", "SOM", "DORS", "VENT", "LIM", "FPC", "DMN")) + xlab(NULL) + ylab("Spearman's rho") + ylim(-.2,0.5)
