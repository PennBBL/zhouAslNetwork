# zhouAslNetwork

Pre: Simple outlier script: asl_outliers.R
relMeanRMSMotion, negativeVoxels, normCoverage
2 SD: 129154 127417 110168 106880 114709 125554 126903
3SD: 110168 106880 114709 125554 126903

coregCrossCorr, coregJaccard, coregDice, coregCoverage
2 SD: 125554, 127305, 106880, 127417

## ASL networks

1. Export xcp and R library paths
```
export XCPEDIR=/home/rciric/xcpAccelerator/xcpEngine
export R_LIBS_USER=$R_LIBS_USER:/data/jux/BBL/applications-from-joy/Rlibraries/3.2:/data/jux/BBL/applications-from-joy/Rlibraries
export R_LIBS=$R_LIBS:/data/jux/BBL/applications-from-joy/Rlibraries/3.2:/data/jux/BBL/applications-from-joy/Rlibraries
```

2. Run design_offlineRecon_9P.dsn through xcp pipeline with n48 offline (minus 3 outliers) and n18 online subjects.
Check that mean perfusion and perfusion images exist in /cbf directory
```
${XCPEDIR}/xcpEngine -d /data/jux/BBL/projects/ASLnetwork/design_offlineRecon_9P.dsn -i /tmp/aslNets -c /data/jux/BBL/projects/ASLnetwork/n48_offlineRecon_aslCohort.csv -o /data/jux/BBL/projects/ASLnetwork -m c -t 2
${XCPEDIR}/xcpEngine -d /data/jux/BBL/projects/ASLnetwork/design_onlineRecon_9P.dsn -i /tmp/aslNets -c /data/jux/BBL/projects/ASLnetwork/n17_onlineRecon_aslCohort1.csv -o /data/jux/BBL/projects/ASLnetwork -m c -t 2
```

3. Collect all fcon connectivity matrix files and average them
```
ls */*/fcon/power264/*_network.txt > network_list.txt
awk '{a[FNR]+=$1;b[FNR]++;}END{for(i=1;i<=FNR;i++)print a[i]/b[i];}' `cat network_list.txt` > average_power264_network.txt
```

4. Generate 14x14 Power node network and compute within-between correlations
_wbNetbyNet.csv is 14x14 adjacency matrix
_wbNetWB.csv is within-between correlations
```
Rscript /home/rciric/xcpAccelerator/xcpEngine/utils/withinBetween.R -m /data/jux/BBL/projects/ASLnetwork/average_power264_network.txt -c /home/rciric/xcpAccelerator/xcpEngine/atlas/power264/power264CommunityAffiliation.1D -o /data/jux/BBL/projects/ASLnetwork/
```

## rs-BOLD networks (documentation in progress)

1. Run design (uses an older version of xcp and R)
```
export XCPEDIR=/data/jux/BBL/applications-from-joy/xcpEngine
export RPATH=/share/apps/R/R-3.1.1/bin/R
${XCPEDIR}/xcpEngine -d /data/joy/BBL/projects/zhouCbfNetworks/zhouCbfNetworksScripts/n101_restbold_20180117_prelim.dsn -m c
```

## Inspecting ASL, BOLD, FA, ODI, and ICVF (documentation in progress)

1. Generate lists of all subjects data in each modality and keep only subject ID (Bash commands not shown in code below)

BOLD
```
cd /data/joy/BBL/projects/zhouCbfNetworks/data/boldNetwork
ls */*/net/SchaeferPNC_200/*network.txt > procBold.txt
```

CBF
```
cd /data/joy/BBL/projects/zhouCbfNetworks/data/boldNetwork/cbfProc
ls */*/fcon/schaefer200/*network.txt > procCbf.txt
```

NODDI
```
cd /data/jux/BBL/projects/multishell_diffusion/processedData/multishellPipelineFall2017
ls */*/tractography/*_matrixts.csv > procNoddi.txt
```

2. Find common subjects among modalities (ceiling will be will be number of NODDI subjects). Relevant scripts:
```
cd /data/joy/BBL/projects/zhouCbfNetworks/zhouCbfNetworksScripts/findCommonSubjects.sh
```

3. Use selectSubjects.sh to write the subjects to include in analysis
```
/data/joy/BBL/projects/zhouCbfNetworks/zhouCbfNetworksScripts/selectSubjects.sh
```

## Run within- and between- modules analysis across edges (documentation in progress)

1. Run within_between_coupling.m Matlab script
```
cd /data/joy/BBL/projects/zhouCbfNetworks/zhouCbfNetworksScripts/within_between_coupling.m
```

2. Visualize within-module coupling as violin plots via violin_plot.R
```
/data/joy/BBL/projects/zhouCbfNetworks/zhouCbfNetworksScripts/violin_plot.R
```