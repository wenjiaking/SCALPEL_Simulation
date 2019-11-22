#========================================================================================
#========================= Simulation 3 =================================================
#========================================================================================

png("Desktop/simulation3/Data/rawframe%04d.png")
for (i in 1:1000) {
  imageVec(Ynew[,i],100,col=gray.colors(500))
}
dev.off()
#ffmpeg -r 30 -i neuronsig%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p neuronfilm2.mp4
mysim3_step0=step0(outputDir = "Desktop/simulation3/myfunSim/",rawdataDir = "./",videoheight=100)
mysim3_step1=step1(mysim3_step0, minSize=25,maxSize=500,maxWidth=30,maxHeight=30,thresholdval=NULL)
mysim3_step2=step2(step1output = mysim3_step1)
mysim3_step3=step3(mysim3_step2)

sim3Scalpel=scalpel("Desktop/simulation3/scalpelSim/","./",videoHeight = 100)
sim3Scalpel=getScalpelStep3("Desktop/simulation3/scalpelSim/")
Ydeltaf=readRDS("Desktop/simulation3/scalpelSim/Step0Data/Ydeltaf_part1.rds")

imageVec(mysim3_step0$deltafY[,420],100,col = gray.colors(100))
imageVec(Ydeltaf[,420],100,col=gray.colors(100))

temporalplot(mysim3_step3$Zhat,ylab = "intensity",title="Estimated Z")

par(mfrow=c(1,2))
spatialplot(sim3Scalpel$Afilter,sim3Scalpel$Zhat,title="Estimated Footprints",border = FALSE,addToPlot = FALSE)
temporalplot(sim3Scalpel$Zhat,title="Estimated Temporal Signals")

plotResults(scalpelOutput = sim3Scalpel,titleA = "Estimated Footprints",titleZ = "Estimated Temporal Signals") 

par(mfrow=c(1,2))
spatialplot(Sim_neuronalInfo2$spatialmatrix,Sim_neuronalInfo2$temporalmatrix,title="Neuron Footprints",border = FALSE,addToPlot = FALSE)
temporalplot(Sim_neuronalInfo2$temporalmatrix,title="Neuron Temporal Signals")
temporalplot(sim3Scalpel$Zhat[c(2,17,6,1,9,10,18,7,13,3,20,19,14,15,11,16,12,8,21,4),],title="Corresponding Estimates")
sim3ScalpelAdj=scalpelStep3(getScalpelStep2("Desktop/simulation3/scalpelSim/"),alpha = 0.2)
spatialplot(sim3ScalpelAdj$Afilter,sim3ScalpelAdj$Zhat,title="Estimated Footprints (alpha=0.2)",border = FALSE,addToPlot = FALSE)
temporalplot(sim3ScalpelAdj$Zhat,title="Temporal Estimates (alpha=0.2)")
#Based on the plots above, it seems that the fake neuron identified by scalpel is among 2,3,4 and 5. 
#We check which one is fake by ploting their responding brightframes
par(mfrow=c(2,2))
plotBrightest(scalpelOutput = sim3Scalpel, AfilterIndex = 2,brightIndex = 1)#819,820,818,35,817
plotBrightest(scalpelOutput = sim3Scalpel, AfilterIndex = 3,brightIndex = 1) #303,424,425,302,304
plotBrightest(scalpelOutput = sim3Scalpel, AfilterIndex = 4,brightIndex = 1)#396
plotBrightest(scalpelOutput = sim3Scalpel, AfilterIndex = 5,brightIndex = 1)#407
which(sim3Scalpel$clusterID==5) #418 420 422 423
imageVec(sim3Scalpel$Azero[,418],100,col = gray.colors(2))
title(main="Preliminary Element 418")
imageVec(sim3Scalpel$Azero[,420],100,col = gray.colors(2))
title(main="Preliminary Element 420")
imageVec(sim3Scalpel$Azero[,422],100,col = gray.colors(2))
title(main="Preliminary Element 422")
imageVec(sim3Scalpel$Azero[,423],100,col = gray.colors(2))
title(main="Preliminary Element 423")
sim3Scalpel$AzeroThreshold[c(418,420,422,423)] #0.0388836
sim3Scalpel$thresholdVec  #0.03888360 0.05781994 0.07675627
mysim3_step1$thresholdval #0.04620095 0.06062502 0.07504909
#cluster 5 is a mis-identified neuron. track its corresponding preliminary elements, and they are 418 420 422 423.
#  all of these preliminary elements are from the minimum threshold 0.0388836, which means the threshod is not large enough.
#  it may be also the reason why the results of my function is better since minimum threshold used for preliminary segmentation is 0.04578128.

#Possible reasons: (1)spatial contamination of the background enlarges the overlapping area of neuron 10 and neuron 20;
#                  (2)the thresholds are so small that the noise pixels fail to be thresholded and consequentily mislead identify connected preliminary elements.
#                  (3)4-connectivity preliminary segmentation is not effective enough
#                  (4)filter the clusters by using larger cluster size.

#check if the background contamination leads to spatial overlapping.
Ydenoise=Ynew-wn
labelVec=c("","","","B3","B7","B4","BG")
partNB=cbind(Sim_neuronalInfo2$spatialmatrix[,c(1,10,20)],Sim_backgroundInfo2$backgroundimage[,c(3,7,4)],teststaticback$staticbackgroundmat[,1])
partNMtemp=rbind(Sim_neuronalInfo2$temporalmatrix[c(1,10,20),],Sim_backgroundInfo2$backgroundfluc[c(3,7,4),])
par(mfrow=c(1,2))
spatialplot(partNB,partNMtemp,title="Neuron & Background Footprints",labelVec = labelVec,border = TRUE,addToPlot =FALSE,pctTransp=0.5)
temporalplot(partNMtemp,title="Neuron & Background Temporal Signals",labelVec = c("N1","N10","N20","B3","B7","B4"))
#Therefore, in this case the local background did not cause overlapping of neuron 10 and neuron 20. 
#Even though N1, N10 and N20 overlapped but they did not spike simultaneourly. 
#Therefore, the fake identification may be caused by insufficient preprocessing.

labelback=paste0("B",1:ncol(Sim_backgroundInfo2$backgroundimage))
spatialplot(Sim_backgroundInfo2$backgroundimage,Sim_backgroundInfo2$backgroundfluc,title="Local Background Footprints",
            labelVec=labelback,border=FALSE,addToPlot = FALSE)
#Comparing the local background locations to the estimated neuronal footprints, 
#  it seems that the estimated neuron location is not accurate and did not remove background corruption.

sim3ScalpelAdj=scalpelStep3(getScalpelStep2("Desktop/simulation3/scalpelSim/"),alpha=0.1)
sum(apply(sim3ScalpelAdj$Zhat,1,max)!=0)
plotResults(sim3ScalpelAdj,titleA = "Estimated Footprints",titleZ = "Estimated Temporal Signals")
sim3ScalpelAdj$lambda #0.10283; lambda*alpha=0.010283(leading to noisy fluctuation);lambda*(1-alpha)=0.092547
sim3Scalpel$lambda #0.02364106; lambda*alpha=0.02127695
#When we use alpha=0.1, the fake neuron can be zeroed out. 
#However, the final tuned soft-threshold lambda*alpha becomes smaller, which means the estimated temporal traces are noisy.

#========================================================================================
#========================= Simulation 4 =================================================
#========================================================================================
#This simulation use the same data but removing the background layer
whitenoise=readRDS("Desktop/Sim_whitenoiseInfo2.rds")
Ysim4=Sim_neuronalInfo2$videomatrix-min(whitenoise)+whitenoise
saveRDS(Ysim4,file = "Y_1.rds")

sim4Scalpel=scalpel("Desktop/simulation4/scalpelSim/","Desktop/simulation4/",videoHeight = 100)
nrow(sim4Scalpel$Zhat)
plotResults(sim4Scalpel,titleA = "Estimated Footprints",titleZ = "Estimated Temporal Signals")

mysim4_step0=step0(outputDir = "Desktop/simulation4/myfunSim/",rawdataDir = "Desktop/simulation4/",videoheight=100)
mysim4_step1=step1(mysim4_step0, minSize=25,maxSize=500,maxWidth=30,maxHeight=30,thresholdval=NULL)
mysim4_step2=step2(step1output = mysim4_step1) #0.0459167688720728, 0.0606326579149862, 0.0753485469578996
mysim4_step3=step3(mysim4_step2)
nrow(mysim4_step3$Zhat)
par(mfrow=c(1,2))
spatialplot(mysim4_step3$Afilter,mysim4_step3$Zhat,title="Estimated Footprints",border = FALSE,addToPlot = FALSE)
temporalplot(mysim4_step3$Zhat,title="Estimated Temporal Signals")
sim4Scalpel$thresholdVec #0.03866044 0.05781640 0.07697236
ncol(sim3Scalpel$Azero)#2767
ncol(sim4Scalpel$Azero)#497
ncol(mysim3_step1$Azero)#2925
ncol(mysim4_step1$Azero)#565
#My preprocessing method will identify more preliminary elements. Next compare preprocessed data:
Ydeltaf4=readRDS("Desktop/simulation4/scalpelSim/Step0Data/Ydeltaf_part1.rds")
par(mfrow=c(1,3))
imageVec(Ydeltaf4[,420],100,col=grDevices::grey(seq(0, 1, length = 256)))
title(main="Scalpel Preprocessed Frame")
imageVec(mysim4_step0$deltafY[,420],100,col=grDevices::grey(seq(0, 1, length = 256)))
title(main="'blur' & 'ksmooth' Preprocessed Frame")
imageVec(Sim_neuronalInfo2$videomatrix[,420],100,col=grDevices::grey(seq(0, 1, length = 256)))
title(main="True Neuron Signal")

###Question: how to compare my temporal estimates with scalpel temporal estimates since their order do not match.
###Question: (1) Add the global background which can be constructed similar to the true data video
###Question: (2) change the spiking shape with longer decaying time (more overlapping spikes of different neurons)
#            and increase spatial density (large population and more overlapping)
#            this mainly used to check the effectiveness of step 3 and try to come up with other reasons for the failure to deal with overlapping cases
###Question: (3)Think about the structure of local background espically its temporal structure
#            may need reference!!
#            Check how does pengcheng generate the simulated data
#========================================================================================
#========================= Simulation 5 =================================================
#========================================================================================
#This experiment used the same true nueronal signal as above simulations,
#   but both the spatially-temporally correalted noise and independent noise are generated by using scalpel method.
#Simulation 5: SNR of correalted noise is 3, and for independent noise SNR=5
#              20 spatially-temporally correlated noise/background components; each component spikes once with width=75
#              random/independent noise is sample form unif(-1,1)
#              both the correlated noise/background and independent noise have range [-1,1]; range of neuronal signal is [0,1]
#              larger SNR means more significant neuronal signals
sim5scalpel=scalpel("Desktop/simulation5/scalpelSim5/","Desktop/simulation5/",videoHeight = 100)
sim5scalpel=getScalpelStep3("Desktop/simulation5/scalpelSim5/",alpha = 0.9,version = 38512)
plotResults(getScalpelStep3("Desktop/simulation5/scalpelSim5/",alpha = 0.9,version = 38512))
nrow(sim5scalpel$Zhat)
ncol(sim5scalpel$Azero) #3173
sum(apply(sim5scalpel$Zhat,1,max)!=0)
sim5scalpel$thresholdVec #0.1777779 0.2171209 0.2564639
par(mfrow=c(1,2))
spatialplot(Sim_neuronalInfo2$spatialmatrix,Sim_neuronalInfo2$temporalmatrix,title="Neuron Footprints",border = FALSE,addToPlot = FALSE)
spatialplot(sim5scalpel$Afilter,sim5scalpel$Zhat,title="Estimated Footprints",border = FALSE,addToPlot = FALSE)
temporalplot(Sim_neuronalInfo2$temporalmatrix,title="Neuron Temporal Signals")
par(mfrow=c(1,1))
temporalplot(sim5scalpel$Zhat,title="Estimated Temporal Signals")

sim5scalpelalpha=scalpelStep3(getScalpelStep2("Desktop/simulation5/scalpelSim5/",version=38512),alpha=0.1)
sum(apply(sim5scalpelalpha$Zhat,1,max)!=0)
par(mfrow=c(2,2))
spatialplot(Sim_neuronalInfo2$spatialmatrix,Sim_neuronalInfo2$temporalmatrix,title="Neuron Footprints",border = FALSE,addToPlot = FALSE)
spatialplot(sim5scalpelalpha$Afilter,sim5scalpelalpha$Zhat,title="Estimated Footprints (alpha=0.1)",border = FALSE,addToPlot = FALSE)
temporalplot(Sim_neuronalInfo2$temporalmatrix,title="Neuron Temporal Signals")
temporalplot(sim5scalpelalpha$Zhat,title="Estimated Temporal Signals (alpha=0.1)")

Ydeltaf5=readRDS("Desktop/simulation5/scalpelSim5/Step0Data/Ydeltaf_part1.rds")
myY5=readRDS("Desktop/simulation5/myfunSim5/Step0Data/deltafY.rds")
diff(range(Ydeltaf5)) #1.291949
diff(range(myY5)) #1.543164

#mysim5_step0=step0(outputDir = "Desktop/simulation5/myfunSim5/",rawdataDir = "Desktop/simulation5/",videoheight = 100)
#mysim5_step1=step1(mysim5_step0)
#mysim5_step2=step2(mysim5_step1)
#mysim5_step3=step3(mysim5_step2)
#par(mfrow=c(1,2))
#spatialplot(mysim5_step3$Afilter,mysim5_step3$Zhat,title="Estimated Footprints",border = FALSE,addToPlot = FALSE)
#temporalplot(mysim5_step3$Zhat,title="Estimated Temporal Signals")
#Ysim5=readRDS("Desktop/simulation5/Y_1.rds")
#imageVec(Ysim5[,420],100)
#imageVec(Ydeltaf5[,420],100,col=grDevices::grey(seq(0, 1, length = 256)))

png("Desktop/simulation5/Data/preprocessedframe%04d.png")
for (i in 1:1000) {
  imageVec(Ydeltaf5[,i],100,col=grDevices::grey(seq(0, 1, length = 256)))
}
dev.off()
#ffmpeg -r 30 -i neuronsig%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p neuronfilm2.mp4
#Seriously fail to remove the background/correlated noise!!
#Based on the preprocessed film, those redundant fake neurons mainly caused by insufficient preprocessing.
#  To be more specific, the neurons are significant enough due to the high SNR, 
#  which means the thresholds used in preliminary segmentation are too small, 
#  and consequently spatially correlated noise/background are mis-identified.

#Try larger thresholdVec
sim5scalpelAdj=scalpelStep1(getScalpelStep0("Desktop/simulation5/scalpelSim5/"),thresholdVec = c(0.2564639,0.5,0.7))
ncol(sim5scalpelAdj$Azero) #1124
ncol(sim5scalpel$Azero) #3173
sim5scalpelAdj2=scalpelStep2(getScalpelStep1("Desktop/simulation5/scalpelSim5/",version=91237))
ncol(sim5scalpelAdj2$A)
sim5scalpelAdj3=scalpelStep3(getScalpelStep2("Desktop/simulation5/scalpelSim5/",version=91237))
sim5scalpelAdj3=getScalpelStep3("Desktop/simulation5/scalpelSim5/",version=91237)
nrow(sim5scalpelAdj3$Zhat) #21
par(mfrow=c(2,2))
spatialplot(Sim_neuronalInfo2$spatialmatrix,Sim_neuronalInfo2$temporalmatrix,title="Neuron Footprints",border = FALSE,addToPlot = FALSE)
spatialplot(sim5scalpelAdj3$Afilter,sim5scalpelAdj3$Zhat,title=paste0(c("Estimated Footprints\n","(larger threshold)")),border = FALSE,addToPlot = FALSE)
temporalplot(Sim_neuronalInfo2$temporalmatrix,title="Neuron Temporal Signals")
temporalplot(sim5scalpelAdj3$Zhat,title=paste0(c("Estimated Temporal Signals\n","(larger threshold)")))

which(sim5scalpelAdj3$clusterID==18) #378 380 382 384 385 387 389 391 393 395 397 399 401 403 405
unique(sim5scalpelAdj3$AzeroThreshold[which(sim5scalpelAdj3$clusterID==18)]) #0.2564639
#The result is much better! There are only one fake neuron identified compared to 40 identified neurons by default thresholds
#The fake neuron corresponds to the lower threshold and is actaully a spataially-temporally correlated noise/background component. 
#  since correlated SNR is larger than independent SNR, the thresholds are still not large enough to filter out some correlated background 
#  there are 20 correlated background components and each component spikes at one time with sin curve
sim5scalpeladjalpha=scalpelStep3(getScalpelStep2("Desktop/simulation5/scalpelSim5/",version=91237),alpha=0.1)
plotResults(sim5scalpelalpha,titleA = "Estimated Footprints (alpha=0.1)",titleZ = "Estimated Temporal Signals (alpha=0.1)")
sum(apply(sim5scalpeladjalpha$Zhat,1,sum)!=0)
#Since the fake neuron is not due to overlapping, tuning alpha may not be helpful.

#========================================================================================
#========================= Simulation 6 =================================================
#========================================================================================
#This experiment used the same true nueronal signal as above simulations,
#   but both the spatially-temporally correalted noise and independent noise are generated by using scalpel method.
#Simulation 6: SNR of correalted noise is 1, and for independent noise SNR=2
#              20 spatially-temporally correlated noise/background components; each component spikes once
#              random/independent noise is sample form unif(-1,1)
#              both the correlated noise/background and independent noise have range [-1,1]; range of neuronal signal is [0,1]
#              larger SNR means more significant neuronal signals
sim6scalpel=scalpel("Desktop/simulation6/scalpelSim6/","Desktop/simulation6/",videoHeight = 100)
sim6scalpel=getScalpelStep3("Desktop/simulation6/scalpelSim6/")
plotResults(sim6scalpel,titleA = "Estimated Footprints",titleZ = "Estimated Temporal Signals")
par(mfrow=c(1,2))
spatialplot(Sim_neuronalInfo2$spatialmatrix,Sim_neuronalInfo2$temporalmatrix,title="Neuron Footprints",border = FALSE,addToPlot = FALSE)
spatialplot(sim6scalpel$Afilter,Sim_neuronalInfo2$temporalmatrix,title="Estimated Footprints",border = FALSE,addToPlot = FALSE)
nrow(sim6scalpel$Zhat)
sum(apply(sim6scalpel$Zhat,1,max)!=0)
sim6scalpel$thresholdVec #0.2083200 0.2487552 0.2891904
table(sim6scalpel$clusterID)
names(sort(table(sim6scalpel$clusterID))) #all mis-identified clusters have miderate size which can not be filtered by setting larger value of minSize.
#This time there are 6 neurons missed but 4 redundant fake neurons are mid-identified (3,5,8,18).
#   which implies that the preprocessing step can not separate the correalted noise/background components from neuronal signals
#   and consequently the thresholds failed to filter away those noise pixels but retain neuronal pixels
#   scalpel may be restriced to a narrow range of SNR, higher(3,5) or lower SNR(1,2) will not be good.
mysim6_step0=step0(outputDir = "Desktop/simulation6/myfunSim6/",rawdataDir = "Desktop/simulation6/",videoheight = 100)
mysim6_step1=step1(mysim6_step0) #threshold: 0.217540621909628, 0.267850521415998, 0.318160420922369
ncol(mysim6_step1$Azero) #266; ncol(sim6scalple$Azero): 162
mysim6_step2=step2(mysim6_step1)
ncol(mysim6_step2$A) #22; ncol(sim6scalpel$A): 18
mysim6_step3=step3(mysim6_step2)
sum(apply(mysim6_step3$Zhat,1, sum)!=0) #22
par(mfrow=c(1,3))
spatialplot(Sim_neuronalInfo2$spatialmatrix,Sim_neuronalInfo2$temporalmatrix,title="Neuron Footprints",border = FALSE,addToPlot = FALSE)
spatialplot(sim6scalpel$Afilter,Sim_neuronalInfo2$temporalmatrix,title="Estimated Footprints (scalpel)",border = FALSE,addToPlot = FALSE)
spatialplot(mysim6_step3$Afilter,mysim6_step3$Zhat,title="Estimated Footprints (myfun)",border = FALSE,addToPlot = FALSE)
labelVec=paste0("N",c(1:5,7,9:14,6,15:18))
temporalplot(Sim_neuronalInfo2$temporalmatrix[c(1:5,7,9:14,6,15:18),],title="Neuron Temporal Signals",labelVec = labelVec)
labelVec=paste0("S",c(1,14,2,10,12,16,4,7,17,9,11,0,13,13,15,6,0))
temporalplot(rbind(sim6scalpel$Zhat[c(1,14,2,10,12,16,4,7,17,9,11),],rep(0,1000),sim6scalpel$Zhat[c(13,13,15,6),],rep(0,1000)),
             title="Estimated Temporal Signals (scalpel)",labelVec = labelVec)
labelVec=paste0("M",c(1,19,3,1,7,20,8,12,18,14,17,10,16,16,15,11,6))
temporalplot(mysim6_step3$Zhat[c(1,19,3,1,7,20,8,12,18,14,17,10,16,16,15,11,6),],title="Estimated Temporal Signals (myfun)",
             labelVec = labelVec)
which(sim6scalpel$clusterID==13)
par(mfrow=c(1,1))
labelVec=c("N6","N15",paste0("F",91:94))
spatialplot(cbind(Sim_neuronalInfo2$spatialmatrix[,c(6,15)],sim6scalpel$Azero[,91],sim6scalpel$Azero[,92],sim6scalpel$Azero[,93],sim6scalpel$Azero[,94]),
            Sim_neuronalInfo2$temporalmatrix,labelVec = labelVec,pctTransp=0.6)
#temporalplot(rbind(sim6scalpel$Zhat[13,],mysim6_step3$Zhat[16,]))
#800+which.max(sim6scalpel$Zhat[13,800:1000]>0) #943
#800+which.max(mysim6_step3$Zhat[16,800:1000]>0) #932
par(mfrow=c(1,2))
imageVec(readRDS("Desktop/simulation6/myfunSim6/Step0Data/deltafY.rds")[,932],100,col=gray.colors(256))
imageVec(readRDS("Desktop/simulation6/scalpelSim6/Step0Data/Ydeltaf_part1.rds")[,943],100,col=gray.colors(256))
sim6scalpel$lambda #0.1882841 
mysim6_step3$lambdaFinal #0.1263048
#Conslusion: (1) With high spatially-temporally correalted random noise/background, 
#                both scalpel and my function will miss some neurons and identify some fake neurons
#            (2) Scalpel missed more neurons which are 6,8,14,18,19,20, and for those ture poisitive, the corresponding temporal estimates are less accurate
#                due to the large tuned lambda (0.1882841), and consequently some spikes are encouraged to 0 and no decaying fluctuations)
#            (3) My function identify 22 neurons in total with only missing 6,8,19,20, but there are more fake identifications.
#                Besides, the temporal estimates are more accurate than scalpel due to small tuned lambda (0.1263048)
#            (4) Both scalpel and my functions seem to merge the 6 and 15 neurons together, which means the temporal estimate of neuron 15 is corrupted by neuron 6,
#                this may be due to the spatial overlap between 6 and 15. Since the neuron 6 is missed out of the dictionary,
#                the regression step needs to attribute the spike of neuron 6 to 15 to reduce bias.
Ysim6=readRDS("Desktop/simulation6/Y_1.rds")
png("Desktop/simulation6/Data/rawframe%04d.png")
for (i in 1:1000) {
  imageVec(Ysim6[,i],100,col=grDevices::grey(seq(0, 1, length = 256)))
}
dev.off()
#ffmpeg -r 30 -i neuronsig%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p neuronfilm2.mp4
Ydeltaf6=readRDS("Desktop/simulation6/scalpelSim6/Step0Data/Ydeltaf_part1.rds")
png("Desktop/simulation6/Data/preprocessedframe%04d.png")
for (i in 1:1000) {
  imageVec(Ydeltaf6[,i],100,col=grDevices::grey(seq(0, 1, length = 256)))
}
dev.off()
myY6=readRDS("Desktop/simulation6/myfunSim6/Step0Data/deltafY.rds")
png("Desktop/simulation6/Data/mypreprocessframe%04d.png")
for (i in 1:1000) {
  imageVec(myY6[,i],100,col=grDevices::grey(seq(0, 1, length = 256)))
}
dev.off()
diff(range(Ydeltaf6)) #0.7101225
diff(range(myY6)) #0.812219
#========================================================================================
#========================= Simulation 7 =================================================
#========================================================================================
#This experiment used the same true nueronal signal as above simulations,
#   but both the spatially-temporally correalted noise and independent noise are generated by using scalpel method.
#Simulation 7: SNR of correalted noise is 1.5, and for independent noise SNR=1
#              20 spatially-temporally correlated noise/background components; each component spikes once
#              random/independent noise is sample form unif(-1,1)
#              both the correlated noise/background and independent noise have range [-1,1]; range of neuronal signal is [0,1]
#              larger SNR means more significant neuronal signals
sim7scalpel=scalpel("Desktop/simulation7/scalpelSim7/","Desktop/simulation7/",videoHeight = 100)
plotResults(sim7scalpel,titleA = "Estimated Footprints",titleZ = "Estimated Temporal Signals")
ncol(sim7scalpel$Azero) #706
sum(apply(sim7scalpel$Zhat,1,max)!=0) #35
par(mfrow=c(1,2))
spatialplot(Sim_neuronalInfo2$spatialmatrix,Sim_neuronalInfo2$temporalmatrix,title="Neuron Footprints",border = FALSE,addToPlot = FALSE)
spatialplot(sim7scalpel$Afilter,Sim_neuronalInfo2$temporalmatrix,title="Estimated Footprints",border = FALSE,addToPlot = FALSE)
temporalplot(Sim_neuronalInfo2$temporalmatrix,title="Neuron Temporal Signals")
temporalplot(sim7scalpel$Zhat[c(2,24,5,1,9,7,13,19,15,25,30,27,21,20,8,22,10,6,34,4),],title = "Corresponding Estimates",labelVec = c(2,24,5,1,9,7,13,19,15,25,30,27,21,20,8,22,10,6,34,4))
temporalplot(sim7scalpel$Zhat,title="Estimated Temporal Signals")
sim7scalpel$thresholdVec #0.1249989 0.1693236 0.2136483
Ythreshold7=readRDS("Desktop/simulation7/scalpelSim7/Step0Data/Ydeltaf_part1.rds")>0.1249989
png("Desktop/simulation7/Data/thresholdframe%04d.png")
for (i in 1:1000) {
  imageVec(Ythreshold7[,i],100,col=grDevices::grey(seq(0, 1, length = 256)))
}
dev.off()
#ffmpeg -r 30 -i neuronsig%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p neuronfilm2.mp4

sim7scalpelAdj=scalpelStep3(getScalpelStep2("Desktop/simulation7/scalpelSim7/"),alpha=0.2)
sum(apply(sim7scalpelAdj$Zhat,1,max)!=0)
par(mfrow=c(1,2))
spatialplot(Sim_neuronalInfo2$spatialmatrix,Sim_neuronalInfo2$temporalmatrix,title="Neuron Footprints",border = FALSE,addToPlot = FALSE)
spatialplot(sim7scalpelAdj$Afilter,sim7scalpelAdj$Zhat,title="Estimated Footprints (alpha=0.2)",border = FALSE,addToPlot = FALSE)
temporalplot(Sim_neuronalInfo2$temporalmatrix,title="Neuron Temporal Signals")
temporalplot(sim7scalpelAdj$Zhat[c(2,24,5,1,9,7,13,19,15,25,30,27,21,20,8,22,10,6,34,4),],
             title="Corresponding Estimates (alpha=0.2)",labelVec = c(2,24,5,1,9,7,13,19,15,25,30,27,21,20,8,22,10,6,34,4))
sim7scalpel$lambda #0.1088601
sim7scalpelAdj$lambda #0.4506187(alpha=0.1) noisy fluctuations 0.3405543(alpha=0.2)
#========================================================================================
#========================= Simulation 8 =================================================
#========================================================================================
#This experiment used the same true nueronal signal as above simulations,
#   but both the spatially-temporally correalted noise and independent noise are generated by using scalpel method.
#Simulation 8: SNR of correalted noise is 1.5, and for independent noise SNR=0.5
#              20 spatially-temporally correlated noise/background components; each component spikes once
#              random/independent noise is sample form unif(-1,1)
#              both the correlated noise/background and independent noise have range [-1,1]; range of neuronal signal is [0,1]
#              larger SNR means more significant neuronal signals
sim8scalpel=scalpel("Desktop/simulation8/scalpelSim8/","Desktop/simulation8/",videoHeight = 100)
plotResults(sim8scalpel,titleA = "Estimated Footprints",titleZ = "Estimated Temporal Signals")
ncol(sim8scalpel$Azero) #180
ncol(sim8scalpel$Afilter) #32
sum(apply(sim8scalpel$Zhat,1,max)!=0)#32
par(mfrow=c(2,2))
spatialplot(Sim_neuronalInfo2$spatialmatrix,Sim_neuronalInfo2$temporalmatrix,title="Neuron Footprints",border = FALSE,addToPlot = FALSE)
spatialplot(sim8scalpel$Afilter,Sim_neuronalInfo2$temporalmatrix,title="Estimated Footprints",border = FALSE,addToPlot = FALSE)
temporalplot(Sim_neuronalInfo2$temporalmatrix,title="Neuron Temporal Signals")
temporalplot(sim8scalpel$Zhat[c(1,28,5,24,11,9,29,12,14,3,7,22,25,18,10,23,21,8,27,2),],title="Corresponding Estimates",labelVec = c(1,28,5,24,11,9,29,12,14,3,7,22,25,18,10,23,21,8,27,2))
sim8scalpel$thresholdVec #0.1027832 0.1432586 0.1837339
sim8scalpelAdj=scalpelStep3(getScalpelStep2("Desktop/simulation8/scalpelSim8/"),alpha=0.1)
sum(apply(sim8scalpelAdj$Zhat,1,max)!=0)
par(mfrow=c(2,2))
spatialplot(Sim_neuronalInfo2$spatialmatrix,Sim_neuronalInfo2$temporalmatrix,title="Neuron Footprints",border = FALSE,addToPlot = FALSE)
spatialplot(sim8scalpelAdj$Afilter,sim8scalpelAdj$Zhat,title="Estimated Footprints (alpha=0.1)",border = FALSE,addToPlot = FALSE)
temporalplot(Sim_neuronalInfo2$temporalmatrix,title="Neuron Temporal Signals")
temporalplot(sim8scalpelAdj$Zhat[c(1,28,5,24,11,9,29,12,14,3,7,22,25,18,10,23,21,8,27,2),],title="Corresponding Estimates (alpha=0.1)",labelVec = c(1,28,5,24,11,9,29,12,14,3,7,22,25,18,10,23,21,8,27,2))
Ythreshold8=readRDS("Desktop/simulation8/scalpelSim8/Step0Data/Ydeltaf_part1.rds")>0.1027832
png("Desktop/simulation8/Data/thresholdframe%04d.png")
for (i in 1:1000) {
  imageVec(Ythreshold8[,i],100,col=grDevices::grey(seq(0, 1, length = 256)))
}
dev.off()
#ffmpeg -r 30 -i neuronsig%04d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p neuronfilm2.mp4

#Conslusion: (1) Since the correlated background is moderate, all neurons are detected but the corresponding temporal estimates are not accurate:
#                the estimate 7,5,18 corresponding to neuron 11,3,14 miss some spike activities
#            (2) Compared to simulation 7 (SNR=1.5,1), higher independent noise (SNR=0.5) renders a small number of preliminary elements (180, compaered to 720 in simulation 7).
#                this may be due to the range of preprocessed data with higher independent noise is narrower than simulation 7.
#            (3) Consequently, the final dictionary of simulation 8 is smaller than simulation 7 (32 and 35 respectively)
#            (4) However, due to the high noise level in simulation 8, the preprocessed step will weaken the neuronal signal, 
#                 which may explain, to some degree, why simulation 8 zeroed out a true neuron by using alpha=0.1 
#                 but simulation 7 successfully zeroed out fake neurons while keep all true neurons.
#            (5) Compared to simulation 5 where both the correlated background and independent noise are low (SNR=3,5), there are much more fake neurons (40).
#                changing alpha to a small value seems not help a lot. In particular, the preliminary dictionary contains 3173 elements,
#                which is due to the inappropriate selection of thresholds in preliminary segmentation step.
#                Using a self-defined larger threshold sequence, we finally identified 21 neurons, and by tracking that fake neuron, 
#                it turns out that this fake cluster contains 15 preliminary elements and all come from the smallest value in the threshold sequence.
#                Also, this fake neuron cannot be zeroed out by changing the value of alpha.
#                In fact, we can discard this cluster before regression by seting the minSize of cluster to be 15, 
#                since this is indeed the smallest cluster, but this is not easy to tune without enough prior knowledge.
#            (6) Compared to simulation 6 (SNR=1,2) where the correlated background is strong, scalpel will miss some true neurons, 
#                which may be due to preprocessing cannot seperate the background and neuronal signal, 
#                and consequently, some neuron pixels are filtered out in preliminary segmentation. 
#                Besides, those fake clusters have moderate size which indicates that they cannot be filtered out by setting larger minSize.
#            (7) One of my suppositions is that if the number of neuron components is much larger than the background components, 
#                the thresholds used for preliminary segmentation may more generally work. Otherwise, the effectiveness may depend on the balance 
#                between the strength of neuronal signal, correlated noise and independent noise. (high independent noise and moderate correlated noise)
#            (8) The fake neurons are more likely from the correlated background components, which implies that the preprocessing step cannot
#                effectively reduce the correlated background. And in both simulation 7 and 8, using a smaller alpha (ie. larger soft-scale)
#                is quite helpful to zero out those fake identification, but this may due to the correlated background components we simulated only spike once.
#                Next, we also need to try whether tuning alpha is still useful when the correlated background spikes several times.
#            (9) Next, we can evaluate scalpel without correlated background contamination, and check its performance on overlapping cases.
#                Especially focus on the effectiveness of group lasso to zero out elements in Afilter.
#========================================================================================
#========================= Simulation 9 =================================================
#========================================================================================
sim9A=readRDS("Desktop/simulation9/data/sim9trueA.rds")
sim9Z=readRDS("Desktop/simulation9/data/sim9trueZ.rds")
Ysim9=sim9A%*%sim9Z+
  readRDS("Desktop/simulation9/data/corNoiseSNR_3.rds")+readRDS("Desktop/simulation9/data/indNoiseSNR_4.rds")
saveRDS(Ysim9,file = "Desktop/simulation9/data/Y_1.rds")
sim9Scalpel=scalpel("Desktop/simulation9/","Desktop/simulation9/data/",videoHeight = 200)
#149
par(mfrow=c(1,2))
spatialplot(sim9A,sim9Z,title="Neuron Footprints",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim9Scalpel$Afilter,sim9Scalpel$Zhat,title="Estimated Footprints (CSNR=3,ISNR=4)",border = TRUE,addToPlot = FALSE,videoheight = 200)
par(mfrow=c(1,1))
imageVec(Ysim9[,56],videoheight = 200)
apply(sim9Z,1,which.max)

#========================================================================================
#========================= Simulation 10 =================================================
#========================================================================================
Ysim10=sim9A%*%sim9Z+
  readRDS("Desktop/simulation9/data/corNoiseSNR_1.5.rds")+readRDS("Desktop/simulation9/data/indNoiseSNR_1.rds")
saveRDS(Ysim10,file = "Desktop/simulation10/data/Y_1.rds")
sim10Scalpel=scalpel("Desktop/simulation10/","Desktop/simulation10/data/",videoHeight = 200,minClusterSize = 5)
ncol(sim10Scalpel$Afilter) #113
#191
par(mfrow=c(1,2))
spatialplot(sim9A,sim9Z,title="Neuron Footprints",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim10Scalpel$Afilter,sim10Scalpel$Zhat,title="Estimated Footprints (CSNR=1.5,ISNR=1)",border = TRUE,addToPlot = FALSE,videoheight = 200)
par(mfrow=c(1,1))
imageVec(Ysim10[,56],videoheight = 200)

#========================================================================================
#========================= Simulation 11 =================================================
#========================================================================================
Ysim11=sim9A%*%sim9Z+
  readRDS("Desktop/simulation9/data/corNoiseSNR_1.5.rds")+readRDS("Desktop/simulation9/data/indNoiseSNR_1.5.rds")
saveRDS(Ysim11,file = "Desktop/simulation11/data/Y_1.rds")
sim11Scalpel=scalpel("Desktop/simulation11/","Desktop/simulation11/data/",videoHeight = 200)
nrow(sim11Scalpel$Zhat)
#173
par(mfrow=c(1,2))
spatialplot(sim9A,sim9Z,title="Neuron Footprints",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim11Scalpel$Afilter,sim11Scalpel$Zhat,title="Estimated Footprints (CSNR=1.5,ISNR=1.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)
par(mfrow=c(1,3))
imageVec(Ysim11[,56],videoheight = 200)
title(main="CSNR=1.5,ISNR=1.5,estn=173")
#========================================================================================
#========================= Simulation 12 =================================================
#========================================================================================
Ysim12=sim9A%*%sim9Z+
  readRDS("Desktop/simulation9/data/corNoiseSNR_1.5.rds")+readRDS("Desktop/simulation9/data/indNoiseSNR_2.rds")
saveRDS(Ysim12,file = "Desktop/simulation12/data/Y_1.rds")
sim12Scalpel=scalpel("Desktop/simulation12/","Desktop/simulation12/data/",videoHeight = 200)
nrow(sim12Scalpel$Zhat) #161
par(mfrow=c(1,2))
spatialplot(sim9A,sim9Z,title="Neuron Footprints",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim12Scalpel$Afilter,sim12Scalpel$Zhat,title="Estimated Footprints (CSNR=1.5,ISNR=2)",border = TRUE,addToPlot = FALSE,videoheight = 200)
par(mfrow=c(1,1))
imageVec(Ysim12[,56],videoheight = 200)
#========================================================================================
#========================= Simulation 13 =================================================
#========================================================================================
Ysim13=sim9A%*%sim9Z+
  readRDS("Desktop/simulation9/data/corNoiseSNR_1.rds")+readRDS("Desktop/simulation9/data/indNoiseSNR_1.5.rds")
saveRDS(Ysim13,file = "Desktop/simulation13/data/Y_1.rds")
sim13Scalpel=scalpel("Desktop/simulation13/","Desktop/simulation13/data/",videoHeight = 200, minClusterSize = 5)
nrow(sim13Scalpel$Zhat) #145
par(mfrow=c(1,2))
spatialplot(sim9A,sim9Z,title="Neuron Footprints",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim13Scalpel$Afilter,sim13Scalpel$Zhat,title="Estimated Footprints (CSNR=1,ISNR=1.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)
par(mfrow=c(1,1))
imageVec(Ysim13[,56],videoheight = 200)
title(main = "CSNR=1,ISNR=1.5,estn=145")
#========================================================================================
#========================= Simulation 14 =================================================
#========================================================================================
Ysim14=sim9A%*%sim9Z+
  readRDS("Desktop/simulation9/data/corNoiseSNR_2.rds")+readRDS("Desktop/simulation9/data/indNoiseSNR_1.5.rds")
saveRDS(Ysim14,file = "Desktop/simulation14/data/Y_1.rds")
sim14Scalpel=scalpel("Desktop/simulation14/","Desktop/simulation14/data/",videoHeight = 200)
nrow(sim14Scalpel$Zhat) #214
par(mfrow=c(1,2))
spatialplot(sim9A,sim9Z,title="Neuron Footprints",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim14Scalpel$Afilter,sim14Scalpel$Zhat,title="Estimated Footprints (CSNR=2,ISNR=1.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)
par(mfrow=c(1,1))
imageVec(Ysim14[,56],videoheight = 200)
title(main = "CSNR=2,ISNR=1.5,estn=214")
#========================================================================================
#========================= Simulation 15 =================================================
#========================================================================================
Ysim15=sim9A%*%sim9Z+
  readRDS("Desktop/simulation9/data/corNoiseSNR_1.5.rds")+readRDS("Desktop/simulation9/data/indNoiseSNR_0.5.rds")
saveRDS(Ysim15,file = "Desktop/simulation15/data/Y_1.rds")
sim15Scalpel=scalpel("Desktop/simulation15/","Desktop/simulation15/data/",videoHeight = 200)
nrow(sim15Scalpel$Zhat) #163
par(mfrow=c(1,2))
spatialplot(sim9A,sim9Z,title="Neuron Footprints",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim15Scalpel$Afilter,sim15Scalpel$Zhat,title="Estimated Footprints (CSNR=1.5,ISNR=0.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)
par(mfrow=c(1,1))
imageVec(Ysim15[,56],videoheight = 200)
title(main = "CSNR=1.5,ISNR=0.5,estn=163")
#========================================================================================
#========================= evaluate (ISNR=1.5) ==========================================
#========================================================================================
eval13=evaluateSims(sim9A,estA = sim13Scalpel$Afilter,estZ = sim13Scalpel$Zhat,
                    cutoff_extra_pixels=0.2,cutoff_matching_pixels=0.5)
eval13$ID.sens #0.91
eval13$ID.prec #0.6275862
which(eval13$m==0) #7 17 28 34 52 71 74 84 91
par(mfrow=c(2,1))
temporalplot(sim9Z[1:3,],title="Temporal Truth")
temporalplot(sim13Scalpel$Zhat[eval13$m[which(eval13$m!=0)],],title="Temporal Estimates (CSNR=1,ISNR=1.5)")
eval13$m[8] #84
#the temporal estimates of 8,82 are not good
spatialplot(sim9A[,which(eval13$m!=0)],sim9Z,title="Spatial Truth",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim13Scalpel$Afilter[,eval13$m[which(eval13$m!=0)]],sim9Z,title="Spatial Estimates (CSNR=1,ISNR=1.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)

###(1.5,1.5)#####
eval11=evaluateSims(sim9A,estA = sim11Scalpel$Afilter,estZ = sim11Scalpel$Zhat,
                    cutoff_extra_pixels=0.2,cutoff_matching_pixels=0.5)
eval11$ID.sens #0.99
eval11$ID.prec #0.5722543
which(eval11$m==0) #23
par(mfrow=c(1,2))
temporalplot(sim9Z[which(eval11$m!=0),],title="Temporal Truth")
temporalplot(sim11Scalpel$Zhat[eval11$m[which(eval11$m!=0)],],title="Temporal Estimates (CSNR=1.5,ISNR=1.5)")
#the temporal estimates are much better
spatialplot(sim9A[,which(eval11$m!=0)],sim9Z,title="Spatial Truth",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim11Scalpel$Afilter[,eval11$m[which(eval11$m!=0)]],sim9Z,title="Spatial Estimates (CSNR=1.5,ISNR=1.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)
###(2,1.5)###
eval14=evaluateSims(sim9A,estA = sim14Scalpel$Afilter,estZ = sim14Scalpel$Zhat,
                    cutoff_extra_pixels=0.2,cutoff_matching_pixels=0.5)
eval14$ID.sens #0.99
eval14$ID.prec #0.4626168
which(eval14$m==0) #23
par(mfrow=c(1,2))
temporalplot(sim9Z[which(eval14$m!=0),],title="Temporal Truth")
temporalplot(sim14Scalpel$Zhat[eval14$m[which(eval14$m!=0)],],title="Temporal Estimates (CSNR=2,ISNR=1.5)")
#the temporal estimates are much better
spatialplot(sim9A[,which(eval14$m!=0)],sim9Z,title="Spatial Truth",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim14Scalpel$Afilter[,eval14$m[which(eval14$m!=0)]],sim9Z,title="Spatial Estimates (CSNR=2,ISNR=1.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)
#========================================================================================
#========================= evaluate (CSNR=1.5) ==========================================
#========================================================================================
###(1.5,1)#####
eval10=evaluateSims(sim9A,estA = Matrix(sim10Scalpel$Afilte),estZ = sim10Scalpel$Zhat,
                    cutoff_extra_pixels=0.2,cutoff_matching_pixels=0.5)
eval10$ID.sens #0.98
eval10$ID.prec #0.513089 #0.8672566
which(eval10$m==0) #23,59
par(mfrow=c(1,2))
temporalplot(sim9Z[which(eval10$m!=0),],title="Temporal Truth")
temporalplot(sim10Scalpel$Zhat[eval10$m[which(eval10$m!=0)],],title="Temporal Estimates (CSNR=1.5,ISNR=1)")
spatialplot(sim9A[,which(eval10$m!=0)],sim9Z,title="Spatial Truth",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim10Scalpel$Afilter[,eval10$m[which(eval10$m!=0)]],sim9Z,title="Spatial Estimates (CSNR=1.5,ISNR=1)",border = TRUE,addToPlot = FALSE,videoheight = 200)

###(1.5,2)#####
eval12=evaluateSims(sim9A,estA = sim12Scalpel$Afilter,estZ = sim12Scalpel$Zhat,
                    cutoff_extra_pixels=0.2,cutoff_matching_pixels=0.5)
eval12$ID.sens #0.98
eval12$ID.prec #0.6086957
which(eval12$m==0) #23,74
par(mfrow=c(1,2))
temporalplot(sim9Z[which(eval12$m!=0),],title="Temporal Truth")
temporalplot(sim12Scalpel$Zhat[eval12$m[which(eval12$m!=0)],],title="Temporal Estimates (CSNR=1.5,ISNR=2)")
spatialplot(sim9A[,which(eval12$m!=0)],sim9Z,title="Spatial Truth",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim12Scalpel$Afilter[,eval12$m[which(eval12$m!=0)]],sim9Z,title="Spatial Estimates (CSNR=1.5,ISNR=2)",border = TRUE,addToPlot = FALSE,videoheight = 200)

###(1.5,0.5)#####
eval15=evaluateSims(sim9A,estA = sim15Scalpel$Afilter,estZ = sim15Scalpel$Zhat,
                    cutoff_extra_pixels=0.2,cutoff_matching_pixels=0.5)
eval15$ID.sens #0.96
eval15$ID.prec #0.5889571
which(eval15$m==0) #19 50 74 94
par(mfrow=c(1,2))
temporalplot(sim9Z[which(eval15$m!=0),],title="Temporal Truth")
temporalplot(sim15Scalpel$Zhat[eval15$m[which(eval15$m!=0)],],title="Temporal Estimates (CSNR=1.5,ISNR=0.5)")
spatialplot(sim9A[,which(eval15$m!=0)],sim9Z,title="Spatial Truth",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim15Scalpel$Afilter[,eval15$m[which(eval15$m!=0)]],sim9Z,title="Spatial Estimates (CSNR=1.5,ISNR=0.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)

par(mfrow=c(2,4))
par(mar=c(10,2,2,2))
imageVec(sim9A%*%sim9Z[,56],videoheight = 200)
title(main = "True neuron in a frame")
imageVec(Ysim13[,56],videoheight = 200)
title(main = "CSNR=1,ISNR=1.5")
imageVec(Ysim11[,56],videoheight = 200)
title(main = "CSNR=1.5,ISNR=1.5")
imageVec(Ysim14[,56],videoheight = 200)
title(main = "CSNR=2,ISNR=1.5")
imageVec(Ysim15[,56],videoheight = 200)
title(main = "CSNR=1.5,ISNR=0.5")
imageVec(Ysim10[,56],videoheight = 200)
title(main = "CSNR=1.5,ISNR=1")
imageVec(Ysim11[,56],videoheight = 200)
title(main = "CSNR=1.5,ISNR=1.5")
imageVec(Ysim12[,56],videoheight = 200)
title(main = "CSNR=1.5,ISNR=2")
x=c(1,1.5,2)
par(mfrow=c(1,1))
par(mar=c(10,4,2,2))
plot(x,c(0.91,0.99,0.99),main = "ISNR=1.5",xlab = "CSNR",ylab = "sensitivity",type = "b",ylim = c(0,1))
plot(x,c(0.63,0.57,0.46),main = "ISNR=1.5",xlab = "CSNR",ylab = "precision",type = "b",ylim = c(0,1))
x=c(0.5,1,1.5,2)
plot(x,c(0.96,0.98,0.99,0.98),main = "CSNR=1.5",xlab = "ISNR",ylab = "sensitivity",type = "b",ylim = c(0,1))
plot(x,c(0.59,0.51,0.57,0.61),main = "CSNR=1.5",xlab = "ISNR",ylab = "precision",type = "b",ylim = c(0,1))


#========================================================================================
#========================= Adjust Simulation 13 =================================================
#========================================================================================
Ysim13=sim9A%*%sim9Z+
  readRDS("Desktop/simulation9/data/corNoiseSNR_1.rds")+readRDS("Desktop/simulation9/data/indNoiseSNR_1.5.rds")
saveRDS(Ysim13,file = "Desktop/simulation13/data/Y_1.rds")
sim13Scalpel=scalpel("Desktop/simulation13/","Desktop/simulation13/data/",videoHeight = 200, minClusterSize = 5)
nrow(sim13Scalpel$Zhat) #145 #116
par(mfrow=c(1,2))
spatialplot(sim9A,sim9Z,title="Neuron Footprints",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim13Scalpel$Afilter,sim13Scalpel$Zhat,title="Estimated Footprints (CSNR=1,ISNR=1.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)
par(mfrow=c(1,1))
imageVec(Ysim13[,56],videoheight = 200)
title(main = "CSNR=1,ISNR=1.5,estn=145")

eval13=evaluateSims(sim9A,estA = sim13Scalpel$Afilter,estZ = sim13Scalpel$Zhat,
                    cutoff_extra_pixels=0.2,cutoff_matching_pixels=0.5)
eval13$ID.sens #0.91 #0.89
eval13$ID.prec #0.6275862 #0.7672414
which(eval13$m==0) #7 17 28 34 52 71 74 84 91 #two more mis neurons: 48, 69
#par(mfrow=c(2,1))
#temporalplot(sim9Z[1:3,],title="Temporal Truth")
#temporalplot(sim13Scalpel$Zhat[eval13$m[which(eval13$m!=0)],],title="Temporal Estimates (CSNR=1,ISNR=1.5)")
#eval13$m[8] #84
#the temporal estimates of 8,82 are not good
#spatialplot(sim9A[,which(eval13$m!=0)],sim9Z,title="Spatial Truth",border = TRUE,addToPlot = FALSE,videoheight = 200)
#spatialplot(sim13Scalpel$Afilter[,eval13$m[which(eval13$m!=0)]],sim9Z,title="Spatial Estimates (CSNR=1,ISNR=1.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)

#========================================================================================
#========================= Adjust Simulation 11 =================================================
#========================================================================================
Ysim11=sim9A%*%sim9Z+
  readRDS("Desktop/simulation9/data/corNoiseSNR_1.5.rds")+readRDS("Desktop/simulation9/data/indNoiseSNR_1.5.rds")
saveRDS(Ysim11,file = "Desktop/simulation11/data/Y_1.rds")
sim11Scalpel=scalpel("Desktop/simulation11/","Desktop/simulation11/data/",videoHeight = 200,minClusterSize = 5)
nrow(sim11Scalpel$Zhat)
#173 #122
par(mfrow=c(1,2))
spatialplot(sim9A,sim9Z,title="Neuron Footprints",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim11Scalpel$Afilter,sim11Scalpel$Zhat,title="Estimated Footprints (CSNR=1.5,ISNR=1.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)
par(mfrow=c(1,3))
imageVec(Ysim11[,56],videoheight = 200)
title(main="CSNR=1.5,ISNR=1.5,estn=173")

eval11=evaluateSims(sim9A,estA = sim11Scalpel$Afilter,estZ = sim11Scalpel$Zhat,
                    cutoff_extra_pixels=0.2,cutoff_matching_pixels=0.5)
eval11$ID.sens #0.99 #0.99
eval11$ID.prec #0.5722543 #0.8114754
which(eval11$m==0) #23
#par(mfrow=c(1,2))
#temporalplot(sim9Z[which(eval11$m!=0),],title="Temporal Truth")
#temporalplot(sim11Scalpel$Zhat[eval11$m[which(eval11$m!=0)],],title="Temporal Estimates (CSNR=1.5,ISNR=1.5)")
#the temporal estimates are much better
#spatialplot(sim9A[,which(eval11$m!=0)],sim9Z,title="Spatial Truth",border = TRUE,addToPlot = FALSE,videoheight = 200)
#spatialplot(sim11Scalpel$Afilter[,eval11$m[which(eval11$m!=0)]],sim9Z,title="Spatial Estimates (CSNR=1.5,ISNR=1.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)


#========================================================================================
#========================= Adjust Simulation 14 =================================================
#========================================================================================
Ysim14=sim9A%*%sim9Z+
  readRDS("Desktop/simulation9/data/corNoiseSNR_2.rds")+readRDS("Desktop/simulation9/data/indNoiseSNR_1.5.rds")
saveRDS(Ysim14,file = "Desktop/simulation14/data/Y_1.rds")
sim14Scalpel=scalpel("Desktop/simulation14/","Desktop/simulation14/data/",videoHeight = 200,minClusterSize = 5)
nrow(sim14Scalpel$Zhat) #214 #121
par(mfrow=c(1,2))
spatialplot(sim9A,sim9Z,title="Neuron Footprints",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim14Scalpel$Afilter,sim14Scalpel$Zhat,title="Estimated Footprints (CSNR=2,ISNR=1.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)
par(mfrow=c(1,1))
imageVec(Ysim14[,56],videoheight = 200)
title(main = "CSNR=2,ISNR=1.5,estn=214")

eval14=evaluateSims(sim9A,estA = sim14Scalpel$Afilter,estZ = sim14Scalpel$Zhat,
                    cutoff_extra_pixels=0.2,cutoff_matching_pixels=0.5)
eval14$ID.sens #0.99
eval14$ID.prec #0.4626168 #0.8181818
which(eval14$m==0) #23
#par(mfrow=c(1,2))
#temporalplot(sim9Z[which(eval14$m!=0),],title="Temporal Truth")
#temporalplot(sim14Scalpel$Zhat[eval14$m[which(eval14$m!=0)],],title="Temporal Estimates (CSNR=2,ISNR=1.5)")
#the temporal estimates are much better
#spatialplot(sim9A[,which(eval14$m!=0)],sim9Z,title="Spatial Truth",border = TRUE,addToPlot = FALSE,videoheight = 200)
#spatialplot(sim14Scalpel$Afilter[,eval14$m[which(eval14$m!=0)]],sim9Z,title="Spatial Estimates (CSNR=2,ISNR=1.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)
#========================================================================================
#========================= Adjust Simulation 15 =================================================
#========================================================================================
Ysim15=sim9A%*%sim9Z+
  readRDS("Desktop/simulation9/data/corNoiseSNR_1.5.rds")+readRDS("Desktop/simulation9/data/indNoiseSNR_0.5.rds")
saveRDS(Ysim15,file = "Desktop/simulation15/data/Y_1.rds")
sim15Scalpel=scalpel("Desktop/simulation15/","Desktop/simulation15/data/",videoHeight = 200,minClusterSize = 5)
nrow(sim15Scalpel$Zhat) #163 #100
par(mfrow=c(1,2))
spatialplot(sim9A,sim9Z,title="Neuron Footprints",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim15Scalpel$Afilter,sim15Scalpel$Zhat,title="Estimated Footprints (CSNR=1.5,ISNR=0.5)",border = TRUE,addToPlot = FALSE,videoheight = 200)
par(mfrow=c(1,1))
imageVec(Ysim15[,56],videoheight = 200)
title(main = "CSNR=1.5,ISNR=0.5,estn=163")
eval15=evaluateSims(sim9A,estA = sim15Scalpel$Afilter,estZ = sim15Scalpel$Zhat,
                    cutoff_extra_pixels=0.2,cutoff_matching_pixels=0.5)
eval15$ID.sens #0.96 #0.95
eval15$ID.prec #0.5889571 #0.95
which(eval15$m==0) #19 50 74 94 #one more missed neuron is 23
#========================================================================================
#========================= Adjust Simulation 10 =================================================
#========================================================================================
Ysim10=sim9A%*%sim9Z+
  readRDS("Desktop/simulation9/data/corNoiseSNR_1.5.rds")+readRDS("Desktop/simulation9/data/indNoiseSNR_1.rds")
saveRDS(Ysim10,file = "Desktop/simulation10/data/Y_1.rds")
sim10Scalpel=scalpel("Desktop/simulation10/","Desktop/simulation10/data/",videoHeight = 200,minClusterSize = 5)
nrow(sim10Scalpel$Zhat) #191 #113
par(mfrow=c(1,2))
spatialplot(sim9A,sim9Z,title="Neuron Footprints",border = TRUE,addToPlot = FALSE,videoheight = 200)
spatialplot(sim10Scalpel$Afilter,sim10Scalpel$Zhat,title="Estimated Footprints (CSNR=1.5,ISNR=1)",border = TRUE,addToPlot = FALSE,videoheight = 200)
eval10=evaluateSims(sim9A,estA = sim10Scalpel$Afilter,estZ = sim10Scalpel$Zhat,
                    cutoff_extra_pixels=0.2,cutoff_matching_pixels=0.5)
eval10$ID.sens #0.98 #0.98
eval10$ID.prec #0.513089 #0.8672566
which(eval10$m==0) #23 59
#========================================================================================
#========================= Adjust Simulation 12 =================================================
#========================================================================================
Ysim12=sim9A%*%sim9Z+
  readRDS("Desktop/simulation9/data/corNoiseSNR_1.5.rds")+readRDS("Desktop/simulation9/data/indNoiseSNR_2.rds")
saveRDS(Ysim12,file = "Desktop/simulation12/data/Y_1.rds")
sim12Scalpel=scalpel("Desktop/simulation12/","Desktop/simulation12/data/",videoHeight = 200,minClusterSize = 5)
nrow(sim12Scalpel$Zhat) #161  #121

eval12=evaluateSims(sim9A,estA = sim12Scalpel$Afilter,estZ = sim12Scalpel$Zhat,
                    cutoff_extra_pixels=0.2,cutoff_matching_pixels=0.5)
eval12$ID.sens #0.98 #0.98
eval12$ID.prec #0.6086957 # 0.8099174
which(eval12$m==0) #23 74
#========================================================================================
#========================= Adjust Simulation 12 =================================================
#========================================================================================
out_1.5_0.5 = getScalpelStep0("Desktop/simulation15/")
Y_1.5_0.5 = getY(out_1.5_0.5)
out_1.5_1 = getScalpelStep0("Desktop/simulation10/")
Y_1.5_1 = getY(out_1.5_1)
out_1.5_2 = getScalpelStep0("Desktop/simulation12/")
Y_1.5_2 = getY(out_1.5_2)
out_1.5_1.5 = getScalpelStep0("Desktop/simulation11/")
Y_1.5_1.5 = getY(out_1.5_1.5)
out_2_1.5 = getScalpelStep0("Desktop/simulation14/")
Y_2_1.5 = getY(out_2_1.5)
out_1_1.5 = getScalpelStep0("Desktop/simulation13/")
Y_1_1.5 = getY(out_1_1.5)

videoHeight=200
frame = 56

par(mfrow=c(2,4))
plotSpatial(A=sim9A[,which(sim9Z[,frame]>0)]>0, videoHeight = videoHeight, number = FALSE, title = "True Neurons in a Single Frame")
plotFrame(scalpelOutput=out_1_1.5, Y=Y_1_1.5, frame = frame, title = "High Spatial Noise (SSCN=1)")
plotFrame(scalpelOutput=out_1.5_1.5, Y=Y_1.5_1.5, frame = frame, title = "Moderate Spatial Noise (SSCN=1.5)")
plotFrame(scalpelOutput=out_2_1.5, Y=Y_2_1.5, frame = frame, title = "Low Spatial Noise (SSCN=2)")
plotFrame(scalpelOutput=out_1.5_0.5, Y=Y_1.5_0.5, frame = frame, title = "Very High Independent Noise (SIN=0.5)")
plotFrame(scalpelOutput=out_1.5_1, Y=Y_1.5_1, frame = frame, title = "High Independent Noise (SIN=1)")
plotFrame(scalpelOutput=out_1.5_1.5, Y=Y_1.5_1.5, frame = frame, title = "Moderate Independent Noise (SIN=1.5)")
plotFrame(scalpelOutput=out_1.5_2, Y=Y_1.5_2, frame = frame, title = "Low Independent Noise (SIN=2)")


rep = 10
folder="Desktop/"
out_1.5_1.5 = getScalpelStep0(paste0(folder, "scalpelResults/rep",rep,"_SNR1.5_1.5/"))
Y_1.5_1.5 = getY(out_1.5_1.5)

out_1.5_1 = getScalpelStep0(paste0(folder, "scalpelResults/rep",rep,"_SNR1.5_1/"))
Y_1.5_1 = getY(out_1.5_1)

out_1.5_0.5 = getScalpelStep0(paste0(folder, "scalpelResults/rep",rep,"_SNR1.5_0.5/"))
Y_1.5_0.5 = getY(out_1.5_0.5)

out_1.5_2 = getScalpelStep0(paste0(folder, "scalpelResults/rep",rep,"_SNR1.5_2/"))
Y_1.5_2 = getY(out_1.5_2)

out_1_1.5 = getScalpelStep0(paste0(folder, "scalpelResults/rep",rep,"_SNR1_1.5/"))
Y_1_1.5 = getY(out_1_1.5)

out_2_1.5 = getScalpelStep0(paste0(folder, "scalpelResults/rep",rep,"_SNR2_1.5/"))
Y_2_1.5 = getY(out_2_1.5)

trueA = readRDS(paste0(folder, "data/trueAZ/rep",rep,"_trueA.rds"))
trueZ = readRDS(paste0(folder, "data/trueAZ/rep",rep,"_trueZ.rds"))
trueY = trueA %*% trueZ

frame = 300
pdf("SCALPELFigure11.pdf", width=11, height=6.2)
par(mfrow=c(2,4))
plotSpatial(A=trueA[,which(trueZ[,frame]>0)]>0, videoHeight = videoHeight, number = FALSE, title = "True Neurons in a Single Frame")
plotFrame(scalpelOutput=out_1_1.5, Y=Y_1_1.5, frame = frame, title = "High Spatial Noise (SSCN=1)")
plotFrame(scalpelOutput=out_1.5_1.5, Y=Y_1.5_1.5, frame = frame, title = "Moderate Spatial Noise (SSCN=1.5)")
plotFrame(scalpelOutput=out_2_1.5, Y=Y_2_1.5, frame = frame, title = "Low Spatial Noise (SSCN=2)")
plotFrame(scalpelOutput=out_1.5_0.5, Y=Y_1.5_0.5, frame = frame, title = "Very High Independent Noise (SIN=0.5)")
plotFrame(scalpelOutput=out_1.5_1, Y=Y_1.5_1, frame = frame, title = "High Independent Noise (SIN=1)")
plotFrame(scalpelOutput=out_1.5_1.5, Y=Y_1.5_1.5, frame = frame, title = "Moderate Independent Noise (SIN=1.5)")
plotFrame(scalpelOutput=out_1.5_2, Y=Y_1.5_2, frame = frame, title = "Low Independent Noise (SIN=2)")
dev.off()
