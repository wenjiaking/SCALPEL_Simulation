library(R.matlab)
library(gstat)
library(sp)
library(scalpel)
library(Matrix)
library(gam)
##global variables
SNR.vec = c(1,1.5,2,0,0,0,0)
SNR_static.vec = c(0,0,0,0.5, 1, 1.5, 2)
videoHeight = 200
cutoff_matching_pixels = 0.5
cutoff_extra_pixels = 0.2

#1. As SIN increases, both the precision and sensitivity decreases, which is not expected
#2. This implies there are more noise pixels and less neuronal pixels
#3. Explore the effect of preprocessing step to independent noise
codeFolder = "Desktop/UCD-COURSES/papers/CNMF-E_PACKAGE/scalpel/SimulationCode/"
#codeFolder = "home/wwj95/scalpel/SimulationCode/"
#folder = "scratch/wenjia/simulation_SSNR/simResults/"
#figurefolder=paste0(folder,"Figurefolder/")
folder="Desktop/"
if (!dir.exists(paste0(folder,"dividedPixels/"))) dir.create(paste0(folder,"dividedPixels/"))
datafolder=paste0(folder,"data/")

###############################################################################################
### Genertate distributions
###############################################################################################

color.vec=c("deeppink","mediumpurple1","seagreen3","blue")
pch.vec=15:18
lty.vec=1:4

noisepixRaw=readRDS(paste0(folder,"dividedPixels/rep1/noisepixRaw.rds"))
spikpixRaw=readRDS(paste0(folder,"dividedPixels/rep1/spikpixRaw.rds"))
pdf(paste0(figurefolder,"rawDistSIN_SSNR.pdf"),width=10, height=5) 
par(mfrow=c(1,2))
plot(density(noisepixRaw[,1]),xlim=c(min(min(noisepixRaw[,1:3]),min(spikpixRaw[,1:3]))*0.9,max(max(noisepixRaw[,1:3]),max(spikpixRaw[,1:3]))*1.1),ylim=c(0,6),col=color.vec[2],xaxt="n",main="SSNR")
lines(density(noisepixRaw[,2]),col=color.vec[3])
lines(density(noisepixRaw[,3]),col=color.vec[4])
lines(density(spikpixRaw[,1]),col=color.vec[2],lty=2)
lines(density(spikpixRaw[,2]),col=color.vec[3],lty=2)
lines(density(spikpixRaw[,3]),col=color.vec[4],lty=2)


plot(density(noisepixRaw[,4]),xlim=c(min(min(noisepixRaw[,4:7]),min(spikpixRaw[,4:7]))*0.9,max(max(noisepixRaw[,4:7]),max(spikpixRaw[,4:7]))*1.1),ylim=c(0,1.1),col=color.vec[1],xaxt="n",main="SIN")
lines(density(noisepixRaw[,5]),col=color.vec[2])
lines(density(noisepixRaw[,6]),col=color.vec[3])
lines(density(noisepixRaw[,7]),col=color.vec[4])
lines(density(spikpixRaw[,4]),col=color.vec[1],lty=2)
lines(density(spikpixRaw[,5]),col=color.vec[2],lty=2)
lines(density(spikpixRaw[,6]),col=color.vec[3],lty=2)
lines(density(spikpixRaw[,7]),col=color.vec[4],lty=2)

dev.off()

####################after filtering
pdf(paste0(figurefolder,"filterDistSIN_SSNR.pdf"),width=10, height=5) 
par(mfrow=c(1,2))
plot(density(noisepixfilter[,1]),xlim=c(min(min(noisepixfilter[,1:3]),min(spikpixfilter[,1:3]))*0.9,max(max(noisepixfilter[,1:3]),max(spikpixfilter[,1:3]))*1.1),ylim=c(0,8),col=color.vec[2],xaxt="n",main="SSNR")
lines(density(noisepixfilter[,2]),col=color.vec[3])
lines(density(noisepixfilter[,3]),col=color.vec[4])
lines(density(spikpixfilter[,1]),col=color.vec[2],lty=2)
lines(density(spikpixfilter[,2]),col=color.vec[3],lty=2)
lines(density(spikpixfilter[,3]),col=color.vec[4],lty=2)


plot(density(noisepixfilter[,4]),xlim=c(min(min(noisepixfilter[,4:7]),min(spikpixfilter[,4:7]))*0.9,max(max(noisepixfilter[,4:7]),max(spikpixfilter[,4:7]))*1.1),ylim=c(0,13),col=color.vec[1],xaxt="n",main="SIN")
lines(density(noisepixfilter[,5]),col=color.vec[2])
lines(density(noisepixfilter[,6]),col=color.vec[3])
lines(density(noisepixfilter[,7]),col=color.vec[4])
lines(density(spikpixfilter[,4]),col=color.vec[1],lty=2)
lines(density(spikpixfilter[,5]),col=color.vec[2],lty=2)
lines(density(spikpixfilter[,6]),col=color.vec[3],lty=2)
lines(density(spikpixfilter[,7]),col=color.vec[4],lty=2)

dev.off()

####################after deltaf over f transform
pdf(paste0(figurefolder,"deltafDistSIN_SSNR.pdf"),width=10, height=5) 
par(mfrow=c(1,2))
plot(density(noisepixDeltaf[,1]),xlim=c(min(min(noisepixDeltaf[,1:3]),min(spikpixDeltaf[,1:3]))*0.9,max(max(noisepixDeltaf[,1:3]),max(spikpixDeltaf[,1:3]))*1.1),ylim=c(0,10),col=color.vec[2],xaxt="n",main="SSNR")
lines(density(noisepixDeltaf[,2]),col=color.vec[3])
lines(density(noisepixDeltaf[,3]),col=color.vec[4])
lines(density(spikpixDeltaf[,1]),col=color.vec[2],lty=2)
lines(density(spikpixDeltaf[,2]),col=color.vec[3],lty=2)
lines(density(spikpixDeltaf[,3]),col=color.vec[4],lty=2)


plot(density(noisepixDeltaf[,4]),xlim=c(min(min(noisepixDeltaf[,4:7]),min(spikpixDeltaf[,4:7]))*0.9,max(max(noisepixDeltaf[,4:7]),max(spikpixDeltaf[,4:7]))*1.1),ylim=c(0,15),col=color.vec[1],xaxt="n",main="SIN")
lines(density(noisepixDeltaf[,5]),col=color.vec[2])
lines(density(noisepixDeltaf[,6]),col=color.vec[3])
lines(density(noisepixDeltaf[,7]),col=color.vec[4])
lines(density(spikpixDeltaf[,4]),col=color.vec[1],lty=2)
lines(density(spikpixDeltaf[,5]),col=color.vec[2],lty=2)
lines(density(spikpixDeltaf[,6]),col=color.vec[3],lty=2)
lines(density(spikpixDeltaf[,7]),col=color.vec[4],lty=2)

dev.off()

#############overall comparison#######
pdf(paste0(figurefolder,"comp4DistSIN_SSNR.pdf"),width=10, height=20) 
par(mfrow=c(4,2))
plot(density(noisepixRaw[,1]),xlim=c(min(min(noisepixRaw[,1:3]),min(spikpixRaw[,1:3]))*0.9,max(max(noisepixRaw[,1:3]),max(spikpixRaw[,1:3]))*1.1),ylim=c(0,6),col=color.vec[2],xaxt="n",main="SSNR")
lines(density(noisepixRaw[,2]),col=color.vec[3])
lines(density(noisepixRaw[,3]),col=color.vec[4])
lines(density(spikpixRaw[,1]),col=color.vec[2],lty=2)
lines(density(spikpixRaw[,2]),col=color.vec[3],lty=2)
lines(density(spikpixRaw[,3]),col=color.vec[4],lty=2)


plot(density(noisepixRaw[,4]),xlim=c(min(min(noisepixRaw[,4:7]),min(spikpixRaw[,4:7]))*0.9,max(max(noisepixRaw[,4:7]),max(spikpixRaw[,4:7]))*1.1),ylim=c(0,1.1),col=color.vec[1],xaxt="n",main="SIN")
lines(density(noisepixRaw[,5]),col=color.vec[2])
lines(density(noisepixRaw[,6]),col=color.vec[3])
lines(density(noisepixRaw[,7]),col=color.vec[4])
lines(density(spikpixRaw[,4]),col=color.vec[1],lty=2)
lines(density(spikpixRaw[,5]),col=color.vec[2],lty=2)
lines(density(spikpixRaw[,6]),col=color.vec[3],lty=2)
lines(density(spikpixRaw[,7]),col=color.vec[4],lty=2)

plot(density(noisepixfilter[,1]),xlim=c(min(min(noisepixfilter[,1:3]),min(spikpixfilter[,1:3]))*0.9,max(max(noisepixfilter[,1:3]),max(spikpixfilter[,1:3]))*1.1),ylim=c(0,8),col=color.vec[2],xaxt="n",main="SSNR")
lines(density(noisepixfilter[,2]),col=color.vec[3])
lines(density(noisepixfilter[,3]),col=color.vec[4])
lines(density(spikpixfilter[,1]),col=color.vec[2],lty=2)
lines(density(spikpixfilter[,2]),col=color.vec[3],lty=2)
lines(density(spikpixfilter[,3]),col=color.vec[4],lty=2)


plot(density(noisepixfilter[,4]),xlim=c(min(min(noisepixfilter[,4:7]),min(spikpixfilter[,4:7]))*0.9,max(max(noisepixfilter[,4:7]),max(spikpixfilter[,4:7]))*1.1),ylim=c(0,13),col=color.vec[1],xaxt="n",main="SIN")
lines(density(noisepixfilter[,5]),col=color.vec[2])
lines(density(noisepixfilter[,6]),col=color.vec[3])
lines(density(noisepixfilter[,7]),col=color.vec[4])
lines(density(spikpixfilter[,4]),col=color.vec[1],lty=2)
lines(density(spikpixfilter[,5]),col=color.vec[2],lty=2)
lines(density(spikpixfilter[,6]),col=color.vec[3],lty=2)
lines(density(spikpixfilter[,7]),col=color.vec[4],lty=2)

plot(density(noisepixBleach[,1]),xlim=c(min(min(noisepixBleach[,1:3]),min(spikpixBleach[,1:3]))*0.9,max(max(noisepixBleach[,1:3]),max(spikpixBleach[,1:3]))*1.1),ylim=c(0,8),col=color.vec[2],xaxt="n",main="SSNR")
lines(density(noisepixBleach[,2]),col=color.vec[3])
lines(density(noisepixBleach[,3]),col=color.vec[4])
lines(density(spikpixBleach[,1]),col=color.vec[2],lty=2)
lines(density(spikpixBleach[,2]),col=color.vec[3],lty=2)
lines(density(spikpixBleach[,3]),col=color.vec[4],lty=2)


plot(density(noisepixBleach[,4]),xlim=c(min(min(noisepixBleach[,4:7]),min(spikpixBleach[,4:7]))*0.9,max(max(noisepixBleach[,4:7]),max(spikpixBleach[,4:7]))*1.1),ylim=c(0,13),col=color.vec[1],xaxt="n",main="SIN")
lines(density(noisepixBleach[,5]),col=color.vec[2])
lines(density(noisepixBleach[,6]),col=color.vec[3])
lines(density(noisepixBleach[,7]),col=color.vec[4])
lines(density(spikpixBleach[,4]),col=color.vec[1],lty=2)
lines(density(spikpixBleach[,5]),col=color.vec[2],lty=2)
lines(density(spikpixBleach[,6]),col=color.vec[3],lty=2)
lines(density(spikpixBleach[,7]),col=color.vec[4],lty=2)

plot(density(noisepixDeltaf[,1]),xlim=c(min(min(noisepixDeltaf[,1:3]),min(spikpixDeltaf[,1:3]))*0.9,max(max(noisepixDeltaf[,1:3]),max(spikpixDeltaf[,1:3]))*1.1),ylim=c(0,10),col=color.vec[2],xaxt="n",main="SSNR")
lines(density(noisepixDeltaf[,2]),col=color.vec[3])
lines(density(noisepixDeltaf[,3]),col=color.vec[4])
lines(density(spikpixDeltaf[,1]),col=color.vec[2],lty=2)
lines(density(spikpixDeltaf[,2]),col=color.vec[3],lty=2)
lines(density(spikpixDeltaf[,3]),col=color.vec[4],lty=2)


plot(density(noisepixDeltaf[,4]),xlim=c(min(min(noisepixDeltaf[,4:7]),min(spikpixDeltaf[,4:7]))*0.9,max(max(noisepixDeltaf[,4:7]),max(spikpixDeltaf[,4:7]))*1.1),ylim=c(0,15),col=color.vec[1],xaxt="n",main="SIN")
lines(density(noisepixDeltaf[,5]),col=color.vec[2])
lines(density(noisepixDeltaf[,6]),col=color.vec[3])
lines(density(noisepixDeltaf[,7]),col=color.vec[4])
lines(density(spikpixDeltaf[,4]),col=color.vec[1],lty=2)
lines(density(spikpixDeltaf[,5]),col=color.vec[2],lty=2)
lines(density(spikpixDeltaf[,6]),col=color.vec[3],lty=2)
lines(density(spikpixDeltaf[,7]),col=color.vec[4],lty=2)

dev.off()
#########################
#dividedfolder=paste0(folder,"dividedPixels/rep1/")
#thresholdedPix=readRDS(paste0(dividedfolder,"thresholdedPix.rds"))
noisePercent=1-thresholdedPix[3,]/thresholdedPix[2,]
pdf(paste0(figurefolder,"quantileCompare.pdf"),width = 10,height = 5)
par(mfrow=c(1,2))
plot(0,type="n",xaxt="n",xlab = "SSRN",ylab="-quantile",ylim=c(0.26,0.27),xlim=c(0.9,2.2))
axis(1,at=SNR.vec[1:3])
points(SNR.vec[1:3],thresholdedPix[1,1:3],col=color.vec[1],type = "b",pch=pch.vec[1])
points(SNR.vec[1:3],noisequant[1:3],col=color.vec[2],type = "b",pch=pch.vec[2])

plot(0,type="n",xaxt="n",xlab = "SIN",ylab="-quantile",ylim=c(0.1,0.11),xlim=c(0.4,2.2))
axis(1,at=SNR_static.vec[4:7])
points(SNR_static.vec[4:7],thresholdedPix[1,4:7],col=color.vec[1],type = "b",pch=pch.vec[1])
points(SNR_static.vec[4:7],noisequant[4:7],col=color.vec[2],type = "b",pch=pch.vec[2])
dev.off()

pdf(paste0(figurefolder,"thresholdedResult.pdf"),width=5,height = 5)
plot(0,type="n",xaxt="n",xlab = "SNR",ylab="percent",ylim=c(17,29),xlim=c(0.4,2.2))
axis(1,at=SNR_static.vec[4:7])
points(SNR.vec[1:3],100*noisePercent[1:3],col=color.vec[2],type = "b",pch=pch.vec[1])
points(SNR_static.vec[4:7],100*noisePercent[4:7],col=color.vec[3],type = "b",pch=pch.vec[2])
dev.off()

nonNeuronpix=matrix(0,nrow=length(nonNeuronpixID)*1000,ncol=length(SNR.vec))
for (i in length(SNR.vec)) {
  SNR=SNR.vec[i]
  SNR_static=SNR_static.vec[i]
  Y=readRDS(paste0(folder,"dividedPixels/rep1/deltafY_SNR",SNR,"_",SNR_static,".rds"))
  nonNeuronpix[,i]=as.vector(Y[nonNeuronpixID,])
}

pdf(paste0(figurefolder,"noisepixVar.pdf"),width=10,height=5)    
par(mfrow=c(1,2))
plot(0,type="n",xaxt="n",xlab = "SSNR",ylab="sd",ylim=c(0.073,0.162),xlim=c(0.9,2.2))
axis(1,at=SNR.vec[1:3])
points(SNR.vec[1:3],noisepixVar[1,1:3],col=color.vec[1],type = "b",pch=pch.vec[1])
points(SNR.vec[1:3],noisepixVar[2,1:3],col=color.vec[2],type = "b",pch=pch.vec[2])
points(SNR.vec[1:3],noisepixVar[3,1:3],col=color.vec[3],type = "b",pch=pch.vec[3])
points(SNR.vec[1:3],noisepixVar[4,1:3],col=color.vec[4],type = "b",pch=pch.vec[4])

plot(0,type="n",xaxt="n",xlab = "SIN",ylab="sd",ylim=c(0.034,1.16),xlim=c(0.4,2.2))
axis(1,at=SNR_static.vec[4:7])
points(SNR_static.vec[4:7],noisepixVar[1,4:7],col=color.vec[1],type = "b",pch=pch.vec[1])
points(SNR_static.vec[4:7],noisepixVar[2,4:7],col=color.vec[2],type = "b",pch=pch.vec[2])
points(SNR_static.vec[4:7],noisepixVar[3,4:7],col=color.vec[3],type = "b",pch=pch.vec[3])
points(SNR_static.vec[4:7],noisepixVar[4,4:7],col=color.vec[4],type = "b",pch=pch.vec[4])

dev.off()

scalpelfolder="Desktop/data/"
test0.5=scalpel(paste0(scalpelfolder,"finalY/rep1_SNR0_0.5/"),outputFolder=paste0(scalpelfolder,"scalpelTest0.5/"),videoHeight = 200)
test2=scalpel(paste0(scalpelfolder,"finalY/rep1_SNR0_2/"),outputFolder=paste0(scalpelfolder,"scalpelTest2/"),videoHeight = 200)

pdf("Desktop/dividedPixels/rep1/thresholdedFrames.pdf",width = 10,height = 10)
par(mfrow=c(2,2))
plotBrightest(test2,AfilterIndex = 35)
plotBrightest(test0.5,AfilterIndex = 53)
plotThresholdedFrame(test2,frame=136,threshold = test2$lowThreshold)
plotThresholdedFrame(test0.5,frame=135,threshold = test0.5$lowThreshold)
dev.off()

testSNR1=scalpel(paste0(scalpelfolder,"finalY/rep1_SNR1_0/"),outputFolder=paste0(scalpelfolder,"scalpelTestSNR1/"),videoHeight = 200)
testSNR2=scalpel(paste0(scalpelfolder,"finalY/rep1_SNR2_0/"),outputFolder=paste0(scalpelfolder,"scalpelTestSNR2/"),videoHeight = 200)
trueA=readRDS(paste0(scalpelfolder,"trueAZ/rep1_trueA.rds"))
pdf("Desktop/dividedPixels/rep1/scalpelresultSSNR.pdf",width = 10,height = 10)
par(mfrow=c(1,3))
plotSpatial(A=trueA,title="true A",videoHeight = 200,border = TRUE)
plotSpatial(scalpelOutput=testSNR1,neuronSet = "Afilter",title="SSNR=1")
plotSpatial(scalpelOutput = testSNR2,neuronSet = "Afilter",title="SSNR=2")
dev.off()
pdf("Desktop/dividedPixels/rep1/scalpelresultSIN.pdf",width = 10,height = 10)
par(mfrow=c(1,3))
plotSpatial(A=trueA,title="true A",videoHeight = 200,border = TRUE)
plotSpatial(scalpelOutput=test0.5,neuronSet = "Afilter",title="SIN=0.5")
plotSpatial(scalpelOutput = test2,neuronSet = "Afilter",title="SIN=2")
dev.off()

pdf("Desktop/dividedPixels/rep1/thresholdedFramesSSNR.pdf",width = 10,height = 10)
par(mfrow=c(2,2))
plotBrightest(testSNR1,AfilterIndex = 48)
plotBrightest(testSNR2,AfilterIndex = 57)
plotThresholdedFrame(testSNR1,frame=136,threshold = testSNR1$lowThreshold)
plotThresholdedFrame(testSNR2,frame=136,threshold = testSNR2$lowThreshold)
dev.off()

pdf("Desktop/dividedPixels/rep1/missingSSNR.pdf",width = 10,height = 10)
par(mfrow=c(2,2))
plotThresholdedFrame(testSNR2,frame=802,threshold = testSNR2$lowThreshold)
plotThresholdedFrame(testSNR1,frame=802,threshold = testSNR1$lowThreshold)
plotThresholdedFrame(testSNR2,frame=185,threshold = testSNR2$lowThreshold)
plotThresholdedFrame(testSNR1,frame=185,threshold = testSNR1$lowThreshold)
dev.off()

pdf("Desktop/dividedPixels/rep1/filterComp1_0.pdf",width = 10,height = 5)
par(mfrow=c(1,3))
plotSpatial(teststep3,neuronSet = "Afilter")
plotSpatial(test0.5step3,neuronSet = "Afilter")
plotSpatial(test4step3,neuronSet = "Afilter")
dev.off()


##TRY LARGER FILTER WIDTH
##FIGURE OUT WHAT HAPPENED IN DELTAF OVER F STEP: WHY THE THRESHOLDEDS ARE SIMILAR BUT MUCH MORE NOISE PIXELS RETAINED?
##WHAT IF SKIP THE DELTAF OVER F STEP??
##TOO MUCH NOISES MAY MIS COMBINE TWO SIMULTANEOUSLY SPIKING NEURONS,GROUPED LASSO DOES NOT WORK WELL?

#==TWO THINGS===
#(1) ABNORMAL PATTERN UNDER INDEPENDENT NOISE CAN BE IMPROVED BY SETTING LARGER FILTER KERNEL?
#(2) WHAT HAPPENED IN DELTAF OVER F SUCH THAT THE NOISE ARE CENTERED AND SCALED, SPIKING PIXELS ARE CENTERED, AND VARIANCE ARE BOOSTED?
#(3) WHY THERE ARE MUCH MORE NOISE PIXELS RETAINED UNDER SIN=2 EVEN THRESHOLD DOES NOT CHANGE MUCH?
