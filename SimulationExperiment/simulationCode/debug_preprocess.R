###############################################################################################
#1. implement the "SCALPEL_step0.R" smoothing helpers
###############################################################################################
#Helper function for imageFastSmooth below
#m is nrow(Y), n is ncol(Y)
imageFastSmoothHelper = function(m, n) {
  M = 2 * m
  N = 2 * n
  xi = -(m - 1):m
  yi = -(n - 1):n
  dd = sqrt((matrix(xi, M, N)^2 + matrix(yi, M, N, byrow = TRUE)^2))
  out = matrix(0.5 * exp(-abs(dd)), nrow = M, ncol = N)
  out2 = matrix(0, M, N)
  out2[m, n] = 1
  W = stats::fft(out)/stats::fft(out2)
  W = W/(M * N)
  
  temp2 = matrix(0, nrow = M, ncol = N)
  temp2[1:m, 1:n] = 1
  temp2 = Re(stats::fft(stats::fft(temp2) * W, inverse = TRUE))[1:m, 1:n]
  return(list(temp2=temp2, W=W))
}

#Function to perform spatial and temporal smoothing
#a simplified version of fields::image.smooth taking into account that we don't have irregularly spaced data, etc.
#takes in Y and prepped, which is the result from calling imageFastSmoothHelper
#returns the same result as fields::image.smooth(Y)$z
imageFastSmooth = function(Y, prepped) {
  m = nrow(Y)
  n = ncol(Y)
  
  temp = matrix(0, nrow = 2*m, ncol = 2*n)
  temp[1:m, 1:n] = Y
  temp = Re(stats::fft(stats::fft(temp) * prepped$W, inverse = TRUE))[1:m, 1:n]
  temp = temp/prepped$temp2
  return(temp)
}
###############################################################################################
#Gaussian Spatially Temporally Filtering
###############################################################################################
filteredY=function(rawY,videoHeight) {
  prepped = imageFastSmoothHelper(m=videoHeight, n=nrow(rawY)/videoHeight)
  for (frame in 1:ncol(rawY)) {
    rawY[,frame] = as.vector(imageFastSmooth(matrix(rawY[,frame], nrow=videoHeight), prepped))
  }
  spatialFilteredY=rawY
  prepped = imageFastSmoothHelper(1, ncol(rawY))
  for (pixel in 1:nrow(rawY)) {
    rawY[pixel,] = as.vector(imageFastSmooth(matrix(rawY[pixel,],nrow=1), prepped))
  }
  return(list(spatialFilteredY=spatialFilteredY,filteredY=rawY))
}

###############################################################################################
### REMOVE BLEACHING AND PERFORM DELTAF OVER F
###############################################################################################
preprocessedY=function(Y) {
  #Remove bleach effect
  bleachVec = apply(Y, 2, stats::median)
  frames = 1:length(bleachVec)
  bleachModel = gam(bleachVec~s(frames,10))
  #nobleachY = t(t(Y) - bleachModel$fitted.values + max(bleachModel$fitted.values))
  nobleachY = t(t(Y) - bleachModel$fitted.values)
  nobleachY=nobleachY+max(bleachModel$fitted.values)
  
  #deltaf over f transformation
  adjFactor = stats::quantile(nobleachY, probs=0.1)
  deltafoverfY = t(apply(nobleachY, 1, function(vec, adjFactor)
    (vec - stats::median(vec))/(stats::median(vec) + adjFactor), adjFactor=adjFactor))
  return(list(nobleachY=nobleachY,deltafoverfY=deltafoverfY))
}


###############################################################################################
### SEPERATE NNOISE PIXELS AND SPIKING PIXELS
###############################################################################################
spikIDlist=vector("list",length(rep.vec))
noiseIDlist=vector("list",length(rep.vec))
for (rep in rep.vec) {
  trueY=readRDS(paste0(datafolder,"trueAZ/rep",rep,"_trueA.rds"))%*%readRDS(paste0(datafolder,"trueAZ/rep",rep,"_trueZ.rds"))
  spikID=which(trueY>0)
  noiseID=which(!(trueY>0))
  spikIDlist[[rep]]=spikID
  noiseIDlist[[rep]]=noiseID
  spikpixRaw=matrix(0,nrow=length(spikID),ncol=length(SNR.vec))
  noisepixRaw=matrix(0,nrow=length(noiseID),ncol=length(SNR.vec))
  spikpixfilter=matrix(0,nrow=length(spikID),ncol=length(SNR.vec))
  noisepixfilter=matrix(0,nrow=length(noiseID),ncol=length(SNR.vec))
  spikpixDeltaf=matrix(0,nrow=length(spikID),ncol=length(SNR.vec))
  noisepixDeltaf=matrix(0,nrow=length(noiseID),ncol=length(SNR.vec))
  thresholdedPix=matrix(0,nrow=3,ncol=length(SNR.vec))
  
  for (noiseIndex in 1:length(SNR.vec)) {
    SNR=SNR.vec[noiseIndex]
    SNR_static = SNR_static.vec[noiseIndex]
    Yraw=readRDS(paste0(datafolder,"finalY/rep",rep,"_SNR",SNR,"_",SNR_static,"/Y_1.rds"))
    Yraw=Yraw-min(Yraw)
    spikpixRaw[,noiseIndex]=Yraw[spikID]
    noisepixRaw[,noiseIndex]=Yraw[noiseID]
    
    Yfiltered=filteredY(Yraw,videoHeight = videoHeight)$filteredY
    spikpixfilter[,noiseIndex]=Yfiltered[spikID]
    noisepixfilter[,noiseIndex]=Yfiltered[noiseID]
    
    deltafY=preprocessedY(Yfiltered)$deltafoverfY
    spikpixDeltaf[,noiseIndex]=deltafY[spikID]
    noisepixDeltaf[,noiseIndex]=deltafY[noiseID]
    
    threshold=-quantile(deltafY,probs = 0.001)
    thresholdedPix[1,noiseIndex]=threshold
    thresholdedPix[2,noiseIndex]=length(which(deltafY>threshold))
    thresholdedPix[3,noiseIndex]=length(which(spikpixDeltaf[,noiseIndex]>threshold))
    
    outputfolder=paste0(folder,"dividedPixels/rep",rep,"/")
    if (!dir.exists(outputfolder)) dir.create(outputfolder)
    saveRDS(spikpixRaw,paste0(outputfolder,"spikpixRaw.rds"))
    saveRDS(noisepixRaw,paste0(outputfolder,"noisepixRaw.rds"))
    saveRDS(spikpixfilter,paste0(outputfolder,"spikpixfilter.rds"))
    saveRDS(noisepixfilter,paste0(outputfolder,"noisepixfilter.rds"))
    saveRDS(spikpixDeltaf,paste0(outputfolder,"spikpixDeltaf.rds"))
    saveRDS(noisepixDeltaf,paste0(outputfolder,"noisepixDeltaf.rds"))
    saveRDS(thresholdedPix,paste0(outputfolder,"thresholdedPix.rds"))
  }
}
saveRDS(spikIDlist,paste0(folder,"dividedPixels/spikID_rep(1-",length(rep.vec),").rds"))
saveRDS(noiseIDlist,paste0(folder,"dividedPixels/noiseID_rep(1-",length(rep.vec),").rds"))
#####################Plot the true footprints #############

plotNeuronOnFrame = function(A, neuron, videoHeight, col, pctTransp=1, border=TRUE) {
  if (border==TRUE) {
    toPlot = as.vector(getBorderMat(mat = matrix(A[,neuron], nrow = videoHeight)))
  } else toPlot = (A[,neuron]==1)
  toPlot[which(toPlot==0)] = NA
  col = transp(col, pctTransp)
  image(z=t(matrix(toPlot, nrow=videoHeight))[,videoHeight:1], zlim=c(0,1),
        axes=FALSE, col=c("white",col), add=TRUE)
}
plotTruth=function(trueA,videoHeight,border=TRUE,pctTransp=0.5) {
  plot.new()
  for (i in 1:ncol(trueA)) {
    plotNeuronOnFrame(A=trueA, neuron=i, videoHeight=videoHeight, col="black", border = border, pctTransp=pctTransp)
  }
}
trueA=readRDS("Desktop/data/trueAZ/rep1_trueA.rds")
pdf("Desktop/truth_rep1.pdf",width=10, height=10)
origMar = graphics::par()$mar
graphics::par(mar=c(0.5,0.5,3,0.5),mfrow=c(1,2))
plotTruth(trueA,videoHeight = 200)
dev.off()













#Using frame 300 as an example, it turns out that the neuronal pixels (811) have larger values, 
#  eg. after preprocessing, the neuron pixles in frame 300 range from -0.09810523 to 0.29903655,
#  compared to the frame ranging from -0.1773774 to 0.2990366.
#The 0.1% quantile of -0.148769 used in threshoding will only keep 73 pixels in the frame, and 70 of them are neuronal pixels
#Thus, with correlated SNR 1.5 and independent SNR 1.5, the propressing works for frame 300

pixelsRem=function(trueY, deltafoverfY) {
  remainPixels=matrix(0,nrow=ncol(trueY),ncol = 4)
  thre=-quantile(deltafoverfY,probs = 0.001)
  for (i in 1:ncol(trueY)) {
    neuPix=which(trueY[,i]>0)
    remainPixels[i,1]=length(neuPix)
    remainPixels[i,2]=sum(deltafoverfY[neuPix,i]>thre)
    remainPixels[i,3]=sum(deltafoverfY[,i]>thre)
    remainPixels[i,4]=remainPixels[i,2]/remainPixels[i,3]
  }
  return(remainPixels)
}
remainPix=pixelsRem(trueY = trueY,deltafoverfY)

#It turns out 147 frame kept the largest number of non-neuronal pixels, 
#  and in 101 frames, remaining neuronal pixels account for less than 50% of kept pixels
par(mfrow=c(1,4))
imageVec(trueY[,147],200)
title(main = "Neuronal signal")
imageVec(Y1_1.5_1.5[,147],200)
title(main = "Raw Data")
imageVec(deltafoverfY[,147],200)
title(main = "deltafoverf Data")
imageVec(deltafoverfY[,147]>(-quantile(deltafoverfY,probs = 0.001)),200)
title(main = "Thresholded Data")


##With weaker correlated noise 
Y1_2_1.5=readRDS("Desktop/data/finalY/rep1_SNR2_1.5/Y_1.rds")
Y1_2_1.5=Y1_2_1.5-min(Y1_2_1.5)
#range(trueY)
range(Y1_2_1.5) #3.648816

Y2=Y1_2_1.5
prepped = imageFastSmoothHelper(m=videoHeight, n=nrow(Y2)/videoHeight)
for (frame in 1:ncol(Y2)) {
  Y2[,frame] = as.vector(imageFastSmooth(matrix(Y2[,frame], nrow=videoHeight), prepped))
}
range(Y2)
#par(mfrow=c(1,3))
#imageVec(Y1_2_1.5[,300],200)
#title(main = "Raw Data")
#imageVec(Y2[,300],200)
#title(main = "Spatially Smoothed Data")

prepped = imageFastSmoothHelper(1, ncol(Y2))
for (pixel in 1:nrow(Y2)) {
  Y2[pixel,] = as.vector(imageFastSmooth(matrix(Y2[pixel,],nrow=1), prepped))
}
#imageVec(Y2[,300],200)
#title(main = "Temporally Smoothed Data")

bleachVec = apply(Y2, 2, stats::median)
frames = 1:length(bleachVec)
bleachModel = gam(bleachVec~s(frames,10))
#nobleachY = t(t(Y) - bleachModel$fitted.values + max(bleachModel$fitted.values))
nobleachY2 = t(t(Y2) - bleachModel$fitted.values)
range(nobleachY2) # -0.5000744  1.3928397
nobleachY2=nobleachY2+max(bleachModel$fitted.values)
range(nobleachY2) #0.6483492 2.5412632
#imageVec(nobleachY2[,300],200)
#title(main = "Bleach Removed Data")

adjFactor = stats::quantile(nobleachY2, probs=0.1)
deltafoverfY2 = t(apply(nobleachY2, 1, function(vec, adjFactor)
  (vec - stats::median(vec))/(stats::median(vec) + adjFactor), adjFactor=adjFactor))
#imageVec(deltafoverfY2[,300],200)
#title(main = "deltafoverf Data")
#hist(deltafoverfY2[,300])
#hist(nobleachY2[,300])
par(mfrow=c(1,4))
imageVec(Y1_1.5_1.5[,147],200)
title(main = "SSCN=1.5")
imageVec(Y1_2_1.5[,147],200)
title(main = "SSCN=2")
imageVec(deltafoverfY[,147],200)
title(main = "SSCN=1.5")
imageVec(deltafoverfY2[,147],200)
title(main = "SSCN=2")
par(mfrow=c(1,3))
imageVec(trueY[,147],200)
title(main = "True neuron")
imageVec(deltafoverfY[,147]>(-quantile(deltafoverfY,probs = 0.001)),200)
title(main = "SSCN=1.5")
imageVec(deltafoverfY2[,147]>(-quantile(deltafoverfY2,probs = 0.001)),200)
title(main = "SSCN=2")

remainPix2=pixelsRem(trueY = trueY,deltafoverfY2)
noiseFrID2=which(remainPix2[,4]<0.8) #209
noiseFrID=which(remainPix[,4]<0.8) #205
noiseFrID2[which.max(remainPix2[noiseFrID2,4])] #612
noiseFrID[which.max(remainPix[noiseFrID,4])] #698

par(mfrow=c(1,3))
imageVec(trueY[,698],200)
title(main = "True neuron: 698")
imageVec(Y1_1.5_1.5[,698],200)
title(main = "Raw Data: 698")
imageVec(deltafoverfY[,698]>(-quantile(deltafoverfY,probs = 0.001)),200)
title(main = "SSCN=1.5,698")
par(mfrow=c(1,3))
imageVec(trueY[,612],200)
title(main = "True neuron: 612")
imageVec(Y1_2_1.5[,612],200)
title(main = "Raw Data: 612")
imageVec(deltafoverfY2[,612]>(-quantile(deltafoverfY2,probs = 0.001)),200)
title(main = "SSCN=2,612")

par(mfrow=c(1,2))
plot(remainPix2[,4],main="SSCN=2 (weaker noise)",pch=20,ylim=c(0.7,1))
plot(remainPix[,4], main="SSCN=1.5",pch=20,ylim=c(0.7,1))

Y1_1_1.5=readRDS("Desktop/data/finalY/rep1_SNR1_1.5/Y_1.rds")
Y1_1_1.5=Y1_1_1.5-min(Y1_1_1.5)
Y1=Y1_1_1.5
prepped = imageFastSmoothHelper(m=videoHeight, n=nrow(Y1)/videoHeight)
for (frame in 1:ncol(Y1)) {
  Y1[,frame] = as.vector(imageFastSmooth(matrix(Y1[,frame], nrow=videoHeight), prepped))
}
range(Y1) #0.559914 3.175733

prepped = imageFastSmoothHelper(1, ncol(Y1))
for (pixel in 1:nrow(Y1)) {
  Y1[pixel,] = as.vector(imageFastSmooth(matrix(Y1[pixel,],nrow=1), prepped))
}
#imageVec(Y2[,300],200)
#title(main = "Temporally Smoothed Data")

bleachVec = apply(Y1, 2, stats::median)
frames = 1:length(bleachVec)
bleachModel = gam(bleachVec~s(frames,10))
#nobleachY = t(t(Y) - bleachModel$fitted.values + max(bleachModel$fitted.values))
nobleachY1 = t(t(Y1) - bleachModel$fitted.values)
range(nobleachY1) # -0.8523546  1.4788639
nobleachY1=nobleachY1+max(bleachModel$fitted.values)
range(nobleachY1) #0.7781192 3.1093377
#imageVec(nobleachY2[,300],200)
#title(main = "Bleach Removed Data")

adjFactor = stats::quantile(nobleachY1, probs=0.1)
deltafoverfY1 = t(apply(nobleachY1, 1, function(vec, adjFactor)
  (vec - stats::median(vec))/(stats::median(vec) + adjFactor), adjFactor=adjFactor))

remainPix1=pixelsRem(trueY = trueY,deltafoverfY1)
length(which(remainPix1[,4]>0.95)) #601
length(which(remainPix[,4]>0.95)) #562
length(which(remainPix2[,4]>0.95)) #419

par(mfrow=c(1,3))
plot(remainPix1[,4],main="SSCN=1 (strong noise)",pch=20,ylim=c(0.8,1))
plot(remainPix[,4], main="SSCN=1.5",pch=20,ylim=c(0.8,1))
plot(remainPix2[,4],main="SSCN=2 (weaker noise)",pch=20,ylim=c(0.8,1))

#The number of frames where remaining neuronal pixels accounts for more than 90% is smaller under weaker noise level (598 vs 667),
#  but more frames where neuronal pixels take up more than 70%, 
#  in other words, between 70% and 100%, weaker noise level shows larger variation than strong noise level.
# One key is that whether less than 10% noise pixels will be identified as elements in the following steps. 
#  If it will, ie only when the neuronal pixels account for more than 90%, the redundant can be ignored in preliminary segmentation,
#     it makes sense that weaker noise level causes more mis-conctructed preliminary elements.
#  One direction is to adjust the way to identify preliminary elements, instead of 4-notion connectivity components,
#   to make the 10% remaining noise pixels neglectable (will not be mis-identified as elements).
#Is this a general situation that as the noise level becomes weaker, the author's assumption does not hold.
##############################
neuPix=which(trueY[,147]>0)
range(deltafoverfY[neuPix,147])
range(deltafoverfY[-neuPix,147])
label=rep(0,nrow(deltafoverfY))
label[neuPix]=1
mat=cbind(deltafoverfY[,147],label)
mat=as.data.frame(mat)
mat$label=as.factor(mat$label)
levels(mat$label)=c("noisePix","neuronPix")
install.packages("sm")
install.packages("XQuartz")
library(sm)

par(mfrow=c(1,3))
plot(density(deltafoverfY1[-neuPix,147]),xlim=range(deltafoverfY1[,147]),main="SSCR=1")
lines(density(deltafoverfY1[neuPix,147]),col="red")
lines(density(deltafoverfY1[,147]),col="blue")

plot(density(deltafoverfY[-neuPix,147]),xlim=range(deltafoverfY[,147]),main="SSCR=1.5")
lines(density(deltafoverfY[neuPix,147]),col="red")
lines(density(deltafoverfY[,147]),col="blue")

plot(density(deltafoverfY2[-neuPix,147]),xlim=range(deltafoverfY2[,147]),main="SSCR=2")
lines(density(deltafoverfY2[neuPix,147]),col="red")
lines(density(deltafoverfY2[,147]),col="blue")



######
dividePixels=function(trueY, Y) {
  noisePixels=c()
  neuronPixels=c()
  for (i in 1:ncol(Y)) {
    neuronID=which(trueY[,i]>0)
    neuronPixels=c(neuronPixels,as.vector(Y[neuronID,i]))
    noisePixels=c(noisePixels,as.vector(Y[-neuronID,i]))
  }
  return(list(neuronPixels=neuronPixels,noisePixels=noisePixels))
}
divided1=dividePixels(trueY,Y=deltafoverfY1)
divided=dividePixels(trueY,Y=deltafoverfY)
divided2=dividePixels(trueY,Y=deltafoverfY2)

pdf("Desktop/PixelsDistribution.pdf", width=10, height=5)
par(mfrow=c(1,3), mar=c(5, 4, 4, 2) + 0.1)

plot(density(divided1$noisePixels),xlim=c(min(deltafoverfY1)*0.9,max(deltafoverfY1)*1.1),col="mediumpurple1",xaxt="n",main="SSNR=1")
axis(1, at=c(quantile(deltafoverfY1,probs = 0.001),-quantile(deltafoverfY1,probs = 0.001),-min(deltafoverfY1)))
lines(density(divided1$neuronPixels),col="seagreen3")
lines(density(deltafoverfY1),lty=3)
abline(v=-quantile(deltafoverfY1,probs = 0.001),lty=2)

plot(density(divided$noisePixels),xlim=c(min(deltafoverfY)*0.9,max(deltafoverfY)*1.1),col="mediumpurple1",xaxt="n",main="SSNR=1.5")
axis(1, at=c(quantile(deltafoverfY,probs = 0.001),-quantile(deltafoverfY,probs = 0.001),-min(deltafoverfY)))
lines(density(divided$neuronPixels),col="seagreen3")
lines(density(deltafoverfY),lty=3)
abline(v=-quantile(deltafoverfY,probs = 0.001),lty=2)

plot(density(divided2$noisePixels),xlim=c(min(deltafoverfY2)*0.9,max(deltafoverfY2)*1.1),col="mediumpurple1",xaxt="n",main="SSNR=2")
axis(1, at=c(quantile(deltafoverfY2,probs = 0.001),-quantile(deltafoverfY2,probs = 0.001),-min(deltafoverfY2)))
lines(density(divided2$neuronPixels),col="seagreen3")
lines(density(deltafoverfY2),lty=3)
abline(v=-quantile(deltafoverfY2,probs = 0.001),lty=2)

dev.off()

######
neuronProp=function(remainPix,propVec=c(0,0.1,0.2,0.3)){
  nnf=rep(0,length(propVec))
  for (i in 1:length(propVec)) {
    nnf[i]=length(which((remainPix[,4]>=(1-propVec[i]))))
  }
  return(nnf)
}

nNoisepix=function(remainPix) {
  nNoise=c(mean(remainPix[,3]-remainPix[,2]),mean(remainPix[,2],na.rm=TRUE))
  return(nNoise)
}

nnf1=neuronProp(remainPix1)
nnf=neuronProp(remainPix)
nnf2=neuronProp(remainPix2)
nPix=rbind(nNoisepix(remainPix1),nNoisepix(remainPix),nNoisepix(remainPix2))

pdf("Desktop/remainPixels.pdf", width=10, height=5)
par(mfrow=c(1,2), mar=c(5, 4, 4, 2) + 0.1)
plot(0, type="n",main="# of frames vs % of retained noise pixels",xaxt="n",xlab="percentage",ylab="nframes",xlim=c(0,31),ylim=c(min(c(min(nnf1),min(nnf),min(nnf2)))*0.9,max(c(max(nnf1),max(nnf),max(nnf2)))*1.1))
axis(1, at=c(0,10,20,30))
points(100*c(0,0.1,0.2,0.3), nnf1, pch=15, col="mediumpurple1", type="b", lty=2)
points(100*c(0,0.1,0.2,0.3), nnf, pch=16, col="seagreen3", type="b", lty=2)
points(100*c(0,0.1,0.2,0.3), nnf2, pch=17, type="b", lty=2)

plot(0, type="n",main="Average number of pixles",xaxt="n",xlab="SSNR",ylab="nPixels",xlim=c(0.9,2.2),ylim=c(30,220))
axis(1, at=c(1,1.5,2))
points(c(1,1.5,2), nPix[,1], pch=15, col="mediumpurple1", type="b", lty=2)
points(c(1,1.5,2), nPix[,2], pch=16, col="seagreen3", type="b", lty=2)

dev.off()
####################explore the same pattern for SIN########
Y1_1.5_0.5=readRDS("Desktop/data/finalY/rep1_SNR1.5_0.5/Y_1.rds")
Y1_1.5_1=readRDS("Desktop/data/finalY/rep1_SNR1.5_1/Y_1.rds")
Y1_1.5_1.5=readRDS("Desktop/data/finalY/rep1_SNR1.5_1.5/Y_1.rds")
Y1_1.5_2=readRDS("Desktop/data/finalY/rep1_SNR1.5_2/Y_1.rds")
Y1_1_1.5=readRDS("Desktop/data/finalY/rep1_SNR1_1.5/Y_1.rds")
Y1_2_1.5=readRDS("Desktop/data/finalY/rep1_SNR2_1.5/Y_1.rds")
preprocess=function(Y1,videoHeight=200) {
  Y1=Y1-min(Y1)
  prepped = imageFastSmoothHelper(m=videoHeight, n=nrow(Y1)/videoHeight)
  for (frame in 1:ncol(Y1)) {
    Y1[,frame] = as.vector(imageFastSmooth(matrix(Y1[,frame], nrow=videoHeight), prepped))
  }
  prepped = imageFastSmoothHelper(1, ncol(Y1))
  for (pixel in 1:nrow(Y1)) {
    Y1[pixel,] = as.vector(imageFastSmooth(matrix(Y1[pixel,],nrow=1), prepped))
  }
  bleachVec = apply(Y1, 2, stats::median)
  frames = 1:length(bleachVec)
  bleachModel = gam(bleachVec~s(frames,10))
  nobleachY1 = t(t(Y1) - bleachModel$fitted.values)
  nobleachY1=nobleachY1+max(bleachModel$fitted.values)
  adjFactor = stats::quantile(nobleachY1, probs=0.1)
  deltafoverfY1 = t(apply(nobleachY1, 1, function(vec, adjFactor)
    (vec - stats::median(vec))/(stats::median(vec) + adjFactor), adjFactor=adjFactor))
  return(deltafoverfY1)
}

deltafoverfY_0.5=preprocess(Y1_1.5_0.5)
deltafoverfY=preprocess(Y1_1.5_1.5)
deltafoverfY_1=preprocess(Y1_1.5_1)
deltafoverfY_2=preprocess(Y1_1.5_2)
deltafoverfY1=preprocess(Y1_1_1.5)
deltafoverfY2=preprocess(Y1_2_1.5)

remainPix_0.5=pixelsRem(trueY = trueY,deltafoverfY_0.5)
remainPix_1=pixelsRem(trueY = trueY,deltafoverfY_1)
remainPix=pixelsRem(trueY = trueY,deltafoverfY)
remainPix_2=pixelsRem(trueY = trueY,deltafoverfY_2)
remainPix1=pixelsRem(trueY = trueY,deltafoverfY1)
remainPix2=pixelsRem(trueY = trueY,deltafoverfY2)

divided_0.5=dividePixels(trueY,Y=deltafoverfY_0.5)
divided_1=dividePixels(trueY,Y=deltafoverfY_1)
divided=dividePixels(trueY,Y=deltafoverfY)
divided_2=dividePixels(trueY,Y=deltafoverfY_2)
divided1=dividePixels(trueY,Y=deltafoverfY1)
divided2=dividePixels(trueY,Y=deltafoverfY2)

pdf("Desktop/PixelsDistribution_SIN.pdf", width=10, height=5)
par(mfrow=c(1,4), mar=c(5, 4, 4, 2) + 0.1)

plot(density(divided_0.5$noisePixels),xlim=c(min(deltafoverfY_0.5)*0.9,max(deltafoverfY_0.5)*1.1),col="mediumpurple1",xaxt="n",main="SIN=0.5")
axis(1, at=c(quantile(deltafoverfY_0.5,probs = 0.001),-quantile(deltafoverfY_0.5,probs = 0.001)))
lines(density(divided_0.5$neuronPixels),col="seagreen3")
lines(density(deltafoverfY_0.5),lty=3)
abline(v=-quantile(deltafoverfY_0.5,probs = 0.001),lty=2)

plot(density(divided_1$noisePixels),xlim=c(min(deltafoverfY_1)*0.9,max(deltafoverfY_1)*1.1),col="mediumpurple1",xaxt="n",main="SIN=1")
axis(1, at=c(quantile(deltafoverfY_1,probs = 0.001),-quantile(deltafoverfY_1,probs = 0.001)))
lines(density(divided_1$neuronPixels),col="seagreen3")
lines(density(deltafoverfY_1),lty=3)
abline(v=-quantile(deltafoverfY_1,probs = 0.001),lty=2)

plot(density(divided$noisePixels),xlim=c(min(deltafoverfY)*0.9,max(deltafoverfY)*1.1),col="mediumpurple1",xaxt="n",main="SIN=1.5")
axis(1, at=c(quantile(deltafoverfY,probs = 0.001),-quantile(deltafoverfY,probs = 0.001)))
lines(density(divided$neuronPixels),col="seagreen3")
lines(density(deltafoverfY),lty=3)
abline(v=-quantile(deltafoverfY,probs = 0.001),lty=2)

plot(density(divided_2$noisePixels),xlim=c(min(deltafoverfY_2)*0.9,max(deltafoverfY_2)*1.1),col="mediumpurple1",xaxt="n",main="SIN=2")
axis(1, at=c(quantile(deltafoverfY_2,probs = 0.001),-quantile(deltafoverfY_2,probs = 0.001)))
lines(density(divided_2$neuronPixels),col="seagreen3")
lines(density(deltafoverfY_2),lty=3)
abline(v=-quantile(deltafoverfY_2,probs = 0.001),lty=2)
dev.off()

nnf_1=neuronProp(remainPix_1)
nnf_0.5=neuronProp(remainPix_0.5)
nnf_2=neuronProp(remainPix_2)
nnf=neuronProp(remainPix)
nPix_SIN=rbind(nNoisepix(remainPix_0.5),nNoisepix(remainPix_1),nNoisepix(remainPix),nNoisepix(remainPix_2))

pdf("Desktop/remainPixels_SIN.pdf", width=10, height=5)
par(mfrow=c(1,2), mar=c(5, 4, 4, 2) + 0.1)
plot(0, type="n",main="# of frames vs % of retained noise pixels",xaxt="n",xlab="percentage",ylab="nframes",xlim=c(0,31),ylim=c(min(c(min(nnf_0.5),min(nnf_1),min(nnf),min(nnf_2)))*0.9,max(c(max(nnf_0.5),max(nnf_1),max(nnf),max(nnf_2)))*1.1))
axis(1, at=c(0,10,20,30))
points(100*c(0,0.1,0.2,0.3), nnf_0.5, pch=18, col="deeppink", type="b", lty=2)
points(100*c(0,0.1,0.2,0.3), nnf_1, pch=15, col="mediumpurple1", type="b", lty=2)
points(100*c(0,0.1,0.2,0.3), nnf, pch=16, col="seagreen3", type="b", lty=2)
points(100*c(0,0.1,0.2,0.3), nnf_2, pch=17, type="b", lty=2)

plot(0, type="n",main="Average number of pixles",xaxt="n",xlab="SSNR",ylab="nPixels",xlim=c(0.4,2.2),ylim=c(32,170))
axis(1, at=c(0.5,1,1.5,2))
points(c(0.5,1,1.5,2), nPix_SIN[,1], pch=15, col="mediumpurple1", type="b", lty=2)
points(c(0.5,1,1.5,2), nPix_SIN[,2], pch=16, col="seagreen3", type="b", lty=2)

dev.off()

pdf("Desktop/SSRN_SIN.pdf", width=10, height=5)
par(mfrow=c(1,2), mar=c(5, 4, 4, 2) + 0.1)
plot(density(divided_0.5$noisePixels),xlim=c(-0.271*0.9,0.665*1.1),col="deeppink",xaxt="n",main="SSRN=1.5")
lines(density(divided_1$noisePixels),col="mediumpurple1")
lines(density(divided$noisePixels),col="seagreen3")
lines(density(divided_2$noisePixels))
lines(density(divided_0.5$neuronPixels),col="deeppink",lty=2)
lines(density(divided_1$neuronPixels),col="mediumpurple1",lty=2)
lines(density(divided$neuronPixels),col="seagreen3",lty=2)
lines(density(divided_2$neuronPixels),lty=2)

plot(density(divided2$noisePixels),xlim=c(-0.283*0.9,0.635*1.1),xaxt="n",main="SIN=1.5")
lines(density(divided$noisePixels),col="seagreen3")
lines(density(divided1$noisePixels),col="mediumpurple1")
lines(density(divided1$neuronPixels),col="mediumpurple1",lty=2)
lines(density(divided$neuronPixels),col="seagreen3",lty=2)
lines(density(divided2$neuronPixels),lty=2)
dev.off()


for (noiseIndex in 1:length(SNR.vec)) {
  SNR=SNR.vec[noiseIndex]
  SNR_static = SNR_static.vec[noiseIndex]
  Yraw=readRDS(paste0(datafolder,"finalY/rep",rep,"_SNR",SNR,"_",SNR_static,"/Y_1.rds"))
  Yraw=Yraw-min(Yraw)
  Yfiltered=filteredY(Yraw,videoHeight = videoHeight)$filteredY
  prepY=preprocessedY(Yfiltered)
  nobleachY=prepY$nobleachY
  spikpixBleach[,noiseIndex]=nobleachY[spikID]
  noisepixBleach[,noiseIndex]=nobleachY[noiseID]
  
  deltafY=prepY$deltafoverfY
  spikpixDeltaf=deltafY[spikID]
  threshold=-quantile(deltafY,probs = 0.001)
  thresholdedPix[1,noiseIndex]=threshold
  thresholdedPix[2,noiseIndex]=length(which(deltafY>threshold))
  thresholdedPix[3,noiseIndex]=length(which(spikpixDeltaf>threshold))
  outputfolder=paste0(folder,"dividedPixels/rep",rep,"/")
  saveRDS(deltafY,paste0(outputfolder,"deltafY_SNR",SNR,"_",SNR_static,".rds"))
 
}
saveRDS(spikpixBleach,paste0(outputfolder,"spikpixBleach.rds"))
saveRDS(noisepixBleach,paste0(outputfolder,"noisepixBleach.rds"))
saveRDS(thresholdedPix,paste0(outputfolder,"thresholdedPix.rds"))

threshold=readRDS("Desktop/dividedPixels/rep1/thresholdedPix.rds")
threshold
