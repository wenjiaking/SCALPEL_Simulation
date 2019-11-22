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
  spikpixBleach=matrix(0,nrow=length(spikID),ncol=length(SNR.vec))
  noisepixBleach=matrix(0,nrow=length(noiseID),ncol=length(SNR.vec))
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
    
    prepY=preprocessedY(Yfiltered)
    nobleachY=prepY$nobleachY
    spikpixBleach[,noiseIndex]=nobleachY[spikID]
    noisepixBleach[,noiseIndex]=nobleachY[noiseID]
    
    deltafY=prepY$deltafoverfY
    spikpixDeltaf[,noiseIndex]=deltafY[spikID]
    noisepixDeltaf[,noiseIndex]=deltafY[noiseID]
    
    threshold=-quantile(deltafY,probs = 0.001)
    thresholdedPix[1,noiseIndex]=threshold
    thresholdedPix[2,noiseIndex]=length(which(deltafY>threshold))
    thresholdedPix[3,noiseIndex]=length(which(spikpixDeltaf[,noiseIndex]>threshold))
    
    outputfolder=paste0(folder,"dividedPixels/rep",rep,"/")
    if (!dir.exists(outputfolder)) dir.create(outputfolder)
    saveRDS(deltafY,paste0(outputfolder,"deltafY_SNR",SNR,"_",SNR_static,".rds"))
  }
}
saveRDS(spikpixRaw,paste0(outputfolder,"spikpixRaw.rds"))
saveRDS(noisepixRaw,paste0(outputfolder,"noisepixRaw.rds"))
saveRDS(spikpixfilter,paste0(outputfolder,"spikpixfilter.rds"))
saveRDS(noisepixfilter,paste0(outputfolder,"noisepixfilter.rds"))
saveRDS(spikpixBleach,paste0(outputfolder,"spikpixBleach.rds"))
saveRDS(noisepixBleach,paste0(outputfolder,"noisepixBleach.rds"))
saveRDS(spikpixDeltaf,paste0(outputfolder,"spikpixDeltaf.rds"))
saveRDS(noisepixDeltaf,paste0(outputfolder,"noisepixDeltaf.rds"))
saveRDS(thresholdedPix,paste0(outputfolder,"thresholdedPix.rds"))
saveRDS(spikIDlist,paste0(folder,"dividedPixels/spikID_rep(1-",length(rep.vec),").rds"))
saveRDS(noiseIDlist,paste0(folder,"dividedPixels/noiseID_rep(1-",length(rep.vec),").rds"))