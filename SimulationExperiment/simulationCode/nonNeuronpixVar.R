###############compare variance##########
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

imageFastSmooth = function(Y, prepped) {
  m = nrow(Y)
  n = ncol(Y)
  
  temp = matrix(0, nrow = 2*m, ncol = 2*n)
  temp[1:m, 1:n] = Y
  temp = Re(stats::fft(stats::fft(temp) * prepped$W, inverse = TRUE))[1:m, 1:n]
  temp = temp/prepped$temp2
  return(temp)
}

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

nonNeuronpixVar=matrix(0,nrow=4,ncol=length(SNR.vec))

for (noiseIndex in 1:length(SNR.vec)) {
  SNR=SNR.vec[noiseIndex]
  SNR_static = SNR_static.vec[noiseIndex]
  Yraw=readRDS(paste0(folder,"data/finalY/rep1_SNR",SNR,"_",SNR_static,"/Y_1.rds"))
  Yraw=Yraw-min(Yraw)
  nonNeuronpixVar[1,noiseIndex]=sd(Yraw[nonNeuronpixID,])
  
  Yfiltered=filteredY(Yraw,videoHeight = videoHeight)$filteredY
  nonNeuronpixVar[2,noiseIndex]=sd(Yfiltered[nonNeuronpixID,])
  
  prepY=preprocessedY(Yfiltered)
  Ynobleach=prepY$nobleachY
  nonNeuronpixVar[3,noiseIndex]=sd(Ynobleach[nonNeuronpixID,])
  
  deltafY=prepY$deltafoverfY
  nonNeuronpixVar[4,noiseIndex]=sd(deltafY[nonNeuronpixID,])
}
saveRDS(nonNeuronpixVar,paste0(folder,"dividedPixels/rep1/nonNeuronpixVar.rds"))