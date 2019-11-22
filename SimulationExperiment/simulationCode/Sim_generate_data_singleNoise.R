#the following script generate the simulated data for all noise scenarios and replicates

###############################################################################################
### Functions
###############################################################################################

#vector giving intensity of spatial noise over time
#vector is length 'nFrames' with the largest value at the element 'peakFrame'
#and the overall pattern of one smooth peak with width 'width' and zero elsewhere
generateTimeTrend = function(nFrames, peakFrame, width) {
  intensity = rep(0, nFrames)
  cosCurve = (cos(seq(-width/2,width/2,length=width)*2*pi/width) + 1)/2
  indices = (peakFrame-round(width/2)):(peakFrame-round(width/2)+width-1)
  replace = match(indices, seq(intensity))
  intensity[replace[!is.na(replace)]] = cosCurve[!is.na(replace)]
  return(intensity)
}

#generate spatially correlated noise for the simulated video
#used the following tutorial when writing this code:
#http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/
generateSpatialNoise = function(videoHeight, psill=0.025, range=25) {
  
  #create structure
  xy = expand.grid(1:videoHeight, 1:videoHeight)
  names(xy) = c("x","y")
  
  #define the gstat object (spatial model)
  #beta is the expected value, psill controls the variability, 
  #range controls the degree of spatial correlation (larger values give 'courser' autocorrelation)
  g.dummy = gstat(formula=z~1, locations=~x+y, dummy=T, beta=0, model=vgm(psill=psill, model="Exp", range=range), nmax=20)
  
  #make a simulation based on the stat object
  yy = predict(g.dummy, newdata=xy, nsim=1)
  #plot it
  #gridded(yy) = ~x+y
  #spplot(yy[1])
  return(yy$sim1)
}

generateCalciumVec = function(peakShape, nFrames=1000, minPeaks = 1, maxPeaks = 3) {
  vec = rep(0, nFrames)
  peakStarts = sample(1:(nFrames-50), sample(minPeaks:maxPeaks, 1))
  if (length(peakStarts)>1) while(min(c(dist(peakStarts)))<=50) {
    peakStarts = sample(1:(nFrames-50), sample(minPeaks:maxPeaks, 1))
    if (length(peakStarts)==1) break
  }
  peakMax = runif(length(peakStarts), 0.8, 1.2)
  for (i in 1:length(peakStarts)) vec[peakStarts[i]:(peakStarts[i]+length(peakShape)-1)] = peakShape * peakMax[i]
  return(vec)
}

#shape of calcium peak and decay, which was generated using Diego's Matlab program
#max of 1
peakShape = c(0.14896831, 0.38829047, 0.68241502, 0.86081054, 0.96901290, 1.00000000, 0.97678037, 0.91378487, 0.85485207,
              0.79972012, 0.74814374, 0.69989368, 0.65475543, 0.61252828, 0.57302450, 0.53606844, 0.50149575, 0.46915275,
              0.43889571, 0.41058999, 0.38410981, 0.35933741, 0.33616264, 0.31448250, 0.29420058, 0.27522669, 0.25747650,
              0.24087108, 0.22533656, 0.21080392, 0.19720855, 0.18448998, 0.17259166, 0.16146070, 0.15104762, 0.14130610,
              0.13219284, 0.12366733, 0.11569165, 0.10823034, 0.10125025, 0.09472031, 0.08861151, 0.08289669, 0.07755043,
              0.07254897, 0.06787007, 0.06349292, 0.05939807, 0.05556731)

###############################################################################################
### End of functions
###############################################################################################

for (rep in rep.vec) {
  
  ###############################################################################################
  ### READ IN TRUE A
  ###############################################################################################
  
  matlabfile = paste0(codeFolder,"trueA/rep",rep,".mat")
  simstuff = readMat(matlabfile)
  trueA = simstuff$Astar
  
  #keep values that are positive and not tiny
  trueA[which(trueA<0.02)] = 0
  #scale so that minimum max neuron value is 1
  trueA = trueA / min(apply(trueA,2,max))
  
  ###############################################################################################
  ### GENERATE TRUE Z
  ###############################################################################################
  
  trueZ = t(replicate(ncol(trueA), generateCalciumVec(peakShape)))
  
  #save trueA and trueZ
  AZfile = paste0(folder,"data/trueAZ/rep",rep,"_")
  saveRDS(trueA, file=paste0(AZfile, "trueA.rds"))
  saveRDS(trueZ, file=paste0(AZfile, "trueZ.rds"))
  
  ###############################################################################################
  ### ADD NOISE TO THE DATA
  ###############################################################################################
  
  for (noiseIndex in 1:3) {
    SNR = SNR.vec[noiseIndex]
    SNR_static = SNR_static.vec[noiseIndex]
    
    set.seed(rep)
    
    #burst of noise every 50 frames
    peaks = seq(50, ncol(trueZ), by = 50)
    #list of spatially correlated noise maps
    noiseMaps = replicate(length(peaks), generateSpatialNoise(videoHeight), simplify = TRUE)
    #matrix of time trends
    timeTrends = sapply(peaks, generateTimeTrend, nFrames = ncol(trueZ), width = 75)
    timeTrends = timeTrends + matrix(rnorm(length(timeTrends), sd=0.05), nrow=nrow(timeTrends))
    spatialNoise = (noiseMaps %*% t(timeTrends))
    spatialNoise = spatialNoise/max(abs(spatialNoise))
    
    #add the spatially correlated noise to the video from Matlab
    newY = (trueA %*% trueZ) + spatialNoise/SNR
    
    ###############################################################################################
    ### SAVE RESULTS
    ###############################################################################################
    
    #save newY (as "Y_1.rds") for analysis with SCALPEL
    YfileFolder = paste0(folder,"data/finalY/rep",rep,"_SNR",SNR,"_",SNR_static,"/")
    if (!dir.exists(YfileFolder)) dir.create(YfileFolder)
    saveRDS(newY, file=paste0(YfileFolder, "Y_1.rds"))
    
  }
  for (noiseIndex in 4:7) {
    SNR = SNR.vec[noiseIndex]
    SNR_static = SNR_static.vec[noiseIndex]
    
    set.seed(rep)
    
    staticNoise = matrix(runif(nrow(trueA)*nrow(trueZ), min=-1, max=1), nrow=nrow(trueA), ncol=ncol(trueZ))
    
    #add the spatially correlated noise to the video from Matlab
    newY = (trueA %*% trueZ)+ staticNoise/SNR_static
    
    ###############################################################################################
    ### SAVE RESULTS
    ###############################################################################################
    
    #save newY (as "Y_1.rds") for analysis with SCALPEL
    YfileFolder = paste0(folder,"data/finalY/rep",rep,"_SNR",SNR,"_",SNR_static,"/")
    if (!dir.exists(YfileFolder)) dir.create(YfileFolder)
    saveRDS(newY, file=paste0(YfileFolder, "Y_1.rds"))
  }
  
}
