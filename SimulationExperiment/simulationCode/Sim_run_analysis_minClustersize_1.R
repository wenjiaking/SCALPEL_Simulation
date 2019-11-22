#runs SCALPEL for the simulated data

for (rep in rep.vec) {
  for (noiseIndex in 1:length(SNR.vec)) {
    SNR = SNR.vec[noiseIndex]
    SNR_static = SNR_static.vec[noiseIndex]
    
    #folder containg the raw data, saved as "Y_1.rds"
    rawDataFolder = paste0(folder, "data/finalY/rep",rep,"_SNR",SNR,"_",SNR_static,"/")
    
    #Create the following folder on your computer or specify an existing folder on your computer
    #in which the results of analyzing the simulated data should be saved
    outputFolder = paste0(folder, "scalpelResults/rep",rep,"_SNR",SNR,"_",SNR_static,"/")
    if (!dir.exists(outputFolder)) dir.create(outputFolder)
    
    #folder for saving estimated A and Z
    folderSaveAZ = paste0(folder, "scalpelResults/estimatedAZ/")
    
    ###############################################################################################
    ### RUN SCALPEL
    ###############################################################################################
    
    simOut = scalpel(outputFolder = outputFolder, rawDataFolder = rawDataFolder, 
                     videoHeight = videoHeight, minClusterSize = 1)
    
    ###############################################################################################
    ### SAVE RESULTS
    ###############################################################################################
    
    #save A hat and Z hat
    AZfile = paste0(folderSaveAZ,"rep",rep,"_SNR",SNR,"_",SNR_static,".RData")
    estA = simOut$Afilter
    estZ = simOut$Zhat
    save(estA, estZ, file=AZfile)
    
    rm(simOut)
  }
}
