## code to evaluate results from simulations
## inputs are true A, estimated A, estimated Z
## outputs are identification sensitivity and precision, and duplicate percentage

for (rep in rep.vec) {
  for (noiseIndex in 1:length(SNR.vec)) {
    SNR = SNR.vec[noiseIndex]
    SNR_static = SNR_static.vec[noiseIndex]
    
    print(paste0("replicate: ",rep))
    
    ###############################################################################################
    ### READ IN DATA TO EVALUATE
    ###############################################################################################
    
    #read in trueA
    trueA = readRDS(paste0(folder,"data/trueAZ/rep",rep,"_trueA.rds"))
    
    #read in estA, estZ
    Rfileresults = paste0(folder,"scalpelResults/estimatedAZ/rep",rep,"_SNR",SNR,"_",SNR_static,".RData")
    load(Rfileresults)
    
    ###############################################################################################
    ### ANALYSIS
    ###############################################################################################
    
    out = evaluateSims(trueA=trueA, estA=estA, estZ=estZ, cutoff_extra_pixels=cutoff_extra_pixels, 
                       cutoff_matching_pixels=cutoff_matching_pixels)
    
    ###############################################################################################
    ### SAVE OUTPUTS
    ###############################################################################################
    
    outfile = paste0(folder,"evalFiles/sensitivityPrecision/scalpel_minClustersize_1_rep",rep,"_SNR",SNR,"_",SNR_static,"_matchcut",cutoff_matching_pixels,"_extracut",cutoff_extra_pixels,".txt")
    #save identification sensitivity/precision, duplicate percentage to a text file
    write.table(data.frame(id.sens=out$ID.sens, id.prec=out$ID.prec, dupl.pct=out$dupl.pct), file = outfile, row.names = F)
    
    outfile = paste0(folder,"evalFiles/matching/scalpel_minClustersize_1_rep",rep,"_SNR",SNR,"_",SNR_static,"_matchcut",cutoff_matching_pixels,"_extracut",cutoff_extra_pixels,".txt")
    #save matching of estimated neurons to true ones
    write.table(out$m, file = outfile, row.names = F, col.names = F)
  }
}