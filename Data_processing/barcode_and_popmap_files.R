
###################
# Creating stacks barcode file and population map
###################

 #Files containing sample names and plate positions as given to the CGRB

  keyfile1 <- read.csv("~/Documents/Columbine/Data/gbs_keyfile_pl1.csv", header = T, as.is = T)
  
  keyfile2 <- read.csv("~/Documents/Columbine/Data/gbs_keyfile_pl2.csv", header = T, as.is = T)

 #Create alpha-numeric sample names for the sake of bioinformatics ease
    keyfile1$Sample <- paste("s", keyfile1$Sample, sep = "")
    keyfile2$Sample <- paste("s", keyfile2$Sample, sep = "")
    
 #Create master list with both plates
    keyfile <- rbind(keyfile1, keyfile2)
    
 ## SANITY CHECK ##
    sum(duplicated(keyfile$Sample))


#### Post-CGRB processing
    
 #Pass-fail files from CGRB 
  # files include barcodes
    
    pf1 <- read.csv("~/Documents/Columbine/Data/gbs0239-passfail.csv", header = T, as.is = T) #Plate 1
    pf2 <- read.csv("~/Documents/Columbine/Data/gbs0240-passfail.csv", header = T, as.is = T) #Plate 2
    
  # alpha-numeric sample names
    
    pf1$Sample <- paste("s", pf1$Sample, sep = "")
    pf2$Sample <- paste("s", pf2$Sample, sep = "")

 # Create master list
  pf <- rbind(pf1, pf2)    

 ## SANITY CHECK ## 
  sum(duplicated(pf$Sample))
  
  all(pf$Sample == keyfile$Sample)
  all(pf1$Sample == keyfile1$Sample)
  all(pf2$Sample == keyfile2$Sample)
  
# Write out barcode files

  barcodes1 <- pf1[,c("Barcode", "Sample")]
  barcodes2 <- pf2[,c("Barcode", "Sample")]
  
  write.table(barcodes1, file = "~/Documents/Columbine/Data/barcodes1_round2.tsv", quote = F, row.names = F, sep = "\t", col.names = F)
  write.table(barcodes2, file = "~/Documents/Columbine/Data/barcodes2_round2.tsv", quote = F, row.names = F, sep = "\t", col.names = F)
  
 
  
  
  
   
#### Creating the population map for stacks
  
  samp_dat <- read.csv("~/Documents/Columbine/Data/sample_data_allplants.csv", header = T, as.is = T)
  
 #alpha-numeric sample name
  
  samp_dat$PLANT_ID <- paste("s", samp_dat$PLANT_ID, sep = "")

 # Dataset for just sequenced samples
  
  seqsamps <- merge(keyfile, samp_dat, by.x = "Sample", by.y = "PLANT_ID")

## SANITY CHECK ##
  
  all(seqsamps$Sample == pf$Sample[order(pf$Sample)])
  
  sum(duplicated(seqsamps))
  sum(duplicated(seqsamps$Sample))

# order by meadow within complex, then write out as popmap
  
 seqsamps <- seqsamps[order(seqsamps$COMPLEX, seqsamps$MEADOW_ID), ]

 write.table(seqsamps, file = "~/Documents/Columbine/Data/sequenced_samples.csv", row.names = F, quote = F, sep = ",")
 
 popmap <- data.frame(seqsamps$Sample, seqsamps$COMPLEX)

 write.table(popmap, file = "~/Documents/Columbine/Data/popmap_round2.tsv", row.names = F, col.names = F, sep = "\t", quote = F)
    
  
  
  
   
  
  
    