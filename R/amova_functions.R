# Functions used in the AMOVA-like analysis used to
# test for altered selection on floral trait genes

require(tidyverse)

# Function to convert the "gt" object of a vcfR object
# to a dosage and loading matrix

allele_load <- function(vcf, sep="\\|"){
  # extract character matrix
  gt <- vcf@gt[,-1]
  ld <- apply(
    gt,
    2,
    FUN = function(v){
      spl <- str_split(v, pattern = sep, n=2, simplify = T)
      spl_num <- apply(spl, 2, as.numeric)
      return(spl_num[,1] + spl_num[,2])
    }
  )
  rownames(ld) <- vcf@fix[,"ID"]
  return(ld)
}

# Function to compute (potentially weighted) (n x L) genetic distance matrix
## Function takes a dosage matrix with n individuals in rows and 
## genotypes at L loci in columns

genDist <- function(dose_mat, W=NULL){
  D <- matrix(data = 0, nrow = dim(dose_mat)[1], ncol = dim(dose_mat)[1])
  if(is.null(W)){
    W <- diag(x=1, nrow = dim(dose_mat)[1], ncol = dim(dose_mat)[1])
  }
  ## Fill distance matrix with squared, weighted distances
  for(i in 1:nrow(D)){
    for(j in 1:nrow(D)){
      g_i <- as.double(dose_mat[i,])
      g_j <- as.double(dose_mat[j,])
      D[i,j] <- (g_i - g_j)%*%W%*%(g_i-g_j)
    }
  }
  return(D)
}


# Function to compute MSD among pops within groups
# based on Excoffier et al. (1992)
#
## The function takes in a distance matrix and a dataframe
## that defines the hypothesized structure of individuals (col 1)
## into populations (col 2) and groups (col 3) to compute the variance
## among populations for each group

sigma_hat_each_grp <- function(dosage, W=NULL){
  # create distance matrix 
  names(dosage)[1:3] <- c("samp", "pop", "grp")
  
  # create numeric pop and grp indicators
  dosage$pop <- factor(dosage$pop, levels = unique(dosage$pop))
  dosage$grp <- factor(dosage$grp, levels = unique(dosage$grp))
  pops_num <- as.numeric(dosage$pop)
  grps_num <- as.numeric(dosage$grp)
  
  # create non-nummeric list of groups
  grps <- unique(dosage$grp)
  
  # create list of row+column ids for grouping individuals into pops
  pop_blocks <- map(
    1:length(unique(pops_num)),
    ~which(pops_num == .x)
  )
  # create list of row/column ids for grouping individuals into groups
  grp_blocks <- map(
    1:length(unique(grps_num)),
    ~which(grps_num==.x)
  )
  # create list of pops in each group
  pops_in_grps <- map(
    1:length(unique(grps)),
    ~ as.numeric(
      unique(
      dosage[
        which(dosage$grp==grps[.x]),
      ]$pop
      )
    )
  )
  
  # create dosage matrix
  dose_mat <- as.matrix(dosage[,-c(1:3)])
  
  # create distance matrix
  if(is.null(W)){
    W <- diag(nrow = ncol(dose_mat), ncol = ncol(dose_mat))
  }
  D <- genDist(dose_mat, W=W)
  
  # SSD within populations
  ssd_eachpop <- map_dbl(
    pop_blocks,
    ~sum(D[.x,.x])/(2*length(.x))
  )
  
  # SSD within groups
  ssd_eachgrp <- map_dbl(
    grp_blocks,
    ~sum(D[.x,.x])/(2*length(.x))
  )
  
  # SSD among pops in each group
  ssd_pop_in_grp <- map_dbl(
    1:length(unique(grps)),
    ~ssd_eachgrp[.x] - 
      sum(ssd_eachpop[pops_in_grps[[.x]]])
  )
  
  # MSD(AP/WG)
  msd_ap_wg <- ssd_pop_in_grp/map_dbl(
    pops_in_grps,
    ~(length(.x)-1)
  )
  
  # n coefficient
  n_igs <- map_dbl(
    pop_blocks,
    ~ length(.x)
  )
  squared_terms <- map_dbl(
    1:length(grps),
    ~ sum(n_igs[pops_in_grps[[.x]]]^2)/sum(n_igs[pops_in_grps[[.x]]])
  )
  ns <- map_dbl(
    1:length(grps),
    ~ (sum(n_igs[pops_in_grps[[.x]]]) - sum(squared_terms[.x]))/
      (length(pops_in_grps[[.x]])-1)
  ) 
  
  # compute sigma_sq parameters
  sigma_sq_1 <- map_dbl(
    1:length(grps),
    ~ sum(ssd_eachpop[pops_in_grps[[.x]]])/
      (sum(n_igs[pops_in_grps[[.x]]]) - length(pops_in_grps[[.x]]))
      
  )
  sigma_sq_2 <- map_dbl(
    1:length(grps),
    ~ (msd_ap_wg[.x] - sigma_sq_1[.x])/ns[.x]
  )
  
  names(sigma_sq_2) <- as.character(grps)
  return(sigma_sq_2)
  
}


# This function combines the previous two functions to
# permute the populations within groups (holding individuals
# fixed within populations) to perform a permutational test of
# significance

perm_group_test <- function(dosage, nperm=100, W=NULL){
  # create list with results
  results <- vector(mode = "list")
  # compute observed values
  results[["obs"]] <- sigma_hat_each_grp(dosage = dosage, W=W) 
  
  # now permute group assignments to pops
  ## reassign pops to groups
  names(dosage)[1:3] <- c("samp", "pop", "grp")
  pops_and_groups <- unique(dosage[,c(2,3)])
  
  # empty matrix of resampled results
  resamp_mat <- matrix(nrow = nperm, ncol = length(unique(dosage$grp)))
  colnames(resamp_mat) <- names(results[["obs"]])
  
  # dosage without samps or groups to merge after shuffle
  dosage_popsonly <- dosage[,-c(1,3)]
  
  #progress bar 
  pb <- txtProgressBar(
    min = 0,
    max = 10,
    initial = 0,
    char = "="
  )
  
  # permute populations
  for(i in 1:nperm){
    pngrps_i <- pops_and_groups
    pngrps_i$grp <- 
      pngrps_i$grp <- pngrps_i$grp[sample(1:nrow(pngrps_i))]
    
    dose_i <- merge(
      pngrps_i,
      dosage_popsonly,
    )
    
    dose_i <- dose_i[
      order(dose_i$grp, dose_i$pop),
    ]
    
    # Add place-holder sample name
    fake_samps <- data.frame(
      samp = paste("s", 1:nrow(dose_i), sep = "")
    )
    dose_i <- cbind(fake_samps, dose_i)
    
    resamp_mat[i,] <- sigma_hat_each_grp(dose_i, W=W)
    setTxtProgressBar(pb, value = round(i/(nperm/10)))
  }
  close(pb)
  
  results[["permutations"]] <- resamp_mat
  
  return(results)
  
}












