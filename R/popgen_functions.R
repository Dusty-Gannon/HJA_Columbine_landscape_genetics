###
# Functions to compute basic population genetic statistics
###

# packages
  require(tidyverse)
  require(hierfstat)

# --------- Observed heterozygosity -------------

## --- Individual level ---

# takes in allele loading matrix and gives
# individual heterozygosity for each sample. 
Ho_ind <- function(L, rows=T){
  if(isTRUE(rows)){
    apply(
      L,
      1,
      function(x){mean(x==1)}
    )
  } else {
    apply(
      2,
      function(x){mean(x==1)}
    )
  }
}


## --- Population level ---
Ho <- function(X){
  
  # get proportion of homozygotes
  onesite_ho <- function(v){
    1-mean(v==0) - mean(v==2)
  }
  
  # estimate this proportion over all loc
  names(X)[1] <- "pop"
  pops <- unique(X[,1])
  
  persite <- map_df(
    pops,
      ~apply(
        subset(X, pop==.x)[,-1],
        2,
        onesite_ho
      )
  )
  
  ho <- apply(persite, 1, mean)
  names(ho) <- pops
  return(ho)
}





# --------- Expected heterozygosity -------------

He <- function(X, bias_correct=T){
  # rename for consistency
  names(X)[1] <- "pop"
  
  # get dimensions
  L <- dim(X)[2]-1
  m <- dim(X)[1]
  
  # function to get maf for each locus
  maf <- function(df){
    v <- apply(df, 2, function(x){sum(x)/(2*length(x))})
    return(
      as.data.frame(
        matrix(data = v, nrow = 1)
      )
    )
  }
  
  # function to center a matrix based on p-hats
  cent_mat <- function(M, phat){
    vlist <- map(
      1:length(phat),
      ~ as.double(M[,.x]) - 2*phat[.x]
    )
    v <- unlist(vlist)
    return(
      matrix(
        data = v,
        nrow = nrow(M),
        ncol = ncol(M)
      )
    )
  }
  
  # get minor allele frequencies by population
  mafs_bypop <- group_by(X,pop) %>%
    group_modify(~maf(.))
  
  # Get pop info
  pops <- unique(X$pop)
  popblocks <- map(pops, ~which(X$pop==.x))
  
  # calculate little h for each locus
  by_loci_h <- apply(
    mafs_bypop[,-1],
    2,
    function(x){
      (1- x^2 - (1-x)^2)
    }
  )
  
  if(isTRUE(bias_correct)){
  
    # center matrices
    Xmcs <- map(
      1:length(pops),
      ~cent_mat(
        as.matrix(X[popblocks[[.x]],-1]),
        as.double(mafs_bypop[.x,-1])
      )
    )
  
    # compute denominators
    denoms <- map_dbl(
      1:length(pops),
      ~2*sum(
        as.double(mafs_bypop[.x,-1])*(1-as.double(mafs_bypop[.x,-1]))
      )
    )
  
    # compute kinship matrices for each pop
    kin_mats <- map(
      1:length(pops),
      ~(Xmcs[[.x]]%*%t(Xmcs[[.x]]))/denoms[.x]
    )
  
    # compute average kinship for pairs of individuals
    phi_bars <- map_dbl(
      kin_mats,
      ~mean(.x)
    )
  
  # compute the gene diversity for each pop
    return(
      apply(by_loci_h, 1, mean)/(1-phi_bars)
    )
  } else{
    
    # sample sizes
    ns <- map_dbl(popblocks, ~length(.x))
    return(
      apply(by_loci_h,1,mean)*(2*ns/(2*ns-1))
    )
    
  }

}




