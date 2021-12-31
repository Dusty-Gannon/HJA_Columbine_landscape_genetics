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




# fill a symmetric matrix based on a vector of the diagonal elements

fillMat_sym <- function(x, rows, cols, upper=T){
  mat <- matrix(data = 0, nrow = rows, ncol = cols)
  if(isTRUE(upper)){
    mat[upper.tri(mat, diag = F)] <- x
  } else{
    mat[lower.tri(mat, diag = F)] <- x
  }
  return(mat + t(mat))
}  
  



# Compute a squared Euclidean distance matrix for a matrix of coordinates
  
gdist2 <- function(mat){
  ones <- rep(1, nrow(mat))
  M <- mat%*%t(mat)
  return(ones%*%t(diag(M)) - 2*M + diag(M)%*%t(ones))
}



