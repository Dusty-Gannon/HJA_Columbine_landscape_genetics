Preparing landscape genetic data by individuals
================
D. G. Gannon
October 2020

## Load relevant data

#### Load genetic data

Load vcf file with SNPs from presumably non-coding regions and convert
genotype matrix into an allelic frequency matrix ![{\\bf
F}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20F%7D "{\\bf F}"),
where ![f\_{ik}\\in
\\{0,0.5,1\\}](https://latex.codecogs.com/png.latex?f_%7Bik%7D%5Cin%20%5C%7B0%2C0.5%2C1%5C%7D
"f_{ik}\\in \\{0,0.5,1\\}") for individual
![i](https://latex.codecogs.com/png.latex?i "i") and locus
![k](https://latex.codecogs.com/png.latex?k "k"): 1 = homozygous for the
randomly selected reference allele; 0.5 = heterozygous; 0 = homozygous
for the alternative allele (Bradburd, Coop, and Ralph 2018).

``` r
  col_vcf <- read.vcfR(here("Data", "AQFO_snps_nc.vcf"))

# Extract genotype matrix
  col_gt <- t(extract.gt(col_vcf, element = "GT", as.numeric = F))
  
# Function to convert to allelic covariance
  allele_freq_random_ref <-  function(x, sep="\\|"){
    spl <- str_split(x, pattern = sep, n=2, simplify = T)
    spl_num <- apply(spl, 2, as.numeric)
    random_ref <- as.numeric(rbernoulli(1,p=0.5))
    apply(spl_num, 1, function(x){mean(x==random_ref)})
  }
  
  col_alfreq <- apply(col_gt,2, allele_freq_random_ref)
  rownames(col_alfreq) <- rownames(col_gt)
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

#### Create allelic scatter matrix ![S](https://latex.codecogs.com/png.latex?S "S").

``` r
  col_alfreq_cent <- col_alfreq - 0.5
  col_alscatter <- (col_alfreq_cent%*%t(col_alfreq_cent))
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

#### Load field data collected on individual plants and meadows

``` r
# convert to dataframe
  col_alfreq_df <- as.data.frame(col_alfreq) 
    col_alfreq_df <- cbind(rownames(col_alfreq_df), col_alfreq_df)
    names(col_alfreq_df)[1] <- "Sample"

# load sample data
  samp_dat <- read.table(here("Data", "sequenced_samples.txt"), 
                              header = T, as.is = T, sep = "\t")

# subset for columns of interest
  samp_dat_sub <- samp_dat[, c("Sample", "COMPLEX", "MEADOW_ID",
                               "MEADOW_AREA", "PROP_FOREST_100",
                               "X_COORD", "Y_COORD", "PL_IN_5M", "COVER",
                               "FLW_STATUS")]
  
 # merge population variables with L matrix
  col_popgen_data <- merge(samp_dat_sub[,c("Sample", "COMPLEX", "MEADOW_ID")], 
                           col_alfreq_df)
  
# sort by complex then meadow
  col_popgen_data <- col_popgen_data[order(col_popgen_data$COMPLEX, 
                                           col_popgen_data$MEADOW_ID), ]
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

## Compute landscape variables of interest

#### Amount of intervening forest

**Create spatial lines object**

The spatial lines object will include a line between all pairs of
individuals.

``` r
# number of meadows
  num_plants <- nrow(col_popgen_data)

# create list to store start-end coords
  node_connects <- vector(mode = "list")

for(i in 1:num_plants){
  if(i < num_plants){
    for(j in (i+1):num_plants){
      coords <- as.matrix(rbind(samp_dat_sub[i, c("X_COORD","Y_COORD")],
                                samp_dat_sub[j, c("X_COORD","Y_COORD")]))
      nm <- paste(samp_dat_sub$Sample[i], samp_dat_sub$Sample[j],
                sep = "-")
      node_connects[[nm]] <- coords
    }
  }
}  
  
#convert to lines
  node_connects_l <- map(node_connects, ~Line(.))
  
# Convert the many lines to a single lines object

  node_connects_ls <- map(1:length(node_connects_l), 
                          ~Lines(node_connects_l[[.]], ID=names(node_connects_l)[.]))
  
# define CRS
  hja_crs <- CRS("+proj=utm +zone=10 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# create one spatial lines object
  node_connects_spl <- SpatialLines(node_connects_ls, proj4string = hja_crs)
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

**Extract the amount of forest between each node**

``` r
  forest <- raster(here("Data","forest.tif"))
  
# extract proportion of forest between the nodes
  prop_for <- extract(forest, node_connects_spl, fun=mean, na.rm=T,
                      df=T)
  prop_for$ID <- names(node_connects_spl)
```

#### Pairwise covariate data

**Create a pairwise list for relevant data**

``` r
  pw_data <- cbind(str_split(prop_for$ID, pattern = "-",
                                       n=2, simplify = T), 
                   prop_for)
  
  names(pw_data) <- c("plant_i", "plant_j", "pair_ij", "prop_forest_ij")
  
# merge in covariates for meadow_i
  covs <- c("Sample", "COMPLEX", "MEADOW_ID", "MEADOW_AREA",
            "PROP_FOREST_100", "PL_IN_5M",
            "COVER", "X_COORD", "Y_COORD")
  pw_data_i <- merge(pw_data, 
                   samp_dat_sub[,colnames(samp_dat_sub) %in% covs],
                   by.x="plant_i", by.y="Sample", sort=F)
  
# rename with "_i" appended to names
  names(pw_data_i)[which(names(pw_data_i) %in% covs)] <- 
    paste(names(pw_data_i)[which(names(pw_data_i) %in% covs)], 
          "i", sep = "_")
  
# merge in covariates for meadow_j
   pw_data_ij <- merge(pw_data_i, 
                       samp_dat_sub[,colnames(samp_dat_sub) %in% covs],
                  by.x="plant_j", by.y="Sample", sort=F)

# rename with "_j" appended to names   
  names(pw_data_ij)[which(names(pw_data_ij) %in% covs)] <- 
    paste(names(pw_data_ij)[which(names(pw_data_ij) %in% covs)], 
          "j", sep = "_")
  
# reorder columns, then sort
  pw_data_ij <- pw_data_ij[,c(2,1,3:ncol(pw_data_ij))]
  pw_data_ij <- pw_data_ij[order(pw_data_ij$MEADOW_ID_i,
                                 pw_data_ij$plant_i,
                                 pw_data_ij$MEADOW_ID_j,
                                 pw_data_ij$plant_j), ]
  
# sanity check
  all.equal(rownames(col_alfreq)[1:191],
            unique(pw_data_ij$plant_i))
  
# geographic distanct
  pw_data_ij <- pw_data_ij %>% 
    mutate(
      dist_ij = sqrt((X_COORD_i - X_COORD_j)^2 +
        (Y_COORD_i - Y_COORD_j)^2)/1000
    )
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

**Compute the geographic distance between all pairs of nodes**

``` r
# create coordinate matrix
  coords <- as.matrix(samp_dat_sub[,names(samp_dat_sub) %in% 
                                     c("X_COORD", "Y_COORD")])
  rownames(coords) <- samp_dat_sub$Sample

# define function to compute squared distance matrix
  gdist2 <- function(mat){
    ones <- rep(1, nrow(mat))
    M <- mat%*%t(mat)
    return(ones%*%t(diag(M)) - 2*M + diag(M)%*%t(ones))
  }
  
# compute distance matrix
  gdist_mat <- sqrt(gdist2(coords))/1000
  rownames(gdist_mat) <- colnames(gdist_mat)
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

#### Stack covariate matrices for covariate array, ![{\\bf X}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20X%7D "{\\bf X}")

The covariate array, ![X\_{(n\\times n \\times
P)}](https://latex.codecogs.com/png.latex?X_%7B%28n%5Ctimes%20n%20%5Ctimes%20P%29%7D
"X_{(n\\times n \\times P)}"), is the 3-dimensional array with
![P](https://latex.codecogs.com/png.latex?P "P"), ![n \\times
n](https://latex.codecogs.com/png.latex?n%20%5Ctimes%20n "n \\times n")
matrix slices, where ![n](https://latex.codecogs.com/png.latex?n "n") is
the number of nodes (samples), that contain covariate data of
![P](https://latex.codecogs.com/png.latex?P "P") covariates. For
example, one slice is the geographic distance matrix between samples

``` r
# for the symmetric matrices (distance and proportion intervening forest),
#  write function to fill in matrix based on vector of values
  fillMat_sym <- function(x, rows, cols, upper=T){
      mat <- matrix(data = 0, nrow = rows, ncol = cols)
      if(isTRUE(upper)){
        mat[upper.tri(mat, diag = F)] <- x
      } else{
        mat[lower.tri(mat, diag = F)] <- x
      }
      return(mat + t(mat))
  }

# character vector of column names
  covs_list <- c("intercept", "geo_dist", "prop_forest_ij", "plant_density_i",
                 "cover_i", "meadow_area_i","forest_100m_i",
                 "plant_density_j", "cover_j", "meadow_area_j",
                 "forest_100m_j")

# create empty array
  X <- array(dim = c(nrow(gdist_mat), ncol(gdist_mat), length(covs_list)),
             dimnames = list(rownames(gdist_mat),
                             colnames(gdist_mat),
                             covs_list))
 
# fill first three slices (already in matrix form)
  X[,,1] <- matrix(1, nrow = dim(X)[1], ncol = dim(X)[2])
  X[,,2] <- gdist_mat
  X[,,3] <- fillMat_sym(pw_data_ij$prop_forest_ij, rows = dim(X)[1],
                        cols = dim(X)[2], upper = F)
  
# Sanity check that the matrix was filled correctly
  rand_row <- sample(1:nrow(pw_data_ij),1)
  row_check <- which(dimnames(X)[[1]] == pw_data_ij$plant_i[rand_row])
  col_check <- which(dimnames(X)[[2]] == pw_data_ij$plant_j[rand_row])
  
  all.equal(X[row_check, col_check,3],
            pw_data_ij$prop_forest_ij[rand_row])
  all.equal(X[row_check, col_check,2],
            pw_data_ij$dist_ij[rand_row])
  
  all.equal(samp_dat_sub$Sample,
            dimnames(X)[[1]])
  
# create remaining slices
  
 # which are covariates pertaining to the "sending" meadow?
  rem_slices <- c("PL_IN_5M", "COVER",
                  "MEADOW_AREA", "PROP_FOREST_100")
  
 # Loop through these, create the matrix
 #    M_k with M_k(ij)=M_k(im) for all j and m
  for(k in 1:length(rem_slices)){
    cov_k <- which(names(samp_dat_sub) == rem_slices[k])
    temp_mat <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
    temp_mat <- sweep(temp_mat,1,
                      samp_dat_sub[,cov_k], FUN = "+")
    X[,,k+3] <- temp_mat
  }
  
# Loop through these, create the matrix
#    M_k with M_k(ij)=M_k(im) for all j and m
  for(k in 1:length(rem_slices)){
    cov_k <- which(names(samp_dat_sub) == rem_slices[k])
    temp_mat <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
    temp_mat <- sweep(temp_mat,2,
                      samp_dat_sub[,cov_k], FUN = "+")
    X[,,k+3+length(rem_slices)] <- temp_mat
  }  
  
# put things on smaller scales
  X[,,"cover_i"] <- X[,,"cover_i"]/100
  X[,,"cover_j"] <- X[,,"cover_j"]/100
  X[,,"meadow_area_i"] <- log(X[,,"meadow_area_i"])
  X[,,"meadow_area_j"] <- log(X[,,"meadow_area_j"])
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

**Save all these data as they are formatted now**

``` r
  save(samp_dat_sub, pw_data_ij, 
       col_popgen_data, X, col_alscatter,
       file = here("Data", "Indiv_lndscp_gen.RData"))
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

**References**

<div id="refs" class="references">

<div id="ref-bradburd2018a">

Bradburd, Gideon S., Graham M. Coop, and Peter L. Ralph. 2018.
“Inferring Continuous and Discrete Population Genetic Structure Across
Space.” *Genetics* 210 (1): 33–52.
<https://doi.org/10.1534/genetics.118.301333>.

</div>

</div>
