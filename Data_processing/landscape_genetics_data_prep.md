Preparing landscape genetic data
================
D. G. Gannon
October 2020

## Load relevant data

#### Load genetic data

Load vcf file with SNPs from presumably non-coding regions and convert
genotype matrix into allelic frequency matrix ![{\\bf
F}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20F%7D "{\\bf F}"),
where ![f\_{ik}\\in
\\{0,0.5,1\\}](https://latex.codecogs.com/png.latex?f_%7Bik%7D%5Cin%20%5C%7B0%2C0.5%2C1%5C%7D
"f_{ik}\\in \\{0,0.5,1\\}") for individual
![i](https://latex.codecogs.com/png.latex?i "i") and locus
![k](https://latex.codecogs.com/png.latex?k "k"): 1 = homozygous for the
randomly selected reference allele; 0.5 = heterozygous; 0 = homozygous
for the alternative allele.

``` r
  col_vcf <- read.vcfR("~/Documents/Columbine/Data/GBS_Data/AQFO_snps_nc.vcf")

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

#### Create allelic covariance matrix ![\\hat{\\boldsymbol \\Omega}](https://latex.codecogs.com/png.latex?%5Chat%7B%5Cboldsymbol%20%5COmega%7D "\\hat{\\boldsymbol \\Omega}").

``` r
  col_alfreq_cent <- col_alfreq - 0.5
  col_alcov <- (col_alfreq_cent%*%t(col_alfreq_cent))/ncol(col_alfreq_cent)
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

#### Load field data collected on individual plants and meadows

``` r
# convert to dataframe
  col_alfreq_df <- as.data.frame(col_alfreq) 
    col_alfreq_df <- cbind(rownames(col_alfreq_df), col_alfreq_df)
    names(col_alfreq_df)[1] <- "Sample"

# load sample data
  samp_dat <- read.table("~/Documents/Columbine/Data/sequenced_samples.txt", header = T, as.is = T, sep = "\t")

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

#### Compute “node” data

We define node locations by averaging the coordinates of plants sampled
in each subpopulation (meadow).

``` r
# find node centers
  node_data <- group_by(samp_dat_sub, MEADOW_ID) %>% 
                  summarise(., n=n(), node_x=mean(X_COORD), node_y=mean(Y_COORD),
                            plant_density=mean(PL_IN_5M), mean_cover=mean(COVER),
                            meadow_area=mean(MEADOW_AREA),
                            forest_100m=mean(PROP_FOREST_100))
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

## Compute landscape variables of interest

#### Amount of intervening forest

**Create spatial lines object**

The spatial lines object will include a line between all pairs of
“nodes”.

``` r
# number of meadows
  num_meads <- nrow(node_data)

# create list to store start-end coords
  node_connects <- vector(mode = "list")

for(i in 1:num_meads){
  if(i < num_meads){
    for(j in (i+1):num_meads){
      coords <- as.matrix(rbind(node_data[i, c("node_x","node_y")],
                                node_data[j, c("node_x","node_y")]))
      nm <- paste(node_data$MEADOW_ID[i], node_data$MEADOW_ID[j],
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
  forest <- raster("~/Documents/Columbine/Data/GIS_Data/forest.tif")
  
# extract proportion of forest between the nodes
  prop_for <- extract(forest, node_connects_spl, fun=mean, na.rm=T, buffer=100,
                      df=T)
  prop_for$ID <- names(node_connects_spl)
```

#### Pairwise covariate data

**Create a pairwise list for relevant data**

``` r
  pw_data <- cbind(str_split(prop_for$ID, pattern = "-",
                                       n=2, simplify = T), 
                   prop_for)
  
  names(pw_data) <- c("meadow_i", "meadow_j", "pair_ij", "prop_forest_ij")
  
# merge in covariates for meadow_i
  covs <- c("MEADOW_ID", "plant_density",
            "mean_cover", "meadow_area",
            "forest_100m")
  pw_data <- merge(pw_data, unique(node_data[,colnames(node_data) %in% 
                                      covs]),
                   by.x="meadow_i", by.y="MEADOW_ID")
  
# rename with "_i" appended to names
  names(pw_data)[which(names(pw_data) %in% covs)] <- 
    paste(names(pw_data)[which(names(pw_data) %in% covs)], 
          "i", sep = "_")
  
# merge in covariates for meadow_j
   pw_data <- merge(pw_data, unique(node_data[,colnames(node_data) %in% 
                                      covs]),
                  by.x="meadow_j", by.y="MEADOW_ID")

# rename with "_j" appended to names   
  names(pw_data)[which(names(pw_data) %in% covs)] <- 
    paste(names(pw_data)[which(names(pw_data) %in% covs)], 
          "j", sep = "_")
  
# reorder columns, then sort
  pw_data <- pw_data[,c(2,1,3:ncol(pw_data))]
  pw_data_byrow <- pw_data[order(pw_data$meadow_i, pw_data$meadow_j), ]
  pw_data_bycol <- pw_data[order(pw_data$meadow_j, pw_data$meadow_i), ]
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

**Compute the geographic distance between all pairs of nodes**

``` r
# create coordinate matrix
  coords <- as.matrix(node_data[,names(node_data)%in%c("node_x", "node_y")])
  rownames(coords) <- node_data$MEADOW_ID

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

The covariate array, ![X\_{(m\\times m \\times
P)}](https://latex.codecogs.com/png.latex?X_%7B%28m%5Ctimes%20m%20%5Ctimes%20P%29%7D
"X_{(m\\times m \\times P)}"), is the 3-dimensional array with ![m
\\times m](https://latex.codecogs.com/png.latex?m%20%5Ctimes%20m
"m \\times m") matrix slices, where
![m](https://latex.codecogs.com/png.latex?m "m") is the number of nodes
(meadows), that contain covariate data of
![P](https://latex.codecogs.com/png.latex?P "P") covariates. For
example, one slice is the geographic distance matrix between nodes.

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
                 "mean_cover_i", "meadow_area_i","forest_100m_i",
                 "plant_density_j", "mean_cover_j", "meadow_area_j",
                 "forest_100m_j")

# create empty array
  X <- array(dim = c(nrow(gdist_mat), ncol(gdist_mat), length(covs_list)),
             dimnames = list(rownames(gdist_mat),
                             colnames(gdist_mat),
                             covs_list))
 
# fill first three slices (already in matrix form)
  X[,,1] <- matrix(1, nrow = dim(X)[1], ncol = dim(X)[2])
  X[,,2] <- gdist_mat
  X[,,3] <- fillMat_sym(pw_data_bycol$prop_forest, rows = dim(X)[1],
                        cols = dim(X)[2], upper = T)
  
# create remaining slices
  
# ensure node_data is sorted correctly
  node_data_sort <- node_data[order(node_data$MEADOW_ID),]
  
# which are covariates pertaining to the "sending" meadow?
  swp_rows <- str_which(covs_list, regex("_i$"))
  
# Loop through these, create the matrix
#    M_k with M_k(ij)=M_k(im) for all j and m
  for(k in swp_rows){
    cov_k <- which(names(node_data_sort) ==
                   str_replace(covs_list[k], "_i", ""))
    temp_mat <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
    temp_mat <- sweep(temp_mat,1,
                      unlist(node_data_sort[,cov_k]), FUN = "+")
    X[,,k] <- temp_mat
  }
  
# which are covariates pertaining to the "receiving" meadow?
  swp_cols <- str_which(covs_list, regex("_j$"))
  
# Loop through these, create the matrix
#    M_k with M_k(ij)=M_k(im) for all j and m
  for(k in swp_cols){
    cov_k <- which(names(node_data_sort) ==
                   str_replace(covs_list[k], "_j", ""))
    temp_mat <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
    temp_mat <- sweep(temp_mat,2,
                      unlist(node_data_sort[,cov_k]), FUN = "+")
    X[,,k] <- temp_mat
  }  
  
# put things on smaller scales
  X[,,"mean_cover_i"] <- X[,,"mean_cover_i"]/100
  X[,,"mean_cover_j"] <- X[,,"mean_cover_j"]/100
  X[,,"meadow_area_i"] <- log(X[,,"meadow_area_i"])
  X[,,"meadow_area_j"] <- log(X[,,"meadow_area_j"])
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

## Compute squared genetic distance matrix ![{\\bf D}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20D%7D "{\\bf D}")

``` r
  col_aldc_sort <- as.matrix(col_popgen_data[,-c(1:3)])
  row.names(col_aldc_sort) <- col_popgen_data$Sample
  
  D <- gdist2(col_aldc_sort)
```

![\~](https://latex.codecogs.com/png.latex?~ "~")

## Define the ![K](https://latex.codecogs.com/png.latex?K "K") matrix

The ![K](https://latex.codecogs.com/png.latex?K "K") matrix maps the
![n](https://latex.codecogs.com/png.latex?n "n") individual samples
(plants) to the ![m](https://latex.codecogs.com/png.latex?m "m") nodes
(meadows), for which we estimate the precision matrix ![{\\bf
Q}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20Q%7D "{\\bf Q}").
Thus, the scale matrix for the Wishart model, ![\\boldsymbol
\\Psi\_{(n\\times
n)}](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5CPsi_%7B%28n%5Ctimes%20n%29%7D
"\\boldsymbol \\Psi_{(n\\times n)}"), becomes

  
![&#10;\\boldsymbol \\Psi = {\\bf K}{\\bf Q}^{-1}{\\bf K}' +
\\tau^2{\\bf
I},&#10;](https://latex.codecogs.com/png.latex?%0A%5Cboldsymbol%20%5CPsi%20%3D%20%7B%5Cbf%20K%7D%7B%5Cbf%20Q%7D%5E%7B-1%7D%7B%5Cbf%20K%7D%27%20%2B%20%5Ctau%5E2%7B%5Cbf%20I%7D%2C%0A
"
\\boldsymbol \\Psi = {\\bf K}{\\bf Q}^{-1}{\\bf K}' + \\tau^2{\\bf I},
")  

where ![{\\bf K}\_{(n\\times
m)}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20K%7D_%7B%28n%5Ctimes%20m%29%7D
"{\\bf K}_{(n\\times m)}") is a matrix of
![n](https://latex.codecogs.com/png.latex?n "n") rows for the
![n](https://latex.codecogs.com/png.latex?n "n") samples and
![m](https://latex.codecogs.com/png.latex?m "m") columns for the
![m](https://latex.codecogs.com/png.latex?m "m") nodes, ![{\\bf
I}\_{(n\\times
n)}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20I%7D_%7B%28n%5Ctimes%20n%29%7D
"{\\bf I}_{(n\\times n)}") is the
rank-![n](https://latex.codecogs.com/png.latex?n "n") identity matrix,
and ![\\tau^2](https://latex.codecogs.com/png.latex?%5Ctau%5E2
"\\tau^2") is the process variation.

``` r
# define empty matrix
  K <- matrix(0, nrow = dim(D)[1], ncol = dim(X)[1])
  colnames(K) <- node_data_sort$MEADOW_ID
  rownames(K) <- col_popgen_data$Sample
  
# fill it in with some ones
  for(i in 1:nrow(node_data_sort)){
    temp.df <- subset(col_popgen_data[,c(1:3)], 
                      MEADOW_ID == node_data_sort$MEADOW_ID[i])
    mappings <- which(col_popgen_data$Sample %in% temp.df$Sample)
    K[mappings,i] <- 1
    
  }
```

**Save all these data as they are formatted now**

``` r
  save(samp_dat_sub, node_data_sort, 
       col_popgen_data, pw_data, X, K, col_alcov,
       file = "~/Documents/Columbine/Data/GBS_Data/spatial_weights_estimation.RData")
```
