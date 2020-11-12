Estimating edge weights connecting HJA columbine genetic networks
================
D. G. Gannon
October 2020

### Data

We have SNP data on
![n=192](https://latex.codecogs.com/png.latex?n%3D192 "n=192")
*Aquilegia formosa* individuals from 25 meadows, which we consider to be
“sub-populations” of the H.J. Andrews *A. formosa* population. We work
with the allelic scatter matrix ![{\\bf
S}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20S%7D "{\\bf S}")
computed from an allelic frequency matrix ![F\_{(n\\times
\\ell)}](https://latex.codecogs.com/png.latex?F_%7B%28n%5Ctimes%20%5Cell%29%7D
"F_{(n\\times \\ell)}"), where
![\\ell](https://latex.codecogs.com/png.latex?%5Cell "\\ell") is the
number of SNP loci sequenced and
![n](https://latex.codecogs.com/png.latex?n "n") is the number of plants
sequenced. Allelic frequencies are coded as ![f\_{ik}
= 0](https://latex.codecogs.com/png.latex?f_%7Bik%7D%20%3D%200
"f_{ik} = 0") if individual ![i](https://latex.codecogs.com/png.latex?i
"i") is homozygous for the randomly selected reference allele at locus
![k](https://latex.codecogs.com/png.latex?k "k"),
![f\_{ik}=0.5](https://latex.codecogs.com/png.latex?f_%7Bik%7D%3D0.5
"f_{ik}=0.5") if individual ![i](https://latex.codecogs.com/png.latex?i
"i") has a heterozygous genotype at site
![k](https://latex.codecogs.com/png.latex?k "k"), and
![f\_{ik}=1](https://latex.codecogs.com/png.latex?f_%7Bik%7D%3D1
"f_{ik}=1") if individual ![i](https://latex.codecogs.com/png.latex?i
"i") is homozygous with two copies of the alternative allele at site
![k](https://latex.codecogs.com/png.latex?k "k"). The data were
previously filtered to include only bi-allelic loci.

We compute ![{\\bf
S}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20S%7D "{\\bf S}") as

  
![&#10;{\\bf S} = \\left({\\bf F} - \\frac{1}{2}{\\bf J}\_{(n\\times
\\ell)}\\right)\\left({\\bf F} - \\frac{1}{2}{\\bf J}\_{(n\\times
\\ell)}\\right)',&#10;](https://latex.codecogs.com/png.latex?%0A%7B%5Cbf%20S%7D%20%3D%20%5Cleft%28%7B%5Cbf%20F%7D%20-%20%5Cfrac%7B1%7D%7B2%7D%7B%5Cbf%20J%7D_%7B%28n%5Ctimes%20%5Cell%29%7D%5Cright%29%5Cleft%28%7B%5Cbf%20F%7D%20-%20%5Cfrac%7B1%7D%7B2%7D%7B%5Cbf%20J%7D_%7B%28n%5Ctimes%20%5Cell%29%7D%5Cright%29%27%2C%0A
"
{\\bf S} = \\left({\\bf F} - \\frac{1}{2}{\\bf J}_{(n\\times \\ell)}\\right)\\left({\\bf F} - \\frac{1}{2}{\\bf J}_{(n\\times \\ell)}\\right)',
")  

where ![{\\bf J}\_{(n\\times
\\ell)}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20J%7D_%7B%28n%5Ctimes%20%5Cell%29%7D
"{\\bf J}_{(n\\times \\ell)}") is an all-ones matrix and
!['](https://latex.codecogs.com/png.latex?%27 "'") denotes the matrix
transpose. This definition is the matrix representation of the allelic
covariance matrix defined in (Bradburd, Coop, and Ralph 2018).

### Model

We assume that ![{\\bf S} \\sim \\text{Wishart}(\\ell,\\ {\\boldsymbol
\\Sigma})](https://latex.codecogs.com/png.latex?%7B%5Cbf%20S%7D%20%5Csim%20%5Ctext%7BWishart%7D%28%5Cell%2C%5C%20%7B%5Cboldsymbol%20%5CSigma%7D%29
"{\\bf S} \\sim \\text{Wishart}(\\ell,\\ {\\boldsymbol \\Sigma})"),
where the number of SNP loci,
![\\ell](https://latex.codecogs.com/png.latex?%5Cell "\\ell"), is the
degrees of freedom parameter and ![\\boldsymbol
\\Sigma](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5CSigma
"\\boldsymbol \\Sigma") is the scale matrix. To model spatial dependence
among individual genotypes (i.e., isolation by distance (Wright 1943) or
isolation by resistance (McRae 2006)), we let

  
![&#10;\\boldsymbol{\\Sigma} = ({\\bf M} - \\rho{\\bf
W})^{-1}.&#10;](https://latex.codecogs.com/png.latex?%0A%5Cboldsymbol%7B%5CSigma%7D%20%3D%20%28%7B%5Cbf%20M%7D%20-%20%5Crho%7B%5Cbf%20W%7D%29%5E%7B-1%7D.%0A
"
\\boldsymbol{\\Sigma} = ({\\bf M} - \\rho{\\bf W})^{-1}.
")  

Above, ![{\\bf W}\_{(n\\times
n)}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20W%7D_%7B%28n%5Ctimes%20n%29%7D
"{\\bf W}_{(n\\times n)}") is the weights matrix. It determines the
degree of connectivity among nodes (defined as plants here) in a spatial
network and is the parameter of interest. The parameter ![{\\bf
M}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20M%7D "{\\bf M}") is
a diagonal matrix with ![m\_{ii} =
\\sum\_{j=1}^nw\_{ij}](https://latex.codecogs.com/png.latex?m_%7Bii%7D%20%3D%20%5Csum_%7Bj%3D1%7D%5Enw_%7Bij%7D
"m_{ii} = \\sum_{j=1}^nw_{ij}") and all off-diagonal elements equal to
zero, and ![\\rho](https://latex.codecogs.com/png.latex?%5Crho "\\rho")
controls the amount of spatial autocorrelation among the genotypes of
nearby individuals. We model the weight between subpopulations
![i](https://latex.codecogs.com/png.latex?i "i") and
![j](https://latex.codecogs.com/png.latex?j "j") as a log-linear
combination of covariates and regression parameters such that

  
![&#10;w\_{ij} = \\exp\\{{\\bf x}\_{ij}'\\boldsymbol
\\beta\\},&#10;](https://latex.codecogs.com/png.latex?%0Aw_%7Bij%7D%20%3D%20%5Cexp%5C%7B%7B%5Cbf%20x%7D_%7Bij%7D%27%5Cboldsymbol%20%5Cbeta%5C%7D%2C%0A
"
w_{ij} = \\exp\\{{\\bf x}_{ij}'\\boldsymbol \\beta\\},
")  
where ![{\\bf
x}\_{ij}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20x%7D_%7Bij%7D
"{\\bf x}_{ij}") is a vector of covariates for individuals
![i](https://latex.codecogs.com/png.latex?i "i") and
![j](https://latex.codecogs.com/png.latex?j "j") and ![\\boldsymbol
\\beta](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Cbeta
"\\boldsymbol \\beta") is a vector of regression parameters. Our
covariates of interest are

  - geographic distance (km)
  - proportion of forested cells that intersect the Euclidean line
    between plants ![i](https://latex.codecogs.com/png.latex?i "i") and
    ![j](https://latex.codecogs.com/png.latex?j "j")
  - interaction between distance and intervening forest
  - percent canopy cover for plant
    ![i](https://latex.codecogs.com/png.latex?i "i")
  - percent canopy cover for plant
    ![j](https://latex.codecogs.com/png.latex?j "j")
  - flowering plant density around plant
    ![i](https://latex.codecogs.com/png.latex?i "i")
    (individuals/25![\\pi
    \\text{m}^2](https://latex.codecogs.com/png.latex?%5Cpi%20%5Ctext%7Bm%7D%5E2
    "\\pi \\text{m}^2"))
  - flowering plant density around plant
    ![j](https://latex.codecogs.com/png.latex?j "j")
    (individuals/25![\\pi
    \\text{m}^2](https://latex.codecogs.com/png.latex?%5Cpi%20%5Ctext%7Bm%7D%5E2
    "\\pi \\text{m}^2"))
  - proportion of landscape that is forested within 100m of plant
    ![i](https://latex.codecogs.com/png.latex?i "i")
  - proportion of landscape that is forested within 100m of plant
    ![j](https://latex.codecogs.com/png.latex?j "j")

### Priors

We assume the following weakly informative priors for the regression
parameters and ![\\rho](https://latex.codecogs.com/png.latex?%5Crho
"\\rho"):

  
![&#10;{\\boldsymbol \\beta} \\sim \\mathcal{N}({\\bf 0}\_P,\\ {\\bf
I}\_{(P\\times
P)}),&#10;](https://latex.codecogs.com/png.latex?%0A%7B%5Cboldsymbol%20%5Cbeta%7D%20%5Csim%20%5Cmathcal%7BN%7D%28%7B%5Cbf%200%7D_P%2C%5C%20%7B%5Cbf%20I%7D_%7B%28P%5Ctimes%20P%29%7D%29%2C%0A
"
{\\boldsymbol \\beta} \\sim \\mathcal{N}({\\bf 0}_P,\\ {\\bf I}_{(P\\times P)}),
")  

and

  
![&#10;\\rho \\sim
\\mathcal{B}(2,2).&#10;](https://latex.codecogs.com/png.latex?%0A%5Crho%20%5Csim%20%5Cmathcal%7BB%7D%282%2C2%29.%0A
"
\\rho \\sim \\mathcal{B}(2,2).
")  

#### Stan Model Code

``` stan

functions{
// function to make weights matrix from regression array X and
// paramater vector v
  matrix make_W(int N, vector v, real[,,] X){
    matrix[N,N] W;
    for(i in 1:N){
      if(i == 1){W[i,i] = 0;}
      else{
        for(j in 1:(i-1)){
          W[i,i] = 0;
          W[i,j] = exp(to_row_vector(X[i,j,])*v);
          W[j,i] = F[i,j];
        }
      }
    }
    return(W);
  }

}

data{

  int<lower=1> N;          //number of individuals sampled
  int<lower=1> L;          //number of SNP loci
  int<lower=1> P;          //number of landscape and node variables
  
  real X[N,N,P];           //design array
  matrix[N,N] S;           //scatter matrix 

}

parameters{

  vector[P] beta;              //regression parameters
  real<lower=0,upper=1> rho;   //spatial dependence
  //real<lower=0> tau;           //

}

transformed parameters{

  matrix[N,N] W;
  vector[N] W_sum;
  matrix[N,N] M;
  cov_matrix[N] Sigma;
  //matrix[N,N] Psi;
  
  W = make_W(N, beta, X);
  
  for(i in 1:N){
    W_sum[i] = sum(W[i,]);
  }
  
  M = diag_matrix(W_sum);
  
  Sigma = inverse((M - (rho*W)));

}

model{
  
//priors
  beta ~ normal(0,1);
  rho ~ beta(2,2);

  
//likelihood
  S ~ wishart(L, Sigma);

}
```

### Loading data and fitting the model

``` r
  load(here("Data", "Indiv_lndscp_gen.RData"))
```

Because we assume the weights matrix ![{\\bf
W}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20W%7D "{\\bf W}") is
symmetric, we average the node-specific covariates. For example, rather
than including an effect for the canopy cover at the location of plant
![i](https://latex.codecogs.com/png.latex?i "i") and one for the canopy
cover at plant ![j](https://latex.codecogs.com/png.latex?j "j") for the
linear predictor of ![\\log
(w\_{ij})](https://latex.codecogs.com/png.latex?%5Clog%20%28w_%7Bij%7D%29
"\\log (w_{ij})"), we include a single effect for the average canopy
cover at plant ![i](https://latex.codecogs.com/png.latex?i "i") and
![j](https://latex.codecogs.com/png.latex?j "j"). We compile the
predictors into a 3-dimensional array, ![{\\bf X}\_{(n\\times n \\times
P)}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20X%7D_%7B%28n%5Ctimes%20n%20%5Ctimes%20P%29%7D
"{\\bf X}_{(n\\times n \\times P)}"), where
![P](https://latex.codecogs.com/png.latex?P "P") is the number of
covariates plus 1 (for the intercept).

``` r
  N <- dim(col_alscatter)[1]
  num_covs <- dim(X_inter)[3]
  num_snps <- dim(col_popgen_data)[2] - 3

  mod_data <- list(N=N,
                   P=num_covs,
                   L=num_snps,
                   X=X_inter,
                   S=col_alscatter)
  
  mfit <- sampling(wishart_edgeweight_model,
                   data=mod_data,
                   iter=1000,
                   chains=2)
```

<div id="refs" class="references">

<div id="ref-bradburd2018">

Bradburd, Gideon S., Graham M. Coop, and Peter L. Ralph. 2018.
“Inferring Continuous and Discrete Population Genetic Structure Across
Space.” *Genetics* 210 (1): 33–52.
<https://doi.org/10.1534/genetics.118.301333>.

</div>

<div id="ref-mcrae2006">

McRae, Brad H. 2006. “Isolation by Resistance.” *Evolution* 60 (8):
1551–61. <https://doi.org/10.1111/j.0014-3820.2006.tb00500.x>.

</div>

<div id="ref-wright1943">

Wright, Sewall. 1943. “Isolation by Distance.” *Genetics* 28 (2):
114–38.

</div>

</div>
