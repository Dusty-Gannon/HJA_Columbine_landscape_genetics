Estimating edge weights connecting HJA columbine genetic networks
================
D. G. Gannon
October 2020

### Data

We have SNP data on
![n=192](https://latex.codecogs.com/png.latex?n%3D192 "n=192")
*Aquilegia formosa* individuals from XX meadows, which we consider to be
“sub-populations” of the H.J. Andrews *A. formosa* population. We work
with a centered, euclidean distance matrix, ![{\\bf D}\_{(n\\times
n)}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20D%7D_%7B%28n%5Ctimes%20n%29%7D
"{\\bf D}_{(n\\times n)}"), computed from a genetic loading matrix
![l\_{(n\\times
\\ell)}](https://latex.codecogs.com/png.latex?l_%7B%28n%5Ctimes%20%5Cell%29%7D
"l_{(n\\times \\ell)}"), where
![\\ell](https://latex.codecogs.com/png.latex?%5Cell "\\ell") is the
number of SNP loci sequenced. Genetic loadings are coded as ![l\_{ik}
= 0](https://latex.codecogs.com/png.latex?l_%7Bik%7D%20%3D%200
"l_{ik} = 0") if individual ![i](https://latex.codecogs.com/png.latex?i
"i") is homozygous for the reference allele at locus
![k](https://latex.codecogs.com/png.latex?k "k"),
![l\_{ik}=1](https://latex.codecogs.com/png.latex?l_%7Bik%7D%3D1
"l_{ik}=1") if individual ![i](https://latex.codecogs.com/png.latex?i
"i") has a heterozygous genotype at site
![k](https://latex.codecogs.com/png.latex?k "k"), and
![l\_{ik}=2](https://latex.codecogs.com/png.latex?l_%7Bik%7D%3D2
"l_{ik}=2") if individual ![i](https://latex.codecogs.com/png.latex?i
"i") is homozygous with two copies of the alternative allele at site
![k](https://latex.codecogs.com/png.latex?k "k"). The data were
previously filtered to include only bi-allelic loci.

Letting ![\\tilde{\\bf
L}](https://latex.codecogs.com/png.latex?%5Ctilde%7B%5Cbf%20L%7D
"\\tilde{\\bf L}") be the column-centered loadings matrix
(i.e. ![\\tilde l\_{ik} = l\_{ik} - \\bar{l}\_{\\cdot
k}](https://latex.codecogs.com/png.latex?%5Ctilde%20l_%7Bik%7D%20%3D%20l_%7Bik%7D%20-%20%5Cbar%7Bl%7D_%7B%5Ccdot%20k%7D
"\\tilde l_{ik} = l_{ik} - \\bar{l}_{\\cdot k}"), where
![\\bar{l}\_{\\cdot
k}](https://latex.codecogs.com/png.latex?%5Cbar%7Bl%7D_%7B%5Ccdot%20k%7D
"\\bar{l}_{\\cdot k}") is the average loading at position
![k](https://latex.codecogs.com/png.latex?k "k")), we compute ![{\\bf
D}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20D%7D "{\\bf D}") as

  
![&#10;{\\bf D} = \\text{diag}({\\bf G}){\\bf 1}' - 2{\\bf G} +
\\text{diag}({\\bf
G}){\\bf 1}',&#10;](https://latex.codecogs.com/png.latex?%0A%7B%5Cbf%20D%7D%20%3D%20%5Ctext%7Bdiag%7D%28%7B%5Cbf%20G%7D%29%7B%5Cbf%201%7D%27%20-%202%7B%5Cbf%20G%7D%20%2B%20%5Ctext%7Bdiag%7D%28%7B%5Cbf%20G%7D%29%7B%5Cbf%201%7D%27%2C%0A
"
{\\bf D} = \\text{diag}({\\bf G}){\\bf 1}' - 2{\\bf G} + \\text{diag}({\\bf G}){\\bf 1}',
")  

where ![{\\bf G} = \\tilde{\\bf L}^\\text{T} \\tilde{\\bf
L}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20G%7D%20%3D%20%5Ctilde%7B%5Cbf%20L%7D%5E%5Ctext%7BT%7D%20%5Ctilde%7B%5Cbf%20L%7D
"{\\bf G} = \\tilde{\\bf L}^\\text{T} \\tilde{\\bf L}") and
![{\\bf 1}](https://latex.codecogs.com/png.latex?%7B%5Cbf%201%7D
"{\\bf 1}") is an ![n](https://latex.codecogs.com/png.latex?n
"n")-dimensional vector of ones.

### Model

We assume that ![-{\\bf D} \\sim \\text{Wishart}(1,\\ 2{\\boldsymbol
\\Sigma})](https://latex.codecogs.com/png.latex?-%7B%5Cbf%20D%7D%20%5Csim%20%5Ctext%7BWishart%7D%281%2C%5C%202%7B%5Cboldsymbol%20%5CSigma%7D%29
"-{\\bf D} \\sim \\text{Wishart}(1,\\ 2{\\boldsymbol \\Sigma})"), where

  
![&#10;\\boldsymbol{\\Sigma} = ({\\bf M} - \\rho{\\bf
W})^{-1}.&#10;](https://latex.codecogs.com/png.latex?%0A%5Cboldsymbol%7B%5CSigma%7D%20%3D%20%28%7B%5Cbf%20M%7D%20-%20%5Crho%7B%5Cbf%20W%7D%29%5E%7B-1%7D.%0A
"
\\boldsymbol{\\Sigma} = ({\\bf M} - \\rho{\\bf W})^{-1}.
")  

Above, ![{\\bf W}\_{(n\\times
n)}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20W%7D_%7B%28n%5Ctimes%20n%29%7D
"{\\bf W}_{(n\\times n)}") is the weights matrix. It determines the
degree of connectivity among nodes (individuals) in a spatial network
and is the parameter of interest. The parameter ![{\\bf
M}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20M%7D "{\\bf M}") is
a diagonal matrix with ![m\_{ii} =
\\sum\_{j=1}^nw\_{ij}](https://latex.codecogs.com/png.latex?m_%7Bii%7D%20%3D%20%5Csum_%7Bj%3D1%7D%5Enw_%7Bij%7D
"m_{ii} = \\sum_{j=1}^nw_{ij}") and all off-diagonal elements equal to
zero, and ![\\rho](https://latex.codecogs.com/png.latex?%5Crho "\\rho")
controls the amount of spatial autocorrelation among the genotypes of
nearby individuals. We model the weight between individuals
![i](https://latex.codecogs.com/png.latex?i "i") and
![j](https://latex.codecogs.com/png.latex?j "j") as a log-linear
combination of covariates and regression parameters such that

  
![&#10;w\_{ij} = \\exp\\{{\\bf x}\_{ij}'\\boldsymbol \\beta +
a\_i\\},&#10;](https://latex.codecogs.com/png.latex?%0Aw_%7Bij%7D%20%3D%20%5Cexp%5C%7B%7B%5Cbf%20x%7D_%7Bij%7D%27%5Cboldsymbol%20%5Cbeta%20%2B%20a_i%5C%7D%2C%0A
"
w_{ij} = \\exp\\{{\\bf x}_{ij}'\\boldsymbol \\beta + a_i\\},
")  
where ![{\\bf
x}\_{ij}](https://latex.codecogs.com/png.latex?%7B%5Cbf%20x%7D_%7Bij%7D
"{\\bf x}_{ij}") is a vector of covariates for individuals
![i](https://latex.codecogs.com/png.latex?i "i") and
![j](https://latex.codecogs.com/png.latex?j "j"), ![\\boldsymbol
\\beta](https://latex.codecogs.com/png.latex?%5Cboldsymbol%20%5Cbeta
"\\boldsymbol \\beta") is a vector of regression parameters, and
![a\_i](https://latex.codecogs.com/png.latex?a_i "a_i") is an effect for
the individual plant. We assume ![a\_1,a\_2,...,a\_n \\sim
\\mathcal{N}(0,\\sigma^2)](https://latex.codecogs.com/png.latex?a_1%2Ca_2%2C...%2Ca_n%20%5Csim%20%5Cmathcal%7BN%7D%280%2C%5Csigma%5E2%29
"a_1,a_2,...,a_n \\sim \\mathcal{N}(0,\\sigma^2)"). Our covariates of
interest are

  - geographic distance (km)
  - proportion of forested cells that intersect the Euclidean line
    between plants ![i](https://latex.codecogs.com/png.latex?i "i") and
    ![j](https://latex.codecogs.com/png.latex?j "j")
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

### Priors

We assume the following weakly informative priors for the regression
parameters, ![\\rho](https://latex.codecogs.com/png.latex?%5Crho
"\\rho"), and ![\\sigma](https://latex.codecogs.com/png.latex?%5Csigma
"\\sigma"):

  
![&#10;{\\boldsymbol \\beta} \\sim \\mathcal{N}({\\bf 0}\_P,\\ {\\bf
I}\_{(P\\times
P)}),&#10;](https://latex.codecogs.com/png.latex?%0A%7B%5Cboldsymbol%20%5Cbeta%7D%20%5Csim%20%5Cmathcal%7BN%7D%28%7B%5Cbf%200%7D_P%2C%5C%20%7B%5Cbf%20I%7D_%7B%28P%5Ctimes%20P%29%7D%29%2C%0A
"
{\\boldsymbol \\beta} \\sim \\mathcal{N}({\\bf 0}_P,\\ {\\bf I}_{(P\\times P)}),
")  

  
![&#10;\\rho \\sim
\\mathcal{U}(0,1),&#10;](https://latex.codecogs.com/png.latex?%0A%5Crho%20%5Csim%20%5Cmathcal%7BU%7D%280%2C1%29%2C%0A
"
\\rho \\sim \\mathcal{U}(0,1),
")  

and

  
![&#10;\\sigma \\sim
\\text{half-Normal}(0,2).&#10;](https://latex.codecogs.com/png.latex?%0A%5Csigma%20%5Csim%20%5Ctext%7Bhalf-Normal%7D%280%2C2%29.%0A
"
\\sigma \\sim \\text{half-Normal}(0,2).
")  

#### Stan Model Code

``` stan

functions{

  matrix make_W(int N, vector v, real[,,] X){
    matrix[N,N] W;
    for(i in 1:N){
      if(i == 1){W[i,i] = 0;}
      else{
        for(j in 1:(i-1)){
          W[i,i] = 0;
          W[i,j] = exp(to_row_vector(X[i,j,])*v);
          W[j,i] = W[i,j];
        }
      }
    }
    return(W);
  }

}

data{

  int<lower=1> N;          //number of individuals sampled
  int<lower=1> L;          //number of SNP loci
  int<lower=1> P;          //number of landscape variables
  
  real X[N,N,P];        //design array
  matrix[N,N] D;           //distance matrix 

}

parameters{

  vector[P] beta;              //regression parameters
  real<lower=0,upper=1> rho;   //spatial dependence
  //real a[N];                   //plant effects
  //real<lower=0> sigma;         //variation of plant effects

}

transformed parameters{

  matrix[N,N] W;
  vector[N] W_sum;
  matrix[N,N] M;
  cov_matrix[N] Sigma;
  
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
  rho ~ uniform(0,1);
  
//likelihood
  D ~ wishart(L, Sigma);


}
```
