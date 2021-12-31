// This Stan file defines a similar model to
// that defined in the Edgeweight_estimation_gendist.Rmd
// document, but does not include the individual meadow 
// effects

// User-defined functions block
functions{

// Function to create contrasts matrix L
  matrix contrasts(int N){
    matrix[N-1,N] L;
    for(i in 1:(N-1)){
      for(j in 1:N){
        if(i==j){L[i,j] = 1;}
        else if(j == (i+1)){L[i,j] = -1;}
        else{L[i,j] = 0;}
      }
    }
    return(L);
  }
  
// function to make weights matrix from regression array X and
// paramater vector v
  matrix make_W(int N, vector v, real[,,] X){
    matrix[N,N] Wmat;
    for(i in 1:N){
      if(i == 1){Wmat[i,i] = 0;}
      else{
        for(j in 1:(i-1)){
          Wmat[i,i] = 0;
          Wmat[i,j] = exp(to_row_vector(X[i,j,])*v);
          Wmat[j,i] = Wmat[i,j];
        }
      }
    }
    return(Wmat);
  }

}


data{

  int<lower=1> N;          //number of individuals sampled
  int<lower=1> K;          //number of SNP loci
  int<lower=1> P;          //number of landscape variables
  
  real X[N,N,P];           //design array with explanatory variables
  matrix[N,K] al_LD;       //allelic load matrix 

}

// Transform allelic load matrix into contrasts on distances
transformed data{

  matrix[N-1,N] L = contrasts(N);
  
  matrix[N-1,N-1] S = 2*(L*al_LD)*((L*al_LD)'); // transformation defined in Petkova et al. (2016)

}


parameters{

  vector[P] beta;                   //regression parameters
  real<lower=0,upper=1> rho;        //spatial dependence
  real<lower=0,upper=1> kappa_std;  // unit-scale degrees of freedom parameter

}


transformed parameters{

  matrix[N,N] W;                   // spatial weights or conductance matrix
  vector[N] W_sum;                 // vector of row sums of W
  matrix[N,N] M;                   // diagonal matrix with W_sum along diagonal
  cov_matrix[N] Sigma;             // inverse of (M-rho*W)
  real<lower=N,upper=K> kappa;     // scaled and shifted degrees of freedom parameter
  
  
  W = make_W(N, beta, X);         // define W based on loglinear model
  
  for(i in 1:N){
    W_sum[i] = sum(W[i,]);
  }
  
  M = diag_matrix(W_sum);
  
  Sigma = inverse((M - (rho*W)));
  
  kappa = (K-N)*kappa_std + N;   // shift and scale degrees of freedom parameter

}


model{
  
//priors
  beta ~ normal(0,5);
  rho ~ beta(20,1);
  kappa_std ~ beta(10, 1);
  

//likelihood
  S ~ wishart(kappa, L*(2*Sigma)*(L'));

}


generated quantities{

  real loglik;         // log-likelihood of each posterior parameter draw given data
  
  loglik = wishart_lpdf(S | kappa, L*(2*Sigma)*(L'));

}

