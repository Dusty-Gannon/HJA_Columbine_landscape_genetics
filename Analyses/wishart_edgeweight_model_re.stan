//

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
  matrix make_W(int N, vector v, real[,,] X, vector a, matrix Z, real s){
    matrix[N,N] Wmat;
    for(i in 1:N){
      if(i == 1){Wmat[i,i] = 0;}
      else{
        for(j in 1:(i-1)){
          Wmat[i,i] = 0;
          Wmat[i,j] = exp(to_row_vector(X[i,j,])*v + 
                           (Z[i,]*a)*s + 
                           (Z[j,]*a)*s);
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
  int<lower=1> P;          //number of landscape and node variables
  int<lower=1> G;          // number of groups (meadows)
  
  real X[N,N,P];           //design array
  matrix[N,G] Z;           // random effect design matrix
  matrix[N,K] al_LD;       //allelic load matrix 

}

// Transform allelic load matrix into contrasts on distances
transformed data{

  matrix[N-1,N] L = contrasts(N);
  
  matrix[N-1,N-1] S = 2*(L*al_LD)*((L*al_LD)');

}


parameters{

  vector[P] beta;              //regression parameters
  real<lower=0,upper=1> rho;   //spatial dependence
  vector[G] alpha_raw;         // random meadow effects
  real<lower=0> tau;           // scale parameter for meadow effects
  real<lower=0,upper=1> kappa_std; // degrees of freedom parameter

}

transformed parameters{

  matrix[N,N] W;                 // spatial weights or conductance matrix
  vector[N] W_sum;               // vector of row sums of W
  matrix[N,N] M;                 // diagonal matrix with W_sum along diagonal
  cov_matrix[N] Sigma;           // inverse of (M-rho*W)
  real<lower=N,upper=K> kappa;   // scaled and shifted degrees of freedom parameter
  
  
  W = make_W(N, beta, X, alpha_raw, Z, tau);  // define W based on loglinear model
  
  for(i in 1:N){
    W_sum[i] = sum(W[i,]);                   // sum of each row
  }
  
  M = diag_matrix(W_sum);                    // diagonal matrix M
  
  Sigma = inverse((M - (rho*W)));            // definition of Sigma
  
  kappa = (K-N)*kappa_std + N;               // shift and scale degrees of freedom parameter

}

model{
  
//priors
  beta ~ normal(0,5);
  rho ~ beta(5,1);
  kappa_std ~ beta(5, 1);
  alpha_raw ~ normal(0,1);
  tau ~ normal(0,2);
  

//likelihood
  S ~ wishart(kappa, L*(2*Sigma)*(L'));

}

generated quantities{

  matrix[N-1,N-1] pred;    // Posterior predictive draws
  real loglik;             // loglikelihood of posterior draws given data
  
  loglik = wishart_lpdf(S | kappa, L*(2*Sigma)*(L'));
  pred = wishart_rng(kappa, L*(2*Sigma)*(L'));

}
