data {
    int<lower=1> Ntotal; // number of observations
    
    int<lower=1> Nlevels1; // number of random intercepts for factor 1 = nlevels of factor 1 
    int<lower=1> NfactorMap1[Ntotal]; // mapping variable to map each observation to a particular level of factor 1
    
    int<lower=1> Nlevels2; // number of random intercepts for factor 2 = nlevels of factor 2 
    int<lower=1> NfactorMap2[Ntotal]; // mapping variable to map each observation to a particular level of factor 2
    
    int<lower=1> Nlevels3; // number of random intercepts for factor 3 = nlevels of factor 3 
    int<lower=1> NfactorMap3[Ntotal]; // mapping variable to map each observation to a particular level of factor 3
    
    int<lower=1> Nlevels4; // number of random intercepts for factor 3 = nlevels of factor 4 
    int<lower=1> NfactorMap4[Ntotal]; // mapping variable to map each observation to a particular level of factor 4
    
    int y[Ntotal]; // response variable binomial distributed
    int Ntrials[Ntotal]; // total trials for binomial variable
}
parameters {
  // parameters to estimate in the model
    real intercept; // constant term
    
    real<lower=0.01> sigmaFactor1; // population level sd for factor 1
    vector[Nlevels1] nCoefFactor1; // number of random jitters, for levels of factor 1
    
    real<lower=0.01> sigmaFactor2; // population level sd for factor 2
    vector[Nlevels2] nCoefFactor2; // number of random jitters, for levels of factor 2
    
    real<lower=0.01> sigmaFactor3; // population level sd for factor 3
    vector[Nlevels3] nCoefFactor3; // number of random jitters, for levels of factor 3
    
    real<lower=0.01> sigmaFactor4; // population level sd for factor 4
    vector[Nlevels4] nCoefFactor4; // number of random jitters, for levels of factor 4
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  mu = intercept + nCoefFactor1[NfactorMap1] + nCoefFactor2[NfactorMap2] + nCoefFactor3[NfactorMap3]
      + nCoefFactor4[NfactorMap4];
  mu = inv_logit(mu);
}
model {
  // using generic priors to start with
  sigmaFactor1 ~ cauchy(0, 1);
  sigmaFactor2 ~ cauchy(0, 1);
  sigmaFactor3 ~ cauchy(0, 1);
  sigmaFactor4 ~ cauchy(0, 1);
  
  intercept ~ cauchy(0, 2);
  
  nCoefFactor1 ~ normal(0, sigmaFactor1);
  nCoefFactor2 ~ normal(0, sigmaFactor2);
  nCoefFactor3 ~ normal(0, sigmaFactor3);
  nCoefFactor4 ~ normal(0, sigmaFactor4);
  
  // likelihood function
  y ~ binomial(Ntrials, mu);
}
