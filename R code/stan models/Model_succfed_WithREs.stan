// remember to open text file and save as .stan file first click 'check on save'

functions{
  
  real P_fed(real a_theta, real a, real b){
    return 1-(exp(b*(1-exp(a*a_theta))/a));
  }
}

data{
  
  int S; // number of data points int because its a count
  int X_survive[S]; // Number of mosquitos that survive and there are S of them. Int because count again
  int N_caught[S]; // number of mosquitos that are caught
  
  
  // data on # successfully fed
  int X_succfed[S]; // need to be ordered so they each correspond to the same thing e.g. X_fed and X_Survive
  
  // geographic
  int nsite;
  int site[S]; // number between 1 and n sites - 1,2,3,...nsite (numeric integer index for each site)
  
  // variables for simulations (credible intervals)
  int S_sim;
  vector[S_sim] theta_sim;
  
}

parameters{
  
  real<lower=0, upper=1> theta[S]; // continuous paramter between 0 and 1 called theta, you have S of these things
  real<lower=0> a;
  real<lower=0> b;
  real phi[nsite];
  real<lower=0> sigma_phi;
}

model{
  
  // survival (to estimate theta)
  
  //likelihood
  for(i in 1:S)
    X_survive[i] ~ binomial(N_caught[i],theta[i]); //theta = probability of surviving, different for each trial
    
  //prior for theta
    theta ~ uniform(0,1); // could be anywhere from 0 to 1
    
  // feeding
  for (i in 1:S)
  X_succfed[i] ~ binomial_logit(N_caught[i], logit(P_fed(theta[i],a,b)) + phi[site[i]]); // prob of feeding depends on theta, RE site
  
  // priors by looking at the data - a (8), b(0.006)
  a ~ normal(8,0.3);
  b ~ normal(0,0.01); // allow a range thats slightly wider than you think it could be - weakly informative prior
  

 // priors on random effects - makes it a bayesian hierarchical model
 phi ~ normal(0,sigma_phi);
 sigma_phi ~ normal(0,1);
 
  
}

// for loglikelihood and credible intervals

generated quantities{
  vector[S_sim] P_fed_sim;
  vector[S] loglikelihood;
  for(i in 1:S_sim){
    P_fed_sim[i]= 1-(exp(b*(1-exp(a*theta_sim[i]))/a));
  }
  for(i in 1:S){
    loglikelihood[i]= (binomial_lpmf(X_survive[i]|N_caught[i],theta[i]) + 
    binomial_logit_lpmf(X_succfed[i] | N_caught[i], logit(P_fed(theta[i],a,b)) + phi[site[i]]));

   }
}
