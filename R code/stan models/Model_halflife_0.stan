// remember to open text file and save as .stan file first click 'check on save'

// functions{
//   
//   real P_halflife(real a_theta, real a, real b){
//     return 1-(exp(b*(1-exp(a*a_theta))/a));
//   }
// }

data{
  
  int S; // number of data points int because its a count
  // int X_surv[S]; // Number of mosquitos that killed and there are S of them. Int because count again
  // int N_caught[S]; // number of mosquitos that are caught
  vector[S] prop_dead; // This is the predictor lp
  
  // data on # successfully fed
  int X_halflife[S]; // mortality for twenty washes
  int N_caught_halflife[S]; // mortality for twenty washes
  
  // geographic
  // int nsite;
  // int site[S]; // number between 1 and n sites - 1,2,3,...nsite (numeric integer index for each site)
  
  // variables for simulations (credible intervals)
  // int S_sim;
  // vector[S_sim] theta_sim;
  
  
}

parameters{
  
  // real<lower=0, upper=1> theta[S]; // continuous paramter between 0 and 1 called theta, we have S of these things
  real a;
  real b;
  // real phi[nsite];
  // real<lower=0> sigma_phi;
}

model{
   real sp[S];
  // mortality (to estimate theta)

    // priors by looking at the data - a (8), b(0.006)
    a ~ normal(0,100);
    b ~ normal(0,100); // allow a range thats slightly wider than you think it could be - weakly informative prior
  
  //likelihood
    for (i in 1:S){
      sp[i] = a  + b * prop_dead[i];
    }
     
  // half-life
    X_halflife ~ binomial_logit(N_caught_halflife, sp); // prob of feeding depends on theta, RE site
    
    
    
    // priors on random effects - makes it a bayesian hierarchical model
    // phi ~ normal(0,sigma_phi);
    // sigma_phi ~ normal(0,1);
    // 
    
}

// for loglikelihood and credible intervals

