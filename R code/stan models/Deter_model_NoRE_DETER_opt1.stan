// remember to open text file and save as .stan file first click 'check on save'

functions{

  real P_inside(real a_theta, real c, real d, real e){ // probability of being caught inside intervention hut
    real P_deter;
    P_deter = e * exp(d * (1 - exp(c * a_theta)) / c);
    return (1-P_deter);
  }
}

data{
  
  int S; // number of data points int because its a count
  int X_survive[S]; // Number of mosquitos that survive and there are S of them. Int because count again
  int N_caught[S]; // number of mosquitos that are caught
  
  // diff in the number caught in control huts vs treated huts - paired up need to make sure they correspond to the same trial
  int X_control[S];
  
  // variables for simulations
  int S_sim;
  vector[S_sim] theta_sim;
  
}

parameters{
  
  real<lower=0, upper=1> theta[S]; // continuous paramter between 0 and 1 called theta, you have S of these things
  
  real<lower=0> N[S]; // vector of non-negative values - continuous. Represents ambient mosquitos.
  real<lower=0.50, upper=1> P_control; // probability of being caught in control hut
  real<lower=0> kappa; // over-dispersion parameter of negative binomial
  
  real<lower=0> c;
  real<lower=0> d;
  real<lower=0, upper=1> e; // to ensure the deterrence probability always has an upper bound of 1

  
}

model{
  
  // survival (to estimate theta)
  
  //likelihood
  for(i in 1:S){
    X_survive[i] ~ binomial(N_caught[i],theta[i]); //theta = probability of surviving, different for each trial
  }
  
  //prior for theta
    theta ~ uniform(0,1); // could be anywhere from 0 to 1
    
 //deterrence
 for (i in 1:S){
 X_control[i] ~ neg_binomial_2(N[i] * P_control, kappa);
 N_caught[i] ~ neg_binomial_2(N[i] * P_inside(theta[i],c,d,e), kappa);
 }
 
  P_control ~ beta(1, 10);
  N ~ lognormal(5, 1);
  kappa ~ cauchy(0, 1);
  
  //priors
  c ~ normal(3, 2);     
  d ~ normal(0.8, 1);
  e ~ normal(0.8, 1);  
} 


generated quantities{
  vector[S_sim] P_deter_sim;
  vector[S_sim] P_inside_sim;
  vector[S_sim] P_ratio;
  vector[S] loglikelihood;
  for(i in 1:S_sim){
    P_deter_sim[i]= e * exp(d * (1 - exp(c * theta_sim[i])) / c);
    P_inside_sim[i] = 1 - P_deter_sim[i];
    P_ratio[i] = (P_control - P_inside_sim[i]) / P_control;
  }
  for(i in 1:S){
    loglikelihood[i]= (binomial_lpmf(X_survive[i]|N_caught[i],theta[i]) + 
    neg_binomial_2_lpmf(X_control[i]|N[i]* P_control,kappa) + 
    neg_binomial_2_lpmf(N_caught[i]|N[i]* P_inside(theta[i],c,d,e),kappa));
      }
}
