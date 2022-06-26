// A1 Log-logistic Model

functions{
  real log_logistic(real x, real b, real a){
    return(1-(1 / (1 + ((1 - x) / b)^(-a))));
  }
}

data{
  int S;
  int N_b[S]; // Number in bioassay
  int N_h[S]; // Number in hut trial
  int X_b[S]; // Number alive in bioassay
  int X_h[S]; // Number alive in hut trial
  
  // Random effects
  int nsite; 
  int site[S]; 
  
  int S_test;
  vector[S_test] theta_b_test;
}

parameters{
  real<lower=0, upper=1> theta_b[S];
  real theta_h_logit[S];
  real<lower=0> b;
  real<lower=0> a;
  real<lower=0> sigma;
  vector[nsite] r_b;
  vector[nsite] r_h;
  real<lower=0> sigma_b;
  real<lower=0> sigma_h;
}

transformed parameters{
 real theta_h_bar_logit[S];
 for(i in 1:S)
  theta_h_bar_logit[i] = logit(log_logistic(theta_b[i], b, a));
}

model{
  // Likelihood for bioassay
  X_b ~ binomial_logit(N_b, to_vector(logit(theta_b)) + r_b[site]);
  
  // Likelihood for huts
  X_h ~ binomial_logit(N_h, to_vector(theta_h_logit) + r_h[site]);
  
  // Priors for shape of log-logistic
  b ~ cauchy(0, 1);
  a ~ normal(0, 2);
  
  // Prior for hut survival
  sigma ~ normal(0, 2);
  for(i in 1:S)
    theta_h_logit[i] ~ normal(theta_h_bar_logit[i], sigma);
    
  // Priors for random effects
  r_b ~ normal(0, sigma_b);
  r_h ~ normal(0, sigma_h);
  sigma_b ~ normal(0, 1);
  sigma_h ~ normal(0, 1);
}

generated quantities{
  vector[S_test] theta_h_test;
  vector[S] logLikelihood;
  for(i in 1:S_test)
    theta_h_test[i] = log_logistic(theta_b_test[i], b, a);
  for(i in 1:S)
    logLikelihood[i] = binomial_logit_lpmf(X_b[i]|N_b[i], logit(theta_b[i]) + r_b[site[i]]) + binomial_logit_lpmf(X_h[i]|N_h[i], theta_h_logit[i] + r_h[site[i]]);

}
