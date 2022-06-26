// bernoulli_logistic transformed data function
data {
  
  int<lower=1> N;                  // rows of data
  
  int<lower=0> n_t[N];             // Total number of mosquitoes entering net huts
  int<lower=0> d_t[N];             // Number mosquites dead net hut

  vector<lower=0>[N] x;       // predictor (mortality in LLIN huts)
  
}

parameters {
  //Consider death. This is the proportion of mosquitoes dying (d_t) in treated huts (n_t)
  real alpha1;
  real alpha2;
  
  //  real<lower=0,upper=10> sigma;
}

model {
  real sp[N];

  alpha1 ~ normal(0,100);
  alpha2 ~ normal(0,100);

  //  study_a ~ normal(0,sigma);
  
  for (n in 1:N) {
    sp[n] = alpha1  + alpha2 * x[n];
  }
  
  d_t ~ binomial_logit(n_t, sp);
}

generated quantities{
  real sp_ppc[100];

    for(t in 1:100){
      sp_ppc[t] = binomial_rng(100, inv_logit(alpha1 + alpha2 * t)) / 100;
  }
}
