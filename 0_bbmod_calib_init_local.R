## binary-binary model for copula calibration - sim array and compile models
## LOCAL version

rm(list=ls())
 
library(rstan)

# set working directory
#wdir <- file.path("~/bayes_cop_pow")

## Define models
#! may not need individual marginal models

# marginal efficacy model
mod_bb_e_code <- "
data {
int<lower=0> N;
matrix[N, 2] x;
int<lower=0, upper=1> y_e[N];
}
parameters {
//params for binary (efficacy) outcome
vector[2] beta_e;
}
model {
// priors
beta_e ~ normal(0,10); 

// marginal for binary (efficacy) outcome
y_e ~ bernoulli(Phi(x*beta_e)); 
}
generated quantities {
// transform to get p_e for placebo and treatment
vector[2] p_e;
p_e = Phi(beta_e);
}
"

# marginal safety model
mod_bb_s_code <- "
data {
int<lower=0> N;
matrix[N, 2] x;
int<lower=0, upper=1> y_s[N];
}
parameters {
//params for binary (safety) outcome
vector[2] beta_s;
}
model {
// priors
beta_s ~ normal(0,10); 

// marginal for binary (safety) outcome
y_s ~ bernoulli(Phi(x*beta_s)); 
}
generated quantities {
// transform to get p_s for placebo and treatment
vector[2] p_s;
p_s = Phi(beta_s);
}
"

# joint models
mod_bb_jnt_code <- "
functions {
#include /bb_normCop_funs.stan
}
data {
int<lower=1> N;
matrix[N, 2] x;
int<lower=0, upper=1> y_e[N];
int<lower=0, upper=1> y_s[N];
int<lower=0> freq[N];
}
parameters {
// params for binary (efficacy) outcome
vector[2] beta_e;

//params for binary (safety) outcome
vector[2] beta_s;

// copula dependence param
vector<lower=-1, upper=1>[2] rho_;  
}
model {
vector[N] pr_e;
vector[N] pr_s;
vector[N] theta;

// priors
// uniform priors (leave beta_e and beta_s unassigned)
//rho ~ uniform(-1,1);

// regularizing priors
// add hyperprior for variance of rho prior??
beta_e ~ normal(0,10);
beta_s ~ normal(0,10); 
rho_ ~ normal(0,0.4); 

// marginal probit for binary (efficacy) outcome
pr_e = Phi(x*beta_e);

// marginal probit for binary (safety) outcome
pr_s = Phi(x*beta_s);

// copula dependence model
theta = x*rho_;

// build log-likelihood
{
  vector[N] loglik;  // vectorize summation
  for (i in 1:N)
  loglik[i] = freq[i]*bi_cop_lp(y_e[i], pr_e[i], y_s[i], pr_s[i], theta[i]);
  target += sum(loglik);
}

}
generated quantities {
vector[2] p_e;
vector[2] p_s;
real p_e_diff;
real p_s_diff;

p_e = Phi(beta_e);
p_s = Phi(beta_s);
p_e_diff = p_e[2]-p_e[1];
p_s_diff = p_s[2]-p_s[1];
}
"


## Compile and save models

# efficacy marginal model
mod_bb_e <- stan_model(model_code = mod_bb_e_code)
saveRDS(mod_bb_e, file = file.path(getwd(),"mod_bb_e_local.rds"))

# safety marginal model
mod_bb_s <- stan_model(model_code = mod_bb_s_code)
saveRDS(mod_bb_s, file = file.path(getwd(),"mod_bb_s_local.rds"))

# joint models
mod_bb_jnt <- stan_model(model_code = mod_bb_jnt_code)
saveRDS(mod_bb_jnt, file = file.path(getwd(),"mod_bb_jnt_local.rds"))
