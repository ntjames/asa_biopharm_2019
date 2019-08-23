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

# joint models with different priors (if time, try diff copula, e.g. Frank)
mod_bb_jnt_code1 <- "
// model 1 - flat priors

functions {
#include /bb_normCop_funs.stan
}
data {
#include /bb_short_dat.stan
}
parameters {
// params for binary efficacy (beta_e) and safety (beta_s) outcomes
vector[2] beta_e; 
vector[2] beta_s;

// copula dependence param
vector<lower=-1, upper=1>[2] rho_;  
}
model {
vector[N] pr_e;
vector[N] pr_s;
vector[N] theta;

// priors (unassigned)
// improper (flat) priors for beta_e and beta_s 
// uniform(-1,1) prior for rho
// rho_ ~ uniform(-1,1);

// marginal probit models for binary efficacy and safety outcomes
pr_e = Phi(x*beta_e);
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
#include /bb_normCop_genquants.stan
}
"



mod_bb_jnt_code2 <- "
// model 2 

functions {
#include /bb_normCop_funs.stan
}
data {
#include /bb_short_dat.stan
}
parameters {
// params for binary efficacy (beta_e) and safety (beta_s) outcomes
vector[2] beta_e; 
vector[2] beta_s;

// copula dependence param
vector<lower=-1, upper=1>[2] rho_;  
}
model {
vector[N] pr_e;
vector[N] pr_s;
vector[N] theta;

// priors
// weakly informative priors for beta_e and beta_s
beta_e ~ normal(0,10);
beta_s ~ normal(0,10);

// weakly informative prior for rho
// add hyperprior for variance of rho prior??
rho_ ~ normal(0,0.4); 

// marginal probit models for binary efficacy and safety outcomes
pr_e = Phi(x*beta_e);
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
#include /bb_normCop_genquants.stan
}
"

mod_bb_jnt_code3 <- "
// model 3

functions {
#include /bb_normCop_funs.stan
}
data {
#include /bb_short_dat.stan
}
parameters {
// params for binary efficacy (beta_e) and safety (beta_s) outcomes
vector[2] beta_e; 
vector[2] beta_s;

// unscaled copula dependence param 
vector<lower=0, upper=1>[2] rho0;  
}
transformed parameters {
// scaled copula dependence param 
vector<lower=-1, upper=1>[2] rho_ = 2*(rho0 - 0.5);
}
model {
vector[N] pr_e;
vector[N] pr_s;
vector[N] theta;

// priors
// weakly informative priors for beta_e and beta_s
beta_e ~ normal(0,10);
beta_s ~ normal(0,10);

// boundary-avoiding beta() prior on unscaled rho0
// see https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
// --> Boundary-avoiding priors for modal estimation
rho0 ~ beta(2,2); 

// marginal probit models for binary efficacy and safety outcomes
pr_e = Phi(x*beta_e);
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
#include /bb_normCop_genquants.stan
}
"


## Compile and save models

if (0){
# efficacy marginal model
mod_bb_e <- stan_model(model_code = mod_bb_e_code)
saveRDS(mod_bb_e, file = file.path(getwd(),"mod_bb_e_local.rds"))

# safety marginal model
mod_bb_s <- stan_model(model_code = mod_bb_s_code)
saveRDS(mod_bb_s, file = file.path(getwd(),"mod_bb_s_local.rds"))
}

# joint models
mod_bb_jnt1 <- stan_model(model_code = mod_bb_jnt_code1)
saveRDS(mod_bb_jnt1, file = file.path(getwd(),"mod_bb_jnt1_local.rds"))

mod_bb_jnt2 <- stan_model(model_code = mod_bb_jnt_code2)
saveRDS(mod_bb_jnt2, file = file.path(getwd(),"mod_bb_jnt2_local.rds"))

mod_bb_jnt3 <- stan_model(model_code = mod_bb_jnt_code3)
saveRDS(mod_bb_jnt3, file = file.path(getwd(),"mod_bb_jnt3_local.rds"))
