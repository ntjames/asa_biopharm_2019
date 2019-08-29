## continuous-binary model for copula calibration - sim array and compile models
rm(list=ls())
library(rstan)

# set working directory
wdir <- file.path("~/bayes_cop_calib")

## Define models
#! may not need individual marginal models

# marginal efficacy model
mod_cb_e_code <- "
data {
int N;
matrix[N, 2] x;
vector[N] y_e;
}
parameters {
// params for continuous (efficacy) outcome
vector[2] beta_e;
real<lower=0> sigma;
}
model {
vector[N] mu;

// priors
beta_e ~ normal(0,10);
sigma ~ inv_gamma(0.001,0.001);

// marginal for continuous (efficacy) outcome
mu = x*beta_e;
y_e ~ normal(mu, sigma);
}
"

# marginal safety model
mod_cb_s_code <- "
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


# joint models with different priors
mod_cb_jnt_code1 <- "
// model 1 - flat priors

functions {
#include /cb_normCop_funs.stan
}
data {
int N;
matrix[N, 2] x;
vector[N] y1;
int<lower=0, upper=1> y2[N];
}
parameters {
// params for continuous (efficacy) outcome
vector[2] beta_e;
vector<lower=0>[2] s;

//params for binary (safety) outcome
vector[2] beta_s;

// copula dependence param
vector<lower=-1, upper=1>[2] omega_;
}
model {
vector[N] mu;
vector[N] sigma;
vector[N] p;
vector[N] theta;

// priors (unassigned)
// improper (flat) priors for beta_e and beta_s
// improper prior for s (std dev)
// uniform(-1,1) prior for omega
// omega_ ~ uniform(-1,1);

// marginal for continuous (efficacy) outcome
mu = x*beta_e;
sigma = x*s;

// marginal for binary (safety) outcome
p = Phi(x*beta_s);

// copula dependence parameter
theta = x*omega_;

// build log-likelihood
{
  vector[N] loglik;  // vectorize summation
  for (i in 1:N)
    loglik[i] = binorm_cop_lp(y1[i], mu[i], sigma[i], y2[i], p[i], theta[i]);;
  target += sum(loglik);
}

}
generated quantities {
#include /cb_normCop_genquants.stan
}
"


mod_cb_jnt_code2 <- "
// model 2

functions {
#include /cb_normCop_funs.stan
}
data {
int N;
matrix[N, 2] x;
vector[N] y1;
int<lower=0, upper=1> y2[N];
}
parameters {
// params for continuous (efficacy) outcome
vector[2] beta_e;
vector<lower=0>[2] s;

//params for binary (safety) outcome
vector[2] beta_s;

// copula dependence param
vector<lower=-1, upper=1>[2] omega_;
}
model {
vector[N] mu;
vector[N] sigma;
vector[N] p;
vector[N] theta;

// weakly informative priors for beta_e and beta_s (and s?)  
beta_e ~ normal(0,10);
beta_s ~ normal(0,10);
s ~ inv_gamma(0.001,0.001);

// weakly informative prior for omega
// add hyperprior for variance of omega prior??
omega_ ~ normal(0,0.4); 

// marginal for continuous (efficacy) outcome
mu = x*beta_e;
sigma = x*s;

// marginal for binary (safety) outcome
p = Phi(x*beta_s);

// copula dependence parameter
theta = x*omega_;

// build log-likelihood
{
  vector[N] loglik;  // vectorize summation
  for (i in 1:N)
    loglik[i] = binorm_cop_lp(y1[i], mu[i], sigma[i], y2[i], p[i], theta[i]);;
  target += sum(loglik);
}

}
generated quantities {
#include /cb_normCop_genquants.stan
}
"


mod_cb_jnt_code3 <- "
// model 3

functions {
#include /cb_normCop_funs.stan
}
data {
int N;
matrix[N, 2] x;
vector[N] y1;
int<lower=0, upper=1> y2[N];
}
parameters {
// params for continuous (efficacy) outcome
vector[2] beta_e;
vector<lower=0>[2] s;

//params for binary (safety) outcome
vector[2] beta_s;

// unscaled copula dependence param 
vector<lower=0, upper=1>[2] omega0;  
}
transformed parameters {
// scaled copula dependence param 
vector<lower=-1, upper=1>[2] omega_ = 2*(omega0 - 0.5);
}
model {
vector[N] mu;
vector[N] sigma;
vector[N] p;
vector[N] theta;

// weakly informative priors for beta_e and beta_s (and s?)  
beta_e ~ normal(0,10);
beta_s ~ normal(0,10);
s ~ inv_gamma(0.001,0.001);

// boundary-avoiding beta() prior on unscaled omega0
// see https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
// --> Boundary-avoiding priors for modal estimation
omega0 ~ beta(2,2); 

// marginal for continuous (efficacy) outcome
mu = x*beta_e;
sigma = x*s;

// marginal for binary (safety) outcome
p = Phi(x*beta_s);

// copula dependence parameter
theta = x*omega_;

// build log-likelihood
{
  vector[N] loglik;  // vectorize summation
  for (i in 1:N)
    loglik[i] = binorm_cop_lp(y1[i], mu[i], sigma[i], y2[i], p[i], theta[i]);;
  target += sum(loglik);
}

}
generated quantities {
#include /cb_normCop_genquants.stan
}
"



## Compile and save models

if (0){
# efficacy marginal model
mod_cb_e <- stan_model(model_code = mod_cb_e_code)
saveRDS(mod_cb_e, file = file.path(wdir,"mod_cb_e.rds"))

# safety marginal model
mod_cb_s <- stan_model(model_code = mod_cb_s_code)
saveRDS(mod_cb_s, file = file.path(wdir,"mod_cb_s.rds"))
}

# joint models
mod_cb_jnt1 <- stan_model(model_code = mod_cb_jnt_code1)
saveRDS(mod_cb_jnt1, file = file.path(wdir,"mod_cb_jnt1.rds"))

mod_cb_jnt2 <- stan_model(model_code = mod_cb_jnt_code2)
saveRDS(mod_cb_jnt2, file = file.path(wdir,"mod_cb_jnt2.rds"))

mod_cb_jnt3 <- stan_model(model_code = mod_cb_jnt_code3)
saveRDS(mod_cb_jnt3, file = file.path(wdir,"mod_cb_jnt3.rds"))
