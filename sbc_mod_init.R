
# do simulation-based calibration (from Carpenter paper) for marginal and joint model
rm(list=ls())
 
library(rstan)

# set working directory
wdir <- file.path("~/bayes_cop_calib")

# marginal efficacy model
mod_bb_e_code <- "
data {
int<lower=0> N;
matrix[N, 2] x;
}
transformed data {
vector[2] beta_e_sim; 
//vector[N] y_sim;
int<lower=0, upper=1> y_sim[N];

beta_e_sim[1] = normal_rng(0,10);
beta_e_sim[2] = normal_rng(0,10);

for (n in 1:N)
  y_sim[n] = bernoulli_rng(Phi(x[n,]*beta_e_sim));
}
parameters {
//params for binary (efficacy) outcome
vector[2] beta_e;
}
model {
// priors
beta_e ~ normal(0,10); 

// marginal for binary (efficacy) outcome
y_sim ~ bernoulli(Phi(x*beta_e)); 
}
generated quantities {
  int<lower=0, upper=1> I_lt_sim[2] 
    = { beta_e_sim[1] < beta_e[1], beta_e_sim[2] < beta_e[2] };
}
"

## Compile and save models

mod_bb_e <- stan_model(model_code = mod_bb_e_code)
saveRDS(mod_bb_e, file = file.path(wdir,"mod_bb_e.rds"))