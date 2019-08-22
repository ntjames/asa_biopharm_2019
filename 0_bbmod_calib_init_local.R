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
mod_bb_jnt_code0 <- "
functions {

// bivariate normal cdf (from Stan manual)
real binormal_cdf(real z1, real z2, real rho){
if (z1!=0 || z2 !=0){
real denom = fabs(rho) < 1.0 ? sqrt((1+rho)*(1+rho)): not_a_number();
real a1 = (z2/z1 - rho)/denom;
real a2 = (z1/z2 - rho)/denom;
real product = z1*z2;
real delt = product < 0 || (product==0 && (z1+z2)<0);
return 0.5*(Phi(z1)+Phi(z2)-delt)-owens_t(z1,a1)-owens_t(z2,a2);
}
return 0.25 + asin(rho)/(2*pi());
}

// bivariate normal copula distribution function 
real pCop_norm(real u1, real u2, real theta){
if (u1==0 || u2==0){ // grounded
return 0;
}
else if (u1==1){ // uniform margin
return u2;
}
else if (u2==1){ // uniform margin
return u1;
} else {
return binormal_cdf(inv_Phi(u1), inv_Phi(u2), theta);
}
}

// copula likelihood for 2 binary margins
real bi_cop_lp(int y1, real p1, int y2, real p2, real theta) {

// pseudo-obs from bernoulli dist
real U1 = bernoulli_cdf(y1, p1); 
real U1_prime = bernoulli_cdf(y1-1, p1);
real U2 = bernoulli_cdf(y2, p2); 
real U2_prime = bernoulli_cdf(y2-1, p2);

// likelihood for joint dist. using differences
// replace with sum() function?
real cP0 = pCop_norm(U1, U2, theta) -
pCop_norm(U1_prime, U2, theta) -
pCop_norm(U1, U2_prime, theta) +
pCop_norm(U1_prime, U2_prime, theta);

// put lower bound on cP to avoid log(0)
// is there a better way to do this??
real cP = cP0 < 1e-200 ? 1e-200: cP0;

// log-likelihood
return log(cP);
}

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
