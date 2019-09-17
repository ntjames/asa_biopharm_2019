## continuous-continuous model for copula calibration - sim array and compile models

rm(list=ls())
library(rstan)

# set working directory
wdir <- file.path("~/bayes_cop_calib")

## Define models
if (0){
#! may not need individual marginal models

# marginal efficacy model
mod_cc_e_code <- "
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
mod_cc_s_code <- "
data {
int N;
matrix[N, 2] x;
vector[N] y_s;
}
parameters {
// params for continuous (safety) outcome
vector[2] beta_s;
real<lower=0> sigma;
}
model {
vector[N] mu;

// priors
beta_s ~ normal(0,10);
sigma ~ inv_gamma(0.001,0.001);

// marginal for continuous (safety) outcome
mu = x*beta_s;
y_s ~ normal(mu, sigma);
}
"
}


# joint models with different priors
mod_cc_jnt_code1 <- "
// model 1 - flat priors

functions {
#include /cc_normCop_funs.stan
}
data {
int<lower=1> N;
matrix[N, 2] x;
real y1[N];
real y2[N];
}
parameters {
// params for margin 1 (efficacy)
vector[2] beta_e;
vector<lower=0>[2] sigma_e;

//params for margin 2 (safety)
vector[2] beta_s;
vector<lower=0>[2] sigma_s;

// copula dependence param
vector<lower=-1, upper=1>[2] omega_;  
}
model {
vector[N] mu_e;
vector[N] sigma_1;
vector[N] mu_s;
vector[N] sigma_2;
vector[N] theta;

// priors (unassigned)
// improper (flat) priors for beta_e and beta_s and sigma_e, sigma_s
// omega_ uniform(-1,1)

// margin 1 model
mu_e = x*beta_e;
sigma_1 = x*sigma_e;

// margin 2 model
mu_s = x*beta_s;
sigma_2 = x*sigma_s;

// copula dependence model
theta = x*omega_;

// build log-likelihood
{
  vector[N] loglik;  // vectorize summation
  for (i in 1:N)
    loglik[i] = bi_cop_lp(y1[i], mu_e[i], sigma_1[i], y2[i], mu_s[i], sigma_2[i], theta[i]);
  target += sum(loglik);
}

}
generated quantities {
#include /cc_normCop_genquants.stan
}
"


mod_cc_jnt_code2 <- "
// model 2

functions {
#include /cc_normCop_funs.stan
}
data {
int<lower=1> N;
matrix[N, 2] x;
real y1[N];
real y2[N];
}
parameters {
// params for margin 1 (efficacy)
vector[2] beta_e;
vector<lower=0>[2] sigma_e;

//params for margin 2 (safety)
vector[2] beta_s;
vector<lower=0>[2] sigma_s;

// copula dependence param
vector<lower=-1, upper=1>[2] omega_;  
}
model {
vector[N] mu_e;
vector[N] sigma_1;
vector[N] mu_s;
vector[N] sigma_2;
vector[N] theta;

// weakly informative priors for beta_e and beta_s (and sigma)  
beta_e ~ normal(0,1000);
beta_s ~ normal(0,1000);
sigma_e ~ inv_gamma(0.001,0.001);
sigma_s ~ inv_gamma(0.001,0.001);

// weakly informative prior for omega
omega_ ~ normal(0,0.4); 

// margin 1 model
mu_e = x*beta_e;
sigma_1 = x*sigma_e;

// margin 2 model
mu_s = x*beta_s;
sigma_2 = x*sigma_s;

// copula dependence model
theta = x*omega_;

// build log-likelihood
{
  vector[N] loglik;  // vectorize summation
  for (i in 1:N)
    loglik[i] = bi_cop_lp(y1[i], mu_e[i], sigma_1[i], y2[i], mu_s[i], sigma_2[i], theta[i]);
  target += sum(loglik);
}

}
generated quantities {
#include /cc_normCop_genquants.stan
}
"


mod_cc_jnt_code3 <- "
// model 3

functions {
#include /cc_normCop_funs.stan
}
data {
int<lower=1> N;
matrix[N, 2] x;
real y1[N];
real y2[N];
}
parameters {
// params for margin 1 (efficacy)
vector[2] beta_e;
vector<lower=0>[2] sigma_e;

//params for margin 2 (safety)
vector[2] beta_s;
vector<lower=0>[2] sigma_s;

// unscaled copula dependence param 
vector<lower=0, upper=1>[2] omega0; 
}
transformed parameters {
// scaled copula dependence param 
vector<lower=-1, upper=1>[2] omega_ = 2*(omega0 - 0.5);
}
model {
vector[N] mu_e;
vector[N] sigma_1;
vector[N] mu_s;
vector[N] sigma_2;
vector[N] theta;

// weakly informative priors for beta_e and beta_s (and sigma)  
beta_e ~ normal(0,1000);
beta_s ~ normal(0,1000);
sigma_e ~ inv_gamma(0.001,0.001);
sigma_s ~ inv_gamma(0.001,0.001);

// boundary-avoiding beta() prior on unscaled omega0
// see https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
// --> Boundary-avoiding priors for modal estimation
omega0 ~ beta(2,2); 

// margin 1 model
mu_e = x*beta_e;
sigma_1 = x*sigma_e;

// margin 2 model
mu_s = x*beta_s;
sigma_2 = x*sigma_s;

// copula dependence model
theta = x*omega_;

// build log-likelihood
{
  vector[N] loglik;  // vectorize summation
  for (i in 1:N)
    loglik[i] = bi_cop_lp(y1[i], mu_e[i], sigma_1[i], y2[i], mu_s[i], sigma_2[i], theta[i]);
  target += sum(loglik);
}

}
generated quantities {
#include /cc_normCop_genquants.stan
}
"




# alternate parameterizations of cc model (not currently used)
if (0){

  mod_cc_jnt_code1 <- "
functions {

// normal copula likelihood for 2 continuous normal margins
real bi_cop_lp(real y1, real mu1, real sigma1, real y2, real mu2, real sigma2, real rho) {
//vector[2] Z;
real Z1;
real Z2;
//vector[2] mu0;
//matrix[2,2] Rho;
real lcP;

// pseudo-obs from marginal normal dists
//real U1 = Phi((y1 - mu1)/sigma1); 
//real U2 = Phi((y2 - mu2)/sigma2); 

// get arguments for multinormal lpdf
//Z = [inv_Phi(U1), inv_Phi(U2)]';

//!! just make Z1 and Z2, no need for vector
Z1 = (y1 - mu1)/sigma1;
Z2 = (y2 - mu2)/sigma2;

//Z = [(y1 - mu1)/sigma1, (y2 - mu2)/sigma2]';
//mu0 = [0, 0]';
//Rho = [[1, rho], [rho, 1]];

// likelihood for joint dist.
// 1st term in likelihood is standard MV normal since using 
// standardized Z and centered at mu0, 2nd arg to multi_normal_lpdf is correlation matrix
//lcP = multi_normal_lpdf(Z|mu0, Rho) - (std_normal_lpdf(Z[1]) + std_normal_lpdf(Z[2])) + normal_lpdf(y1|mu1, sigma1) + normal_lpdf(y2|mu2, sigma2);
//lcP = multi_normal_lpdf(Z|mu0, Rho) - (std_normal_lpdf(Z[1]) + std_normal_lpdf(Z[2])) + log_sum_exp(normal_lpdf(y1|mu1, sigma1), normal_lpdf(y2|mu2, sigma2));

//from Joe text p. 163 & 226
//something's not quite right w/ copula density - maybe I copied down wrong?
//lcP = -0.5*log1m(rho^2) - 0.5*(Z1^2+Z2^2- 2*rho*Z1*Z2)*(1-rho^2) + 0.5*(Z1^2+Z2^2) + normal_lpdf(y1|mu1, sigma1) + normal_lpdf(y2|mu2, sigma2);

// use copula density from Meyer - bivariate normal copula paper
//lcP = -0.5*log1m(rho^2) + ((2*rho*Z1*Z2 - rho^2*(Z1^2+Z2^2))/(2*(1-rho^2))) + normal_lpdf(y1|mu1, sigma1) + normal_lpdf(y2|mu2, sigma2);

// put lower bound on cP to avoid log(0)
//real cP = cP0 < 1e-200 ? 1e-200: cP0;

// log-likelihood
return lcP;
}

}
data {
int<lower=1> N;
matrix[N, 2] x;
real y1[N];
real y2[N];
}
parameters {
// params for margin 1 (efficacy)
vector[2] beta_e;
vector<lower=0>[2] sigma_e;

//params for margin 2 (safety)
vector[2] beta_s;
vector<lower=0>[2] sigma_s;

// copula dependence param
vector<lower=-1, upper=1>[2] omega_;  
}
model {
vector[N] mu_e;
vector[N] sigma_1;
vector[N] mu_s;
vector[N] sigma_2;
vector[N] theta;

// priors (unassigned)
// improper (flat) priors for beta_e and beta_s and sigma_e, sigma_s
// dependence 
// uniform(-1,1) prior for omega
// omega_ ~ uniform(-1,1);

// weakly informative priors for sigma_e and sigma_s (std dev)
//sigma_e ~ inv_gamma(0.001,0.001);
//sigma_s ~ inv_gamma(0.001,0.001);

// margin 1 model
mu_e = x*beta_e;
sigma_1 = x*sigma_e;

// margin 2 model
mu_s = x*beta_s;
sigma_2 = x*sigma_s;

// copula dependence model
theta = x*omega_;

// build log-likelihood
{
  vector[N] loglik;  // vectorize summation
  for (i in 1:N)
    loglik[i] = bi_cop_lp(y1[i], mu_e[i], sigma_1[i], y2[i], mu_s[i], sigma_2[i], theta[i]);
  target += sum(loglik);
}

}
"
  
  
# Stan code for 2 cont. normal margins with normal copula - 
# modify to avoid Phi() and inv_Phi() transform 
mod_cc_jnt_code2 <- "
functions {

// normal copula likelihood for 2 continuous normal margins
real bi_cop_lp(real y1, real mu1, real sigma1, real y2, real mu2, real sigma2, 
matrix Rho) {
vector[2] Z;
vector[2] mu0;
real lcP;

// pseudo-obs from marginal normal dists
//real U1 = Phi((y1 - mu1)/sigma1); 
//real U2 = Phi((y2 - mu2)/sigma2); 

// get arguments for multinormal lpdf
//Z = [inv_Phi(U1), inv_Phi(U2)]';

//Z1; 
//Z2;
Z = [(y1 - mu1)/sigma1, (y2 - mu2)/sigma2]';
mu0 = [0, 0]';

// likelihood for joint dist.
// 1st term in likelihood is standard MV normal since using 
// standardized Z and centered at mu0, 2nd arg to multi_normal_lpdf is correlation matrix
lcP = multi_normal_lpdf(Z|mu0, Rho) - (std_normal_lpdf(Z[1]) + std_normal_lpdf(Z[2])) + normal_lpdf(y1|mu1, sigma1) + normal_lpdf(y2|mu2, sigma2);

// put lower bound on cP to avoid log(0)
//real cP = cP0 < 1e-200 ? 1e-200: cP0;

// log-likelihood
return lcP;
}

}
data {
int<lower=1> N;
real y_1[N];
real y_2[N];
}
parameters {
// params for margin 1
real beta_1;
real<lower=0> sigma_1;

//params for margin 2
real beta_2;
real<lower=0> sigma_2;

// copula dependence param
corr_matrix[2] omega;  
}
model {

// priors
// margin 1
beta_1 ~ normal(0,10);
sigma_1 ~ exponential(1);

// margin 2
beta_2 ~ normal(0,10); 
sigma_2 ~ exponential(1);

// dependence 
// better to work with cholesky factor of omega? fix this
omega ~ lkj_corr(1);

// margin 1 model
//mu1 = beta1;

// margin 2 model
//mu2 = beta2;

// copula dependence model
//theta = x*omega;

// build log-likelihood
{
  vector[N] loglik;  // vectorize summation
  for (i in 1:N)
  loglik[i] = bi_cop_lp(y_1[i], beta_1, sigma_1, y_2[i], beta_2, sigma_2, omega);
  target += sum(loglik);
}

}
"

# Stan code for 2 cont. normal margins with normal copula - 
# use SUR model w/ cholesky decomp. from Stan user guide (v2.18) sec 2.15
mod_cc_jnt_code_sur <- "
data {
  int<lower=1> K; // number of separate outcomes
  int<lower=1> J; // number of predictors
  int<lower=0> N; // number of observations
  vector[J] x[N]; // N x J covariate matrix
  vector[K] y[N]; // K x J observation matrix 
}
parameters {
  matrix[K, J] beta; 
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] L_sigma;
}
model {
  vector[K] mu[N];
  matrix[K, K] L_Sigma;
  
  for (n in 1:N)
    mu[n] = beta * x[n];
    
  L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
  
  to_vector(beta) ~ normal(0,10);
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ cauchy(0, 2.5);
  
  y ~ multi_normal_cholesky(mu, L_Sigma);
}
generated quantities {
//matrix[K, K] L_Sigma;
//matrix[K, K] L_Omega;
//matrix[K, K] Sigma;
matrix[K, K] Omega;

//Sigma = L_Sigma * L_Sigma';
Omega = L_Omega * L_Omega';
}
"

}

## Compile and save models

if (0){
# efficacy marginal model
mod_cc_e <- stan_model(model_code = mod_cc_e_code)
saveRDS(mod_cc_e, file = file.path(getwd(),"mod_cc_e.rds"))

# safety marginal model
mod_cc_s <- stan_model(model_code = mod_cc_s_code)
saveRDS(mod_cc_s, file = file.path(getwd(),"mod_cc_s.rds"))
}

# joint models
mod_cc_jnt1 <- stan_model(model_code = mod_cc_jnt_code1)
saveRDS(mod_cc_jnt1, file = file.path(wdir,"mod_cc_jnt1.rds"))

mod_cc_jnt2 <- stan_model(model_code = mod_cc_jnt_code2)
saveRDS(mod_cc_jnt2, file = file.path(wdir,"mod_cc_jnt2.rds"))

mod_cc_jnt3 <- stan_model(model_code = mod_cc_jnt_code3)
saveRDS(mod_cc_jnt3, file = file.path(wdir,"mod_cc_jnt3.rds"))

