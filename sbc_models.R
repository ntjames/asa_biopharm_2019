
# do simulation-based calibration for marginal and joint model
# (from Carpenter 2019 Simulation-Based Calibration with Stan and RStan) 
# no data --> draw from prior, 


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


# precompile model
mod_bb_e <-stan_model(model_code = mod_bb_e_code)


# @param model: precompile Stan model
# @param data: data for model (defaults to empty list)
# @return size of the generated quatity array I_lt_sim
num_monitored_params <- function(model, data=list()){
  fit <- sampling(model, data=data, iter=1, chains=1,
                  warmup=0, refresh=0, seed=4321)
  fit@par_dims$I_lt_sim
}


if (0){
nsamps<-20
trt1<-sort(rep(c(1,0),nsamps/2))
trt2<-1-trt1
mod_data_bb <- list(N=nsamps, 
                    x=data.frame(cbind(trt1,trt2)))


num_monitored_params(mod_bb_e, mod_data_bb)


fit<-sampling(mod_bb_e, data=mod_data_bb, chains=1, iter=4000,
         control=list(adapt_delta=0.99))

ranks <- matrix(nrow=100,ncol=2)

lt_sim <- rstan::extract(fit)$I_lt_sim

for (i in 1:2)
  ranks[1, i] <- sum(lt_sim[, i]) +1
}


sbc <- function(model, data=list(),
                sbc_sims=1000, stan_sims=999,
                init_thin = 4, max_thin = 64,
                target_n_eff = 0.8 * stan_sims, 
                refresh = 0) {
  num_params <- num_monitored_params(model, data)
  ranks <- matrix(nrow = sbc_sims, ncol = num_params)
  thins <- rep(NA, sbc_sims)
  for (n in 1:sbc_sims){
    n_eff <- 0
    thin <- init_thin
    while (TRUE) {
      fit <- sampling(model, 
                      data = data,
                      chains = 1,
                      iter = 2 * thin * stan_sims,
                      thin = thin,
                      control = list(adapt_delta=0.99),
                      refresh = refresh)
      fit_summary <- summary(fit, pars = c("lp__"), probs=c())$summary
      n_eff <- fit_summary["lp__", "n_eff"]
      if (n_eff >= target_n_eff || (2*thin)> max_thin) break;
      thin <- 2* thin
    }
    thins[n] <- thin
    lt_sim <- rstan::extract(fit)$I_lt_sim
    for (i in 1::num_params)
      ranks[n, i] <- sum(lt_sim[, i]) +1
  }
  
  list(rank=ranks, thin = thins)
}


sbc(mod_bb_e, data=mod_data_bb, sbc_sims=5, refresh=2000)

# parallelize so can be used on cluster

sbc_par <- function(model, data=list(), sbc_sim_num=1,
                stan_sims=999,
                init_thin = 4, max_thin = 64,
                target_n_eff = 0.8 * stan_sims, 
                refresh = 0) {
  
  num_params <- num_monitored_params(model, data)
  ranks <- matrix(nrow = 1, ncol = num_params)
 # thins <- rep(NA, sbc_sims)
 # for (n in 1:sbc_sims){
    n_eff <- 0
    thin <- init_thin
    while (TRUE) {
      fit <- sampling(model, 
                      data = data,
                      chains = 1,
                      iter = 2 * thin * stan_sims,
                      thin = thin,
                      control = list(adapt_delta=0.99),
                      refresh = refresh)
      fit_summary <- summary(fit, pars = c("lp__"), probs=c())$summary
      n_eff <- fit_summary["lp__", "n_eff"]
      if (n_eff >= target_n_eff || (2*thin)> max_thin) break;
      thin <- 2* thin
    }
   # thins[n] <- thin
    lt_sim <- rstan::extract(fit)$I_lt_sim
    for (i in 1:num_params) { ranks[1, i] <- sum(lt_sim[, i]) +1 }
  list(rank=ranks, thin = thin)
}


sbc_sims<-replicate(2, sbc_par(mod_bb_e, data=mod_data_bb, refresh=2000), simplify=FALSE)

sbc_sims<-sapply(1:2, sbc_par(mod_bb_e, data=mod_data_bb, refresh=2000), simplify=FALSE)


library(parallel)
mycores <- detectCores() - 1
clust <- makeCluster(mycores, type="FORK")

sbc_sims_par<-parSapplyLB(cl = clust, 1:3, function(x) 
                          sbc_par(mod_bb_e, data=mod_data_bb, refresh=2000), simplify=FALSE)

stopCluster(clust)


set.seed(134)
pn<-4
probs<-pnorm(rnorm(pn,0,10))
n <- 2e5
mat<-matrix(rbinom(n,1,prob=probs), nrow=n/pn, byrow = TRUE)
colMeans(mat)
probs
