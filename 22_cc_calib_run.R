## continuous-continuous model for copula calibration - simulate data and sample from posterior 
rm(list=ls())

libs <- c("copula", "magrittr", "rstan", "sfsmisc")
invisible(lapply(libs, library, character.only = TRUE))

sessionInfo()
head(Sys.cpuinfo(),19)

# set working directory
wdir <- file.path("~/bayes_cop_calib")

# load sim array
cc_calib_simarray <- readRDS(file.path(wdir,"cc_calib_simarray.rds"))

# select parameters for given sim id using SLURM_ARRAY_TASK_ID
arg <- commandArgs(trailing=TRUE)
s_id <- as.integer(arg[1])

sim_parm <- cc_calib_simarray[cc_calib_simarray$sim_id == s_id,]

### Simulate data ###

##' @title Make bivariate samples given 2 marg. Normal means and vars and Normal copula parameter
##' @param n number of samples
##' @param mu_1 mean for 1st marginal dist. 
##' @param var_1 var for 1st marginal dist. 
##' @param mu_2 mean for 1st marginal dist. 
##' @param var_2 var for 1st marginal dist. 
##' @param theta Normal copula dependence parameter (Pearson correlation)
##' @return n x 2 matrix
mk_samps<-function(n, mu_1, var_1, mu_2, var_2, theta){
nc <- normalCopula(theta)
dist <- mvdc(nc, margins = c("norm", "norm"),
             paramMargins = list(list(mean = mu_1, sd=sqrt(var_1)), 
                                 list(mean = mu_2, sd=sqrt(var_2))) )

# n random draws from multivariate dist.
rMvdc(n, dist)   
}

# set seed
set.seed(sim_parm$dat_seed)

# number of samples per arm
n <- sim_parm$n

## placebo group
pbo_samps <- mk_samps(n, sim_parm$mu_e1_tr,sim_parm$sigma2_e1, 
                      sim_parm$mu_s1_tr,sim_parm$sigma2_s1, 
                      sim_parm$theta_1_tr)

## treatment group
trt_samps <- mk_samps(n, sim_parm$mu_e2_tr,sim_parm$sigma2_e2, 
                      sim_parm$mu_s2_tr,sim_parm$sigma2_s2, 
                      sim_parm$theta_2_tr)


#combine placebo and treatment data
dat_cc <- rbind(pbo_samps,trt_samps) %>% 
          cbind(sort(rep(c(0,1),n)),
              sort(rep(c(0,1),n),decreasing=TRUE),
              sort(rep(c(0,1),n))) %>% 
          as.data.frame() 

names(dat_cc) <- c("efficacy","safety","treatment","trt1","trt2")

# format data for Stan
mod_data_cc <- list(N=nrow(dat_cc), 
                    x=dat_cc[,c("trt1","trt2")], 
                    y1=dat_cc$efficacy, 
                    y2=dat_cc$safety)


if (0){ # SUR parameterization
mod_data_cc2 <- list(
  K = 2, 
  J = 1,
  N = nrow(dat2),
  x = dat2[,3,drop=F],
  y = dat2[,1:2]
)

#int<lower=1> K; // number of separate outcomes
#int<lower=1> J; // number of predictors
#int<lower=0> N; // number of observations
#vector[J] x[N]; // N x J covariate matrix
#vector[K] y[N]; // K x J observation matrix 
}


## run Stan models

# MCMC parameters
options(mc.cores = parallel::detectCores())
n_chains <- 2
n_warmup <- 5000
n_iter <- n_warmup + 2500

samp_seed <- sim_parm$samp_seed

## Load pre-compiled models from 21_ccmod_calib_init.R ##
#mod_cc_e <- readRDS(file.path(wdir,"mod_cc_e.rds"))
#mod_cc_s <- readRDS(file.path(wdir,"mod_cc_s.rds"))
#mod_cc_jnt <- readRDS(file.path(wdir,"mod_cc_jnt.rds"))

# define which model should be used
# 1 - flat priors; 2 - ...
mod_num <- sim_parm$mod_num

# for ACCRE
mod_path<-paste0("mod_cc_jnt",mod_num,".rds")

mod_cc_jnt <- readRDS(file.path(wdir,mod_path))

## Sample from compiled models ##

# efficacy marginal model
if (0){
fit_cc_e <- sampling(mod_cc_e, data=mod_data_cc, seed=samp_seed, 
                 iter=n_iter1, warmup=n_warmup1, chains=n_chains,
                 control = list(adapt_delta = 0.95))

assign(paste0("summ_cc_e_",s_id), summary(fit_cc_e)$summary)

# safety marginal model
fit_cc_s <- sampling(mod_cc_s, data=mod_data_cc, seed=samp_seed,
                 iter=n_iter1, warmup=n_warmup1, chains=n_chains,
                 control = list(adapt_delta = 0.95))

assign(paste0("summ_cc_s_",s_id), summary(fit_cc_s)$summary)

}

## joint copula model
# efficacy marginal model MLE for initialization
mle_cc_e <- summary(lm(efficacy~trt1+trt2-1, data=dat_cc))
  
# safety marginal model MLE for initialization
mle_cc_s <- summary(lm(safety~trt1+trt2-1, data=dat_cc))
  
#initalize margins at jittered MLE estimate
init_list0 <- list(beta_e=mle_cc_e$coefficients[,1], 
                   beta_s=mle_cc_s$coefficients[,1],
                   sigma_e=rep(mle_cc_e$sigma,2),
                   sigma_s=rep(mle_cc_s$sigma,2))
init_list <- lapply(1:n_chains, function(x) lapply(init_list0 ,jitter))

fit_cc_jnt <- sampling(mod_cc_jnt, data=mod_data_cc, seed=samp_seed,
                             iter=n_iter, chains=n_chains, warmup=n_warmup,
                             init=init_list, control = list(adapt_delta = 0.8))

# get summary of posterior
assign(paste0("summ_cc_jnt_",s_id), summary(fit_cc_jnt)$summary)

# get number of divergences 
#!! add other diagnostic measures??
n_divs <- do.call(c, lapply(1:n_chains, function(x)
        get_sampler_params(fit_cc_jnt, inc_warmup=FALSE)[[x]][,'divergent__'])) %>%
        sum()

#get_sampler_params(fit_cc_jnt, inc_warmup=FALSE)[[1]][,'divergent__']

# keep summaries and samples, dataset (dat_cc), sim parameters
assign(paste0("dat_cc_",s_id),dat_cc)
assign(paste0("sim_parm_",s_id),sim_parm)
assign(paste0("n_divs_",s_id),n_divs)

# save samples
savelist <- c("summ_cc_jnt", "dat_cc", "sim_parm", "n_divs")


save(list=paste0(savelist,"_",s_id), 
     file=file.path(wdir,"ccsims", paste0("cc_sim_",s_id,".RData")))

