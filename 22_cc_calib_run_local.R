## continuous-continuous model for copula calibration - simulate data and sample from posterior 
## LOCAL version
rm(list=ls())

libs <- c("copula", "magrittr", "rstan", "sfsmisc")
invisible(lapply(libs, library, character.only = TRUE))

sessionInfo()
head(Sys.cpuinfo(),19)

# set working directory
#wdir <- file.path("~/bayes_cop_pow")
wdir<-file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/Biopharm2019/bayes_cop_calib_code")

# load sim array
#cb_calib_simarray <- readRDS(file.path(wdir,"cb_calib_simarray.rds"))

#! remove and use above for ACCRE version
cb_calib_simarray <- readRDS(file.path(wdir,"cb_calib_simarray_local.rds"))

# select parameters for given sim id using SLURM_ARRAY_TASK_ID
arg <- commandArgs(trailing=TRUE)
s_id <- as.integer(arg[1])

#! remove and use above for ACCRE version
# s_id <- 7

sim_parm <- cb_calib_simarray[cb_calib_simarray$sim_id == s_id,]

### Simulate data ###
# NB: can estimate polyserial correlation in samples with polycor::polyserial

##' @title Make bivariate samples given marg. Normal mean and binomial prob and Normal copula parameter
##' @param n number of samples
##' @param p_1 probability of success for 1st marginal dist. 
##' @param p_2 probability of success for 2nd marginal dist. '
##' @param theta Normal copula dependence parameter
##' @return n x 2 matrix
mk_samps<-function(n, mu_1, var_1, p_2, theta){
nc <- normalCopula(theta)
dist <- mvdc(nc, margins = c("norm", "binom"),
             paramMargins = list(list(mean = mu_1, sd=sqrt(var_1)), 
                                 list(size = 1, prob = p_2)) )

# n random draws from multivariate dist.
rMvdc(n, dist)   
}

# set seed
set.seed(sim_parm$dat_seed)

# number of samples per arm
n <- sim_parm$n

## placebo group
pbo_samps <- mk_samps(n, sim_parm$mu_e1_tr,sim_parm$sigma2_e1, sim_parm$p_s1_tr, sim_parm$theta_1_tr)

## treatment group
trt_samps <- mk_samps(n, sim_parm$mu_e2_tr,sim_parm$sigma2_e2, sim_parm$p_s2_tr, sim_parm$theta_2_tr)


#combine placebo and treatment data
dat_cb <- rbind(pbo_samps,trt_samps) %>% 
          cbind(sort(rep(c(0,1),n)),
              sort(rep(c(0,1),n),decreasing=TRUE),
              sort(rep(c(0,1),n))) %>% 
          as.data.frame() 

names(dat_cb) <- c("efficacy","safety","treatment","trt1","trt2")

# format data for Stan
mod_data_cb <- list(N=nrow(dat_cb), 
                    x=dat_cb[,c("trt1","trt2")], 
                    y1=dat_cb$efficacy, 
                    y2=dat_cb$safety)

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


## run Stan models

# MCMC parameters
options(mc.cores = parallel::detectCores())
# n_chains <- 4
# n_warmup <- 5000
n_chains <- 4
n_warmup <- 5000
n_iter <- n_warmup+1250

samp_seed <- sim_parm$samp_seed

## Load pre-compiled models from 11_cbmod_calib_init.R ##
#mod_cb_e <- readRDS(file.path(wdir,"mod_cb_e.rds"))
#mod_cb_s <- readRDS(file.path(wdir,"mod_cb_s.rds"))
#mod_cb_jnt <- readRDS(file.path(wdir,"mod_cb_jnt.rds"))

#! remove and use above for accre version
#mod_cb_e <- readRDS(file.path(wdir,"mod_cb_e_local.rds"))
#mod_cb_s <- readRDS(file.path(wdir,"mod_cb_s_local.rds"))

# define which model should be used
# 1 - flat priors; 2 - ...
mod_num <- sim_parm$mod_num

# for ACCRE
# mod_path<-paste0("mod_cb_jnt",mod_num,".rds")

#! remove and use above for accre version
mod_path <- paste0("mod_cb_jnt",mod_num,"_local.rds")

mod_cb_jnt <- readRDS(file.path(wdir,mod_path))

## Sample from compiled models ##

# efficacy marginal model
if (0){
fit_cb_e <- sampling(mod_cb_e, data=mod_data_cb, seed=samp_seed, 
                 iter=n_iter1, warmup=n_warmup1, chains=n_chains,
                 control = list(adapt_delta = 0.95))

assign(paste0("summ_cb_e_",s_id), summary(fit_cb_e)$summary)

# safety marginal model
fit_cb_s <- sampling(mod_cb_s, data=mod_data_cb, seed=samp_seed,
                 iter=n_iter1, warmup=n_warmup1, chains=n_chains,
                 control = list(adapt_delta = 0.95))

assign(paste0("summ_cb_s_",s_id), summary(fit_cb_s)$summary)

}

## joint copula model
# efficacy marginal model MLE for initialization
mle_cb_e <- summary(lm(efficacy~trt1+trt2-1, data=dat_cb))
  
# safety marginal model MLE for initialization
mle_cb_s <- glm(safety~trt1+trt2-1, data=dat_cb, family=binomial(link="probit"))

#initalize margins at jittered MLE estimate??
init_list0 <- list(beta_e=mle_cb_e$coefficients[,1], beta_s=mle_cb_s$coefficients,
                   s=rep(mle_cb_e$sigma,2))
init_list <- lapply(1:n_chains, function(x) lapply(init_list0 ,jitter))





fit_cc_jnt <- sampling(mod_cc_jnt, data=mod_data_cc, seed=samp_seed,
                             iter=n_iter, chains=n_chains, warmup=n_warmup,
                             init=init_list, control = list(adapt_delta = 0.8))


fit_cc_jnt2 <- sampling(mod_cc_jnt2, data=mod_data_cc2, seed=samp_seed,
                       iter=n_iter, chains=n_chains, warmup=n_warmup,
                       init=init_list, control = list(adapt_delta = 0.8))


# get summary of posterior
assign(paste0("summ_cb_jnt_",s_id), summary(fit_cb_jnt)$summary)

# get number of divergences 
#!! add other diagnostic measures??
n_divs <- do.call(c, lapply(1:n_chains, function(x)
        get_sampler_params(fit_cb_jnt, inc_warmup=FALSE)[[x]][,'divergent__'])) %>%
        sum()

#get_sampler_params(fit_cb_jnt, inc_warmup=FALSE)[[1]][,'divergent__']

# check pairs plot
# pairs(fit_cb_jnt,pars=c("beta_e[1]","beta_e[2]",
#                        "beta_s[1]","beta_s[2]",
#                        "rho_[1]","rho_[2]","lp__"))



# keep summaries and samples, dataset (dat_cb), sim parameters
assign(paste0("dat_cb_",s_id),dat_cb)
assign(paste0("sim_parm_",s_id),sim_parm)
assign(paste0("n_divs_",s_id),n_divs)

# save samples
#savelist <- c("summ_cb_e", "samp_cb_e", "summ_cb_s", "samp_cb_s", 
#"summ_cb_jnt", "samp_cb_jnt", "dat_cb", "sim_parm")

# save separate marginal models
#savelist <- c("summ_cb_e", "summ_cb_s","summ_cb_jnt", "dat_cb_short", "sim_parm")

savelist <- c("summ_cb_jnt", "dat_cb", "sim_parm", "n_divs")

#save(list=paste0(savelist,"_",s_id), 
#     file=file.path(wdir,"cbsims", paste0("cb_sim_",s_id,".RData")))


save(list=paste0(savelist,"_",s_id), 
     file=file.path(wdir,"sims_local","cbsims", paste0("cb_sim_",s_id,".RData")))

