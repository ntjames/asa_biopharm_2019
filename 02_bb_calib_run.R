## binary-binary model for copula calibration - simulate data and sample from posterior 
rm(list=ls())

libs <- c("copula", "magrittr", "rstan", "sfsmisc")
invisible(lapply(libs, library, character.only = TRUE))

sessionInfo()
head(Sys.cpuinfo(),19)

# set working directory
wdir <- file.path("~/bayes_cop_calib")

# load sim array
bb_calib_simarray <- readRDS(file.path(wdir,"bb_calib_simarray.rds"))

# select parameters for given sim id using SLURM_ARRAY_TASK_ID
arg <- commandArgs(trailing=TRUE)
s_id <- as.integer(arg[1])
sim_parm <- bb_calib_simarray[bb_calib_simarray$sim_id == s_id,]

### Simulate data ###
# NB: can estimate tetrachoric correlation in samples with polycor::polychor

##' @title Make bivariate samples given marg. binomial probs and Normal copula parameter
##' @param n number of samples
##' @param p_1 probability of success for 1st marginal dist. 
##' @param p_2 probability of success for 2nd marginal dist. '
##' @param rho Normal copula dependence parameter
##' @return n x 2 matrix
mk_samps<-function(n, p_1, p_2, rho){
nc <- normalCopula(rho)
dist <- mvdc(nc, margins = c("binom", "binom"),
             paramMargins = list(list(size = 1, prob = p_1), 
                                 list(size = 1, prob = p_2)) )

# n random draws from multivariate dist.
rMvdc(n, dist)   
}

# set seed
set.seed(sim_parm$dat_seed)

# number of samples per arm
n <- sim_parm$n

## placebo group
pbo_samps <- mk_samps(n, sim_parm$p_e1_tr, sim_parm$p_s1_tr, sim_parm$rho_1_tr)

## treatment group
trt_samps <- mk_samps(n, sim_parm$p_e2_tr, sim_parm$p_s2_tr, sim_parm$rho_2_tr)

#combine placebo and treatment data
dat_bb <- rbind(pbo_samps,trt_samps) %>% 
          cbind(sort(rep(c(0,1),n)),
              sort(rep(c(0,1),n),decreasing=TRUE),
              sort(rep(c(0,1),n))) %>% 
          as.data.frame() 

names(dat_bb) <- c("efficacy","safety","treatment","trt1","trt2")

dat_bb_short <- plyr::count(dat_bb,vars=names(dat_bb))

# format data for Stan
if (0){ # not currently used
mod_data_bb <- list(N=nrow(dat_bb), 
                    x=dat_bb[,c("trt1","trt2")], 
                    y_e=dat_bb$efficacy, 
                    y_s=dat_bb$safety)
}

mod_data_bb_short <- list(N=nrow(dat_bb_short),
                          x=dat_bb_short[,c("trt1","trt2")],
                          y_e=dat_bb_short$efficacy,
                          y_s=dat_bb_short$safety,
                          freq=dat_bb_short$freq)

## run Stan models

# MCMC parameters
options(mc.cores = parallel::detectCores())
n_chains <- 2
n_warmup <- 5000
n_iter <- n_warmup + 2500

samp_seed <- sim_parm$samp_seed

## Load pre-compiled models from 0_bbmod_calib_init.R ##

# define which model should be used
# 1 - flat priors; 2 - ...
mod_num <- sim_parm$mod_num
mod_path <- paste0("mod_bb_jnt",mod_num,".rds")
mod_bb_jnt <- readRDS(file.path(wdir,mod_path))

## Sample from compiled model ##

## efficacy marginal model
if (0){
fit_bb_e <- sampling(mod_bb_e, data=mod_data_bb, seed=samp_seed, 
                 iter=n_iter1, warmup=n_warmup1, chains=n_chains,
                 control = list(adapt_delta = 0.95))

assign(paste0("summ_bb_e_",s_id), summary(fit_bb_e)$summary)

# this keeps actual samples, but may not be needed
# assign(paste0("samp_bb_e_",s_id), as.matrix(fit_bb_e, pars=c("p_e")))
}

## safety marginal model
if (0){
fit_bb_s <- sampling(mod_bb_s, data=mod_data_bb, seed=samp_seed,
                 iter=n_iter1, warmup=n_warmup1, chains=n_chains,
                 control = list(adapt_delta = 0.95))

assign(paste0("summ_bb_s_",s_id), summary(fit_bb_s)$summary)

# this keeps actual samples, but may not be needed
#assign(paste0("samp_bb_s_",s_id), as.matrix(fit_bb_s, pars=c("p_s")))

# check pairs plot
#pairs(fit_bb_s,pars=c("beta_s[1]","beta_s[2]","lp__"))
}

## joint copula model 
# efficacy marginal model MLE for initialization
mle_bb_e <- glm(efficacy~trt1+trt2-1, data=dat_bb, family=binomial(link="probit"))

# safety marginal model MLE for initialization
mle_bb_s <- glm(safety~trt1+trt2-1, data=dat_bb, family=binomial(link="probit"))

#initalize margins at jittered MLE estimate??
init_list0 <- list(beta_e=mle_bb_e$coefficients, beta_s=mle_bb_s$coefficients)
init_list <- lapply(1:n_chains, function(x) lapply(init_list0 ,jitter, amount=1.5))

fit_bb_jnt <- sampling(mod_bb_jnt, data=mod_data_bb_short, seed=samp_seed,
                       chains=n_chains, iter=n_iter, warmup=n_warmup,
                       init=init_list, control = list(adapt_delta = 0.99))

# get summary of posterior
assign(paste0("summ_bb_jnt_",s_id), summary(fit_bb_jnt)$summary)

# get number of divergences 
#!! add other diagnostic measures??
n_divs <- do.call(c, lapply(1:n_chains, function(x)
        get_sampler_params(fit_bb_jnt, inc_warmup=FALSE)[[x]][,'divergent__'])) %>%
        sum()

# keep summaries and samples, dataset (dat_bb_short), sim parameters
assign(paste0("dat_bb_short_",s_id),dat_bb_short)
assign(paste0("sim_parm_",s_id),sim_parm)
assign(paste0("n_divs_",s_id),n_divs)

savelist <- c("summ_bb_jnt", "dat_bb_short", "sim_parm", "n_divs")

save(list=paste0(savelist,"_",s_id), 
     file=file.path(wdir,"bbsims", paste0("bb_sim_",s_id,".RData")))

if (0){
  # store results in scratch
  sdir <- file.path("/gpfs23/scratch/jamesnt")
  
  # keep 3 fits, dataset (dat_bb), sim parameters
  save(list=paste0(savelist,"_",s_id), 
       file=file.path(sdir,"bbsims", paste0("bb_sim_",s_id,".RData")))
}