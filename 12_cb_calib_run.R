## continuous-binary model for copula calibration - simulate data and sample from posterior 
rm(list=ls())

libs <- c("copula", "magrittr", "rstan", "sfsmisc")
invisible(lapply(libs, library, character.only = TRUE))

sessionInfo()
head(Sys.cpuinfo(),19)

# set working directory
wdir <- file.path("~/bayes_cop_calib")

# load sim array
cb_calib_simarray <- readRDS(file.path(wdir,"cb_calib_simarray.rds"))

# select parameters for given sim id using SLURM_ARRAY_TASK_ID
arg <- commandArgs(trailing=TRUE)
s_id <- as.integer(arg[1])
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

## run Stan models

# MCMC parameters
options(mc.cores = parallel::detectCores())
n_chains <- 2
n_warmup <- 5000
n_iter <- n_warmup + 2500

samp_seed <- sim_parm$samp_seed

## Load pre-compiled models from 11_cbmod_calib_init.R ##


# define which model should be used
# 1 - flat priors; 2 - ...
mod_num <- sim_parm$mod_num
mod_path <- paste0("mod_cb_jnt",mod_num,".rds")
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
init_list <- lapply(1:n_chains, function(x) lapply(init_list0, jitter))

fit_cb_jnt <- sampling(mod_cb_jnt, data=mod_data_cb, seed=samp_seed,
                       iter=n_iter, chains=n_chains, warmup=n_warmup,
                       init=init_list, control = list(adapt_delta = 0.8))


# get summary of posterior
assign(paste0("summ_cb_jnt_",s_id), summary(fit_cb_jnt)$summary)

# get number of divergences 
#!! add other diagnostic measures??
n_divs <- do.call(c, lapply(1:n_chains, function(x)
        get_sampler_params(fit_cb_jnt, inc_warmup=FALSE)[[x]][,'divergent__'])) %>%
        sum()

# keep summaries and samples, dataset (dat_cb), sim parameters
assign(paste0("dat_cb_",s_id),dat_cb)
assign(paste0("sim_parm_",s_id),sim_parm)
assign(paste0("n_divs_",s_id),n_divs)

savelist <- c("summ_cb_jnt", "dat_cb", "sim_parm", "n_divs")

save(list=paste0(savelist,"_",s_id), 
     file=file.path(wdir,"cbsims", paste0("cb_sim_",s_id,".RData")))

if (0){
  # store results in scratch
  sdir <- file.path("/gpfs23/scratch/jamesnt")
  
  # keep 3 fits, dataset (dat_cb), sim parameters
  save(list=paste0(savelist,"_",s_id), 
       file=file.path(sdir,"cbsims", paste0("cb_sim_",s_id,".RData")))
}