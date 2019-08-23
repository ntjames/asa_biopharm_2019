## binary-binary model for copula calibration - simulate data and sample from posterior 
## LOCAL version

libs <- c("copula", "magrittr", "rstan", "sfsmisc")
invisible(lapply(libs, library, character.only = TRUE))

sessionInfo()
head(Sys.cpuinfo(),19)

# set working directory
#wdir <- file.path("~/bayes_cop_pow")
wdir<-file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/Biopharm2019/bayes_cop_calib_code")

# load sim array
#bb_calib_simarray <- readRDS(file.path(wdir,"bb_calib_simarray.rds"))

#! remove and use above for accre version
bb_calib_simarray <- readRDS(file.path(wdir,"bb_calib_simarray_local.rds"))


#select parameters for given sim id using SLURM_ARRAY_TASK_ID
arg <- commandArgs(trailing=TRUE)
s_id <- as.integer(arg[1])

#! remove and use above for accre version
#s_id <- 18

sim_params <- bb_calib_simarray[bb_calib_simarray$sim_id == s_id,]

### Simulate data ###

# function to make n_samps samples given marg. binomial probs and Normal copula rho
mk_samps<-function(n_samps, p_e, p_s, rh){
nc <- normalCopula(rh)
dist <- mvdc(nc, margins = c("binom", "binom"),
             paramMargins = list(list(size = 1, prob = p_e), 
                                 list(size = 1, prob = p_s)) )

samps <- rMvdc(n_samps, dist)   
return(samps)
}

# set seed
set.seed(sim_params$dat_seed)

# number of samples per arm
n <- sim_params$n

## placebo group
#p_e1 <- 0.2 # prob effective
#p_s1 <- 0.1 # prob AE
#rho_1 <- 0.1 # tetrachoric correlation, can estimate with polycor::polychor
pbo_samps <- mk_samps(n, sim_params$p_e1_tr, sim_params$p_s1_tr, sim_params$rho_1_tr)


#! normal copula
#! nc_p <- normalCopula(sim_params$rho1_tr)
#! pbo_dist <- mvdc(nc_p, margins = c("binom", "binom"),
#!                  paramMargins = list(list(size = 1, prob = sim_params$p_e1_tr), 
#!                                      list(size = 1, prob = sim_params$p_s1_tr)) )
#! 
#! pbo_samps <- rMvdc(n, pbo_dist)

## treatment group
#p_e2 <- sim_params$p_e2 # prob effective
#p_s2 <- sim_params$p_s2 # prob AE
#rho_2 <- sim_params$rho_2 # tetrachoric corr

trt_samps <- mk_samps(n, sim_params$p_e2_tr, sim_params$p_s2_tr, sim_params$rho_2_tr)

#! normal copula
#! nc_t <- normalCopula(sim_params$rho2_tr)
#! trt_dist <- mvdc(nc_t, margins = c("binom", "binom"),
#!                  paramMargins = list(list(size = 1, prob = sim_params$p_e2_tr), 
#!                                      list(size = 1, prob = sim_params$p_s2_tr)) )
#! 
#! trt_samps <- rMvdc(n, trt_dist)

#combine placebo and treatment data
dat_bb <- rbind(pbo_samps,trt_samps) %>% cbind(sort(rep(c(0,1),n)),
                                               sort(rep(c(0,1),n),decreasing=TRUE),
                                               sort(rep(c(0,1),n))) %>% as.data.frame() 
names(dat_bb) <- c("efficacy","safety","treatment","trt1","trt2")

dat_bb_short <- plyr::count(dat_bb,vars=names(dat_bb))

# format data for stan
mod_data_bb <- list(N=nrow(dat_bb), 
                    x=dat_bb[,c("trt1","trt2")], 
                    y_e=dat_bb$efficacy, 
                    y_s=dat_bb$safety)

mod_data_bb_short <- list(N=nrow(dat_bb_short),
                          x=dat_bb_short[,c("trt1","trt2")],
                          y_e=dat_bb_short$efficacy,
                          y_s=dat_bb_short$safety,
                          freq=dat_bb_short$freq)

## run Stan models

# MCMC parameters
options(mc.cores = parallel::detectCores())
n_chains <- 4
#n_warmup1 <- 3000
#n_iter1 <- n_warmup1+2500
n_warmup2 <- 5000
n_iter2 <- n_warmup2+1250

samp_seed <- sim_params$samp_seed

# Load pre-compiled models from bbmod_init.R
#mod_bb_e <- readRDS(file.path(wdir,"mod_bb_e.rds"))
#mod_bb_s <- readRDS(file.path(wdir,"mod_bb_s.rds"))
#mod_bb_jnt <- readRDS(file.path(wdir,"mod_bb_jnt.rds"))

#! remove and use above for accre version
#mod_bb_e <- readRDS(file.path(wdir,"mod_bb_e_local.rds"))
#mod_bb_s <- readRDS(file.path(wdir,"mod_bb_s_local.rds"))

# define which model should be used
# 1 - uniform; 2 - ...
mod_num<-sim_params$mod_num

# for ACCRE
# mod_path<-paste0("mod_bb_jnt",mod_num,".rds")

#! remove and use above for accre version
mod_path<-paste0("mod_bb_jnt",mod_num,"_local.rds")

mod_bb_jnt <- readRDS(file.path(wdir,mod_path))

## Sample from compiled models


# efficacy marginal model
if (0){
fit_bb_e <- sampling(mod_bb_e, data=mod_data_bb, seed=samp_seed, 
                 iter=n_iter1, warmup=n_warmup1, chains=n_chains,
                 control = list(adapt_delta = 0.95))

assign(paste0("summ_bb_e_",s_id), summary(fit_bb_e)$summary)

# this keeps actual samples, but may not be needed
# assign(paste0("samp_bb_e_",s_id), as.matrix(fit_bb_e, pars=c("p_e")))
}

# safety marginal model
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
                   iter=n_iter2, chains=n_chains, warmup=n_warmup2,
                   init=init_list, control = list(adapt_delta = 0.99))


#str(rstan::extract(fit_bb_jnt))

assign(paste0("summ_bb_jnt_",s_id), summary(fit_bb_jnt)$summary)

# save divergences 
# add other diagnostic measures??
divs <- do.call(c, lapply(1:n_chains, function(x)
        get_sampler_params(fit_bb_jnt, inc_warmup=FALSE)[[x]][,'divergent__']))

n_divs <- sum(divs)

#get_sampler_params(fit_bb_jnt, inc_warmup=FALSE)[[1]][,'divergent__']

# this keeps actual samples, but may not be needed
#assign(paste0("samp_bb_jnt_",s_id), as.matrix(fit_bb_jnt, pars=c("omega","p_e","p_s")))

# check pairs plot
 # pairs(fit_bb_jnt,pars=c("beta_e[1]","beta_e[2]",
 #                        "beta_s[1]","beta_s[2]",
 #                        "rho_[1]","rho_[2]","lp__"))


if (0){
# store results in scratch
sdir <- file.path("/gpfs23/scratch/jamesnt")

# keep 3 fits, dataset (dat_bb), sim parameters
save(fit_bb_e, fit_bb_s, fit_bb_jnt, dat_bb, sim_params, 
     file=file.path(sdir,"bbsims", paste0("bb_sim_",s_id,".RData")))
}


# keep summaries and samples, dataset (dat_bb_short), sim parameters
assign(paste0("dat_bb_short_",s_id),dat_bb_short)
assign(paste0("sim_params_",s_id),sim_params)
assign(paste0("n_divs_",s_id),n_divs)

# save samples
#savelist <- c("summ_bb_e", "samp_bb_e", "summ_bb_s", "samp_bb_s", 
#"summ_bb_jnt", "samp_bb_jnt", "dat_bb", "sim_params")

# save separate marginal models
#savelist <- c("summ_bb_e", "summ_bb_s","summ_bb_jnt", "dat_bb_short", "sim_params")

savelist <- c("summ_bb_jnt", "dat_bb_short", "sim_params", "n_divs")

#save(list=paste0(savelist,"_",s_id), 
#     file=file.path(wdir,"bbsims", paste0("bb_sim_",s_id,".RData")))


save(list=paste0(savelist,"_",s_id), 
     file=file.path(wdir,"sims_local","bbsims", paste0("bb_sim_",s_id,".RData")))

