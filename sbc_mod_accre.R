
# do simulation-based calibration (from Carpenter paper) for marginal and joint model
rm(list=ls())

libs<-c("rstan", "sfsmisc")
invisible(lapply(libs, library, character.only = TRUE))

sessionInfo()
head(Sys.cpuinfo(),19)

#select parameters for given sim id using SLURM_ARRAY_TASK_ID
arg <- commandArgs(trailing=TRUE)
s_id <- as.integer(arg[1])

# set working directory
wdir <- file.path("~/bayes_cop_calib")

# @param model: precompile Stan model
# @param data: data for model (defaults to empty list)
# @return size of the generated quatity array I_lt_sim
num_monitored_params <- function(model, data=list()){
  fit <- sampling(model, data=data, iter=1, chains=1,
                  warmup=0, refresh=0, seed=4321)
  fit@par_dims$I_lt_sim
}

sbc_par <- function(model, data=list(),
                stan_sims=999,
                init_thin = 4, max_thin = 64,
                target_n_eff = 0.8 * stan_sims, 
                refresh = 0) {
  
  num_params <- num_monitored_params(model, data)
  ranks <- matrix(nrow = 1, ncol = num_params)
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
    lt_sim <- rstan::extract(fit)$I_lt_sim
    for (i in 1:num_params) { ranks[1, i] <- sum(lt_sim[, i]) +1 }
  list(rank=ranks, thin = thin)
}


# load precompiled model
mod_bb_e <- readRDS(file.path(wdir,"mod_bb_e.rds"))

# sample size and trt (move to sep file)
nsamps<-200
trt1<-sort(rep(c(1,0),nsamps/2))
trt2<-1-trt1
mod_data_bb <- list(N=nsamps, 
                    x=data.frame(cbind(trt1,trt2)))

sbc_bb_e <- sbc_par(mod_bb_e, data=mod_data_bb, refresh=4000)

savelist <- c("sbc_bb_e")

save(list=savelist, file=file.path(wdir,"sbcsims","bb_e", paste0("sbc_bb_e_",s_id,".RData")))
