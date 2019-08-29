## binary-binary model for copula calibration - sim array and compile models
rm(list=ls())
 
# set working directory
wdir <- file.path("~/bayes_cop_calib")

##' @title Calculate Beta distribution parameters given mean and variance
##' @param mu mean of Beta distribution
##' @param var variance of Beta distribution
##' @return alpha and beta parameters
##' @references https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
estBetaParams <- function(mu, var) {
  # check if var >= mu(1-mu)
  if (var>=mu*(1-mu)) stop("must have var < mu*(1-mu)")
  
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(c(alpha = alpha, beta = beta))
}

##' @title Make random draws from (scaled) Beta dist. given means and variance var
##' @param mn_vec a vector of means of (scaled) Beta distribution
##' @param var scalar variance of (scaled) Beta distribution
##' @param lb lower bound of scaled Beta dist (default=0)
##' @param ub upper bound of scaled Beta dist (default=1)
##' @return alpha and beta parameters
drawBeta <- function(mn_vec, var, lb=0, ub=1){
  # scale from (lb,ub) to (0,1)
  mn_vec_sc <- (mn_vec - lb)/(ub - lb) 
  var_sc <- var/(ub-lb)^2
  
  betaparms <- lapply(mn_vec_sc, function(x) estBetaParams(x,var_sc))
  sapply(betaparms, function(x) rbeta(1,x[[1]],x[[2]])*(ub - lb) + lb)
} 

## make sim params array
set.seed(261484)
nreps <- 50 #! up for actual sims
nmods <- 3 # number of models (w/ diff priors) for each scenario

# placebo group parameters
# same for all
p_e1 <- 0.2 # prob effective in placebo
p_s1 <- 0.1 # prob AE in placebo
rho_1 <- 0.1 # tetrachoric corr in placebo

ns <- c(50, 100, 200) # sample size (n per arm) 

rho_2s <- c(0.1, 0.35, 0.6) # tetrachoric correlation in treatment group

p_e2s <-p_e1 + c(0, 0.3, 0.6) # prob effective in trt

p_s2s <-p_s1 + c(0, 0.3, 0.6) # prob AE in trt

# build param scenarios data frame
a1 <- expand.grid(ns,rho_1,rho_2s,p_e1,p_e2s,p_s1,p_s2s)
names(a1) <- c("n","rho_1","rho_2","p_e1","p_e2","p_s1","p_s2")

# replicate each scenario nreps times
a2 <- do.call("rbind", replicate(nreps, a1, simplify = FALSE))

# reorder 
a3 <- a2[order(a2$n,a2$rho_2,a2$p_e2,a2$p_s2),]

## add draws from true dist.
p_draws0 <- sapply(1:nrow(a3), function(x) drawBeta(a3[x,4:7],var=0.15^2))
p_draws <- t(p_draws0)
colnames(p_draws) <- paste0(rownames(p_draws0),"_tr")

rho_draws0 <- sapply(1:nrow(a3), function(x) drawBeta(a3[x,2:3],var=0.05^2,lb=-1,ub=1))
rho_draws <- t(rho_draws0)
colnames(rho_draws) <- paste0(rownames(rho_draws0),"_tr")

a4 <- cbind(a3,p_draws,rho_draws)

# add rep id and scenario id
a4$rep_id <- rep(1:nreps,nrow(a1))
a4$scn_id <- sort(rep(1:nrow(a1),nreps))

# randomize order
a5 <- a4[sample(1:nrow(a4)),]

# add seeds 
a5$dat_seed <- sample.int(1e7,size=nrow(a5))
a5$samp_seed <- sample.int(1e7,size=nrow(a5))

# duplicate scenarios for different models
bb_calib_simarray <- do.call("rbind", replicate(nmods, a5, simplify = FALSE))
bb_calib_simarray$mod_num<-sort(rep(1:nmods,nrow(a5)))

# add simulation id
bb_calib_simarray$sim_id <- 1:nrow(bb_calib_simarray)

# save simulation array
saveRDS(bb_calib_simarray, file = file.path(wdir,"bb_calib_simarray.rds"))
