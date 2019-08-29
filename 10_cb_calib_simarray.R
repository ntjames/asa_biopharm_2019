## continuous-binary model for copula calibration - sim array and compile models
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

##' @title Make random draws from Normal dist. given means and variance var
##' @param mn_vec a vector of means of Normal distribution
##' @param var_vec a vector of variances of Normal distribution
##' @return mu parameter
drawNorm <- function(mn_vec, var_vec){
  mapply(function(x,y) rnorm(1,x,y), mn_vec, sqrt(var_vec))
} 


if (0){ 
  # check bounds for poly-serial correlation given true mean eff. and prop. safety
  library(polycor)
  
  n<-1e2
  polyserial(rnorm(n,150,100),rbinom(n,1,0.6))
  
  rho_bounds<-function(p,q){
    r<-c(max(0,p+q-1),min(p,q))
    s1<-1-p-q+r[1]
    s2<-1-p-q+r[2]
    
    lb1 <- 6*r[1]*s1*(r[1]+s1)-1
    ub1 <- 1-6*(p-r[1])*(q-r[1])*(1-r[1]-s1)
    
    lb2 <- 6*r[2]*s2*(r[2]+s2)-1
    ub2 <- 1-6*(p-r[2])*(q-r[2])*(1-r[2]-s2)
    
    return(list(c(lb1,ub1),c(lb2,ub2)))
  }
  
  rho_bounds(0.6,0.8)
  
  library(copula)
  nc_p <- normalCopula(0.4)
  pbo_dist <- mvdc(nc_p, margins = c("norm", "binom"),
                   paramMargins = list(list(mean = -150, sd = 100), 
                                       list(size = 1, prob = 0.8)) )
  pbo_samps <- rMvdc(1e5, pbo_dist)
  
  rowMeans(t(pbo_samps))
  cor(pbo_samps[,1],pbo_samps[,2])
  polyserial(pbo_samps[,1],pbo_samps[,2])
  
  table(pbo_samps[,1], pbo_samps[,2])
}


## make sim params array
set.seed(634884)
nreps <- 50
nmods <- 3 # number of models (w/ diff priors) for each scenario

# placebo group parameters
# same for all
mu_e1 <- -150 # mean efficacy in placebo
sigma2_e1 <- 100^2 # var efficacy in placebo
p_s1 <- 0.1 # prob AE in placebo
theta_1 <- 0.1 # poly-serial corr in placebo

ns <- c(50, 100, 200) # sample size (n per arm) 

theta_2s <- c(0.1, 0.35, 0.6) # poly-serial correlation in treatment group

mu_e2s <-mu_e1 + c(0, 100, 150) # effectiveness in trt
sigma2_e2s <- 100^2

p_s2s <-p_s1 + c(0, 0.3, 0.6) # prob AE in trt

# build param scenarios data frame
a1 <- expand.grid(ns,theta_1,theta_2s, mu_e1, sigma2_e1, mu_e2s, sigma2_e2s, p_s1, p_s2s)
names(a1) <- c("n","theta_1","theta_2","mu_e1","sigma2_e1","mu_e2","sigma2_e2","p_s1","p_s2")

# replicate each scenario nreps times
a2 <- do.call("rbind", replicate(nreps, a1, simplify = FALSE))

# reorder 
a3 <- a2[order(a2$n,a2$theta_2,a2$mu_e2,a2$p_s2),]

## add draws from true dist.
p_draws0 <- sapply(1:nrow(a3), function(x) drawBeta(a3[x,8:9],var=0.15^2))
p_draws <- t(p_draws0)
colnames(p_draws) <- paste0(rownames(p_draws0),"_tr")

mu_draws0 <- sapply(1:nrow(a3), function(x) drawNorm(a3[x,c(4,6)], a3[x,c(5,7)]))
mu_draws <- t(mu_draws0)
colnames(mu_draws) <- paste0(rownames(mu_draws0),"_tr")

theta_draws0 <- sapply(1:nrow(a3), function(x) drawBeta(a3[x,2:3],var=0.05^2,lb=-1,ub=1))
theta_draws <- t(theta_draws0)
colnames(theta_draws) <- paste0(rownames(theta_draws0),"_tr")

a4 <- cbind(a3,mu_draws,p_draws,theta_draws)

# add rep id and scenario id
a4$rep_id <- rep(1:nreps,nrow(a1))
a4$scn_id <- sort(rep(1:nrow(a1),nreps))

# randomize order
a5 <- a4[sample(1:nrow(a4)),]

# add seeds 
a5$dat_seed <- sample.int(1e7,size=nrow(a5))
a5$samp_seed <- sample.int(1e7,size=nrow(a5))

# duplicate scenarios for different models
cb_calib_simarray <- do.call("rbind", replicate(nmods, a5, simplify = FALSE))
cb_calib_simarray$mod_num<-sort(rep(1:nmods,nrow(a5)))

# add simulation id
cb_calib_simarray$sim_id <- 1:nrow(cb_calib_simarray)

# save simulation array
saveRDS(cb_calib_simarray, file = file.path(wdir,"cb_calib_simarray.rds"))
