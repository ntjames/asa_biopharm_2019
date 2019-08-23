## binary-binary model for copula calibration - sim array and compile models
rm(list=ls())
 
# set working directory
# ! for ACCRE
# wdir <- file.path("~/bayes_cop_calib")

## function to calculate beta dist. parameters given mean proportion and var from
# https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
estBetaParams <- function(p, var) {
  # check if var >= p(1-p)
  if (var>=p*(1-p)) stop("must have var < p*(1-p)")
  
  alpha <- ((1 - p) / var - 1 / p) * p ^ 2
  beta <- alpha * (1 / p - 1)
  return(c(alpha = alpha, beta = beta))
}

#beta_parms<-estBetaParams(0.5,0.02)
#hist(rbeta(1000,beta_parms[1],beta_parms[2]))

## function to make random draw from (scaled) beta dist. 
## with means vec and variance var
## lb and ub are lower and upper bounds of scaled beta draws
drawBeta <- function(mn_vec, var=0.05^2, lb=0, ub=1){
  # scale from (lb,ub) to (0,1)
  mn_vec_sc <- (mn_vec - lb)/(ub - lb) 
  var_sc <- var/(ub-lb)^2
  
  betaparms <- lapply(mn_vec_sc, function(x) estBetaParams(x,var_sc))
  sapply(betaparms, function(x) rbeta(1,x[[1]],x[[2]])*(ub - lb) + lb)
} 


## function to make random draw from beta distribution with means pvec and variance v
# for prob. eff and prob AE
# beta_draw0 <-function(p_vec, v=0.02){
#   betaparms <- lapply(p_vec, function(x) estBetaParams(x,v))
#   sapply(betaparms, function(x) rbeta(1,x[[1]],x[[2]]))
# }

#beta_draw(a3[1,4:7])

#! scratch
if (0){
  beta_draw0<-function(p, v=0.02){
    p_par <- estBetaParams(p,v)
    rbeta(1,p_par$alpha,p_par$beta)
  }
  
  beta_draw0(p_e1)
  
  beta_draw1<-function(p, v=0.02){
    p1_par <- estBetaParams(p[1],v)
    p2_par <- estBetaParams(p[2],v)
    p3_par <- estBetaParams(p[3],v)
    p4_par <- estBetaParams(p[4],v)
    c(rbeta(1,p1_par[[1]],p1_par[[2]]),
      rbeta(1,p2_par[[1]],p2_par[[2]]),
      rbeta(1,p3_par[[1]],p3_par[[2]]),
      rbeta(1,p4_par[[1]],p4_par[[2]]))
  }
  
  foo<-lapply(a3[1,4:7],function(x) estBetaParams(x,v=0.02))
  
  sapply(foo,function(x) rbeta(1,x[[1]],x[[2]]))
  
  p_e1_par <- estBetaParams(p_e1,0.02)
  p_e1_dr <- rbeta(1,p_e1_par$alpha,p_e1_par$beta)
}

## function to make random draw from scaled beta dist. with means rho_vec and variance v
# for tetrachoric corr (rho) between outcomes
# omega_draw <- function(rho_vec, v=0.05^2){
#   # scale from (a,c) to (0,1)
#   a <- -1
#   c <- 1
#   rho_vec_sc <- (rho_vec - a)/(c- a) 
#   v_sc <- v/(c-a)^2
#   
#   betaparms <- lapply(rho_vec_sc, function(x) estBetaParams(x,v_sc))
#   sapply(betaparms, function(x) rbeta(1,x[[1]],x[[2]])*(c - a) + a)
# } 
# 

# verify correct scaling for beta dist.
if (0){
  #draw between (a,c) centered at true rho 
  a<- -1
  c<- 1
  
  rho <- 0.75
  v <- 0.05^2
  rho_sc <- (rho - a)/(c- a) # scale from (-1,1) to (0,1); (p + 1)/(2)
  v_sc <- v/(c-a)^2
  pars<-estBetaParams(rho_sc, v_sc)
  vals<-rbeta(5000,pars[1],pars[2])*(c - a) + a
  
  mean(vals)
  sd(vals)
  hist(vals)
}

#! scratch
if (0){ 
  # bounds for tetrachoric correlation given true proportions?
  library(polycor)
  
  n<-1e5
  polychor(rbinom(n,1,0.01),rbinom(n,1,0.99))
  
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
  nc_p <- normalCopula(0.8)
  pbo_dist <- mvdc(nc_p, margins = c("binom", "binom"),
                   paramMargins = list(list(size = 1, prob = 0.1), 
                                       list(size = 1, prob = 0.8)) )
  pbo_samps <- rMvdc(1e5, pbo_dist)
  
  rowMeans(t(pbo_samps))
  cor(pbo_samps[,1],pbo_samps[,2])
  polychor(pbo_samps[,1],pbo_samps[,2])
  
  table(pbo_samps[,1], pbo_samps[,2])
}

## make sim params array
set.seed(341884)
nreps <- 3 #! up for actual sims
nmods <- 3 # number of models (w/ diff priors) for each scenario

# placebo group parameters
# same for all
p_e1 <- 0.2 # prob effective in placebo
p_s1 <- 0.1 # prob AE in placebo
rho_1 <- 0.1 # tetrachoric corr in placebo

#! ns <- c(50, 100, 200, 400) # sample size (n per arm) 
ns <- c(150) # sample size (n per arm) 

#! rho_2s <- c(0.1, 0.35, 0.6) # tetrachoric correlation in treatment group
rho_2s <- c(0.1, 0.35) # tetrachoric correlation in treatment group

#! p_e2s <-p_e1 + c(0, 0.3, 0.6) # prob effective in trt
p_e2s <-p_e1 + c(0.3) # prob effective in trt

#! p_s2s <-p_s1 + c(0, 0.3, 0.6) # prob AE in trt
p_s2s <-p_s1 + c(0.3) # prob AE in trt

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
#saveRDS(bb_calib_simarray, file = file.path(wdir,"bb_calib_simarray.rds"))
saveRDS(bb_calib_simarray, file = file.path(getwd(),"bb_calib_simarray_local.rds"))

