## binary-binary model for copula calibration - sim array and compile models
rm(list=ls())
 
# set working directory
# ! for ACCRE
# wdir <- file.path("~/bayes_cop_calib")

## make sim array
set.seed(241884)
nreps <- 10 #! up for actual sims

# same for all
p_e1 <- 0.2 # prob effective in placebo
p_s1 <- 0.1 # prob AE in placebo
rho_1 <- 0.1 # tetrachoric corr in placebo

# function to calculate beta dist. given mean proportion and var from
# https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
estBetaParams <- function(p, var) {
  # error check if var >= p(1-p)
  if (var>=p*(1-p)) stop("must have var < p*(1-p)")
  
  alpha <- ((1 - p) / var - 1 / p) * p ^ 2
  beta <- alpha * (1 / p - 1)
  return(c(alpha = alpha, beta = beta))
}

#beta_parms<-estBetaParams(0.5,0.02)
#hist(rbeta(1000,beta_parms[1],beta_parms[2]))


#! ns <- c(50, 100, 200, 400) # sample size (n per arm) 
ns <- c(150,250) # sample size (n per arm) 

#! rho_2s <- c(0.1, 0.35, 0.6) # tetrachoric correlation in treatment group
rho_2s <- c(0.35) # tetrachoric correlation in treatment group

#! p_e2s <-p_e1 + c(0, 0.3, 0.6) # prob effective in trt
p_e2s <-p_e1 + c(0.3) # prob effective in trt

#! p_s2s <-p_s1 + c(0, 0.3, 0.6) # prob AE in trt
p_s2s <-p_s1 + c(0.3) # prob AE in trt


#expand.grid(rho_2s,p_e2s,p_s2s)
a1 <- expand.grid(ns,rho_1,rho_2s,p_e1,p_e2s,p_s1,p_s2s)
names(a1) <- c("n","rho_1","rho_2","p_e1","p_e2","p_s1","p_s2")

# replicate each scenario nreps times
a2 <- do.call("rbind", replicate(nreps, a1, simplify = FALSE))

# reorder 
a3<-a2[order(a2$n,a2$rho_2,a2$p_e2,a2$p_s2),]

# add draws from true dist.

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

beta_draw<-function(p_vec, v=0.02){
  betaparms<-lapply(p_vec, function(x) estBetaParams(x,v))
  sapply(betaparms, function(x) rbeta(1,x[[1]],x[[2]]))
}

#beta_draw(a3[1,4:7])

p_draws0<-sapply(1:nrow(a3), function(x) beta_draw(a3[x,4:7]))
p_draws<-t(p_draws0)

colnames(p_draws)<-paste0(rownames(p_draws0),"_tr")

a4<-cbind(a3,p_draws)

# add rep id and scenario id
a4$rep_id<-rep(1:nreps,nrow(a1))
a4$scn_id<-sort(rep(1:nrow(a1),nreps))

# randomize order
bb_calib_simarray <- a4[sample(1:nrow(a4)),]

# add seeds and simulation id
bb_calib_simarray$dat_seed <- sample.int(1e7,size=nrow(bb_calib_simarray))
bb_calib_simarray$samp_seed <- sample.int(1e7,size=nrow(bb_calib_simarray))
bb_calib_simarray$sim_id <- 1:nrow(bb_calib_simarray)



# save simulation array
#saveRDS(bb_calib_simarray, file = file.path(wdir,"bb_calib_simarray.rds"))
saveRDS(bb_calib_simarray, file = file.path(getwd(),"bb_calib_simarray_local.rds"))

