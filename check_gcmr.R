# try to fit model using gcmr
# NB! this can't account for different correlation in placebo & treatment grp

library(gcmr)
library(copula)
library(dplyr)
library(tidyr)
library(magrittr)

# example
# data(epilepsy)
# gcmr(counts ~ offset(log(time)) + visit + trt + visit:trt, data = epilepsy, 
#      subset = (id != 49), marginal = negbin.marg, cormat = cluster.cormat(id, "ar1"), 
#      options=gcmr.options(seed=123, nrep=100 ))

# for Linux
#bb_calib_simarray <- readRDS(file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/Biopharm2019/bayes_cop_calib_code/bb_calib_simarray_local.rds"))

# for Windows
bb_calib_simarray <- readRDS(file.path("C:/Users/nj115/Dropbox/njames/school/PhD/misc/conferences/Biopharm2019/bayes_cop_calib_code/bb_calib_simarray_local.rds"))


### Simulate data ###
s_id <- 1
sim_params <- bb_calib_simarray[bb_calib_simarray$sim_id == s_id,]

# set seed
set.seed(sim_params$dat_seed)

# number of samples per arm
n <- sim_params$n

# tetrachoric correlation, can estimate with polycor::polychor
# only have one rho since gcmr can't adjust correlation based on covariates. !! check this
rho_fixed <- 0.2 

# normal copula
nc <- normalCopula(rho_fixed)

## placebo group
pbo_dist <- mvdc(nc, margins = c("binom", "binom"),
                 paramMargins = list(list(size = 1, prob = sim_params$p_e1_tr), 
                                     list(size = 1, prob = sim_params$p_s1_tr)) )

pbo_samps <- rMvdc(n, pbo_dist)

## treatment group
trt_dist <- mvdc(nc, margins = c("binom", "binom"),
                 paramMargins = list(list(size = 1, prob = sim_params$p_e2_tr), 
                                     list(size = 1, prob = sim_params$p_s2_tr)) )

trt_samps <- rMvdc(n, trt_dist)

#combine placebo and treatment data
dat_bb <- rbind(pbo_samps,trt_samps) %>% cbind(sort(rep(c(0,1),n)),
                                               sort(rep(c(0,1),n),decreasing=TRUE),
                                               sort(rep(c(0,1),n))) %>% as.data.frame() 
names(dat_bb) <- c("efficacy","safety","treatment","trt1","trt2")

# format data for gcmr
gcmr_dat_bb <- dat_bb %>% select(efficacy,safety,treatment) %>% 
  mutate(id=1:nrow(dat_bb)) %>%
  gather("type", value="outcome", efficacy, safety) %>%
  arrange(id)

# fit with gcmr
gcmr.fit <- gcmr(outcome~treatment*type, data = gcmr_dat_bb, 
          marginal = binomial.marg(link="probit"), 
          cormat = cluster.cormat(id, "unstructured"), 
          options=gcmr.options(seed=123, nrep=100))

# get probability estimates
probs<-pnorm( c(
              coef(gcmr.fit)[1],
              sum( coef(gcmr.fit)[1:2] ),
              sum( coef(gcmr.fit)[c(1,3)] ),
              sum( coef(gcmr.fit)[1:4] ) ))

names(probs)<-c("p_e1","p_e2","p_s1","p_s2")


eff, placebo 1
0, 0, 0 

eff, trt 1+2
0, 1, 0

saf, placebo 1+3
1, 0, 0 

saf, trt 1+2+3+4

1, 1, 1