## binary-binary model plots for copula calibration
rm(list=ls())

libs <- c("copula", "magrittr", "rstan", "sfsmisc", "tidyverse")
invisible(lapply(libs, library, character.only = TRUE))

# set working directory
#wdir <- file.path("~/bayes_cop_pow")
wdir<-file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/Biopharm2019/bayes_cop_calib_code")

# load binary-binary simulation data
load(file.path(wdir,"sims","bb_simulations.RData"))


# load different simulations scenarios 
# !!! need to use bb_calib_simarray (from ACCRE) not bb_calib_simarray_local
bb_calib_simarray <- readRDS(file.path(wdir,"bb_calib_simarray.rds"))


## plot

# efficacy risk difference
bb_sim_dat %>% filter(scn_id==17 & param=="p_e_diff" & mod=="jnt") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) + 
  geom_point() + geom_abline(slope=1,intercept=0, alpha=0.5) +
  geom_smooth(se=FALSE)

# safety risk difference
bb_sim_dat %>% filter(scn_id==17 & param=="p_s_diff" & mod=="jnt") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) + 
  geom_point() + geom_abline(slope=1,intercept=0, alpha=0.5) +
  geom_smooth(se=FALSE)

# placebo group correlation
bb_sim_dat %>% filter(scn_id==17 & param=="rho_1" & mod=="jnt") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) + 
  geom_point() + geom_abline(slope=1,intercept=0, alpha=0.5) +
  geom_smooth(se=FALSE)

# treatment group correlation
bb_sim_dat %>% filter(scn_id==17 & param=="rho_2" & mod=="jnt") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) + 
  geom_point() + geom_abline(slope=1,intercept=0, alpha=0.5) +
  geom_smooth(se=FALSE)

# bb_sim_dat %>% filter(scn_id==18 & param=="p_s1" & mod=="jnt") %>%
#   ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) + 
#   geom_point() + geom_abline(slope=1,intercept=0, alpha=0.5) +
#   geom_smooth(se=FALSE)

# bb_sim_dat %>% filter(scn_id==18 & param=="p_e1" & mod=="jnt") %>%
#   ggplot(aes(x=draw_par, y=mean, col=factor(mod_num))) + 
#   geom_point() + geom_abline(slope=1, intercept=0, alpha=0.5) +
#   geom_smooth(se=FALSE)


# beter to show risk difference, e.g. p_e2 -p_e1 for joint model and marginal model?  
  
  
if (0){
    x<-runif(50)
    y<-runif(50)
    
    data.frame(x,y) %>% ggplot(aes(x,y) )+geom_point()
    
    # show marginal and joint - 
    # trying to convey that marginal params can be calibrated but joint not calibrated
    
    nc_p <- normalCopula(0.3)
    dist <- mvdc(nc_p, margins = c("unif", "unif"),
                 paramMargins = list(list(min=0, max=1), 
                                     list(min=0, max=1)) )
    
    samps <- rMvdc(1000, dist) %>% data.frame()
    
    samps %>% ggplot(aes(X1,X2))+geom_point()
    #+geom_density2d()
    
    ggplot(samps, aes(X1))+geom_histogram(binwidth=0.1, boundary=0)
    
    ggplot(samps, aes(X2))+geom_histogram(binwidth=0.1, boundary=0)
    
    
    #select parameters for given sim id using SLURM_ARRAY_TASK_ID
    #arg <- commandArgs(trailing=TRUE)
    #s_id <- as.integer(arg[1])
  }
