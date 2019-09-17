## continuous-binary model postprocessing for copula calibration - postprocess
rm(list=ls())

libs <- c("copula", "magrittr", "rstan", "sfsmisc", "tidyverse")
invisible(lapply(libs, library, character.only = TRUE))

# set working directory
wdir<-file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/Biopharm2019/bayes_cop_calib_code")

# plot directory
pdir<-file.path(wdir,"plots")

# load continuous-binary simulation data
load(file.path(wdir,"sims", "cb_simulations.RData"))

# load different simulations scenarios 
# use cb_calib_simarray (from ACCRE) not cb_calib_simarray_local
cb_calib_simarray <- readRDS(file.path(wdir,"cb_calib_simarray.rds"))

cb_calib_simarray_trim<-cb_calib_simarray %>% select(n, theta_2, mu_e2, p_s2, rep_id, scn_id, sim_id)

cb_sim_merged<-merge(cb_sim_dat, cb_calib_simarray_trim, by=c("scn_id","rep_id","sim_id"))


# sims
# (3 true p_e2) x (3 true p_s2) x (3 trt group corr) x (50 reps)
#  (3 samp size) x (3 models) 

## make all 9 (3 samp size, 3 mod) 3x3 plots for each of the 4 params 

# make calibration plot of drawn parameter vs. posterior mean
mkplot <- function(nn, mm){
  dat<-subset(cb_sim_merged,  n==nn & mod_num==mm)
  
  add_abline <- function() {geom_abline(slope=1,intercept=0, alpha=0.3)}
  add_smooth <- function() {stat_smooth(method="loess", size=0.7, se=FALSE)}
  add_facet <- function() {facet_grid(rows=vars(mu_e2), cols=vars(p_s2),
                                      labeller = "label_both")}
  no_legend <- function() {theme(legend.position = "none")}
  
  d1 <-dat %>% filter(param=="mu_e_diff")
  d2 <-dat %>% filter(param=="p_s_diff")
  d3 <-dat %>% filter(param=="theta_1")
  d4 <-dat %>% filter(param=="theta_2")
  
  p1<-ggplot(d1,aes(x=draw_par, y=mean, col=factor(theta_2))) +
    add_abline() + add_smooth() + add_facet() + no_legend() + ggtitle("mu_e_diff")
  
  p2<-ggplot(d2,aes(x=draw_par, y=mean, col=factor(theta_2))) +
    add_abline() + add_smooth() + add_facet() + ggtitle("p_s_diff")
  
  p3<-ggplot(d3,aes(x=draw_par, y=mean, col=factor(theta_2))) +
    add_abline() + add_smooth() + add_facet() + no_legend() + ylim(-1,1) + ggtitle("theta_1")
  
  p4<-ggplot(d4,aes(x=draw_par, y=mean, col=factor(theta_2))) +
    add_abline() + add_smooth() + add_facet() + ylim(-1,1) + ggtitle("theta_2")
  
  title <- ggdraw() + 
    draw_label(paste0("continuous-binary model ",mm," (n=",nn,")"),
               fontface = 'bold', x = 0, hjust = 0) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))
  
  pg<-plot_grid(p1,p2,p3,p4, rel_widths = c(1,1.5))
  
  plot_grid(title, pg, ncol = 1,
            # rel_heights values control vertical title margins
            rel_heights = c(0.1, 1))
}

# make histogram of divergences
mkplot2 <- function(nn, mm){
  dat<-subset(cb_sim_merged,  n==nn & mod_num==mm)
  s1<-dat %>% pull(sim_id) %>% unique()
  
  cb_parm_div_dat %>% filter(sim_id %in% s1) %>% 
    ggplot(aes(x=n_div,fill=factor(theta_2))) +
    geom_histogram(bins=40)+
    facet_grid(rows=vars(mu_e2),cols=vars(p_s2),labeller="label_both")+
    ggtitle(paste0("divergences continuous-binary model ",mm," (n=",nn,")"))
  
}

mkplot(50,1)
mkplot2(50,1)


ns<-c(50,100,200)
mods<-1:3
all<-expand.grid(ns,mods)

for (i in 1:nrow(all)){
  title <- paste0("cb_mod",all[i,2],"_n",all[i,1],".pdf")
  p <- mkplot(all[i,1],all[i,2])
  p2 <- mkplot2(all[i,1],all[i,2])
  
  pdf(file.path(pdir,title),width=11,height=8.5)
  print(p)
  print(p2)
  dev.off()
}



if (0){
# for n=200, model 1, 3*3 table of mu_e_diff, p_s_diff, theta_1, theta_2
# with rows faceted by true mu_e2, cols by true p_s2

cb_sim_merged %>% filter(n==200, mod_num==1, param=="mu_e_diff") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(theta_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(mu_e2),cols=vars(p_s2))

cb_sim_merged %>% filter(n==200, mod_num==1, param=="p_s_diff") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(theta_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(mu_e2),cols=vars(p_s2))

cb_sim_merged %>% filter(n==200, mod_num==1, param=="theta_1") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(theta_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(mu_e2),cols=vars(p_s2))

cb_sim_merged %>% filter(n==200, mod_num==1, param=="theta_2") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(theta_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(mu_e2),cols=vars(p_s2))

## add divergence info
s1<-cb_sim_merged %>% filter(n==200, mod_num==1) %>% pull(sim_id) %>% unique()
cb_parm_div_dat %>% filter(sim_id %in% s1) %>% 
  ggplot(aes(x=n_div,fill=factor(theta_2))) +
  geom_histogram()+
  facet_grid(rows=vars(mu_e2),cols=vars(p_s2))
}

if (0){
  
  ## plot
  cb_sim_dat %>% filter(scn_id==2 & param=="mu_e1" & mod=="jnt") %>%
    ggplot(aes(x=draw_par, y=mean, col=factor(mod_num))) + 
    geom_point() + geom_abline(slope=1, intercept=0, alpha=0.5) +
    geom_smooth(se=FALSE)
  
  cb_sim_dat %>% filter(scn_id==2 & param=="mu_e_diff" & mod=="jnt") %>%
    ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) + 
    geom_point() + geom_abline(slope=1,intercept=0, alpha=0.5) +
    geom_smooth(se=FALSE)
  
  cb_sim_dat %>% filter(scn_id==2 & param=="theta_1" & mod=="jnt") %>%
    ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) + 
    geom_point() + geom_abline(slope=1,intercept=0, alpha=0.5) +
    geom_smooth(se=FALSE)
  
  cb_sim_dat %>% filter(scn_id==2 & param=="p_s1" & mod=="jnt") %>%
    ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) + 
    geom_point() + geom_abline(slope=1,intercept=0, alpha=0.5) +
    geom_smooth(se=FALSE)
  
  # beter to show risk difference, e.g. p_e2 -p_e1 for joint model and marginal model?  


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
