## continuous-continuous model postprocessing for copula calibration - postprocess
rm(list=ls())

libs <- c("copula", "magrittr", "rstan", "sfsmisc", "tidyverse")
invisible(lapply(libs, library, character.only = TRUE))

# set working directory
wdir<-file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/Biopharm2019/bayes_cop_calib_code")

# plot directory
pdir<-file.path(wdir,"plots")

# load continuous-continuous simulation data
load(file.path(wdir,"sims", "cc_simulations.RData"))

# load different simulations scenarios 
# use cc_calib_simarray (from ACCRE) not cc_calib_simarray_local
cc_calib_simarray <- readRDS(file.path(wdir,"cc_calib_simarray.rds"))

cc_calib_simarray_trim<-cc_calib_simarray %>% select(n, theta_2, mu_e2, mu_s2, rep_id, scn_id, sim_id)

cc_sim_merged<-merge(cc_sim_dat, cc_calib_simarray_trim, by=c("scn_id","rep_id","sim_id"))



# sims
# (3 true p_e2) x (3 true p_s2) x (3 trt group corr) x (50 reps)
#  (3 samp size) x (3 models) 

## make all 9 (3 samp size, 3 mod) 3x3 plots for each of the 4 params 

# make calibration plot of drawn parameter vs. posterior mean
mkplot <- function(nn, mm){
  dat<-subset(cc_sim_merged,  n==nn & mod_num==mm)
  
  add_abline <- function() {geom_abline(slope=1,intercept=0, alpha=0.3)}
  add_smooth <- function() {stat_smooth(method="loess", size=0.7, se=FALSE)}
  add_facet <- function() {facet_grid(rows=vars(mu_e2), cols=vars(mu_s2),
                                      labeller = "label_both")}
  no_legend <- function() {theme(legend.position = "none")}
  
  d1 <-dat %>% filter(param=="mu_e_diff")
  d2 <-dat %>% filter(param=="mu_s_diff")
  d3 <-dat %>% filter(param=="theta_1")
  d4 <-dat %>% filter(param=="theta_2")
  
  p1<-ggplot(d1,aes(x=draw_par, y=mean, col=factor(theta_2))) +
    add_abline() + add_smooth() + add_facet() + no_legend() + ggtitle("mu_e_diff")
  
  p2<-ggplot(d2,aes(x=draw_par, y=mean, col=factor(theta_2))) +
    add_abline() + add_smooth() + add_facet() + ggtitle("mu_s_diff")
  
  p3<-ggplot(d3,aes(x=draw_par, y=mean, col=factor(theta_2))) +
    add_abline() + add_smooth() + add_facet() + no_legend() + ylim(-1,1) + ggtitle("theta_1")
  
  p4<-ggplot(d4,aes(x=draw_par, y=mean, col=factor(theta_2))) +
    add_abline() + add_smooth() + add_facet() + ylim(-1,1) + ggtitle("theta_2")
  
  title <- ggdraw() + 
    draw_label(paste0("continuous-continuous model ",mm," (n=",nn,")"),
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
  dat<-subset(cc_sim_merged,  n==nn & mod_num==mm)
  s1<-dat %>% pull(sim_id) %>% unique()
  
  cc_parm_div_dat %>% filter(sim_id %in% s1) %>% 
    ggplot(aes(x=n_div,fill=factor(theta_2))) +
    geom_histogram(bins=40)+
    facet_grid(rows=vars(mu_e2),cols=vars(mu_s2),labeller="label_both")+
    ggtitle(paste0("divergences continuous-continuous model ",mm," (n=",nn,")"))
  
}

mkplot(50,1)
mkplot2(50,1)


ns<-c(50,100,200)
mods<-1:3
all<-expand.grid(ns,mods)

for (i in 1:nrow(all)){
title <- paste0("cc_mod",all[i,2],"_n",all[i,1],".pdf")
p <- mkplot(all[i,1],all[i,2])
p2 <- mkplot2(all[i,1],all[i,2])
  
pdf(file.path(pdir,title),width=11,height=8.5)
  print(p)
  print(p2)
dev.off()
}






if (0){
# for n=200, model 1, 3*3 table of mu_e_diff, mu_s_diff, theta_1, theta_2
# with rows faceted by true mu_e2, cols by true mu_s2

cc_sim_merged %>% filter(n==200, mod_num==1, param=="mu_e_diff") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(theta_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(mu_e2),cols=vars(mu_s2))

cc_sim_merged %>% filter(n==200, mod_num==1, param=="mu_s_diff") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(theta_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(mu_e2),cols=vars(mu_s2))

cc_sim_merged %>% filter(n==200, mod_num==1, param=="theta_1") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(theta_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(mu_e2),cols=vars(mu_s2))

cc_sim_merged %>% filter(n==200, mod_num==1, param=="theta_2") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(theta_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(mu_e2),cols=vars(mu_s2))
}
