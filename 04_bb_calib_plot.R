## binary-binary model plots for copula calibration
rm(list=ls())

libs <- c("magrittr", "cowplot", "tidyverse")
invisible(lapply(libs, library, character.only = TRUE))

# set working directory
wdir<-file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/Biopharm2019/bayes_cop_calib_code")

# plot directory
pdir<-file.path(wdir,"plots")


# load binary-binary simulation data
load(file.path(wdir,"sims","bb_simulations.RData"))

# load different simulations scenarios 
# use bb_calib_simarray (from ACCRE) not bb_calib_simarray_local
bb_calib_simarray <- readRDS(file.path(wdir,"bb_calib_simarray.rds"))

bb_calib_simarray_trim<-bb_calib_simarray %>% select(n, rho_2, p_e2, p_s2, rep_id, scn_id, sim_id)

bb_sim_merged<-merge(bb_sim_dat, bb_calib_simarray_trim, by=c("scn_id","rep_id","sim_id"))


# sims
# (3 true p_e2) x (3 true p_s2) x (3 trt group corr) x (50 reps)
#  (3 samp size) x (3 models) 

## make all 9 (3 samp size, 3 mod) 3x3 plots for each of the 4 params 

# make calibration plot of drawn parameter vs. posterior mean
mkplot <- function(nn, mm){
  dat<-subset(bb_sim_merged,  n==nn & mod_num==mm)
  
  add_abline <- function() {geom_abline(slope=1,intercept=0, alpha=0.3)}
  add_smooth <- function() {stat_smooth(method="loess", size=0.7, se=FALSE)}
  add_facet <- function() {facet_grid(rows=vars(p_e2), cols=vars(p_s2),
                                      labeller = "label_both")}
  no_legend <- function() {theme(legend.position = "none")}
  
  d1 <-dat %>% filter(param=="p_e_diff")
  d2 <-dat %>% filter(param=="p_s_diff")
  d3 <-dat %>% filter(param=="rho_1")
  d4 <-dat %>% filter(param=="rho_2")
  
  p1<-ggplot(d1,aes(x=draw_par, y=mean, col=factor(rho_2))) +
      add_abline() + add_smooth() + add_facet() + no_legend() + ggtitle("p_e_diff")
  
  p2<-ggplot(d2,aes(x=draw_par, y=mean, col=factor(rho_2))) +
    add_abline() + add_smooth() + add_facet() + ggtitle("p_s_diff")
  
  p3<-ggplot(d3,aes(x=draw_par, y=mean, col=factor(rho_2))) +
    add_abline() + add_smooth() + add_facet() + no_legend() + ylim(-1,1) + ggtitle("rho_1")
  
  p4<-ggplot(d4,aes(x=draw_par, y=mean, col=factor(rho_2))) +
    add_abline() + add_smooth() + add_facet() + ylim(-1,1) + ggtitle("rho_2")
  
  title <- ggdraw() + 
    draw_label(paste0("binary-binary model ",mm," (n=",nn,")"),
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
  dat<-subset(bb_sim_merged,  n==nn & mod_num==mm)
  s1<-dat %>% pull(sim_id) %>% unique()
  
  bb_parm_div_dat %>% filter(sim_id %in% s1) %>% 
    ggplot(aes(x=n_div,fill=factor(rho_2))) +
    geom_histogram(bins=40)+
    facet_grid(rows=vars(p_e2),cols=vars(p_s2),labeller="label_both")+
    ggtitle(paste0("divergences binary-binary model ",mm," (n=",nn,")"))
  
}
  
ns<-c(50,100,200)
mods<-1:3
all<-expand.grid(ns,mods)

for (i in 1:nrow(all)){
title <- paste0("bb_mod",all[i,2],"_n",all[i,1],".pdf")
p <- mkplot(all[i,1],all[i,2])
p2 <- mkplot2(all[i,1],all[i,2])

pdf(file.path(pdir,title),width=11,height=8.5)
  print(p)
  print(p2)
dev.off()
}



if(0){
  
  t1 <-bb_sim_merged %>% filter(n==50, mod_num==1)
  t2 <-bb_sim_merged %>% filter(n==100, mod_num==1)
  
  t1 %>% filter(param=="p_e_diff") %>%
    ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
    geom_abline(slope=1,intercept=0, alpha=0.3) +
    stat_smooth(method="loess", size=0.7, se=FALSE)+
    facet_grid(rows=vars(p_e2),cols=vars(p_s2))
  
  t2 %>% filter(param=="p_e_diff") %>%
    ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
    geom_abline(slope=1,intercept=0, alpha=0.3) +
    stat_smooth(method="loess", size=0.7, se=FALSE)+
    facet_grid(rows=vars(p_e2),cols=vars(p_s2))

pdf(file.path(pdir,"bb_mod1_n50.pdf"))
mkplot(50,1)
dev.off()

mkplot(50,2)
mkplot(50,3)

mkplot(100,1)
mkplot(100,2)
mkplot(100,3)

mkplot(200,1)
mkplot(200,2)
mkplot(200,3)

# for n=200, model 1, 3x3 table of p_e_diff, p_s_diff, rho_1, rho_2
# with rows faceted by true p_e2, cols by true p_s2


bb_sim_merged %>% filter(n==200, mod_num==1, param=="p_e_diff") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(p_e2),cols=vars(p_s2))

bb_sim_merged %>% filter(n==200, mod_num==1, param=="p_s_diff") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(p_e2),cols=vars(p_s2))

bb_sim_merged %>% filter(n==200, mod_num==1, param=="rho_1") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(p_e2),cols=vars(p_s2))

bb_sim_merged %>% filter(n==200, mod_num==1, param=="rho_2") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(p_e2),cols=vars(p_s2))+ ylim(-0.75,1)

## add divergence info
s1<-bb_sim_merged %>% filter(n==200, mod_num==1) %>% pull(sim_id) %>% unique()
bb_parm_div_dat %>% filter(sim_id %in% s1) %>% 
  ggplot(aes(x=n_div,fill=factor(rho_2))) +
  geom_histogram()+
  facet_grid(rows=vars(p_e2),cols=vars(p_s2))

# n=200 model 2

bb_sim_merged %>% filter(n==200, mod_num==2, param=="p_e_diff") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(p_e2),cols=vars(p_s2))

bb_sim_merged %>% filter(n==200, mod_num==2, param=="p_s_diff") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(p_e2),cols=vars(p_s2))

bb_sim_merged %>% filter(n==200, mod_num==2, param=="rho_1") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(p_e2),cols=vars(p_s2))

bb_sim_merged %>% filter(n==200, mod_num==2, param=="rho_2") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(p_e2),cols=vars(p_s2))+ ylim(-0.75,1)




# n=200 model 3

bb_sim_merged %>% filter(n==200, mod_num==3, param=="p_e_diff") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(p_e2),cols=vars(p_s2))

bb_sim_merged %>% filter(n==200, mod_num==3, param=="p_s_diff") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(p_e2),cols=vars(p_s2))

bb_sim_merged %>% filter(n==200, mod_num==3, param=="rho_1") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(p_e2),cols=vars(p_s2))

bb_sim_merged %>% filter(n==200, mod_num==3, param=="rho_2") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE) +
  facet_grid(rows=vars(p_e2),cols=vars(p_s2)) + ylim(-0.75,1)
}


if (0){

bb_sim_merged %>% filter(n==100, mod_num==1, param=="p_e_diff") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(rho_2))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE)+
  facet_grid(rows=vars(p_e2),cols=vars(p_s2))
p1_ids<-bb_calib_simarray %>% filter(n==100, rho_2==0.1) %>% pull(sim_id)

#bar <- bb_calib_simarray %>% filter(n==100, rho_2==0.1, mod_num==1)

# plot by 3 x 3 matrix p_e_diff
p1_parms <- c("p_e_diff")
#p2_parms <- c("rho_1","rho_2")

bb_sim_dat %>% filter(sim_id %in% p1_ids,param =="p_e_diff") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE) +
 # stat_smooth(aes(x=draw_par, y=X2.5.), size=0.4, method="loess", se=FALSE) +
 #  stat_smooth(aes(x=draw_par, y=X97.5.), size=0.4, method="loess", se=FALSE) +
  facet_grid(rows=vars(true_par_mn),cols=vars(param))
  

bb_sim_dat %>% filter(sim_id %in% p1_ids,param =="p_s_diff") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE) +
#  stat_smooth(aes(x=draw_par, y=X2.5.), size=0.4, method="loess", se=FALSE) +
#  stat_smooth(aes(x=draw_par, y=X97.5.), size=0.4, method="loess", se=FALSE) +
  facet_grid(rows=vars(true_par_mn),cols=vars(param))


p2_ids<-bb_calib_simarray %>% mutate(p_e_diff=p_e2-p_e1,p_s_diff=p_s2-p_s1) %>% 
  filter(n==100, p_e_diff==0, p_s_diff==0) %>% pull(sim_id)


bb_calib_simarray %>% mutate(p_e_diff=p_e2-p_e1,p_s_diff=p_s2-p_s1) %>% 
  filter(n==100, p_e_diff==0, p_s_diff==0) %>% pull(scn_id) %>% table()

bb_sim_dat %>% filter(sim_id %in% p2_ids, param =="rho_2") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) +
  geom_abline(slope=1,intercept=0, alpha=0.3) +
  stat_smooth(method="loess", size=0.7, se=FALSE) +
  # stat_smooth(aes(x=draw_par, y=X2.5.), size=0.4, method="loess", se=FALSE) +
  #  stat_smooth(aes(x=draw_par, y=X97.5.), size=0.4, method="loess", se=FALSE) +
  facet_grid(rows=vars(true_par_mn),cols=vars(param))

#bb_sim_dat %>% filter(sim_id %in% p1_ids,param %in% p1_parms) %>%
#  ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) +
#  geom_abline(slope=1,intercept=0, alpha=0.5) + geom_smooth(se=FALSE) +  
#  facet_grid(rows=vars(true_par_mn),cols=vars(param))

# build plot object for rendering 
gg1 <- ggplot_build(g1)

df2 <- data.frame(x =gg1$data[[2]]$x,
                  ymin = gg1$data[[2]]$y,
                  ymax = gg1$data[[3]]$y) 

g1 +
  geom_ribbon(data=df2, aes(x=x, ymin=ymin, ymax=ymin ),fill="grey",
              inherit.aes = FALSE, alpha = 0.8)


  
  geom_smooth(aes(ymin=X2.5.), se=FALSE) +
  geom_smooth(aes(ymin=X97.5.), se=FALSE)
  
  
  +  
  facet_grid(rows=vars(true_par_mn),cols=vars(param))

#foo<-bb_sim_dat %>% filter(sim_id %in% p1_ids,param %in% p1_parms)

## plots

# efficacy risk difference
bb_sim_dat %>% filter(scn_id==17 & param=="p_e_diff") %>%
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


}

if (0){
  
  # bb_sim_dat %>% filter(scn_id==18 & param=="p_s1" & mod=="jnt") %>%
  #   ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) + 
  #   geom_point() + geom_abline(slope=1,intercept=0, alpha=0.5) +
  #   geom_smooth(se=FALSE)
  
  # bb_sim_dat %>% filter(scn_id==18 & param=="p_e1" & mod=="jnt") %>%
  #   ggplot(aes(x=draw_par, y=mean, col=factor(mod_num))) + 
  #   geom_point() + geom_abline(slope=1, intercept=0, alpha=0.5) +
  #   geom_smooth(se=FALSE)
  
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
