## binary-binary model postprocessing for copula calibration - postprocess
## LOCAL version

libs <- c("copula", "magrittr", "rstan", "sfsmisc", "tidyverse")
invisible(lapply(libs, library, character.only = TRUE))

# set working directory
#wdir <- file.path("~/bayes_cop_pow")
wdir<-file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/Biopharm2019/bayes_cop_calib_code")

# load sim array
#bb_sim_array <- readRDS(file.path(wdir,"bb_calib_simarray.rds"))

#! remove and use above for accre version
#bb_calib_simarray <- readRDS(file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/Biopharm2019/bayes_cop_calib_code/bb_calib_simarray_local.rds"))

#! remove for accre
#!s_id <- 4 
#!load(file.path(wdir,"sims_local","bbsims", paste0("bb_sim_",s_id,".RData")))


# Function to extract & format summaries
bb_summ_sim0 <- function(s_id){
  
  eff_summ <-get(paste0("summ_bb_e_",s_id)) %>% data.frame() %>% 
    rownames_to_column(var = "param") %>%
    filter(str_detect(param, "p_e")) %>% mutate(mod="eff")
  
  safety_summ<-get(paste0("summ_bb_s_",s_id)) %>% data.frame() %>% 
    rownames_to_column(var = "param") %>%
    filter(str_detect(param, "p_s")) %>% mutate(mod="saf")
  
  jnt_summ<-get(paste0("summ_bb_jnt_",s_id)) %>% data.frame() %>% 
    rownames_to_column(var = "param") %>%
    filter(str_detect(param, "rho|p_e|p_s")) %>% mutate(mod="jnt")
  
  bind_rows(eff_summ,safety_summ,jnt_summ) %>% mutate(sim_id=s_id)
  
}

bb_summ_sim <- function(s_id){
  
  jnt_summ<-get(paste0("summ_bb_jnt_",s_id)) %>% data.frame() %>% 
    rownames_to_column(var = "param") %>%
    filter(str_detect(param, "rho_|p_e|p_s")) %>% 
    mutate(mod="jnt", sim_id=s_id)
  
  return(jnt_summ)
  
}

# Function to batch read-in and processing of data
# batches are deleted after loading
doBatch <- function(all_sim_ids, batchsize){
  
  sim_ids_split <- split(all_sim_ids, ceiling(seq_along(all_sim_ids)/batchsize))
  batches <- length(sim_ids_split)
  
  for (i in 1:batches){
    
    sims_in_batch <- sim_ids_split[[i]] 
    
    for (j in seq_along(sims_in_batch)){
      # load data
      s_id <- sims_in_batch[j]
      load(file.path(wdir,"sims_local","bbsims", paste0("bb_sim_",s_id,".RData")), .GlobalEnv)
    }
    
    # print(sims_in_batch) # for debugging
    
    # summ_bb_sim for batch - summaries
      assign(paste0("bb_summ_dat0_",i),
             data.frame( bind_rows( lapply(sims_in_batch, bb_summ_sim)) ) ,
             envir=.GlobalEnv)
      
    # combine params for batch
      sim_parms_vec <- paste0("sim_params_",sims_in_batch) 
      
      assign(paste0("bb_params_dat0_",i),
        bind_rows(sapply(sim_parms_vec, get, simplify = FALSE)),
        envir=.GlobalEnv)
      
    # remove loaded datasets, samples, summaries and sim params for batch
    rm(list=ls(envir=.GlobalEnv)[grep("dat_bb_short_|summ_bb_|sim_params_",
                                      ls(envir=.GlobalEnv))], envir=.GlobalEnv)
  }
  
}



# simulation file names
bb_sim_files <- dir(file.path(wdir,"sims_local","bbsims"), pattern=".RData")

#!bb_sim_files_num0 <- substr(bb_sim_files,8,50)
#!bb_sim_files_num <- as.numeric(do.call(rbind,strsplit(bb_sim_files_num0,".R"))[,1])

bb_sim_files_num <- substr(bb_sim_files,8,50) %>% strsplit(".R") %>% 
  do.call(rbind,.) %>% `[`(,1) %>% as.numeric()

#!cb_sim_array_ord_run <- subset(cb_sim_array_ord, sim_id %in% cb_sim_files_num)
#!rm(cb_sim_array, cb_sim_array_ord, cb_sim_files_num0, nreps, nscn, sc_id)

# load and process sims in batches
bsize <- 6
doBatch(bb_sim_files_num, bsize)

bb_summ_dat_vec <- ls()[grep("bb_summ_dat0_",ls())]
bb_summ_dat <- bind_rows(sapply(bb_summ_dat_vec, get,simplify = FALSE)) %>% 
  mutate(param=str_replace_all(param,"\\[|\\]",""))

bb_params_dat_vec <- ls()[grep("bb_params_dat0_",ls())]
bb_params_dat <- bind_rows(sapply(bb_params_dat_vec, get, simplify = FALSE)) %>% 
  select(matches("p_e|p_s|rho|_id|mod_num")) %>%
  mutate(p_e_diff = p_e2 - p_e1, p_s_diff = p_s2 - p_s1,
    p_e_diff_tr = p_e2_tr-p_e1_tr, p_s_diff_tr = p_s2_tr - p_s1_tr) %>% 
  gather(p_e1, p_e2, p_e1_tr, p_e2_tr,
         p_s1, p_s2, p_s1_tr, p_s2_tr,
         rho_1, rho_1_tr, rho_2, rho_2_tr,
         p_e_diff, p_s_diff, p_e_diff_tr, p_s_diff_tr, 
         key="var", value = "val") %>%
  separate(var,sep="_t",into=c("param","suffix"), fill="right") %>%
  spread(suffix,val) %>% rename(draw_par=r, true_par_mn=`<NA>`)

rm(list=c(bb_summ_dat_vec,bb_params_dat_vec))

#merge bb_params_dat and bb_summ_dat by sim_id and param

bb_sim_dat <-left_join(bb_summ_dat,bb_params_dat,by=c("sim_id","param"))



## plot
bb_sim_dat %>% filter(scn_id==2 & param=="p_e1" & mod=="jnt") %>%
  ggplot(aes(x=draw_par, y=mean, col=factor(mod_num))) + 
  geom_point() + geom_abline(slope=1, intercept=0, alpha=0.5) +
  geom_smooth(se=FALSE)

bb_sim_dat %>% filter(scn_id==2 & param=="p_e_diff" & mod=="jnt") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) + 
  geom_point() + geom_abline(slope=1,intercept=0, alpha=0.5) +
  geom_smooth(se=FALSE)

bb_sim_dat %>% filter(scn_id==2 & param=="rho_1" & mod=="jnt") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) + 
  geom_point() + geom_abline(slope=1,intercept=0, alpha=0.5) +
  geom_smooth(se=FALSE)

bb_sim_dat %>% filter(scn_id==2 & param=="p_s1" & mod=="jnt") %>%
  ggplot(aes(x=draw_par,y=mean,col=factor(mod_num))) + 
  geom_point() + geom_abline(slope=1,intercept=0, alpha=0.5) +
  geom_smooth(se=FALSE)

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
