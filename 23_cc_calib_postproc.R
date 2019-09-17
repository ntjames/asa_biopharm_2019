## continuous-continuous model postprocessing for copula calibration - postprocess
rm(list=ls())

libs <- c("copula", "magrittr", "rstan", "sfsmisc", "tidyverse")
invisible(lapply(libs, library, character.only = TRUE))

# set working directory
wdir<-file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/Biopharm2019/bayes_cop_calib_code")

# load sim array
#cc_sim_array <- readRDS(file.path(wdir,"cc_calib_simarray.rds"))

if (0){
# summary function for separate marginal mods in addition to joint mod
cc_summ_sim0 <- function(s_id){
  
  eff_summ <-get(paste0("summ_cc_e_",s_id)) %>% data.frame() %>%
    rownames_to_column(var = "param") %>%
    filter(str_detect(param, "p_e")) %>% mutate(mod="eff")
  
  safety_summ<-get(paste0("summ_cc_s_",s_id)) %>% data.frame() %>%
    rownames_to_column(var = "param") %>%
    filter(str_detect(param, "p_s")) %>% mutate(mod="saf")
  
  jnt_summ<-get(paste0("summ_cc_jnt_",s_id)) %>% data.frame() %>%
    rownames_to_column(var = "param") %>%
    filter(str_detect(param, "rho|p_e|p_s")) %>% mutate(mod="jnt")
  
  bind_rows(eff_summ,safety_summ,jnt_summ) %>% mutate(sim_id=s_id)
  
}

}

##' @title Extract and format posterior summaries
##' @param s_id simulation id
##' @return data frame with posterior parameter summaries 
cc_summ_sim <- function(s_id){
  
  jnt_summ <- get(paste0("summ_cc_jnt_",s_id)) %>% 
    data.frame() %>% 
    rownames_to_column(var = "param") %>%
    filter(str_detect(param, "theta|mu_e|mu_s")) %>% 
    mutate(mod="jnt", sim_id=s_id)
  
  return(jnt_summ)
}

# 3, 10, 14
#! load(file.path(wdir,"sims","ccsims", paste0("cc_sim_",3,".RData")), .GlobalEnv)
#! cc_summ_sim(3)

##' @title Batch read-in and processing of simulation data, 
##' @description Loaded files in a batch are removed from workspace after extracting summaries
##' @param sim_ids vector containing all simulation ids to be processed
##' @param batchsize number of simulations to load in one batch
##' @return ...
doBatch <- function(sim_ids, batchsize){
  
  sim_ids_split <- split(sim_ids, ceiling(seq_along(sim_ids)/batchsize))
  batches <- length(sim_ids_split)
  
  for (i in 1:batches){
    
    sims_in_batch <- sim_ids_split[[i]] 
    
    for (j in seq_along(sims_in_batch)){
      # load data
      s_id <- sims_in_batch[j]
      # file path to simulation results
      load(file.path(wdir,"sims","ccsims", paste0("cc_sim_",s_id,".RData")), .GlobalEnv)
    }
    
    # print(sims_in_batch) # for debugging
    
    # summ_cc_sim for batch - summaries
      assign(paste0("cc_summ_dat0_",i),
             data.frame( bind_rows( lapply(sims_in_batch, cc_summ_sim)) ) ,
             envir=.GlobalEnv)
      
    # combine params for batch
      sim_parms_vec <- paste0("sim_parm_",sims_in_batch) 
      
      assign(paste0("cc_params_dat0_",i),
        bind_rows(sapply(sim_parms_vec, get, simplify = FALSE)),
        envir=.GlobalEnv)
      
    # remove loaded datasets, samples, summaries and sim params for batch
    rm(list=ls(envir=.GlobalEnv)[grep("dat_cc_|summ_cc_|sim_parm_",
                                      ls(envir=.GlobalEnv))], envir=.GlobalEnv)
  }
  
}


# simulation file names
cc_sim_files <- dir(file.path(wdir,"sims","ccsims"), pattern=".RData")

#!cc_sim_files_num0 <- substr(cc_sim_files,8,50)
#!cc_sim_files_num <- as.numeric(do.call(rbind,strsplit(cc_sim_files_num0,".R"))[,1])

cc_sim_files_num <- substr(cc_sim_files,8,50) %>% strsplit(".R") %>% 
  do.call(rbind,.) %>% `[`(,1) %>% as.numeric() %>% sort()

# load and process sims in batches
bsize <- 600
doBatch(cc_sim_files_num, bsize)

cc_summ_dat_vec <- ls()[grep("cc_summ_dat0_",ls())]
cc_summ_dat <- bind_rows(sapply(cc_summ_dat_vec, get, simplify = FALSE)) %>% 
  mutate(param=str_replace_all(param,"\\[|\\]",""))

cc_params_dat_vec <- ls()[grep("cc_params_dat0_",ls())]
cc_params_dat0 <-bind_rows(sapply(cc_params_dat_vec, get, simplify = FALSE))
cc_params_dat <- cc_params_dat0 %>% 
  select(matches("mu_e|mu_s|theta|_id|mod_num")) %>%
  mutate(mu_e_diff = mu_e2 - mu_e1, mu_s_diff = mu_s2 - mu_s1,
    mu_e_diff_tr = mu_e2_tr-mu_e1_tr, mu_s_diff_tr = mu_s2_tr - mu_s1_tr) %>% 
  gather(mu_e1, mu_e2, mu_e1_tr, mu_e2_tr,
         mu_s1, mu_s2, mu_s1_tr, mu_s2_tr,
         theta_1, theta_1_tr, theta_2, theta_2_tr,
         mu_e_diff, mu_s_diff, mu_e_diff_tr, mu_s_diff_tr, 
         key="var", value = "val") %>%
  separate(var, sep="_t", into=c("param","suffix"), fill="right") %>%
  spread(suffix,val) %>% 
  rename(draw_par=r, true_par_mn=`<NA>`)

# merge cc_params_dat and cc_summ_dat by sim_id and param
cc_sim_dat <-left_join(cc_summ_dat,cc_params_dat,by=c("sim_id","param"))

rm(list=c(cc_summ_dat_vec,cc_params_dat_vec))

cc_ndivs_vec <- ls()[grep("n_divs_",ls())]

cc_ndivs_dat<-bind_rows(sapply(cc_ndivs_vec, get, simplify = FALSE)) %>% 
  gather(cc_ndivs_vec, key="id",value="n_div") %>%
  separate(id, sep="_", into=c("a","b","sim_id")) %>%
  mutate(sim_id=as.numeric(sim_id)) %>% select(sim_id, n_div)

cc_parm_div_dat<-left_join(cc_params_dat0,cc_ndivs_dat,by=c("sim_id"))



rm(list=c(cc_ndivs_vec))

# save continuous-continuous simulation data
save(list=c("cc_sim_dat","cc_parm_div_dat"), 
     file=file.path(wdir,"sims", "cc_simulations.RData"))

rm(cc_ndivs_dat,cc_params_dat0,cc_params_dat, cc_summ_dat)
