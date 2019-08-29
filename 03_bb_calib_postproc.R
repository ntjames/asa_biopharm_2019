## binary-binary model postprocessing for copula calibration - postprocess
rm(list=ls())

libs <- c("copula", "magrittr", "rstan", "sfsmisc", "tidyverse")
invisible(lapply(libs, library, character.only = TRUE))

# set working directory
#wdir <- file.path("~/bayes_cop_pow")
wdir<-file.path("/home/nathan/Dropbox/njames/school/PhD/misc/conferences/Biopharm2019/bayes_cop_calib_code")

if (0){
# summary function for separate marginal mods in addition to joint mod
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

}

##' @title Extract and format posterior summaries
##' @param s_id simulation id
##' @return data frame with posterior parameter summaries 
bb_summ_sim <- function(s_id){
  
  jnt_summ <- get(paste0("summ_bb_jnt_",s_id)) %>% 
    data.frame() %>% 
    rownames_to_column(var = "param") %>%
    filter(str_detect(param, "rho_|p_e|p_s")) %>% 
    mutate(mod="jnt", sim_id=s_id)
  
  return(jnt_summ)
}

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
      load(file.path(wdir,"sims","bbsims", paste0("bb_sim_",s_id,".RData")), .GlobalEnv)
    }
    
    # print(sims_in_batch) # for debugging
    
    # summ_bb_sim for batch - summaries
      assign(paste0("bb_summ_dat0_",i),
             data.frame( bind_rows( lapply(sims_in_batch, bb_summ_sim)) ) ,
             envir=.GlobalEnv)
      
    # combine params for batch
      sim_parms_vec <- paste0("sim_parm_",sims_in_batch) 
      
      assign(paste0("bb_params_dat0_",i),
        bind_rows(sapply(sim_parms_vec, get, simplify = FALSE)),
        envir=.GlobalEnv)
      
    # remove loaded datasets, samples, summaries and sim params for batch
    rm(list=ls(envir=.GlobalEnv)[grep("dat_bb_short_|summ_bb_|sim_parm_",
                                      ls(envir=.GlobalEnv))], envir=.GlobalEnv)
  }
  
}


# simulation file names
bb_sim_files <- dir(file.path(wdir,"sims","bbsims"), pattern=".RData")

#!bb_sim_files_num0 <- substr(bb_sim_files,8,50)
#!bb_sim_files_num <- as.numeric(do.call(rbind,strsplit(bb_sim_files_num0,".R"))[,1])

bb_sim_files_num <- substr(bb_sim_files,8,50) %>% strsplit(".R") %>% 
  do.call(rbind,.) %>% `[`(,1) %>% as.numeric() %>% sort()

# load and process sims in batches
bsize <- 500
doBatch(bb_sim_files_num, bsize)

bb_summ_dat_vec <- ls()[grep("bb_summ_dat0_",ls())]
bb_summ_dat <- bind_rows(sapply(bb_summ_dat_vec, get, simplify = FALSE)) %>% 
  mutate(param=str_replace_all(param,"\\[|\\]",""))

bb_params_dat_vec <- ls()[grep("bb_params_dat0_",ls())]
bb_params_dat0 <-bind_rows(sapply(bb_params_dat_vec, get, simplify = FALSE))
bb_params_dat <- bb_params_dat0 %>% 
  select(matches("p_e|p_s|rho|_id|mod_num")) %>%
  mutate(p_e_diff = p_e2 - p_e1, p_s_diff = p_s2 - p_s1,
    p_e_diff_tr = p_e2_tr-p_e1_tr, p_s_diff_tr = p_s2_tr - p_s1_tr) %>% 
  gather(p_e1, p_e2, p_e1_tr, p_e2_tr,
         p_s1, p_s2, p_s1_tr, p_s2_tr,
         rho_1, rho_1_tr, rho_2, rho_2_tr,
         p_e_diff, p_s_diff, p_e_diff_tr, p_s_diff_tr, 
         key="var", value = "val") %>%
  separate(var, sep="_t", into=c("param","suffix"), fill="right") %>%
  spread(suffix,val) %>% 
  rename(draw_par=r, true_par_mn=`<NA>`)

rm(list=c(bb_summ_dat_vec,bb_params_dat_vec))

# merge bb_params_dat and bb_summ_dat by sim_id and param
bb_sim_dat <-left_join(bb_summ_dat,bb_params_dat,by=c("sim_id","param"))


bb_ndivs_vec <- ls()[grep("n_divs_",ls())]

bb_ndivs_dat<-bind_rows(sapply(bb_ndivs_vec, get, simplify = FALSE)) %>% 
  gather(bb_ndivs_vec, key="id",value="n_div") %>%
  separate(id, sep="_", into=c("a","b","sim_id")) %>%
  mutate(sim_id=as.numeric(sim_id)) %>% select(sim_id, n_div)

bb_parm_div_dat<-left_join(bb_params_dat0,bb_ndivs_dat,by=c("sim_id"))

rm(list=c(bb_ndivs_vec))

# save binary-binary simulation data
save(list=c("bb_sim_dat","bb_parm_div_dat"), 
     file=file.path(wdir,"sims", "bb_simulations.RData"))

rm(bb_ndivs_dat,bb_params_dat0,bb_params_dat, bb_summ_dat)
