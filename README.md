# asa_biopharm_2019
Code for "Bayesian Joint Benefit-Risk Copula Model Calibration" poster presented at ASA Biopharm 2019

For each simulation scenario (binary-binary, continuous-binary, continuous-continuous) there are 5 files:

*0_xx_calib_simarray.R - make simulation array with all parameter & model combinations  

*1_xxmod_calib_init.R - compile Stan models (this also calls in the appropriate *.stan files)  

*2_xx_calib_run.R - create data based on simulation parameters, call & run compiled Stan model, and save output; designed to be called by an embarrasingly parallel array job  

*3_xx_calib_postproc.R - process all the simulation output files into one dataset  

*4_xx_calib_plot.R - create plots summarizing the simulation results  

