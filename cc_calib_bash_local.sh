#!/bin/bash

# loop over sim_ids in cb_calib_simarray
for i in {1..18}; do
    Rscript --no-save 22_cc_calib_run_local.R $i
done
