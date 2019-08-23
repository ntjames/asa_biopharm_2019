#!/bin/bash

# loop over sim_ids in bb_calib_simarray
for i in {1..20}; do
    Rscript --no-save 1_bb_calib_run_local.R $i
done
