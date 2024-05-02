#!/bin/bash
PROJ_DIR="/uscms_data/d1/$(whoami)/CMSSW_10_6_20/src/llp_analyzer"
EOS_DIR="/eos/uscms/store/user/$(whoami)/llp"

# Define arrays for run years and muon types
runs=("2022E" "2022F" "2022G")

mkdir -p jobs;

# Loop over run years
for run in "${runs[@]}"; do
    # Loop over muon types
    input_dir=/eos/uscms/store/user/amalbert/MDSTriggerEff/skim_Run${run}_moreJobs
    output_dir="${EOS_DIR}/skim_hadd_Run${run}"
    job_dir=${PROJ_DIR}/"jobs_2022/hadd_skimmed_Run${run}"
    mkdir -p "${job_dir}"
    
    # Run Python script with appropriate arguments
    python3 submitJob_LPC_llpMerge_hadd_skimmed.py --input "${input_dir}" --output "${output_dir}" --job "${job_dir}" --n 20
done