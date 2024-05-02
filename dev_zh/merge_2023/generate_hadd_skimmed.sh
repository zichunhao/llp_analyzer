#!/bin/bash

PROJ_DIR="/uscms_data/d1/$(whoami)/CMSSW_10_6_20/src/llp_analyzer"
EOS_DIR="/eos/uscms/store/user/$(whoami)/llp"

set -xe

runs=("Muon0_Run2023B" "Muon1_Run2023B" "Muon0_Run2023C" "Muon1_Run2023C" "Muon0_Run2023D" "Muon1_Run2023D")
for run in "${runs[@]}"; do
    # Loop over muon types
    input_dir="${EOS_DIR}/skim_${run}"
    output_dir="${EOS_DIR}/skim_hadd_${run}"
    job_dir="${PROJ_DIR}/jobs_2023/skim_hadd_${run}"

    mkdir -p "${output_dir}"
    mkdir -p "${job_dir}"
    
    # Run Python script with appropriate arguments
    python3 submitJob_LPC_llpMerge_hadd_skimmed.py \
    --input "${input_dir}" \
    --output "${output_dir}" \
    --job "${job_dir}" \
    --n 20
done