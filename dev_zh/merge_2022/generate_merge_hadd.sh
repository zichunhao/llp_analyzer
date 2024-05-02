#!/bin/bash

# Define input and output directories
EOS_DIR="/eos/uscms/store/user/$(whoami)/llp"
OUT_DIR="${EOS_DIR}/merge_Run2022"
SRC_DIR="/uscms_data/d1/$(whoami)/CMSSW_10_6_20/src/llp_analyzer"

mkdir -p ${OUT_DIR}

# Define arrays for run years and muon types
runs=("2022E" "2022F" "2022G")

mkdir -p ${output_base}

for run in "${runs[@]}"; do
    input_dir="${EOS_DIR}/merge_split_Run${run}"
    output_file="${OUT_DIR}/merge_displacedJetMuon_ntupler_Run${run}.root"
    job_dir="${SRC_DIR}/jobs_2022/merge_hadd_Run${run}"
    mkdir -p ${job_dir}

    python3 ${SRC_DIR}/submitJob_LPC_llpMerge_hadd_merged.py --input "${input_dir}" --output "${output_file}" --job "${job_dir}"
done