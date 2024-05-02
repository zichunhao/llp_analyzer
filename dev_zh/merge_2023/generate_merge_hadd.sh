#!/bin/bash

# Define input and output directories
input_base="/eos/uscms/store/user/$(whoami)/llp"
output_base="/eos/uscms/store/user/$(whoami)/llp/merge_Run2023"
mkdir -p "${output_base}"
SRC_DIR="/uscms_data/d1/$(whoami)/CMSSW_10_6_20/src/llp_analyzer"

# Define arrays for run years and muon types
datasets=("0" "1")

# 2023C
for data_num in "${datasets[@]}"; do
    for version_num in "1" "2" "3" "4"; do
        input_dir="${input_base}/merge_split_Muon${data_num}_Run2023C_v${version_num}"
        output_file="${output_base}/merge_displacedJetMuon_ntupler_Muon${data_num}_Run2023C_v${version_num}.root"
        job_dir=${SRC_DIR}/jobs_2023/"merge_hadd_Muon${data_num}_Run2023C_v${version_num}"
        mkdir -p "${job_dir}"

        python3 ${SRC_DIR}/submitJob_LPC_llpMerge_hadd_merged.py --input "${input_dir}" --output "${output_file}" --job "${job_dir}"
    done
done

# 2023D
for data_num in "${datasets[@]}"; do
    for version_num in "1" "2"; do
        input_dir="${input_base}/merge_split_Muon${data_num}_Run2023D_v${version_num}"
        output_file="${output_base}/merge_displacedJetMuon_ntupler_Muon${data_num}_Run2023D_v${version_num}.root"
        job_dir=${SRC_DIR}/jobs_2023/"merge_hadd_Muon${data_num}_Run2023D_v${version_num}"
        mkdir -p "${job_dir}"

        python3 ${SRC_DIR}/submitJob_LPC_llpMerge_hadd_merged.py --input "${input_dir}" --output "${output_file}" --job "${job_dir}"
    done
done