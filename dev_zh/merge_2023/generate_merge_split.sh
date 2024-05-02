set -xe

current_dir=$(pwd)

PROJ_DIR="/uscms_data/d1/$(whoami)/CMSSW_10_6_20/src/llp_analyzer"
EOS_DIR="/eos/uscms/store/user/$(whoami)/llp"
MUON_NUMS=("0" "1")

cd ${PROJ_DIR};

# 2023B
for muon_num in "${MUON_NUMS[@]}"; do
    ntupler_dir=${EOS_DIR}/skim_hadd_Muon${muon_num}_Run2023B;
    job_dir=${PROJ_DIR}/jobs_2023/merge_split_Muon${muon_num}_Run2023B;
    output_dir=${EOS_DIR}/merge_split_Muon${muon_num}_Run2023B;
    nanoAOD_path=${PROJ_DIR}/nanoAOD2023/Muon${muon_num}_2023B.list;
    python3 ${PROJ_DIR}/submitJob_LPC_llpMerge_split.py \
    --exe ./MergeNtuples \
    --dir-ntupler ${ntupler_dir} \
    --NanoAOD ${nanoAOD_path} \
    --dir-job ${job_dir} \
    --dir-output ${output_dir} \
    --n 50 --copy-root --memory 2048;
    mkdir -p ${output_dir};
done


# 2023C
for muon_num in "${MUON_NUMS[@]}"; do
    for version_num in "1" "2" "3" "4"; do
        ntupler_dir=${EOS_DIR}/skim_hadd_Muon${muon_num}_Run2023C;
        job_dir=${PROJ_DIR}/jobs_2023/merge_split_Muon${muon_num}_Run2023C_v${version_num};
        output_dir=${EOS_DIR}/merge_split_Muon${muon_num}_Run2023C_v${version_num};
        nanoAOD_path=${PROJ_DIR}/nanoAOD2023/Muon${muon_num}_2023C_v${version_num}.list;
        python3 ${PROJ_DIR}/submitJob_LPC_llpMerge_split.py \
        --exe ./MergeNtuples \
        --dir-ntupler ${ntupler_dir} \
        --NanoAOD ${nanoAOD_path} \
        --dir-job ${job_dir} \
        --dir-output ${output_dir} \
        --n 50 --copy-root --memory 2048;
        mkdir -p ${output_dir};
    done
done

# 2023D
for muon_num in "${MUON_NUMS[@]}"; do
    for version_num in "1" "2"; do
        ntupler_dir=${EOS_DIR}/skim_hadd_Muon${muon_num}_Run2023D;
        job_dir=${PROJ_DIR}/jobs_2023/merge_split_Muon${muon_num}_Run2023D_v${version_num};
        output_dir=${EOS_DIR}/merge_split_Muon${muon_num}_Run2023D_v${version_num};
        nanoAOD_path=${PROJ_DIR}/nanoAOD2023/Muon${muon_num}_2023D_v${version_num}.list;
        python3 ${PROJ_DIR}/submitJob_LPC_llpMerge_split.py \
        --exe ./MergeNtuples \
        --dir-ntupler ${ntupler_dir} \
        --NanoAOD ${nanoAOD_path} \
        --dir-job ${job_dir} \
        --dir-output ${output_dir} \
        --n 50 --copy-root --memory 2048;
        mkdir -p ${output_dir};
    done
done

cd $current_dir;