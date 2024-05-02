# set -xe
current_dir=$(pwd)
PROJ_DIR="/uscms_data/d1/$(whoami)/CMSSW_10_6_20/src/llp_analyzer"
EOS_DIR="/eos/uscms/store/user/$(whoami)/llp"
runs=("2022E" "2022F" "2022G")

cd ${PROJ_DIR};

for run in "${runs[@]}"; do
    job_dir=${PROJ_DIR}/jobs_2022/merge_split_Run${run};
    out_dir=${EOS_DIR}/merge_split_Run${run};
    mkdir -p ${job_dir};
    mkdir -p ${out_dir}

    python3 ${PROJ_DIR}/submitJob_LPC_llpMerge_split.py \
    --exe ./MergeNtuples \
    --dir-ntupler ${EOS_DIR}/skim_hadd_Run${run} \
    --NanoAOD ${PROJ_DIR}/lists/Run3_nanoAOD/Muon_${run}.txt \
    --dir-job ${job_dir} \
    --dir-output ${out_dir} \
    --n 50 --copy-root --memory 2048;
done

cd $current_dir