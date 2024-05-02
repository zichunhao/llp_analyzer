CURRENT_DIR=$(pwd);
JOB_DIR="/uscms_data/d1/$(whoami)/CMSSW_10_6_20/src/llp_analyzer/jobs_2023";
datasets=("0" "1")

set -xe;

# rm -rf "${JOB_DIR}"/merge_hadd_Muon*_Run2023*_v*;

# 2023B
for data_num in "${datasets[@]}"; do
    cd "${JOB_DIR}"/merge_hadd_Muon${data_num}_Run2023B; condor_submit runjob.jdl;
done

# 2023C
for data_num in "${datasets[@]}"; do
    for version_num in "1" "2" "3" "4"; do
        cd "${JOB_DIR}"/merge_hadd_Muon${data_num}_Run2023C_v${version_num}; condor_submit runjob.jdl;
    done
done

# 2023D
for data_num in "${datasets[@]}"; do
    for version_num in "1" "2"; do
        cd "${JOB_DIR}"/merge_hadd_Muon${data_num}_Run2023D_v${version_num}; condor_submit runjob.jdl;
    done
done