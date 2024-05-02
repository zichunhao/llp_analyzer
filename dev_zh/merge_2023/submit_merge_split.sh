CURRENT_DIR=$(pwd)
JOBS_DIR="/uscms_data/d1/$(whoami)/CMSSW_10_6_20/src/llp_analyzer/jobs_2023"
EOS_DIR="/eos/uscms/store/user/$(whoami)/llp"
set -xe

# rm -rf "${EOS_DIR}"/merge_split_Muon*_Run2023*_v*;

cd "${JOBS_DIR}"/merge_split_Muon0_Run2023B; condor_submit runjob.jdl;
cd "${JOBS_DIR}"/merge_split_Muon1_Run2023B; condor_submit runjob.jdl;

for version_num in "1" "2" "3" "4"; do
    cd "${JOBS_DIR}"/merge_split_Muon0_Run2023C_v${version_num}; condor_submit runjob.jdl;
    cd "${JOBS_DIR}"/merge_split_Muon1_Run2023C_v${version_num}; condor_submit runjob.jdl;
done

for version_num in "1" "2"; do
    cd "${JOBS_DIR}"/merge_split_Muon0_Run2023D_v${version_num}; condor_submit runjob.jdl;
    cd "${JOBS_DIR}"/merge_split_Muon1_Run2023D_v${version_num}; condor_submit runjob.jdl;
done

cd $CURRENT_DIR;
