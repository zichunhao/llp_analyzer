CURRENT_DIR=$(pwd)
SRC_DIR="/uscms_data/d1/$(whoami)/CMSSW_10_6_20/src/llp_analyzer"
EOS_DIR="/eos/uscms/store/user/$(whoami)/llp"
set -xe;
rm -rf "${EOS_DIR}"/hadd_skimmed_2023*;
"${SRC_DIR}"/dev_zh/merge/generate_hadd_skimmed.sh;
cd "${SRC_DIR}"/hadd_skimmed_2023C_Muon0; condor_submit runjob.jdl;
cd "${SRC_DIR}"/hadd_skimmed_2023C_Muon1; condor_submit runjob.jdl;
cd "${SRC_DIR}"/hadd_skimmed_2023D_Muon0; condor_submit runjob.jdl;
cd "${SRC_DIR}"/hadd_skimmed_2023D_Muon1; condor_submit runjob.jdl;
cd "${CURRENT_DIR}";