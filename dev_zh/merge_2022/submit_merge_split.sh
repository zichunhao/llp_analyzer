curr_dir=$(pwd)
runs=("2022E" "2022F" "2022G")
PROJ_DIR="/uscms_data/d1/$(whoami)/CMSSW_10_6_20/src/llp_analyzer"

for run in "${runs[@]}"; do
    cd ${PROJ_DIR}/jobs_2022/merge_split_Run${run}
    condor_submit runjob.jdl;
done

cd $curr_dir

