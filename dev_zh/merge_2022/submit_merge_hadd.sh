curr_dir=$(pwd)
PROJ_DIR="/uscms_data/d1/$(whoami)/CMSSW_10_6_20/src/llp_analyze"
JOB_DIR="${PROJ_DIR}/jobs_2022"
runs=("2022E" "2022F" "2022G")

set -xe
for run in "${runs[@]}"; do
    cd ${JOB_DIR}/merge_hadd_Run${run}
    condor_submit runjob.jdl;
done

cd $curr_dir

