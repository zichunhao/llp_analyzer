curr_dir=$(pwd)
runs=("Run2022E" "Run2022F" "Run2022G")

for run in "${runs[@]}"; do
    cd /uscms_data/d1/$(whoami)/CMSSW_10_6_20/src/llp_analyzer_alex/jobs/hadd_skimmed_1_${run}
    condor_submit runjob.jdl;
done

cd $curr_dir

