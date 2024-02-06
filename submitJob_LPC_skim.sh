# NOTE: before you submit (every time) you need to have a valid proxy:
# voms-proxy-init --valid 192:00 -voms cms

# check flag --submit
if [[ $1 == "--submit" ]]; then
    submit=true
else
    submit=false
fi

echo "Before you submit (every time) you need to have a valid proxy"
echo "voms-proxy-init --valid 192:00 -voms cms"
# Check the proxy. If it is expired, renew it.
echo "Type 'y' or 'yes' if you checked the proxy and it is valid."
read -r response
if [[ $response != "y" && $response != "yes" ]]; then
    echo "Please check the proxy and run the script again."
    exit
fi

username=$(whoami)

exe="SkimNtuple"
file_path="/uscms_data/d1/${username}/CMSSW_10_6_20/src/llp_analyzer/lists/displacedJetMuonNtuple/V1p19/Data2023/MuonHitsOnly"
working_dir="/uscms_data/d1/${username}/CMSSW_10_6_20/src/llp_analyzer"
cd ${working_dir}

# set -xe
# run through era in (B, C, D)
for era in B C D; do
    python submitJob_LPC_skim.py --exe ${exe} -i ${file_path}/"Muon0_Run2023${era}-RAW-v1.txt" --njobs 50  -o skim_Run2023${era}_Muon0 --dryRun
    python submitJob_LPC_skim.py --exe ${exe} -i ${file_path}/"Muon1_Run2023${era}-RAW-v1.txt" --njobs 50  -o skim_Run2023${era}_Muon1 --dryRun
    if [[ $submit == true ]]; then
        cd ${working_dir}/skim_Run2023${era}_Muon0
        condor_submit runjob.jdl
        cd ${working_dir}/skim_Run2023${era}_Muon1
        condor_submit runjob.jdl
    fi
done