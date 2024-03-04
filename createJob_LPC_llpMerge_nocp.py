from pathlib import Path

RUN_Vs = [
    # ('B', 0),
    # ('B', 1),
    ("C", 0),
    ("C", 1),
    ("D", 0),
    ("D", 1),
]

DIR_OUTPUT = Path("/uscms_data/d1/zhao1/CMSSW_10_6_20/src/llp_analyzer/")
N_QUEUES = 100

for run_num, muon_num in RUN_Vs:
    sub_dir_output = DIR_OUTPUT / f"merge_Run2023{run_num}_Muon{muon_num}"
    sub_dir_output.mkdir(parents=True, exist_ok=True)
    script_sh = f"""
#!/bin/bash
date
MAINDIR=`pwd`
ls
voms-proxy-info --all
#CMSSW from scratch (only need for root)
export CWD=${{PWD}}
export PATH=${{PATH}}:/cvmfs/cms.cern.ch/common
export SCRAM_ARCH=slc7_amd64_gcc700
scramv1 project CMSSW CMSSW_10_6_20
cp MergeNtuples CMSSW_10_6_20/src
cd CMSSW_10_6_20/src
eval `scramv1 runtime -sh` # cmsenv
echo "Untar JEC:" 
echo "After Untar: "
echo "Inside $MAINDIR:"
ls -lah
echo "Running job..."
# jobs
PATH_NTUPLES="root://cmseos.fnal.gov//store/user/zhao1/llp/skim_Run2023{run_num}_Muon{muon_num}/displacedJetMuon_ntupler.root"
PATH_NANOAODS="${{MAINDIR}}/Muon{muon_num}_2023{run_num}.txt"
./MergeNtuples ${{PATH_NTUPLES}} ${{PATH_NANOAODS}} merge_displacedJetMuon_ntupler_2023{run_num}_Muon{muon_num}.root ""
echo "Inside $MAINDIR:"
ls -lah
echo "coping to eos: "+xrdcp -f merge_displacedJetMuon_ntupler_2023{run_num}_Muon{muon_num}.root root://cmseos.fnal.gov//store/user/zhao1/llp
xrdcp -f merge_displacedJetMuon_ntupler_2023{run_num}_Muon{muon_num}.root root://cmseos.fnal.gov//store/user/zhao1/llp
echo "DELETING..."
rm -rf CMSSW_10_6_20
rm -rf *.pdf *.C core*
cd $MAINDIR  
echo "remove output local file"
rm -rf *.root 
ls
date
""".strip()

    script_jdl = f"""
universe = vanilla
Executable = runjob.sh
Should_Transfer_Files = YES 
WhenToTransferOutput = ON_EXIT_OR_EVICT
Transfer_Input_Files = runjob.sh,/uscms_data/d3/zhao1/CMSSW_10_6_20/src/llp_analyzer/MergeNtuples,/uscms_data/d3/zhao1/CMSSW_10_6_20/src/llp_analyzer/convertListMerge.py,/uscms/home/zhao1/nobackup/CMSSW_10_6_20/src/llp_analyzer/lists/Run3_nanoAOD/Muon{muon_num}_2023{run_num}.txt
Output = runjob.$(Process).$(Cluster).stdout
Error  = runjob.$(Process).$(Cluster).stdout
Log    = runjob.$(Process).$(Cluster).log
Arguments = $(Process) {N_QUEUES}
Queue {N_QUEUES}
""".strip()

    (sub_dir_output / "runjob.sh").write_text(script_sh)
    (sub_dir_output / "runjob.jdl").write_text(script_jdl)

    print(f"Scripts for Run2023{run_num}_Muon{muon_num} saved to {sub_dir_output}")
