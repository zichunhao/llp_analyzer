#!/bin/sh
#dir=/mnt/hadoop/store/group/phys_susy/razor/Run2Analysis/

if [ "$1" == "" ] || [ "$2" == "" ]; then
    echo "LLPRun_LPC <list of input files> <analyzer name> <options>"
else
   echo "Creating tarball with analyzer"
   tar cvf input.tar ${CMSSW_BASE}/src/llp_analyzer/bin/Run$2
   echo "Copying JEC and event weights"
   cp ${CMSSW_BASE}/src/llp_analyzer/data/JEC.tar.gz .
   cp ${CMSSW_BASE}/src/llp_analyzer/bin/Run$2 .
   #cp ${CMSSW_BASE}/src/llp_analyzer/data/JEC/Autumn18_RunsABCD_V19_DATA.tar.gz .
   #cp ${CMSSW_BASE}/src/llp_analyzer/data/JEC/Autumn18_V19_MC.tar.gz .
   cp ${CMSSW_BASE}/src/llp_analyzer/data/ScaleFactors/2017/Muon* .

   cp ${CMSSW_BASE}/src/llp_analyzer/data/ScaleFactors/METTriggers_SF.root .
   cp ${CMSSW_BASE}/src/llp_analyzer/data/PileupWeights/PileupReweight_MC_Fall17_ZToMuMu_NNPDF31_13TeV-powheg_M_50_120.root .
   cp ${CMSSW_BASE}/src/llp_analyzer/data/PileupWeights/PileupReweight_MC_Fall17_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root .
   cp ${CMSSW_BASE}/src/llp_analyzer/data/PileupWeights/PileupReweight_MC_Fall18_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root .
   cp ${CMSSW_BASE}/src/llp_analyzer/data/PileupWeights/PileupReweight_MC_Summer16_ggH_HToSSTobbbb_MH-125_TuneCUETP8M1_13TeV-powheg-pythia8.root .
   cp ${CMSSW_BASE}/src/llp_analyzer/data/PileupWeights/PileupReweight_source18_target17.root .
   cp ${CMSSW_BASE}/src/llp_analyzer/data/PileupWeights/PileupReweight_source18_target16.root .
   cp ${CMSSW_BASE}/src/llp_analyzer/data/HiggsPtWeights/ggH_HiggsPtReweight_NNLOPS.root .
   tar vxzf JEC.tar.gz
   echo $1 $2
    ./Run$2 $1 ${@:3}
    rm -rf Run$2
    rm -rf PileupReweight_MC_Fall17_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root
    rm -rf PileupReweight_MC_Fall18_ggH_HToSSTobbbb_MH-125_TuneCP5_13TeV-powheg-pythia8.root
    rm -rf PileupReweight_MC_Summer16_ggH_HToSSTobbbb_MH-125_TuneCUETP8M1_13TeV-powheg-pythia8.root
    rm -rf PileupReweight_MC_Fall17_ZToMuMu_NNPDF31_13TeV-powheg_M_50_120.root
    rm -rf ggH_HiggsPtReweight_NNLOPS.root
    rm -rf MuonEfficienciesAndSF_RunBtoF_Nov17Nov2017.root
    rm -rf MuonRunBCDEF_*.2017.root
    rm -rf JEC.tar.gz
    rm -rf JEC
    rm -rf METTriggers_SF.root
    #rm Autumn18_RunsABCD_V19_DATA.tar.gz
    #rm Autumn18_V19_MC.tar.gz
fi
