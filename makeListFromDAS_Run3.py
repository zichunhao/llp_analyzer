#!/usr/bin/python

import os
import datetime
import time
import subprocess
import glob
import sys
import json


datasets = {

####################################################################################################
#2022 Datasets
####################################################################################################   
#'JetMET_2022C':'/JetMET/Run2022C-PromptNanoAODv10_v1-v1/NANOAOD',
#'JetMET_2022D-v1':'/JetMET/Run2022D-PromptNanoAODv10_v1-v1/NANOAOD',
#'JetMET_2022D-v2':'/JetMET/Run2022D-PromptNanoAODv10_v2-v1/NANOAOD',
#'JetMET_2022E':'/JetMET/Run2022E-PromptNanoAODv10_v1-v3/NANOAOD',
#'JetMET_2022F':'/JetMET/Run2022F-PromptNanoAODv10_v1-v2/NANOAOD',
#'JetMET_2022G':'/JetMET/Run2022G-PromptNanoAODv10_v1-v1/NANOAOD',
#'Muon_2022E':'/Muon/Run2022E-PromptNanoAODv10_v1-v3/NANOAOD',
#'Muon_2022F':'/Muon/Run2022F-PromptNanoAODv10_v1-v2/NANOAOD',
#'Muon_2022G':'/Muon/Run2022G-PromptNanoAODv10_v1-v1/NANOAOD',
#"Muon0_2023B":"/Muon0/Run2023B-PromptNanoAODv11p9_v1-v2/NANOAOD",
"Muon0_2023C_v1":"/Muon0/Run2023C-PromptNanoAODv11p9_v1-v1/NANOAOD",
"Muon0_2023C_v2":"/Muon0/Run2023C-22Sep2023_v2-v1/NANOAOD",  ## no Prompt for -v2
"Muon0_2023C_v3":"/Muon0/Run2023C-PromptNanoAODv12_v3-v1/NANOAOD", 
"Muon0_2023C_v4":"/Muon0/Run2023C-PromptNanoAODv12_v4-v1/NANOAOD", 
"Muon1_2023C_v1":"/Muon1/Run2023C-PromptNanoAODv11p9_v1-v1/NANOAOD",
"Muon1_2023C_v2":"/Muon1/Run2023C-22Sep2023_v2-v1/NANOAOD",  ## no Prompt for -v2
"Muon1_2023C_v3":"/Muon1/Run2023C-PromptNanoAODv12_v3-v1/NANOAOD", 
"Muon1_2023C_v4":"/Muon1/Run2023C-PromptNanoAODv12_v4-v1/NANOAOD", 
"Muon0_2023D_v1":"/Muon0/Run2023D-PromptReco-v1/NANOAOD",
"Muon0_2023D_v2":"/Muon0/Run2023D-PromptReco-v2/NANOAOD",
"Muon1_2023D_v1":"/Muon1/Run2023D-PromptReco-v1/NANOAOD",
"Muon1_2023D_v2":"/Muon1/Run2023D-PromptReco-v2/NANOAOD",
#"Muon1_2023B":"/Muon1/Run2023B-PromptNanoAODv11p9_v1-v2/NANOAOD",
#"Muon1_2023D":"/Muon1/Run2023D-PromptReco-v2/NANOAOD"

}

# if (len(sys.argv) -1 < 1):
#     print "Error. Not enough arguments provided.\n"
#     print "Usage: python printFilesInGivenBlocks.py\n"
#     exit()

#datasetName = sys.argv[1]


for processName in datasets.keys():

    outputFile = open(processName+".list","w")
    print(processName)
    #command = "dasgoclient -query=\"file dataset=" + datasets[processName] + " instance=prod/phys03 \" -json > tmpOutput.json"
    command = "dasgoclient -query=\"file dataset=" + datasets[processName] + " \" -json > tmpOutput.json"
    print(command)
    os.system(command)

    jsonFile = open("tmpOutput.json","r")
    data = json.load(jsonFile)

    for p in data:
        blockName = p['file'][0]['block.name']
        fileName = p['file'][0]['name']
        outputFile.write("root://cmsxrootd.fnal.gov/"+fileName+"\n")
