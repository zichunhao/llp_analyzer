#include "TTree.h"
#include <vector>

void matchMuonEventsWithNanoAODs(*TTree muonEvents, *TTree nanoAODs, vector<*TTree> *matchedMuonEventIdx, *Tree mergedTree)
{
    int numMuonEvents = muonEvents->GetEntries();
    int numNanoAODs = nanoAODs->GetEntries();
    
    for (int i = 0; i < numNanoAODs; i++)
    {
        nanoAODs->GetEntry(i);
        for (int j = 0; j < numMuonEvents; j++)
        {
            if (matchedMuonEventIdx->contains(j))
            {
                // already matched -> skip
                continue;
            }

            muonEvents->GetEntry(j);
            
            // Match muon event number and run number with nanoAODs
            if ((muonEvents->event == nanoAODs->event) && (muonEvents->run == nanoAODs->run))
            {
                matchedMuonEventIdx->push_back(j);
                // TODO: Merge muonEvents and nanoAODs

            }

        }
    }
}