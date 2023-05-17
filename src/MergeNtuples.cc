#include <fstream>
#include <sstream>
#include <iterator>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TKey.h"
#include <assert.h>
#include <TRandom3.h>
#include "TTreeFormula.h"
#include <iostream>
#include "RazorAnalyzer.h"
#include "HNLMuonSystemTree.h"
#include "RazorHelper.h"
using namespace std;


std::string ParseCommandLine( int argc, char* argv[], std::string opt )
{
  for (int i = 1; i < argc; i++ )
    {
      std::string tmp( argv[i] );
      if ( tmp.find( opt ) != std::string::npos )
        {
          if ( tmp.find( "=" )  != std::string::npos ) return tmp.substr( tmp.find_last_of("=") + 1 );
	  if ( tmp.find( "--" ) != std::string::npos ) return "yes";
	}
    }
  
  return "";
};

struct leptons
{
  TLorentzVector lepton;
  int pdgId;
  float dZ;
  // bool passLooseId;
  // bool passMediumId;
  bool passId;
  bool passVetoId;
  bool passLooseIso;
  bool passTightIso;
  bool passVTightIso;
  bool passVVTightIso;

};


struct jets
{
  TLorentzVector jet;
  float time;
  bool passId;
  // bool passLooseId;
  // bool passMediumId;
  // bool passTightId;
  bool isCSVL;
  int ecalNRechits;
  float jetCISV;
  float jetCMVA;
  float ecalRechitE;
  float jetChargedEMEnergyFraction;
  float jetNeutralEMEnergyFraction;
  float jetChargedHadronEnergyFraction;
  float jetNeutralHadronEnergyFraction;
  bool jetPassMuFrac;
  float jetPtJESUp;
  float jetPtJESDown;
  float jetEJESUp;
  float jetEJESDown;
  float JecUnc;

  float electronEnergyFraction;
  float neutralEmEnergyFraction;
  float chargedHadronEnergyFraction;
  float neutralHadronEnergyFraction;
  float muonEnergyFraction;

};

//lepton highest pt comparator
struct largest_pt
{
  inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt;

//jet highest pt comparator
struct largest_pt_jet
{
  inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_jet;



//get list of files to open, add normalization branch to the tree in each file
int main(int argc, char* argv[]) {

    //parse input list to get names of ROOT files
    if(argc < 4){
        cerr << "usage MergeNtuple [inputfile1] [inputfile2] [outputfile] [isData]" << endl;
        return -1;
    }
    string inputfilename1(argv[1]);
    string inputfilename2(argv[2]);
    string outputfilename(argv[3]);
    string analysisTag(argv[4]); //analysisTag

    const float ELE_MASS = 0.000511;
    const float MU_MASS  = 0.105658;
    const float Z_MASS   = 91.2;
  
    if (analysisTag == ""){
      analysisTag = "Razor2016_80X";
  
    }

    //create output file
    TFile *outputFile = new TFile(outputfilename.c_str(), "RECREATE");

    //open files
    //TFile *inputFile1 = TFile::Open(inputfilename1.c_str(), "READ");
    HNLMuonSystemTree *MuonSystem = new HNLMuonSystemTree;
    MuonSystem->LoadTree(inputfilename1.c_str());
    if (!MuonSystem) cout << "Input Tree not found in file " << inputfilename1 << "\n";
    assert(MuonSystem);


    //TFile *inputFile2 = TFile::Open(inputfilename2.c_str(), "READ");
    //TTree *inputTree2 = (TTree*)inputFile2->Get("tree");
    //if (!inputTree2) cout << "Input Tree not found in file " << inputfilename2 << "\n";
    //assert(inputTree2);

    //build the TChain
    //tree name is set give the structure in the first root file, see while loop below
    TChain* theChain = new TChain();
    string curFileName;

    ifstream inputFile(inputfilename2.c_str());
    int NFilesLoaded = 0;
    if ( !inputFile ){
      cerr << "Error: input file not found!" << endl;
      return -1;
    }

    while ( getline(inputFile, curFileName) )      {
	if ( NFilesLoaded == 0 )	  {
	    /*
	      checks root file structure and add first file
	    */
		std::cout << "[INFO]: loading file: " << curFileName.c_str() <<std::endl;
	    TFile* f_0 = TFile::Open( curFileName.c_str() );
	    f_0->ls();	
		theChain->SetName("ntuples/llp");
		std::cout << "[INFO]: setting name of tchain = ntuples/llp"<< std::endl;
	    theChain->Add( curFileName.c_str() );
        std::cout << "Loaded  " << theChain->GetEntries() << " events\n";
	    delete f_0;
	  }	else	  {
	    //Addind remaining files after file structure is decided
	    theChain->Add( curFileName.c_str() );
	  }
        NFilesLoaded++;
      }     
    std::cout << "Loaded Total of " << NFilesLoaded << " files\n";
    std::cout << "Loaded Total of " << theChain->GetEntries() << " events\n";
    if ( theChain == NULL ) return -1;
 

    uint NEventsTree1 = MuonSystem->tree_->GetEntries(); 
    uint NEventsTree2 = theChain->GetEntries();
    //NEventsTree1 = 10000;

    //*****************************************************************************************
    //Make map of event number in tree 1 to event number in tree 2
    //*****************************************************************************************
    int numberOfNonMatchedEvents = 0;
    std::map<uint,uint> EventIndexToEventIndexMap;

    //Read AOD input with analyzer
    RazorAnalyzer ra  =  RazorAnalyzer(theChain);
    ra.EnableAll();

    //loop over tree2
    std::vector<std::pair<uint,uint> > eventList2;
    std::vector<bool> matchedevent;
    cout << "building tree 2 index map \n";
    for (uint m=0; m < NEventsTree2;m++) { 
      if ( m % 10000 == 0) cout << "Event " << m << "\n";
      ra.fChain->GetEntry(m);    
      std::pair<uint,uint> p (ra.runNum, ra.eventNum);
      //cout << " runNum = " <<p.first<<" LS ="<<ra.lumiNum<< " evtNum = "<< p.second << "  \n";
      eventList2.push_back(p);
    }
    
    cout << "Looping over tree 1:  \n";
    cout << "Total Input Entries: " << NEventsTree1 << "\n";
    //loop over tree1
    for (uint n=0;n<NEventsTree1;n++) { 
      if ( n % 10000 == 0) cout << "Event " << n << "\n";
      MuonSystem->tree_->GetEntry(n);    
      //cout << " runNum = " <<MuonSystem->runNum<<" LS = "<< MuonSystem->lumiSec<<" evtNum = "<< MuonSystem->evtNum << "  \n";

      //loop over tree2
      bool matchFound = false;
      for (uint m=0; m < eventList2.size();m++) 
      { 
          	if (eventList2[m].first != MuonSystem->runNum) continue;
          	if (eventList2[m].second != MuonSystem->evtNum) continue;
          	//if matched, then add entry to the map
          	EventIndexToEventIndexMap[n] = m;
          	//cout << "Match : " << n << " --> " << m << "\n";
          	matchFound = true;
          	break;
      }            

      if (!matchFound) 
      { 
          numberOfNonMatchedEvents++;
          matchedevent.push_back(false);
        //	cout << "Event " << MuonSystem->runNum << ":" << MuonSystem->evtNum <<" Not Matched\n";
      }
      else{
          matchedevent.push_back(true);
      }
      
    }
 
    cout << "Matched events    = "<< NEventsTree1-numberOfNonMatchedEvents << " / " << NEventsTree1 <<" \n";
    cout << "Un-Matched events = "<< numberOfNonMatchedEvents << " / " << NEventsTree1 <<" \n";

    //*****************************************************************************************
    //Produce Output Tree
    //*****************************************************************************************       
    cout << "Merging tree " << MuonSystem->tree_->GetName() << endl;

    cout << "Events in the ntuple: " << MuonSystem->tree_->GetEntries() << endl;
    //create output tree
    outputFile->cd();
    
    // clone tree with 0 entries, copy all the branch addresses only
    TTree *outputTree = MuonSystem->tree_->CloneTree(0);

    RazorHelper *helper = 0;
    bool isData = true;
    helper = new RazorHelper(analysisTag, isData, false);
  
    //loop over Tree1 and add all the branches from tree2
    for (uint n=0;n<NEventsTree1;n++) { 
      if (n%10000==0) cout << "Processed Event " << n << "\n";

      //Check if found a match
      if(matchedevent[n] == false) continue;

      //Get entries
      MuonSystem->tree_->GetEntry(n);
      //theChain->GetEntry(EventIndexToEventIndexMap[n]);
      ra.fChain->GetEntry(EventIndexToEventIndexMap[n]);

      //Copy branches from big nTuple 
      MuonSystem->nRpc = ra.nRpc;
      MuonSystem->rho = ra.fixedGridRhoFastjetAll;
      MuonSystem->met = ra.metType1Pt;
      MuonSystem->metPhi = ra.metType1Phi;
      MuonSystem->metJESUp = MuonSystem->met;
      MuonSystem->metJESDown = MuonSystem->met;

      //cout << "Extract MET" << endl;
      MuonSystem->metSF = helper->getMetTriggerSF(MuonSystem->met);

      std::pair<double,double> corrected_met;
      if (analysisTag=="Razor2016_07Aug2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(ra.metType1Pt, ra.metType1Phi, ra.runNum, 2016, !isData, ra.nPV);
      else if (analysisTag=="Razor2017_17Nov2017Rereco") corrected_met = helper->METXYCorr_Met_MetPhi(ra.metType1Pt, ra.metType1Phi, ra.runNum, 2017, !isData, ra.nPV);
      else if (analysisTag=="Razor2018_17SeptEarlyReReco") corrected_met = helper->METXYCorr_Met_MetPhi(ra.metType1Pt, ra.metType1Phi, ra.runNum, 2018, !isData, ra.nPV);
  
      MuonSystem->metXYCorr = corrected_met.first;
      MuonSystem->metPhiXYCorr = corrected_met.second;
  
      MuonSystem->metEENoise = MuonSystem->metXYCorr;
      MuonSystem->metPhiEENoise = MuonSystem->metPhiXYCorr;

      for(int i = 0; i < NTriggersMAX; i++){
        MuonSystem->HLTDecision[i] = ra.HLTDecision[i];
      }
    // flags
    MuonSystem->Flag_HBHENoiseFilter = ra.Flag_HBHENoiseFilter;
    MuonSystem->Flag_HBHEIsoNoiseFilter = ra.Flag_HBHEIsoNoiseFilter;
    MuonSystem->Flag_BadPFMuonFilter = ra.Flag_BadPFMuonFilter;
    MuonSystem->Flag_globalSuperTightHalo2016Filter = ra.Flag_globalSuperTightHalo2016Filter;
    MuonSystem->Flag_goodVertices = ra.Flag_goodVertices;
    MuonSystem->Flag_ecalBadCalibFilter = ra.Flag_ecalBadCalibFilter;
    MuonSystem->Flag_BadChargedCandidateFilter = ra.Flag_BadChargedCandidateFilter;
    MuonSystem->Flag_eeBadScFilter = ra.Flag_eeBadScFilter;

    MuonSystem->Flag2_HBHENoiseFilter = ra.Flag2_HBHENoiseFilter;
    MuonSystem->Flag2_HBHEIsoNoiseFilter = ra.Flag2_HBHEIsoNoiseFilter;
    MuonSystem->Flag2_BadPFMuonFilter = ra.Flag2_BadPFMuonFilter;
    MuonSystem->Flag2_globalSuperTightHalo2016Filter = ra.Flag2_globalSuperTightHalo2016Filter;
    MuonSystem->Flag2_globalTightHalo2016Filter = ra.Flag2_globalTightHalo2016Filter;
    MuonSystem->Flag2_BadChargedCandidateFilter =ra.Flag2_BadChargedCandidateFilter;
    // Flag2_goodVertices = Flag2_goodVertices;
    MuonSystem->Flag2_EcalDeadCellTriggerPrimitiveFilter = ra.Flag2_EcalDeadCellTriggerPrimitiveFilter;
    MuonSystem->Flag2_ecalBadCalibFilter = ra.Flag2_ecalBadCalibFilter;
    MuonSystem->Flag2_eeBadScFilter = ra.Flag2_eeBadScFilter;
    MuonSystem->Flag2_all = ra.Flag2_HBHENoiseFilter && ra.Flag2_HBHEIsoNoiseFilter && ra.Flag2_BadPFMuonFilter && ra.Flag2_globalSuperTightHalo2016Filter && ra.Flag2_EcalDeadCellTriggerPrimitiveFilter;
    if (isData) MuonSystem->Flag2_all = MuonSystem->Flag2_all && ra.Flag2_eeBadScFilter;

    if (analysisTag!="Razor2016_07Aug2017Rereco")
      {
        MuonSystem->Flag2_all = MuonSystem->Flag2_all && ra.Flag2_ecalBadCalibFilter;
      }
    //*************************************************************************
    //Start Object Selection
    //*************************************************************************

    std::vector<leptons> Leptons;
    //-------------------------------
    //Muons
    //-------------------------------
    for( int i = 0; i < ra.nMuons; i++ ) {

      if(!ra.isMuonPOGTightMuon(i)) continue;
      if(ra.muonPt[i] < 25) continue;
      if(fabs(ra.muonEta[i]) > 2.4) continue;

      //remove overlaps
      bool overlap = false;
      for(auto& lep : Leptons)
        {
          if (ra.deltaR(ra.muonEta[i],ra.muonPhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
      if(overlap) continue;

      leptons tmpMuon;
      tmpMuon.lepton.SetPtEtaPhiM(ra.muonPt[i],ra.muonEta[i], ra.muonPhi[i], MU_MASS);
      tmpMuon.pdgId = 13 * -1 * ra.muonCharge[i];
      tmpMuon.dZ = ra.muon_dZ[i];
      tmpMuon.passId = ra.isMuonPOGTightMuon(i);
      float muonIso = (ra.muon_chargedIso[i] + fmax(0.0,  ra.muon_photonIso[i] + ra.muon_neutralHadIso[i] - 0.5*ra.muon_pileupIso[i])) / ra.muonPt[i];

      tmpMuon.passLooseIso = muonIso<0.25;
      tmpMuon.passTightIso = muonIso<0.15;
      tmpMuon.passVTightIso = muonIso<0.10;
      tmpMuon.passVVTightIso = muonIso<0.05;

      tmpMuon.passVetoId = false;
      Leptons.push_back(tmpMuon);
    }

    //-------------------------------
    //Electrons
    //-------------------------------
    for( int i = 0; i < ra.nElectrons; i++ ) {

      if (!ra.isEGammaPOGTightElectron(i, true, false, true, "vid")) continue;

      //if(elePt[i] < 35) continue;
      if(ra.elePt[i] < 30) continue;
      if(fabs(ra.eleEta[i]) > 2.5) continue;

      //remove overlaps
      bool overlap = false;
      for(auto& lep : Leptons)
        {
          if (ra.deltaR(ra.eleEta[i],ra.elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
      if(overlap) continue;
      leptons tmpElectron;
      tmpElectron.lepton.SetPtEtaPhiM(ra.elePt[i],ra.eleEta[i], ra.elePhi[i], ELE_MASS);
      tmpElectron.pdgId = 11 * -1 * ra.eleCharge[i];
      tmpElectron.dZ = ra.ele_dZ[i];
      tmpElectron.passId = ra.isEGammaPOGTightElectron(i, true, true, true, "Summer16");
      tmpElectron.passVetoId = ra.isEGammaPOGVetoElectron(i, true, true, true, "Summer16");
      Leptons.push_back(tmpElectron);
    }

    sort(Leptons.begin(), Leptons.end(), my_largest_pt);

    MuonSystem->nLeptons = 0;
    for (const auto& tmp : Leptons ) {
      MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
      MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
      MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
      MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
      MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
      MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
      MuonSystem->lepPassId[MuonSystem->nLeptons] = tmp.passId;
      MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
      MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
      MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
      MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
      MuonSystem->nLeptons++;
    }
    for(int i = 0; i < ra.nMuons; i++) {
      if (fabs(ra.muonEta[i]>3.0)) continue;
      float muonIso = (ra.muon_chargedIso[i] + fmax(0.0,  ra.muon_photonIso[i] + ra.muon_neutralHadIso[i] - 0.5*ra.muon_pileupIso[i])) / ra.muonPt[i];
        MuonSystem->muonPt[i]   = ra.muonPt[i];
        MuonSystem->muonE[i]   = ra.muonE[i];
        MuonSystem->muonEta[i]  = ra.muonEta[i];
        MuonSystem->muonPhi[i]  = ra.muonPhi[i];
        MuonSystem->muonIso[i]  = muonIso;
        MuonSystem->muonIsGlobal[i]  = ra.muon_isGlobal[i];
        MuonSystem->muonTightId[i]  = ra.isMuonPOGTightMuon(i);
        MuonSystem->muonLooseId[i]  = ra.isMuonPOGLooseMuon(i);
    }

    MuonSystem->nMuons = ra.nMuons;

    std::vector<jets> Jets;
    float MetXCorr_JESUp = 0.;
    float MetYCorr_JESUp = 0.;
    float MetXCorr_JESDown = 0.;
    float MetYCorr_JESDown = 0.;
    float MetXCorr_HEM = 0.;
    float MetYCorr_HEM = 0.;
    float MetXCorr_EENoise = 0.;
    float MetYCorr_EENoise = 0.;

    for(int i = 0; i < ra.nJets; i++) {

      //------------------------------------------------------------
      //exclude selected muons and electrons from the jet collection
      //------------------------------------------------------------
      double deltaR = -1;
      for(auto& lep : Leptons){
        double thisDR = ra.deltaR(ra.jetEta[i],ra.jetPhi[i],lep.lepton.Eta(),lep.lepton.Phi());
        if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
      }
      if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

      //------------------------------------------------------------
      //Apply Jet Energy and Resolution Corrections
      //------------------------------------------------------------
      // double JEC = JetEnergyCorrectionFactor(ra.jetPt[i], ra.jetEta[i], ra.jetPhi[i], jetE[i],
      // fixedGridRhoFastjetAll, jetJetArea[i] , JetCorrector);
      // cout<<"before JEC"<<endl;
      // double JEC = JetEnergyCorrectionFactor(ra.jetPt[i], ra.jetEta[i], ra.jetPhi[i], jetE[i],
      //  				   fixedGridRhoFastjetAll, jetJetArea[i],
      //  				   runNum,
      //  				   JetCorrectorIOV,JetCorrector);
      double JEC = 1.0;
      double jetCorrPt = ra.jetPt[i]*JEC;
      double jetCorrE = ra.jetE[i]*JEC;
      TLorentzVector thisJet = ra.makeTLorentzVector( jetCorrPt, ra.jetEta[i], ra.jetPhi[i], jetCorrE );

      if (thisJet.Eta()>-3.0 && thisJet.Eta()<-1.3 && thisJet.Phi() >-1.57 && thisJet.Phi() <-0.87 && analysisTag == "Razor2018_17SeptEarlyReReco")
        {
          MetXCorr_HEM += thisJet.Px();
          MetYCorr_HEM += thisJet.Py();
        }
      if (fabs(thisJet.Eta())> 2.65 && fabs(thisJet.Eta())<3.139 && thisJet.Pt() < 50  && analysisTag == "Razor2017_17Nov2017Rereco")
        {
          MetXCorr_EENoise += thisJet.Px();
          MetYCorr_EENoise += thisJet.Py();
          if (ra.eventNum==123969624)
            {
              cout<<i<<", "<<thisJet.Eta()<<","<<ra.jetEta[i]<<endl;
              cout<<MetXCorr_EENoise<<", "<<MetYCorr_EENoise<<", "<<thisJet.Px()<<", "<<thisJet.Py()<<","<<thisJet.Pt()<<","<<ra.jetPt[i]<<endl;
            }
        }
      if (fabs(thisJet.Eta())> 2.25 && fabs(thisJet.Eta())<3.0 && thisJet.Pt() > 100 && (analysisTag == "Razor2016_07Aug2017Rereco" || analysisTag == "Razor2017_17Nov2017Rereco"))
        {
          MuonSystem->EE_prefiring = false;
        }


      if (fabs(thisJet.Eta()) >= 3.0)continue;

      jets tmpJet;
      tmpJet.jet    = thisJet;
      tmpJet.passId = ra.jetPassIDTight[i];
      tmpJet.jetCISV = ra.jetCISV[i];
      tmpJet.jetCMVA = ra.jetCMVA[i];

      // calculate jet energy scale uncertainty
      double unc = helper->getJecUnc( jetCorrPt, ra.jetEta[i], ra.runNum ); //use run=999 as default
      tmpJet.jetPtJESUp = jetCorrPt*(1+unc);
      tmpJet.jetPtJESDown = jetCorrPt*(1-unc);
      tmpJet.jetEJESUp = jetCorrE*(1+unc);
      tmpJet.jetEJESDown = jetCorrE*(1-unc);
      tmpJet.JecUnc = unc;
      TLorentzVector thisJetJESUp = ra.makeTLorentzVector(tmpJet.jetPtJESUp, ra.jetEta[i], ra.jetPhi[i], tmpJet.jetEJESUp);
      TLorentzVector thisJetJESDown = ra.makeTLorentzVector(tmpJet.jetPtJESDown, ra.jetEta[i], ra.jetPhi[i], tmpJet.jetEJESDown);
      if (tmpJet.jetPtJESUp > 10) {
        MetXCorr_JESUp += -1 * (thisJetJESUp.Px() - thisJet.Px());
        MetYCorr_JESUp += -1 * (thisJetJESUp.Py() - thisJet.Py());
      }
      if (tmpJet.jetPtJESDown > 10) {
        MetXCorr_JESDown += -1 * (thisJetJESDown.Px() - thisJet.Px());
        MetYCorr_JESDown += -1 * (thisJetJESDown.Py() - thisJet.Py());
      }
      // if ( !jetPassIDLoose[i] ) continue;
      if(!ra.isPFTightJet(i, true,analysisTag))continue;
      if( thisJet.Pt() < 20 ) continue;//According to the April 1st 2015 AN
      Jets.push_back(tmpJet);
    }
    MuonSystem->nJets = 0;

    // if( Jets.size() < 2 ) continue;
    // if (triggered) trig_lepId_dijet->Fill(1);
    sort(Jets.begin(), Jets.end(), my_largest_pt_jet);
    if (Jets.size()>0) {
      MuonSystem->jetMet_dPhi = ra.deltaPhi(ra.jetPhi[0],ra.metType1Phi);
    }
    else{
      MuonSystem->jetMet_dPhi = -999.;
    }
    double jetMet_dPhiMin_temp = 999.;
    double jetMet_dPhiMin4_temp = 999.;
    int nbjet =0;
    for ( auto &tmp : Jets ) {
      // if(tmp.jet.Pt()<50)continue;
      if(abs(tmp.jet.Eta())>2.4)continue;
      MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
      MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
      MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
      MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
      MuonSystem->jetTime[MuonSystem->nJets] = tmp.time;
      MuonSystem->jetPtJESUp[MuonSystem->nJets] = tmp.jetPtJESUp;
      MuonSystem->jetPtJESDown[MuonSystem->nJets] = tmp.jetPtJESDown;
      MuonSystem->jetEJESUp[MuonSystem->nJets] = tmp.jetEJESUp;
      MuonSystem->jetEJESDown[MuonSystem->nJets] = tmp.jetEJESDown;
      MuonSystem->JecUnc[MuonSystem->nJets] = tmp.JecUnc;
      MuonSystem->jetCISV[MuonSystem->nJets] = tmp.jetCISV;
      MuonSystem->jetCMVA[MuonSystem->nJets] = tmp.jetCMVA;
      if (jetMet_dPhiMin4_temp > abs(ra.deltaPhi(tmp.jet.Phi(),ra.metType1Phi)) && MuonSystem->nJets < 4) {
        jetMet_dPhiMin4_temp = abs(ra.deltaPhi(tmp.jet.Phi(),ra.metType1Phi));

      }
      if (jetMet_dPhiMin_temp > abs(ra.deltaPhi(tmp.jet.Phi(),ra.metType1Phi))) {
        if (tmp.jet.Pt()>30 && abs(tmp.jet.Eta())<2.4) {
          jetMet_dPhiMin_temp = abs(ra.deltaPhi(tmp.jet.Phi(),ra.metType1Phi));
        }
      }
      MuonSystem->jetTightPassId[MuonSystem->nJets] = tmp.passId;
      MuonSystem->nJets++;
    }

    //cout<< "Select MET " << endl;

    MuonSystem-> jetMet_dPhiMin = jetMet_dPhiMin_temp;
    MuonSystem-> jetMet_dPhiMin4 = jetMet_dPhiMin4_temp;
    // TLorentzVector PFMET = makeTLorentzVectorPtEtaPhiM(metType1Pt, 0, metType1Phi, 0);
    TLorentzVector PFMET = ra.makeTLorentzVectorPtEtaPhiM(MuonSystem->metXYCorr, 0, MuonSystem->metPhiXYCorr, 0);

    //JES up
    float PFMetXJESUp   = PFMET.Px() + MetXCorr_JESUp;
    float PFMetYJESUp   = PFMET.Py() + MetYCorr_JESUp;
    MuonSystem->metJESUp    = sqrt( pow(PFMetXJESUp,2) + pow(PFMetYJESUp,2) );
    MuonSystem->metPhiJESUp    = atan(PFMetYJESUp/PFMetXJESUp);
    if  (PFMetXJESUp < 0.0) MuonSystem->metPhiJESUp = ra.deltaPhi(TMath::Pi() + MuonSystem->metPhiJESUp,0.0);
    MuonSystem->metJESUpSF = helper->getMetTriggerSF(MuonSystem->metJESUp)/MuonSystem->metSF;

    //JES down
    float PFMetXJESDown   = PFMET.Px() + MetXCorr_JESDown;
    float PFMetYJESDown   = PFMET.Py() + MetYCorr_JESDown;
    MuonSystem->metJESDown    =  sqrt( pow(PFMetXJESDown,2) + pow(PFMetYJESDown,2) );
    MuonSystem->metPhiJESDown    = atan(PFMetYJESDown/PFMetXJESDown);
    if  (PFMetXJESUp < 0.0) MuonSystem->metPhiJESDown = ra.deltaPhi(TMath::Pi() + MuonSystem->metPhiJESDown,0.0);
    MuonSystem->metJESDownSF = helper->getMetTriggerSF(MuonSystem->metJESDown)/MuonSystem->metSF;

    //HEM
    float PFMetXHEM   = PFMET.Px() + MetXCorr_HEM;
    float PFMetYHEM   = PFMET.Py() + MetYCorr_HEM;
    MuonSystem->metHEM    = sqrt( pow(PFMetXHEM,2) + pow(PFMetYHEM,2) );
    MuonSystem->metPhiHEM    = atan(PFMetYHEM/PFMetXHEM);
    if  (PFMetXHEM < 0.0) MuonSystem->metPhiHEM = ra.deltaPhi(TMath::Pi() + MuonSystem->metPhiHEM,0.0);


    //EENoise
    float PFMetXEENoise   = PFMET.Px() + MetXCorr_EENoise;
    float PFMetYEENoise   = PFMET.Py() + MetYCorr_EENoise;
    MuonSystem->metEENoise    = sqrt( pow(PFMetXEENoise,2) + pow(PFMetYEENoise,2) );
    MuonSystem->metPhiEENoise    = atan(PFMetYEENoise/PFMetXEENoise);
    if  (PFMetXEENoise < 0.0) MuonSystem->metPhiEENoise = ra.deltaPhi(TMath::Pi() + MuonSystem->metPhiEENoise,0.0);


    // Do cluster Matchings
    for(int iCls = 0; iCls < MuonSystem->nCscRechitClusters3; iCls++) {

      MuonSystem->cscRechitCluster3JetVetoEta[iCls]     = -999.0;
      MuonSystem->cscRechitCluster3JetVetoPhi[iCls]     = -999.0;
      MuonSystem->cscRechitCluster3JetVetoPt[iCls]     = -999.0;
      MuonSystem->cscRechitCluster3JetVetoE[iCls]      = -999.0;
      MuonSystem->cscRechitCluster3MuonVetoEta[iCls]    = -999.0;
      MuonSystem->cscRechitCluster3MuonVetoPhi[iCls]    = -999.0;
      MuonSystem->cscRechitCluster3MuonVetoPt[iCls]    = -999.0;
      MuonSystem->cscRechitCluster3MuonVetoE[iCls]     = -999.0;

      // jet veto
      for(int i = 0; i < ra.nJets; i++) {
        if (fabs(ra.jetEta[i]>3.0)) continue;
        if (ra.deltaR(ra.jetEta[i], ra.jetPhi[i], MuonSystem->cscRechitCluster3Eta[iCls],MuonSystem->cscRechitCluster3Phi[iCls]) < 0.4
            && ra.jetPt[i] > MuonSystem->cscRechitCluster3JetVetoPt[iCls] ) {
          MuonSystem->cscRechitCluster3JetVetoPt[iCls]  = ra.jetPt[i];
          MuonSystem->cscRechitCluster3JetVetoEta[iCls]  = ra.jetEta[i];
          MuonSystem->cscRechitCluster3JetVetoPhi[iCls]  = ra.jetPhi[i];

        }
        if (ra.deltaR(ra.jetEta[i], ra.jetPhi[i], MuonSystem->cscRechitCluster3Eta[iCls], MuonSystem->cscRechitCluster3Phi[iCls]) < 0.4
            && ra.jetE[i] > MuonSystem->cscRechitCluster3JetVetoE[iCls] ) {
          MuonSystem->cscRechitCluster3JetVetoE[iCls]  = ra.jetE[i];
        }
      }
      // muon vetos
      for(int i = 0; i < ra.nMuons; i++) {
        if (fabs(ra.muonEta[i]>3.0)) continue;
        float muonIso = (ra.muon_chargedIso[i] + fmax(0.0,  ra.muon_photonIso[i] + ra.muon_neutralHadIso[i] - 0.5*ra.muon_pileupIso[i])) / ra.muonPt[i];
        if (ra.deltaR(ra.muonEta[i], ra.muonPhi[i], MuonSystem->cscRechitCluster3Eta[iCls],MuonSystem->cscRechitCluster3Phi[iCls]) < 0.8
            && ra.muonPt[i] > MuonSystem->cscRechitCluster3MuonVetoPt[iCls] ) {
          MuonSystem->cscRechitCluster3MuonVetoPt[iCls]  = ra.muonPt[i];
          MuonSystem->cscRechitCluster3MuonVetoE[iCls]  = ra.muonE[i];
          MuonSystem->cscRechitCluster3MuonVetoPhi[iCls]  = ra.muonPhi[i];
          MuonSystem->cscRechitCluster3MuonVetoEta[iCls]  = ra.muonEta[i];
          MuonSystem->cscRechitCluster3MuonVetoLooseIso[iCls]  = muonIso<0.25;
          MuonSystem->cscRechitCluster3MuonVetoTightIso[iCls]  = muonIso<0.15;
          MuonSystem->cscRechitCluster3MuonVetoVTightIso[iCls]  = muonIso<0.10;
          MuonSystem->cscRechitCluster3MuonVetoVVTightIso[iCls]  = muonIso<0.05;
          MuonSystem->cscRechitCluster3MuonVetoTightId[iCls]  = ra.isMuonPOGTightMuon(i);
          MuonSystem->cscRechitCluster3MuonVetoLooseId[iCls]  = ra.isMuonPOGLooseMuon(i);
        }
      }

      // match DT segments
      for (int i = 0; i < ra.nDtSeg; i++) {
        if (ra.deltaR(ra.dtSegEta[i], ra.dtSegPhi[i], MuonSystem->cscRechitCluster3Eta[iCls],MuonSystem->cscRechitCluster3Phi[iCls]) < 0.4 )
          {
            MuonSystem->cscRechitCluster3_match_dtSeg_0p4[iCls] ++;
            if (ra.dtSegStation[i] == 1) MuonSystem->cscRechitCluster3_match_MB1Seg_0p4[iCls] ++;
          }

        if (ra.deltaR(ra.dtSegEta[i], ra.dtSegPhi[i], MuonSystem->cscRechitCluster3Eta[iCls],MuonSystem->cscRechitCluster3Phi[iCls]) < 0.6 )
          {
            MuonSystem->cscRechitCluster3_match_dtSeg_0p6[iCls] ++;
            if (ra.dtSegStation[i] == 1) MuonSystem->cscRechitCluster3_match_MB1Seg_0p6[iCls] ++;
          }

      }
      //match to RPC hits in RE1/2
      for (int i = 0; i < ra.nRpc; i++) {
        float rpcR = sqrt(ra.rpcX[i]*ra.rpcX[i] + ra.rpcY[i]*ra.rpcY[i]);
        if (ra.deltaR(ra.rpcEta[i], ra.rpcPhi[i], MuonSystem->cscRechitCluster3Eta[iCls],MuonSystem->cscRechitCluster3Phi[iCls]) < 0.4 )
          {
            if (rpcR < 461.0 && rpcR > 275 && abs(ra.rpcZ[i]) > 663 && abs(ra.rpcZ[i]) < 730)
              {
                MuonSystem->cscRechitCluster3_match_RE12_0p4[iCls] ++;
              }
            if (rpcR < 470 && rpcR > 380 && abs(ra.rpcZ[i]) < 661)
              {
                MuonSystem->cscRechitCluster3_match_RB1_0p4[iCls] ++;
              }

          }
        if (ra.deltaR(ra.rpcEta[i], ra.rpcPhi[i], MuonSystem->cscRechitCluster3Eta[iCls],MuonSystem->cscRechitCluster3Phi[iCls]) < 0.6 )
          {
            if (rpcR < 461.0 && rpcR > 275 && abs(ra.rpcZ[i]) > 663 && abs(ra.rpcZ[i]) < 730)
              {
                MuonSystem->cscRechitCluster3_match_RE12_0p6[iCls] ++;
              }
            if (rpcR < 470 && rpcR > 380 && abs(ra.rpcZ[i]) < 661)
              {
                MuonSystem->cscRechitCluster3_match_RB1_0p6[iCls] ++;
              }
          }
      }
      MuonSystem->cscRechitCluster3Met_dPhi[iCls] =  ra.deltaPhi(MuonSystem->cscRechitCluster3Phi[iCls],MuonSystem->metPhi);
      MuonSystem->cscRechitCluster3MetXYCorr_dPhi[iCls] =  ra.deltaPhi(MuonSystem->cscRechitCluster3Phi[iCls],MuonSystem->metPhiXYCorr);
      MuonSystem->cscRechitCluster3MetEENoise_dPhi[iCls] =  ra.deltaPhi(MuonSystem->cscRechitCluster3Phi[iCls],MuonSystem->metPhiEENoise);
      MuonSystem->cscRechitCluster3MetHEM_dPhi[iCls] =  ra.deltaPhi(MuonSystem->cscRechitCluster3Phi[iCls],MuonSystem->metPhiHEM);

    }// end CSC cluster loop

    MuonSystem->nRpc = ra.nRpc;


    // Matchings for DT clusters 
    for(int iCls = 0; iCls < MuonSystem->nDtRechitClusters; iCls++) {
        //Jet veto/ muon veto
        MuonSystem->dtRechitClusterJetVetoEta[iCls] = -999.9;
        MuonSystem->dtRechitClusterJetVetoPhi[iCls] = -999.9;
        MuonSystem->dtRechitClusterJetVetoPt[iCls] = -999.9;
        MuonSystem->dtRechitClusterJetVetoE[iCls]  = -999.9;
        MuonSystem->dtRechitClusterMuonVetoEta[iCls]= -999.9;
        MuonSystem->dtRechitClusterMuonVetoPhi[iCls]= -999.9;
        MuonSystem->dtRechitClusterMuonVetoPt[iCls]= -999.9;
        MuonSystem->dtRechitClusterMuonVetoE[iCls] = -999.9;


        // jet veto
        for(int i = 0; i < ra.nJets; i++)
        {
          if (fabs(ra.jetEta[i]>3.0)) continue;
          if (ra.deltaR(ra.jetEta[i], ra.jetPhi[i], MuonSystem->dtRechitClusterEta[iCls],MuonSystem->dtRechitClusterPhi[iCls]) < 0.4 && ra.jetPt[i] > MuonSystem->dtRechitClusterJetVetoPt[iCls] ) {
            MuonSystem->dtRechitClusterJetVetoPt[iCls]  = ra.jetPt[i];
            MuonSystem->dtRechitClusterJetVetoEta[iCls]  = ra.jetEta[i];
            MuonSystem->dtRechitClusterJetVetoPhi[iCls]  = ra.jetPhi[i];
          }
          if (ra.deltaR(ra.jetEta[i], ra.jetPhi[i], MuonSystem->dtRechitClusterEta[iCls], MuonSystem->dtRechitClusterPhi[iCls]) < 0.4 && ra.jetE[i] > MuonSystem->dtRechitClusterJetVetoE[iCls] ) {
            MuonSystem->dtRechitClusterJetVetoE[iCls]  = ra.jetE[i];
          }
        }


        for(int i = 0; i < ra.nMuons; i++)
        {
          if (fabs(ra.muonEta[i]>3.0)) continue;
          if (ra.deltaR(ra.muonEta[i], ra.muonPhi[i], MuonSystem->dtRechitClusterEta[iCls],MuonSystem->dtRechitClusterPhi[iCls]) < 0.8 && ra.muonPt[i] > MuonSystem->dtRechitClusterMuonVetoPt[iCls] ) {
            MuonSystem->dtRechitClusterMuonVetoPt[iCls]  = ra.muonPt[i];
            MuonSystem->dtRechitClusterMuonVetoE[iCls]  = ra.muonE[i];
            MuonSystem->dtRechitClusterMuonVetoPhi[iCls]  = ra.muonPhi[i];
            MuonSystem->dtRechitClusterMuonVetoEta[iCls]  = ra.muonEta[i];
            MuonSystem->dtRechitClusterMuonVetoGlobal[iCls]  = ra.muon_isGlobal[i];
            MuonSystem->dtRechitClusterMuonVetoLooseId[iCls]  = ra.isMuonPOGLooseMuon(i);
          }
        }

        std::vector<int> dtRechitCluster_match_rpcBx;

        //match to RPC hits with dPhi<0.5 and same wheel in DT
        for (int i = 0; i < ra.nRpc; i++) {
          float rpcR = sqrt(ra.rpcX[i]*ra.rpcX[i] + ra.rpcY[i]*ra.rpcY[i]);
          if (ra.rpcRegion[i]!=0) continue;
          if (abs(ra.deltaPhi(ra.rpcPhi[i], MuonSystem->dtRechitClusterPhi[iCls])) < 0.5 )
          {
            if (ra.rpcRing[i] == MuonSystem->dtRechitClusterWheel[iCls])
            {
              dtRechitCluster_match_rpcBx.push_back(ra.rpcBx[i]);
              MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[iCls]++;
              if (rpcR < 470 && rpcR > 380 && abs(ra.rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_dPhi0p5[iCls] ++;

            }
          }
          if(ra.deltaR(ra.rpcEta[i], ra.rpcPhi[i], MuonSystem->dtRechitClusterEta[iCls], MuonSystem->dtRechitClusterPhi[iCls]) < 0.4 )
          {
            if (rpcR < 470 && rpcR > 380 && abs(ra.rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_0p4[iCls] ++;
          }

        }
        int max_occurence = 0;
        int max_bx = -999;
        for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++)
        {
          int counter = 0;
          for(unsigned int j = 0; j < dtRechitCluster_match_rpcBx.size(); j ++)
          {
            if (dtRechitCluster_match_rpcBx[j] == dtRechitCluster_match_rpcBx[l]) counter++;
          }
          if (counter>max_occurence)
          {
            max_occurence = counter;
            max_bx = dtRechitCluster_match_rpcBx[l];
          }
        }
          MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[iCls] = max_bx;

        MuonSystem->dtRechitClusterMetEENoise_dPhi[iCls] =  ra.deltaPhi(MuonSystem->dtRechitClusterPhi[iCls],MuonSystem->metPhiEENoise);
        if (MuonSystem->nLeptons > 0)
        {
          MuonSystem->dtRechitClusterLep_dPhi[iCls] =  ra.deltaPhi(MuonSystem->dtRechitClusterPhi[iCls],MuonSystem->lepPhi[0]);

        }

    }

      //Fill out tree if conditions are met
      outputTree->Fill();		
    }
    //save information
    outputFile->cd();
    outputTree->Write();
    cout << "Closing output file." << endl;
    outputFile->Close();
    delete outputFile;
}
