#include "RazorHelper.h"
#include "SusyLLPTree.h"
#include "assert.h"
#include "TTree.h"

#define N_MAX_LEPTONS 100
#define N_MAX_JETS 100
#define NTriggersMAX 601 //Number of trigger in the .dat file

// Constructor
SusyLLPTree::SusyLLPTree()
{
  InitVariables();
};

SusyLLPTree::~SusyLLPTree()
{
  if (f_) f_->Close();
};

//Initialize
void SusyLLPTree::InitVariables()
{ 
    //evnt info
    runNum=0; lumiSec=0; evtNum=0; category=0;
    npv=0; npu=0; rho=-1; weight=-1;
    met=-1; metPhi=-1;

    //leptons
    nLeptons = 0;
    for( int i = 0; i < N_MAX_LEPTONS; i++ )
    {
      lepE[i]      = -999.;
      lepPt[i]     = -999.;
      lepEta[i]    = -999.;
      lepPhi[i]    = -999.;
      lepPdgId[i]  = -999;
      lepDZ[i]     = -999.;
      // lepLoosePassId[i] = false;
      // lepMediumPassId[i] = false;
      // lepTightPassId[i] = false;
      lepPassId[i] = false;
    }

    //Z-candidate
    ZMass = -999.; ZPt = -999.; ZEta = -999.; ZPhi = -999.;
    MT = -999.;
    ZleptonIndex1 = -999; ZleptonIndex2 = -999;

    //jets
    nJets = 0;
    nCaloJets = 0;
    for( int i = 0; i < N_MAX_JETS; i++ )
    {
      jetE[i]      = -999.;
      jetEt[i] = -999.;
      jetPt[i]     = -999.;
      jetEta[i]    = -999.;
      jetPhi[i]    = -999.;
      jetTime[i]   = -999.;
      // jetLoosePassId[i] = false;
      jetPassId[i] = false;
      matched[i] = false;
      jet_energy_frac[i] = -999.;
      jet_sig_et1[i] = -999.;
      jet_sig_et2[i] = -999.;
      ecalNRechits[i] = -999;
      ecalRechitE[i] = -999.;
      jetChargedEMEnergyFraction[i] = -999.;
      jetNeutralEMEnergyFraction[i] = -999.;
      jetChargedHadronEnergyFraction[i] = -999.;
      jetNeutralHadronEnergyFraction[i] = -999.;
      jetGammaMax_ET[i] = -999.;
      jetMinDeltaRPVTracks[i] = -999.;
      jetPtAllPVTracks[i] = -999.;

      calojetE[i] = -999.;
      calojetEt[i] = -999.;
      calojetPt[i] = -999.;
      calojetEta[i] = -999.;
      calojetPhi[i] = -999.;
      calojetTime[i] = -999.;
      calojetPassId[i] = -999.;
      calojetGammaMax_ET[i] = -999.;
      calojetMinDeltaRPVTracks[i] = -999.;
      calojetPtAllPVTracks[i] = -999.;
      calojetNRechits[i] = -999.;
      calojetRechitE[i] = -999.;
    }

    for(int i = 0; i <NTriggersMAX; i++){
      HLTDecision[i] = false;
    }

};

//SetBranchAddress
void SusyLLPTree::InitTree()
{
    assert(tree_);
    InitVariables();

    //event info
    tree_->SetBranchAddress("runNum",      &runNum);
    tree_->SetBranchAddress("lumiSec",     &lumiSec);
    tree_->SetBranchAddress("evtNum",      &evtNum);
    tree_->SetBranchAddress("category",    &category);
    tree_->SetBranchAddress("npv",         &npv);
    tree_->SetBranchAddress("npu",         &npu);
    tree_->SetBranchAddress("weight",      &weight);
    tree_->SetBranchAddress("rho",         &rho);
    tree_->SetBranchAddress("met",         &met);
    tree_->SetBranchAddress("metPhi",      &metPhi);

    //Leptons
    tree_->SetBranchAddress("nLeptons",    &nLeptons);
    tree_->SetBranchAddress("lepE",        lepE);
    tree_->SetBranchAddress("lepPt",       lepPt);
    tree_->SetBranchAddress("lepEta",      lepEta);
    tree_->SetBranchAddress("lepPhi",      lepPhi);
    tree_->SetBranchAddress("lepPdgId",  lepPdgId);
    tree_->SetBranchAddress("lepDZ",     lepDZ);
    // tree_->SetBranchAddress("lepLoosePassId", lepLoosePassId);
    // tree_->SetBranchAddress("lepMediumPassId", lepMediumPassId);
    // tree_->SetBranchAddress("lepTightPassId", lepTightPassId);
    tree_->SetBranchAddress("lepPassId", lepPassId);

    //Z-candidate
    tree_->SetBranchAddress("ZMass",       &ZMass);
    tree_->SetBranchAddress("ZPt",         &ZPt);
    tree_->SetBranchAddress("ZEta",        &ZEta);
    tree_->SetBranchAddress("ZPhi",        &ZPhi);
    tree_->SetBranchAddress("ZleptonIndex1", &ZleptonIndex1);
    tree_->SetBranchAddress("ZleptonIndex2", &ZleptonIndex2);
    tree_->SetBranchAddress("MT", &MT);

    //jets
    tree_->SetBranchAddress("nJets",     &nJets);
    tree_->SetBranchAddress("jetE",      jetE);
    tree_->SetBranchAddress("jetEt",      jetEt);

    tree_->SetBranchAddress("jetPt",     jetPt);
    tree_->SetBranchAddress("jetEta",    jetEta);
    tree_->SetBranchAddress("jetPhi",    jetPhi);
    tree_->SetBranchAddress("jetTime",   jetTime);
    tree_->SetBranchAddress("jetPassId", jetPassId);
    tree_->SetBranchAddress("matched", matched);
    tree_->SetBranchAddress("jet_energy_frac", jet_energy_frac);
    tree_->SetBranchAddress("jet_sig_et1", jet_sig_et1);
    tree_->SetBranchAddress("jet_sig_et2", jet_sig_et2);

    tree_->SetBranchAddress("ecalNRechits",   ecalNRechits);
    tree_->SetBranchAddress("ecalRechitE", ecalRechitE);
    tree_->SetBranchAddress("jetChargedEMEnergyFraction", jetChargedEMEnergyFraction);
    tree_->SetBranchAddress("jetNeutralEMEnergyFraction", jetNeutralEMEnergyFraction);
    tree_->SetBranchAddress("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction);
    tree_->SetBranchAddress("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction);
    tree_->SetBranchAddress("jetGammaMax_ET", jetGammaMax_ET);
    tree_->SetBranchAddress("jetMinDeltaRPVTracks", jetMinDeltaRPVTracks);
    tree_->SetBranchAddress("jetPtAllPVTracks", jetPtAllPVTracks);

    //calojets
    tree_->SetBranchAddress("nCaloJets",     &nCaloJets);
    tree_->SetBranchAddress("calojetE",      calojetE);
    tree_->SetBranchAddress("calojetEt",      calojetEt);
    tree_->SetBranchAddress("calojetPt",     calojetPt);
    tree_->SetBranchAddress("calojetEta",    calojetEta);
    tree_->SetBranchAddress("calojetPhi",    calojetPhi);
    tree_->SetBranchAddress("calojetTime",   calojetTime);
    tree_->SetBranchAddress("calojetPassId", calojetPassId);
    tree_->SetBranchAddress("calojetRechitE",   calojetRechitE);
    // tree_->SetBranchAddress("calojetRechitT_rms", calojetRechitT_rms);
    // tree_->SetBranchAddress("calojet_EMEnergyFraction", calojet_EMEnergyFraction);
    // tree_->SetBranchAddress("calojet_HadronicEnergyFraction", calojet_HadronicEnergyFraction);

    tree_->SetBranchAddress("calojetGammaMax_ET", calojetGammaMax_ET);
    tree_->SetBranchAddress("calojetMinDeltaRPVTracks", calojetMinDeltaRPVTracks);
    tree_->SetBranchAddress("calojetPtAllPVTracks", calojetPtAllPVTracks);
    // tree_->SetBranchAddress("jetLoosePassId", jetLoosePassId);
    // tree_->SetBranchAddress("jetTightPassId", jetTightPassId);
    // triggers
    tree_->SetBranchAddress("HLTDecision",   HLTDecision);
/*
  // gLLP
  tree_->SetBranchAddress("gLLP_eta",    gLLP_eta);
  tree_->SetBranchAddress("gLLP_phi",    gLLP_phi);
  tree_->SetBranchAddress("gLLP_csc",    gLLP_csc);
  tree_->SetBranchAddress("gLLP_ctau",    gLLP_ctau);

  tree_->SetBranchAddress("gLLP_decay_vertex_r",    gLLP_decay_vertex_r);
  tree_->SetBranchAddress("gLLP_decay_vertex_x",    gLLP_decay_vertex_x);
  tree_->SetBranchAddress("gLLP_decay_vertex_y",    gLLP_decay_vertex_y);
  tree_->SetBranchAddress("gLLP_decay_vertex_z",    gLLP_decay_vertex_z);
*/
};

//LoadTree
void SusyLLPTree::LoadTree(const char* file)
{
    f_ = TFile::Open(file);
    assert(f_);
    tree_ = dynamic_cast<TTree*>(f_->Get("SusyLLPTree"));
    InitTree();
    assert(tree_);
};

//CreateTree
void SusyLLPTree::CreateTree()
{
    //tree
    tree_ = new TTree("SusyLLPTree","SusyLLPTree");
    f_ = 0;

    //event info
    tree_->Branch("runNum",      &runNum,     "runNum/i");      // event run number
    tree_->Branch("lumiSec",     &lumiSec,    "lumiSec/i");     // event lumi section
    tree_->Branch("evtNum",      &evtNum,     "evtNum/i");      // event number
    tree_->Branch("category",    &category,   "category/i");    // dilepton category
    tree_->Branch("npv",         &npv,        "npv/i");         // number of primary vertices
    tree_->Branch("npu",         &npu,        "npu/i");         // number of in-time PU events (MC)
    tree_->Branch("weight",      &weight,     "weight/F");
    tree_->Branch("rho",         &rho,        "rho/F");
    tree_->Branch("met",         &met,        "met/F");         // MET
    tree_->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)

    //leptons
    tree_->Branch("nLeptons",  &nLeptons, "nLeptons/I");
    tree_->Branch("lepE",      lepE,      "lepE[nLeptons]/F");
    tree_->Branch("lepPt",     lepPt,     "lepPt[nLeptons]/F");
    tree_->Branch("lepEta",    lepEta,    "lepEta[nLeptons]/F");
    tree_->Branch("lepPhi",    lepPhi,    "lepPhi[nLeptons]/F");
    tree_->Branch("lepPdgId",  lepPdgId,  "lepPdgId[nLeptons]/I");
    tree_->Branch("lepDZ",     lepDZ,     "lepDZ[nLeptons]/F");
    tree_->Branch("lepPassId", lepPassId, "lepPassId[nLeptons]/O");
    // tree_->Branch("lepLoosePassId", lepLoosePassId, "lepLoosePassId[nLeptons]/O");
    // tree_->Branch("lepMediumPassId", lepMediumPassId, "lepMediumPassId[nLeptons]/O");
    // tree_->Branch("lepTightPassId", lepTightPassId, "lepTightPassId[nLeptons]/O");

    //Z-candidate
    tree_->Branch("MT",      &MT,  "MT/F");
    tree_->Branch("ZMass",      &ZMass,  "ZMass/F");
    tree_->Branch("ZPt",        &ZPt,    "ZPt/F");
    tree_->Branch("ZEta",       &ZEta,   "ZEta/F");
    tree_->Branch("ZPhi",       &ZPhi,   "ZPhi/F");
    tree_->Branch("ZleptonIndex1", &ZleptonIndex1, "ZleptonIndex1/I");
    tree_->Branch("ZleptonIndex2", &ZleptonIndex2, "ZleptonIndex2/I");

    //jets
    tree_->Branch("nJets",     &nJets,    "nJets/I");
    tree_->Branch("jetE",      jetE,      "jetE[nJets]/F");
    tree_->Branch("jetEt",      jetEt,      "jetEt[nJets]/F");
    tree_->Branch("jetPt",     jetPt,     "jetPt[nJets]/F");
    tree_->Branch("jetEta",    jetEta,    "jetEta[nJets]/F");
    tree_->Branch("jetPhi",    jetPhi,    "jetPhi[nJets]/F");
    tree_->Branch("jetTime",   jetTime,   "jetTime[nJets]/F");
    tree_->Branch("jetPassId", jetPassId, "jetPassId[nJets]/O");
    tree_->Branch("matched", matched, "matched[nJets]/O");
    tree_->Branch("jet_energy_frac", jet_energy_frac, "jet_energy_frac[nJets]/F");
    tree_->Branch("jet_sig_et1", jet_sig_et1, "jet_sig_et1[nJets]/F");
    tree_->Branch("jet_sig_et2", jet_sig_et2, "jet_sig_et2[nJets]/F");


    tree_->Branch("jetGammaMax_ET", jetGammaMax_ET, "jetGammaMax_ET[nJets]/F");
    tree_->Branch("jetMinDeltaRPVTracks", jetMinDeltaRPVTracks, "jetMinDeltaRPVTracks[nJets]/F");
    tree_->Branch("jetPtAllPVTracks", jetPtAllPVTracks, "jetPtAllPVTracks[nJets]/F");
    tree_->Branch("ecalNRechits",   ecalNRechits,   "ecalNRechits[nJets]/F");
    tree_->Branch("ecalRechitE", ecalRechitE, "ecalRechitE[nJets]/F");
    tree_->Branch("jetChargedEMEnergyFraction",   jetChargedEMEnergyFraction,   "jetChargedEMEnergyFraction[nJets]/F");
    tree_->Branch("jetNeutralEMEnergyFraction",   jetNeutralEMEnergyFraction,   "jetNeutralEMEnergyFraction[nJets]/F");
    tree_->Branch("jetChargedHadronEnergyFraction",   jetChargedHadronEnergyFraction,   "jetChargedHadronEnergyFraction[nJets]/F");
    tree_->Branch("jetNeutralHadronEnergyFraction",   jetNeutralHadronEnergyFraction,   "jetNeutralHadronEnergyFraction[nJets]/F");
    // tree_->Branch("jetLoosePassId", jetLoosePassId, "jetLoosePassId[nJets]/O");
    // tree_->Branch("jetTightPassId", jetTightPassId, "jetTightPassId[nJets]/O");


    //calojet
    tree_->Branch("nCaloJets",     &nCaloJets,    "nCaloJets/I");
    tree_->Branch("calojetE",      calojetE,      "calojetE[nCaloJets]/F");
    tree_->Branch("calojetEt",      calojetEt,      "calojetEt[nCaloJets]/F");
    tree_->Branch("calojetPt",     calojetPt,     "calojetPt[nCaloJets]/F");
    tree_->Branch("calojetEta",    calojetEta,    "calojetEta[nCaloJets]/F");
    tree_->Branch("calojetPhi",    calojetPhi,    "calojetPhi[nCaloJets]/F");
    tree_->Branch("calojetTime",   calojetTime,   "calojetTime[nCaloJets]/F");
    tree_->Branch("calojetPassId", calojetPassId, "calojetPassId[nCaloJets]/O");
    tree_->Branch("calojetGammaMax_ET", calojetGammaMax_ET, "calojetGammaMax_ET[nCaloJets]/F");
    tree_->Branch("calojetMinDeltaRPVTracks", calojetMinDeltaRPVTracks, "calojetMinDeltaRPVTracks[nCaloJets]/F");
    tree_->Branch("calojetPtAllPVTracks", calojetPtAllPVTracks, "calojetPtAllPVTracks[nCaloJets]/F");
    tree_->Branch("calojetNRechits",   calojetNRechits,   "calojetNRechits[nCaloJets]/F");
    tree_->Branch("calojetRechitE", calojetRechitE, "calojetRechitE[nCaloJets]/F");

    //HLT
    tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[601]/O"); //hardcoded
/*
  //gLLP branches
  tree_->Branch("gLLP_eta",          gLLP_eta,          "gLLP_eta[2]/F");
  tree_->Branch("gLLP_phi",          gLLP_phi,          "gLLP_phi[2]/F");
  tree_->Branch("gLLP_ctau",          gLLP_ctau,          "gLLP_ctau[2]/F");

  tree_->Branch("gLLP_decay_vertex_r",          gLLP_decay_vertex_r,          "gLLP_decay_vertex_r[2]/F");
  tree_->Branch("gLLP_decay_vertex_x",          gLLP_decay_vertex_x,          "gLLP_decay_vertex_x[2]/F");
  tree_->Branch("gLLP_decay_vertex_y",          gLLP_decay_vertex_y,          "gLLP_decay_vertex_y[2]/F");
  tree_->Branch("gLLP_decay_vertex_z",          gLLP_decay_vertex_z,          "gLLP_decay_vertex_z[2]/F");
*/
};
