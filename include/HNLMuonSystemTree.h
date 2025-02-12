#ifndef HNLMuonSystemTree_H
#define HNLMuonSystemTree_H

#define N_MAX_LEPTONS 200
#define N_MAX_JETS 100
#define N_MAX_CSC 20
#define N_MAX_CSCRECHITS 1000
#define N_MAX_DTRECHITS 100
#define NTriggersMAX 1201 // Number of trigger in the .dat file
#define N_MAX_GPARTICLES 5000

#include <iostream>
#include <string>
#include <sys/stat.h>
#include "assert.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TRandom.h"
#include "FactorizedJetCorrector.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
#include "SimpleJetResolution.h"
#include "BTagCalibrationStandalone.h"
#include "TTree.h"
#include "DBSCAN.h"

#include "RazorAnalyzer.h"

#include "RazorHelper.h"

class HNLMuonSystemTree
{

public:
  HNLMuonSystemTree();
  ~HNLMuonSystemTree();

  TTree *tree_;
  TFile *f_;

  UInt_t  runNum, lumiSec;
  //ULong64_t evtNum;
  UInt_t evtNum; ULong64_t evtNumLong;
  UInt_t  category;
  UInt_t  npv, npu;
  float rho, weight;
  float pileupWeight, pileupWeightUp, pileupWeightDown;
  float scaleWeights[9];
  int ZCategory;
  float metSF;
  float lepOverallSF;
  float sf_facScaleUp, sf_facScaleDown, sf_renScaleUp, sf_renScaleDown, sf_facRenScaleUp, sf_facRenScaleDown;
  bool EE_prefiring;
  float met, metPhi, HT, jetMet_dPhi, jetMet_dPhiMin,jetMet_dPhiMin4, metJESUp, metPhiJESUp, metJESDown, metPhiJESDown, metJESDownSF, metJESUpSF, metNoMu, metEENoise, metPhiEENoise, metHEM, metPhiHEM, metHEMXYCorr, metPhiHEMXYCorr;
  float metXYCorr, metPhiXYCorr, metEENoiseXYCorr, metPhiEENoiseXYCorr;
  bool Flag_HBHENoiseFilter, Flag_HBHEIsoNoiseFilter, Flag_BadPFMuonFilter, Flag_globalSuperTightHalo2016Filter,
  Flag_CSCTightHaloFilter, Flag_BadChargedCandidateFilter, Flag_eeBadScFilter, Flag_goodVertices, Flag_ecalBadCalibFilter, Flag_all;
  int mH, mX, ctau;

  bool Flag2_HBHENoiseFilter, Flag2_HBHEIsoNoiseFilter, Flag2_BadPFMuonFilter, Flag2_globalSuperTightHalo2016Filter,
  Flag2_globalTightHalo2016Filter, Flag2_BadChargedCandidateFilter, Flag2_EcalDeadCellTriggerPrimitiveFilter,
  Flag2_ecalBadCalibFilter, Flag2_eeBadScFilter, Flag2_all;

  float gWPt;
  int gLepId;
  float gLepPt, gLepPhi, gLepEta, gLepE;
  int nGenParticle;
  int gParticleId[N_MAX_GPARTICLES];
  int gParticleStatus[N_MAX_GPARTICLES];
  int gParticleMotherId[N_MAX_GPARTICLES];
  int gParticleMotherIndex[N_MAX_GPARTICLES];
  float gParticlePt[N_MAX_GPARTICLES];
  float gParticleEta[N_MAX_GPARTICLES];
  float gParticlePhi[N_MAX_GPARTICLES];
  float gParticleE[N_MAX_GPARTICLES];
  float gParticleProdVertexX[N_MAX_GPARTICLES];
  float gParticleProdVertexY[N_MAX_GPARTICLES];
  float gParticleProdVertexZ[N_MAX_GPARTICLES];

  float genMetPtTrue;
  float genMetPhiTrue;
  float genMetPtCalo;
  float genMetPhiCalo;

  //
   int nGenJets;
   float genJetE[N_MAX_GPARTICLES];
   float genJetPt[N_MAX_GPARTICLES];
   float genJetEta[N_MAX_GPARTICLES];
   float genJetPhi[N_MAX_GPARTICLES];
  //
  //
  //
   float genJetMET[N_MAX_GPARTICLES];

  //dt


  int           nDTRechits;
  float         dtPosTpeak;
  float         dtNegTpeak;
  int           nDtRings;
  int           nDTPositiveYRechits;
  int           nDTNegativeYRechits;
  int           nDTRechitsChamberMinus12;
  int           nDTRechitsChamberMinus11;
  int           nDTRechitsChamber10;
  int           nDTRechitsChamberPlus11;
  int           nDTRechitsChamberPlus12;
  int           nDTRechitsChamberMinus22;
  int           nDTRechitsChamberMinus21;
  int           nDTRechitsChamber20;
  int           nDTRechitsChamberPlus21;
  int           nDTRechitsChamberPlus22;
  int           nDTRechitsChamberMinus32;
  int           nDTRechitsChamberMinus31;
  int           nDTRechitsChamber30;
  int           nDTRechitsChamberPlus31;
  int           nDTRechitsChamberPlus32;
  int           nDTRechitsChamberMinus42;
  int           nDTRechitsChamberMinus41;
  int           nDTRechitsChamber40;
  int           nDTRechitsChamberPlus41;
  int           nDTRechitsChamberPlus42;
  int           nDTRechitsStation1;
  int           nDTRechitsStation2;
  int           nDTRechitsStation3;
  int           nDTRechitsStation4;
  int           nDTRechitsWheelMinus2;
  int           nDTRechitsWheelMinus1;
  int           nDTRechitsWheel0;
  int           nDTRechitsWheelPlus1;
  int           nDTRechitsWheelPlus2;
  int           nDtStations25;
  int           nDtWheels25;

  // rpc
  int           nRpc;
  float         rpcX[N_MAX_CSCRECHITS];   //[npc]
  float         rpcY[N_MAX_CSCRECHITS];   //[npc]
  float         rpcZ[N_MAX_CSCRECHITS];   //[npc]
  float         rpcPhi[N_MAX_CSCRECHITS];   //[npc]
  float         rpcEta[N_MAX_CSCRECHITS];   //[npc]
  int           rpcBx[N_MAX_CSCRECHITS];   //[npc]
  int           rpcStation[N_MAX_CSCRECHITS];   //[npc]
  int           rpcRing[N_MAX_CSCRECHITS];   //[npc]
  int           rpcSector[N_MAX_CSCRECHITS];   //[npc]
  int           rpcLayer[N_MAX_CSCRECHITS];   //[npc]


  // DT

  int nDtRechitClusters;


  float     dtRechitClusterX[N_MAX_DTRECHITS];
  float     dtRechitClusterY[N_MAX_DTRECHITS];
  float     dtRechitClusterZ[N_MAX_DTRECHITS];
  int       dtRechitClusterWheel[N_MAX_DTRECHITS];
  float     dtRechitClusterEta[N_MAX_DTRECHITS];
  float     dtRechitClusterPhi[N_MAX_DTRECHITS];
  int       dtRechitClusterSize[N_MAX_DTRECHITS];

  int       dtRechitClusterNSegmentStation1[N_MAX_DTRECHITS];
  int       dtRechitClusterNSegmentStation2[N_MAX_DTRECHITS];
  int       dtRechitClusterNSegmentStation3[N_MAX_DTRECHITS];
  int       dtRechitClusterNSegmentStation4[N_MAX_DTRECHITS];

  int       dtRechitClusterMaxChamber[N_MAX_DTRECHITS];
  float     dtRechitClusterMaxChamberRatio[N_MAX_DTRECHITS];
  int       dtRechitClusterNChamber[N_MAX_DTRECHITS];
  int       dtRechitClusterMaxStation[N_MAX_DTRECHITS];
  float     dtRechitClusterMaxStationRatio[N_MAX_DTRECHITS];
  int       dtRechitClusterNStation[N_MAX_DTRECHITS];
  int       dtRechitClusterNStation10[N_MAX_DTRECHITS];
  float       dtRechitClusterAvgStation10[N_MAX_DTRECHITS];

  float     dtRechitClusterJetVetoPt[N_MAX_DTRECHITS];
  float     dtRechitClusterJetVetoE[N_MAX_DTRECHITS];
  float     dtRechitClusterJetVetoPhi[N_MAX_DTRECHITS];
  float     dtRechitClusterJetVetoEta[N_MAX_DTRECHITS];

  float     dtRechitClusterMuonVetoPt[N_MAX_DTRECHITS];
  float     dtRechitClusterMuonVetoE[N_MAX_DTRECHITS];
  float     dtRechitClusterMuonVetoPhi[N_MAX_DTRECHITS];
  float     dtRechitClusterMuonVetoEta[N_MAX_DTRECHITS];
  bool      dtRechitClusterMuonVetoGlobal[N_MAX_DTRECHITS];
  bool      dtRechitClusterMuonVetoLooseId[N_MAX_DTRECHITS];


  float     dtRechitCluster_match_gParticle_deltaR[N_MAX_DTRECHITS];
  int       dtRechitCluster_match_gParticle_Id[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gParticle_Pt[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gParticle_Eta[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gParticle_Phi[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gParticle_E[N_MAX_DTRECHITS];
  int       dtRechitCluster_match_gParticle_Status[N_MAX_DTRECHITS];
  int       dtRechitCluster_match_gParticle_MotherId[N_MAX_DTRECHITS];

  bool      dtRechitCluster_match_gLLP[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gLLP_minDeltaR[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gLLP_eta[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gLLP_phi[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gLLP_decay_x[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gLLP_decay_y[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gLLP_decay_z[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gLLP_ctau[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gLLP_beta[N_MAX_DTRECHITS];
  bool      dtRechitCluster_match_gLLP_csc[N_MAX_DTRECHITS];
  bool      dtRechitCluster_match_gLLP_dt[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gLLP_e[N_MAX_DTRECHITS];
  float     dtRechitCluster_match_gLLP_pt[N_MAX_DTRECHITS];
  float     dtRechitClusterLep_dPhi[N_MAX_DTRECHITS];


  int       dtRechitCluster_match_MB1hits_0p5[N_MAX_DTRECHITS];
  int       dtRechitCluster_match_MB1hits_0p4[N_MAX_DTRECHITS];
  int       dtRechitCluster_match_MB1hits_cosmics_plus[N_MAX_DTRECHITS];
  int       dtRechitCluster_match_MB1hits_cosmics_minus[N_MAX_DTRECHITS];
  int       dtRechitCluster_match_RPCBx_dPhi0p5[N_MAX_DTRECHITS];
  int       dtRechitCluster_match_RB1_0p4[N_MAX_DTRECHITS];
  int       dtRechitCluster_match_RB1_dPhi0p5[N_MAX_DTRECHITS];
  int       dtRechitCluster_match_RPChits_dPhi0p5[N_MAX_DTRECHITS];
  float     dtRechitClusterMetEENoise_dPhi[N_MAX_DTRECHITS];








  //csc
  int           nCscRechits;
  int           nEarlyCscRechits;
  int           nLateCscRechits;
  int           nEarly2CscRechits;
  int           nLate2CscRechits;
  int           nCscRings;
  int           nCscPositiveYRechits;
  int           nCscNegativeYRechits;
  float         cscPosTpeak;
  float         cscNegTpeak;
  int           nCscRechitsChamberPlus11;
  int           nCscRechitsChamberPlus12;
  int           nCscRechitsChamberPlus13;
  int           nCscRechitsChamberPlus21;
  int           nCscRechitsChamberPlus22;
  int           nCscRechitsChamberPlus31;
  int           nCscRechitsChamberPlus32;
  int           nCscRechitsChamberPlus41;
  int           nCscRechitsChamberPlus42;

  int           nCscRechitsChamberMinus11;
  int           nCscRechitsChamberMinus12;
  int           nCscRechitsChamberMinus13;
  int           nCscRechitsChamberMinus21;
  int           nCscRechitsChamberMinus22;
  int           nCscRechitsChamberMinus31;
  int           nCscRechitsChamberMinus32;
  int           nCscRechitsChamberMinus41;
  int           nCscRechitsChamberMinus42;



  int           nCscRechitClusters3;

  int           cscRechitCluster3_match_Me1112_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_Me1112_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_Me1112_0p8[N_MAX_CSC];
  int           cscRechitCluster3_match_Me11_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_Me11_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_Me11_0p8[N_MAX_CSC];
  int           cscRechitCluster3_match_Me12_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_Me12_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_Me12_0p8[N_MAX_CSC];
  int           cscRechitCluster3_match_cscRechits_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_cscSeg_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_ME11Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_ME12Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_cscSeg_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_ME11Seg_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_ME12Seg_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_dtRechits_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_dtRechits_phi0p2[N_MAX_CSC];
  int           cscRechitCluster3_match_MB1_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_dtRechits_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_MB1_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_dtSeg_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_MB1Seg_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_dtSeg_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_MB1Seg_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_RB1_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_RE12_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_RB1_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_RE12_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_highEta_0p4[N_MAX_CSC];
  int           cscRechitCluster3_match_highEta_0p6[N_MAX_CSC];
  int           cscRechitCluster3_match_highEta_0p8[N_MAX_CSC];
  float           cscRechitCluster3_match_cluster_dR[N_MAX_CSC];
  int           cscRechitCluster3_match_cluster_index[N_MAX_CSC];
  bool          cscRechitCluster3_match_gParticle[N_MAX_CSC];
  float         cscRechitCluster3_match_gParticle_minDeltaR[N_MAX_CSC];
  int           cscRechitCluster3_match_gParticle_index[N_MAX_CSC];
  int           cscRechitCluster3_match_gParticle_id[N_MAX_CSC];
  float         cscRechitCluster3_match_gParticle_eta[N_MAX_CSC];
  float         cscRechitCluster3_match_gParticle_phi[N_MAX_CSC];
  float         cscRechitCluster3_match_gParticle_E[N_MAX_CSC];
  float         cscRechitCluster3_match_gParticle_pt[N_MAX_CSC];
  int           cscRechitCluster3_match_gParticle_MotherId[N_MAX_CSC];

  bool          cscRechitCluster3_match_gLLP[N_MAX_CSC];
  int           cscRechitCluster3_match_gLLP_index[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_minDeltaR[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_eta[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_phi[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_decay_r[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_decay_x[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_decay_y[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_decay_z[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_ctau[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_beta[N_MAX_CSC];
  bool         cscRechitCluster3_match_gLLP_csc[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_e[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_pt[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_lepdPhi[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_daughter0_deltaR[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_daughter1_deltaR[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_daughter2_deltaR[N_MAX_CSC];
  float         cscRechitCluster3_match_gLLP_daughter3_deltaR[N_MAX_CSC];

  float         cscRechitCluster3X[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3Y[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3Z[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3Time[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3TimeTotal[N_MAX_CSC];
  float         cscRechitCluster3TimeWeighted[N_MAX_CSC];
  float         cscRechitCluster3TimeWire[N_MAX_CSC];
  float         cscRechitCluster3TimeWirePruned[N_MAX_CSC];

  float         cscRechitCluster3TimeSpread[N_MAX_CSC];
  float         cscRechitCluster3TimeSpreadWeighted[N_MAX_CSC];
  float         cscRechitCluster3TimeSpreadWeightedAll[N_MAX_CSC];
  float         cscRechitCluster3TimeWireSpread[N_MAX_CSC];
  float         cscRechitCluster3TimeTotalSpread[N_MAX_CSC];
  float         cscRechitCluster3TimeTotalSpreadPruned[N_MAX_CSC];

  float         cscRechitCluster3GenMuonDeltaR[N_MAX_CSC];
  float         cscRechitCluster3MajorAxis[N_MAX_CSC];
  float         cscRechitCluster3MinorAxis[N_MAX_CSC];
  float         cscRechitCluster3XSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3YSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3ZSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3XYSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3RSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3DeltaRSpread[N_MAX_CSC];


  float         cscRechitCluster3EtaPhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3EtaSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3PhiSpread[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3Eta[N_MAX_CSC];   //[nCsc]
  float         cscRechitCluster3Phi[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster3Size[N_MAX_CSC];
  float         cscRechitCluster3Me11Ratio[N_MAX_CSC];
  float         cscRechitCluster3Me12Ratio[N_MAX_CSC];


  float         cscRechitCluster3MaxStationRatio[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster3MaxStation[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster3NStation[N_MAX_CSC];
  int           cscRechitCluster3NStation5[N_MAX_CSC];
  int           cscRechitCluster3NStation10[N_MAX_CSC];
  int           cscRechitCluster3NStation10perc[N_MAX_CSC];
  float          cscRechitCluster3AvgStation[N_MAX_CSC];
  float          cscRechitCluster3AvgStation5[N_MAX_CSC];
  float          cscRechitCluster3AvgStation10[N_MAX_CSC];
  float          cscRechitCluster3AvgStation10perc[N_MAX_CSC];
  float         cscRechitCluster3MaxChamberRatio[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster3MaxChamber[N_MAX_CSC];   //[nCsc]
  int           cscRechitCluster3NChamber[N_MAX_CSC];

  float         cscRechitCluster3JetVetoElectronEnergyFraction[N_MAX_CSC];
  float         cscRechitCluster3JetVetoPhotonEnergyFraction[N_MAX_CSC];
  float         cscRechitCluster3JetVetoChargedHadronEnergyFraction[N_MAX_CSC];
  float         cscRechitCluster3JetVetoNeutralHadronEnergyFraction[N_MAX_CSC];
  float         cscRechitCluster3JetVetoMuonEnergyFraction[N_MAX_CSC];
  float         cscRechitCluster3JetVetoEta[N_MAX_CSC];
  float         cscRechitCluster3JetVetoPhi[N_MAX_CSC];
float         cscRechitCluster3JetVetoPt[N_MAX_CSC];
  float         cscRechitCluster3JetVetoE[N_MAX_CSC];
  float         cscRechitCluster3GenJetVetoPt[N_MAX_CSC];
  float         cscRechitCluster3GenJetVetoE[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoPt[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoE[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoPhi[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoEta[N_MAX_CSC];
  float         cscRechitCluster3JetVetoPt_0p6[N_MAX_CSC];
  float         cscRechitCluster3JetVetoPt_0p8[N_MAX_CSC];
  float         cscRechitCluster3JetVetoE_0p6[N_MAX_CSC];
  float         cscRechitCluster3JetVetoE_0p8[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoPt_0p6[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoPt_0p8[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoE_0p6[N_MAX_CSC];
  float         cscRechitCluster3MuonVetoE_0p8[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1[N_MAX_CSC];
  bool         cscRechitCluster3ZLep2[N_MAX_CSC];
  int         cscRechitCluster3ZLep1Id[N_MAX_CSC];
  int         cscRechitCluster3ZLep2Id[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1TightId[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1LooseIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1TightIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1VTightIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1VVTightIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep1Tag[N_MAX_CSC];

  bool         cscRechitCluster3ZLep2TightId[N_MAX_CSC];
  bool         cscRechitCluster3ZLep2LooseIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep2TightIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep2VTightIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep2VVTightIso[N_MAX_CSC];
  bool         cscRechitCluster3ZLep2Tag[N_MAX_CSC];
  bool          cscRechitCluster3MuonVetoLooseIso[N_MAX_CSC];
  bool          cscRechitCluster3MuonVetoTightIso[N_MAX_CSC];
  bool          cscRechitCluster3MuonVetoVTightIso[N_MAX_CSC];
  bool          cscRechitCluster3MuonVetoVVTightIso[N_MAX_CSC];
  bool          cscRechitCluster3MuonVetoTightId[N_MAX_CSC];
  bool          cscRechitCluster3MuonVetoLooseId[N_MAX_CSC];

  bool          cscRechitCluster3MuonVetoIso[N_MAX_CSC];
  float         cscRechitCluster3IsoMuonVetoPt[N_MAX_CSC];
  float         cscRechitCluster3IsoMuonVetoE[N_MAX_CSC];
  float         cscRechitCluster3IsoMuonVetoPhi[N_MAX_CSC];
  float         cscRechitCluster3IsoMuonVetoEta[N_MAX_CSC];
  float         cscRechitCluster3GenMuonVetoPt[N_MAX_CSC];
  float         cscRechitCluster3GenMuonVetoE[N_MAX_CSC];
  float         cscRechitCluster3GenMuonVetoProdX[N_MAX_CSC];
  float         cscRechitCluster3GenMuonVetoProdY[N_MAX_CSC];
  float         cscRechitCluster3GenMuonVetoProdZ[N_MAX_CSC];
  int         cscRechitCluster3GenMuonVetoLLPIndex[N_MAX_CSC];
  float         cscRechitCluster3GenMuonVetoLLPDist[N_MAX_CSC];
  int           cscRechitCluster3MuonVetoType[N_MAX_CSC];


  int           cscRechitCluster3NRechitChamberPlus11[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus12[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus13[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus21[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus22[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus31[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus32[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus41[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberPlus42[N_MAX_CSC];

  int           cscRechitCluster3NRechitChamberMinus11[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus12[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus13[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus21[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus22[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus31[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus32[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus41[N_MAX_CSC];
  int           cscRechitCluster3NRechitChamberMinus42[N_MAX_CSC];
  float         cscRechitCluster3Met_dPhi[N_MAX_CSC];
  float         cscRechitCluster3MetXYCorr_dPhi[N_MAX_CSC];

  float     cscRechitCluster3MetHEM_dPhi[N_MAX_CSC];
  float     cscRechitCluster3MetHEMXYCorr_dPhi[N_MAX_CSC];
  float     cscRechitCluster3MetEENoise_dPhi[N_MAX_CSC];
  float     cscRechitCluster3MetEENoiseXYCorr_dPhi[N_MAX_CSC];
  float     cscRechitCluster3MetJesUp_dPhi[N_MAX_CSC];
  float     cscRechitCluster3MetJesDown_dPhi[N_MAX_CSC];


  //gLLP
  float gLLP_eta[1];
  float gLLP_phi[1];
  float gLLP_csc[1];
  float gLLP_dt[1];
  float gLLP_beta[1];
  float gLLP_e[1];
  float gLLP_pt[1];
  float gLLP_lepdPhi[1];
  int gLLP_multiplicity[1];



  float gLLP_ctau[1];
  float gLLP_decay_vertex_r[1];
  float gLLP_decay_vertex_x[1];
  float gLLP_decay_vertex_y[1];
  float gLLP_decay_vertex_z[1];

  float gLLP_daughter_deltaR[1];
  float gLLP_daughter_pt[3];
  int   gLLP_daughter_id[3];
  float gLLP_daughter_eta[3];
  float gLLP_daughter_phi[3];
  float gLLP_daughter_e[3];
  float gLLP_daughter_mass[3];
  float gHiggsPt;
  float gHiggsEta;
  float gHiggsPhi;
  float gHiggsE;

  //leptons
  int nMuons;
  float muonPt[N_MAX_LEPTONS];
  float muonE[N_MAX_LEPTONS];
  float muonEta[N_MAX_LEPTONS];
  float muonPhi[N_MAX_LEPTONS];
  float muonIso[N_MAX_LEPTONS];
  bool  muonIsGlobal[N_MAX_LEPTONS];
  bool  muonTightId[N_MAX_LEPTONS];
  bool  muonLooseId[N_MAX_LEPTONS];
  int nLeptons;
  float lepE[N_MAX_LEPTONS];
  float lepPt[N_MAX_LEPTONS];
  float lepEta[N_MAX_LEPTONS];
  float lepPhi[N_MAX_LEPTONS];
  int  lepPdgId[N_MAX_LEPTONS];
  float lepDZ[N_MAX_LEPTONS];
  float lepTriggerSF[N_MAX_LEPTONS];
  float lepTightIdSF[N_MAX_LEPTONS];
  float lepLooseIdSF[N_MAX_LEPTONS];
  float lepTightIsoSF[N_MAX_LEPTONS];
  float lepLooseIsoSF[N_MAX_LEPTONS];
  float lepTriggerMCEfficiency[N_MAX_LEPTONS];
  float lepTightIdMCEfficiency[N_MAX_LEPTONS];
  float lepLooseIdMCEfficiency[N_MAX_LEPTONS];
  float lepTightIsoMCEfficiency[N_MAX_LEPTONS];
  float lepLooseIsoMCEfficiency[N_MAX_LEPTONS];
  float lepEff[N_MAX_LEPTONS];
  float lepSF[N_MAX_LEPTONS];
  bool lepTag[N_MAX_LEPTONS];

  bool lepLoosePassId[N_MAX_LEPTONS];
  bool lepMediumPassId[N_MAX_LEPTONS];
  bool lepTightPassId[N_MAX_LEPTONS];
  bool lepPassVetoId[N_MAX_LEPTONS];
  bool lepFromZ[N_MAX_LEPTONS];
  bool lepPassId[N_MAX_LEPTONS];
  bool lepPassLooseIso[N_MAX_LEPTONS];
  bool lepPassTightIso[N_MAX_LEPTONS];
  bool lepPassVTightIso[N_MAX_LEPTONS];
  bool lepPassVVTightIso[N_MAX_LEPTONS];

  //Z-candidate
  float MT;
  float ZMass1;
  float ZMass;
  float ZPt;
  float ZEta;
  float ZPhi;
  int ZleptonIndex1;
  int ZleptonIndex2;
  //jets
  int nJets;
  float jetE[N_MAX_JETS];
  float jetPt[N_MAX_JETS];
  float jetEta[N_MAX_JETS];
  float jetPhi[N_MAX_JETS];
  float jetTime[N_MAX_JETS];
  float jetPtJESUp[N_MAX_JETS];
  float jetPtJESDown[N_MAX_JETS];
  float jetEJESUp[N_MAX_JETS];
  float jetEJESDown[N_MAX_JETS];
  float JecUnc[N_MAX_JETS];
  float jetCISV[N_MAX_JETS]; 
  float jetCMVA[N_MAX_JETS]; 


  float ecalNRechits[N_MAX_JETS];
  float ecalRechitE[N_MAX_JETS];
  float jetElectronEnergyFraction[N_MAX_JETS];
  float jetPhotonEnergyFraction[N_MAX_JETS];
  float jetChargedHadronEnergyFraction[N_MAX_JETS];
  float jetNeutralHadronEnergyFraction[N_MAX_JETS];

  float jetMuonEnergyFraction[N_MAX_JETS];


  bool jetPassMuFrac[N_MAX_JETS];
  float jet_match_genJet_minDeltaR[N_MAX_JETS];
  int jet_match_genJet_index[N_MAX_JETS];
  float jet_match_genJet_pt[N_MAX_JETS];

  // bool jetLoosePassId[N_MAX_JETS];
  bool jetPassId[N_MAX_JETS];
  bool jetTightPassId[N_MAX_JETS];
  bool HLTDecision[NTriggersMAX];
  bool SingleMuonTrigger;
  bool SingleEleTrigger;
  bool SingleLepTrigger;
  UInt_t wzevtNum,trig, trig_lepId, trig_lepId_dijet; //number of events that pass each criteria




  void InitVariables();
  void InitTree();
  void LoadTree(const char* file);
  void CreateTree();

};
#endif
