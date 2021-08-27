#include "RazorHelper.h"
#include "HNLMuonSystemTree.h"
#include "assert.h"
#include "TTree.h"
#include "DBSCAN.h"

// Constructor
HNLMuonSystemTree::HNLMuonSystemTree()
{
  InitVariables();
};
HNLMuonSystemTree::~HNLMuonSystemTree()
{
  if (f_) f_->Close();
};
void HNLMuonSystemTree::InitVariables()
{
  runNum=0; lumiSec=0; evtNum=0; category=0;
  npv=0; npu=0;
  pileupWeight = 0; pileupWeightUp = 0; pileupWeightDown = 0;
  lepOverallSF = 1.0;
  metSF = 1.0;

  ZCategory = 3;
  sf_facScaleUp = 0; sf_facScaleDown = 0; sf_renScaleUp = 0; sf_renScaleDown = 0; sf_facRenScaleUp = 0; sf_facRenScaleDown = 0;

  for (int i = 0; i < 9; i++)
  {
    scaleWeights[i] = 0.0;
  }
  weight=-1.0;rho=-1;
  met=-1; metPhi=-1; metXYCorr=-1; metPhiXYCorr=-1; HT = 0.0; jetMet_dPhi = -999.;jetMet_dPhiMin = 999.;jetMet_dPhiMin4 = 999.;
  metNoMu = -1.;metPhiJESUp = -999.; metPhiJESDown = -999.;  metHEM=-999.; metPhiHEM = -999.;metHEMXYCorr=-999.; metPhiHEMXYCorr = -999.;
  EE_prefiring = true;
  Flag_HBHENoiseFilter = false; Flag_HBHEIsoNoiseFilter = false; Flag_BadPFMuonFilter = false; Flag_CSCTightHaloFilter = false; Flag_goodVertices = false;
  Flag_ecalBadCalibFilter = false; Flag_all = false; Flag_globalSuperTightHalo2016Filter = false; Flag_BadChargedCandidateFilter = false; Flag_eeBadScFilter = false;

  Flag2_HBHENoiseFilter = false; Flag2_HBHEIsoNoiseFilter = false; Flag2_BadPFMuonFilter = false; Flag2_globalSuperTightHalo2016Filter = false;
  Flag2_globalTightHalo2016Filter = false; Flag2_BadChargedCandidateFilter = false; Flag2_EcalDeadCellTriggerPrimitiveFilter = false; Flag2_ecalBadCalibFilter = false;
  Flag2_eeBadScFilter = false;
  Flag2_all = false;
  mH = 0; mX = 0; ctau = 0;



  metJESUp = -999.;metJESDown = -999.; metEENoise = -999.;metPhiEENoise = -999.;metEENoiseXYCorr = -999.;metPhiEENoiseXYCorr = -999.;
  gWPt = 0.0;
  gLepId = 0;
  gLepPt = 0.; gLepPhi = 0.; gLepEta = 0.; gLepE = 0.;
  gHiggsPt = 0.; gHiggsPhi = 0.; gHiggsEta = 0.; gHiggsE = 0.;
  //CSC
  nCscRechitClusters3 = 0;
  nCscRechits = 0;
  nEarlyCscRechits = 0;
  nLateCscRechits = 0;
  nEarly2CscRechits = 0;
  nLate2CscRechits = 0;
  nCscPositiveYRechits = 0;
  nCscNegativeYRechits = 0;
  cscPosTpeak = 0.0;
  cscNegTpeak = 0.0;
  nDtRings = 0;
  nCscRings = 0;


  nCscRechitsChamberPlus11 = 0;
  nCscRechitsChamberPlus12 = 0;
  nCscRechitsChamberPlus13 = 0;
  nCscRechitsChamberPlus21 = 0;
  nCscRechitsChamberPlus22 = 0;
  nCscRechitsChamberPlus31 = 0;
  nCscRechitsChamberPlus32 = 0;
  nCscRechitsChamberPlus41 = 0;
  nCscRechitsChamberPlus42 = 0;
  nCscRechitsChamberMinus11 = 0;
  nCscRechitsChamberMinus12 = 0;
  nCscRechitsChamberMinus13 = 0;
  nCscRechitsChamberMinus21 = 0;
  nCscRechitsChamberMinus22 = 0;
  nCscRechitsChamberMinus31 = 0;
  nCscRechitsChamberMinus32 = 0;
  nCscRechitsChamberMinus41 = 0;
  nCscRechitsChamberMinus42 = 0;

  nDTRechits = 0;
  nDTNegativeYRechits  = 0;
  nDTPositiveYRechits = 0;
  nDTRechitsChamberMinus12 = 0;
  nDTRechitsChamberMinus11 = 0;
  nDTRechitsChamber10 = 0;
  nDTRechitsChamberPlus11 = 0;
  nDTRechitsChamberPlus12 = 0;
  nDTRechitsChamberMinus22 = 0;
  nDTRechitsChamberMinus21 = 0;
  nDTRechitsChamber20 = 0;
  nDTRechitsChamberPlus21 = 0;
  nDTRechitsChamberPlus22 = 0;
  nDTRechitsChamberMinus32 = 0;
  nDTRechitsChamberMinus31 = 0;
  nDTRechitsChamber30 = 0;
  nDTRechitsChamberPlus31 = 0;
  nDTRechitsChamberPlus32 = 0;
  nDTRechitsChamberMinus42 = 0;
  nDTRechitsChamberMinus41 = 0;
  nDTRechitsChamber40 = 0;
  nDTRechitsChamberPlus41 = 0;
  nDTRechitsChamberPlus42 = 0;
  nDtStations25 = 0;
  nDtWheels25 = 0;
  nDTRechitsStation1 = 0;
  nDTRechitsStation2 = 0;
  nDTRechitsStation3 = 0;
  nDTRechitsStation4 = 0;

  nDTRechitsWheelMinus2 = 0;
  nDTRechitsWheelMinus1 = 0;
  nDTRechitsWheel0 = 0;
  nDTRechitsWheelPlus1 = 0;
  nDTRechitsWheelPlus2 = 0;

  for( int i = 0; i < N_MAX_CSC; i++ )
  {
      cscRechitCluster3_match_Me1112_0p4[i] = 0;
      cscRechitCluster3_match_Me1112_0p6[i] = 0;
      cscRechitCluster3_match_Me1112_0p8[i] = 0;
      cscRechitCluster3_match_Me11_0p4[i] = 0;
      cscRechitCluster3_match_Me11_0p6[i] = 0;
      cscRechitCluster3_match_Me11_0p8[i] = 0;
      cscRechitCluster3_match_Me12_0p4[i] = 0;
      cscRechitCluster3_match_Me12_0p6[i] = 0;
      cscRechitCluster3_match_Me12_0p8[i] = 0;

      cscRechitCluster3_match_cscRechits_0p4[i] = 0;

      cscRechitCluster3_match_cscSeg_0p4[i] = 0;
      cscRechitCluster3_match_ME11Seg_0p4[i] = 0;
      cscRechitCluster3_match_ME12Seg_0p4[i] = 0;
      cscRechitCluster3_match_cscSeg_0p6[i] = 0;
      cscRechitCluster3_match_ME11Seg_0p6[i] = 0;
      cscRechitCluster3_match_ME12Seg_0p6[i] = 0;

      cscRechitCluster3_match_dtRechits_0p4[i] = 0;
      cscRechitCluster3_match_dtRechits_0p6[i] = 0;
      cscRechitCluster3_match_dtRechits_phi0p2[i] = 0;
      cscRechitCluster3_match_MB1_0p4[i] = 0;
      cscRechitCluster3_match_MB1_0p6[i] = 0;
      cscRechitCluster3_match_dtSeg_0p4[i] = 0;
      cscRechitCluster3_match_dtSeg_0p6[i] = 0;
      cscRechitCluster3_match_MB1Seg_0p4[i] = 0;
      cscRechitCluster3_match_MB1Seg_0p6[i] = 0;
      cscRechitCluster3_match_RB1_0p4[i] = 0;
      cscRechitCluster3_match_RE12_0p4[i] = 0;
      cscRechitCluster3_match_RB1_0p6[i] = 0;
      cscRechitCluster3_match_RE12_0p6[i] = 0;
      cscRechitCluster3_match_highEta_0p4[i] = 0;
      cscRechitCluster3_match_highEta_0p6[i] = 0;
      cscRechitCluster3_match_highEta_0p8[i] = 0;
      cscRechitCluster3_match_cluster_dR[i] = 999.;
      cscRechitCluster3_match_cluster_index[i] = 999;


      cscRechitCluster3_match_gParticle[i] = false;
      cscRechitCluster3_match_gParticle_minDeltaR[i] = -999.;
      cscRechitCluster3_match_gParticle_index[i] = -999;
      cscRechitCluster3_match_gParticle_id[i] = -999;
      cscRechitCluster3_match_gParticle_eta[i] = -999.;
      cscRechitCluster3_match_gParticle_phi[i] = -999.;
      cscRechitCluster3_match_gParticle_E[i] = -999.;
      cscRechitCluster3_match_gParticle_pt[i] = -999.;
      cscRechitCluster3_match_gParticle_MotherId[i]  = -999;

        cscRechitCluster3_match_gLLP[i] = false;
        cscRechitCluster3_match_gLLP_minDeltaR[i] = 999;
        cscRechitCluster3_match_gLLP_index[i] = 999;
        cscRechitCluster3_match_gLLP_eta[i] = 999.;
        cscRechitCluster3_match_gLLP_phi[i] = 999.;
        cscRechitCluster3_match_gLLP_decay_r[i] = 999.;
        cscRechitCluster3_match_gLLP_decay_x[i] = 999.;
        cscRechitCluster3_match_gLLP_decay_y[i] = 999.;
        cscRechitCluster3_match_gLLP_decay_z[i] = 999.;
        cscRechitCluster3_match_gLLP_ctau[i] = 999.;
        cscRechitCluster3_match_gLLP_beta[i] = 999.;
        cscRechitCluster3_match_gLLP_csc[i] = false;
        cscRechitCluster3_match_gLLP_e[i] = 999.;
        cscRechitCluster3_match_gLLP_pt[i] = 999.;
        cscRechitCluster3_match_gLLP_lepdPhi[i] = 999.;
        cscRechitCluster3_match_gLLP_daughter0_deltaR[i] = 999.0;
        cscRechitCluster3_match_gLLP_daughter1_deltaR[i] = 999.0;
        cscRechitCluster3_match_gLLP_daughter2_deltaR[i] = 999.0;
        cscRechitCluster3_match_gLLP_daughter3_deltaR[i] = 999.0;

        cscRechitCluster3Size[i] = -999;
        cscRechitCluster3X[i] = -999.;
        cscRechitCluster3Y[i] = -999.;
        cscRechitCluster3Z[i] = -999.;
        cscRechitCluster3Time[i] = -999.;
        cscRechitCluster3TimeTotal[i] = -999.;
        cscRechitCluster3TimeWire[i] = -999.;
        cscRechitCluster3TimeWirePruned[i] = -999.;


        cscRechitCluster3GenMuonDeltaR[i] = 999.;
        cscRechitCluster3TimeSpread[i] = -999.;
        cscRechitCluster3TimeTotalSpread[i] = -999.;
        cscRechitCluster3TimeTotalSpreadPruned[i] = -999.;
        cscRechitCluster3TimeWireSpread[i] = -999.;

        cscRechitCluster3MajorAxis[i] = -999.;
        cscRechitCluster3MinorAxis[i] = -999.;
        cscRechitCluster3XSpread[i] = -999.;
        cscRechitCluster3YSpread[i] = -999.;
        cscRechitCluster3ZSpread[i] = -999.;
        cscRechitCluster3EtaPhiSpread[i] = -999.;
        cscRechitCluster3XYSpread[i] = -999.;
        cscRechitCluster3RSpread[i] = -999.;
        cscRechitCluster3DeltaRSpread[i] =999.;

        cscRechitCluster3EtaSpread[i] = -999.;
        cscRechitCluster3PhiSpread[i] = -999.;
        cscRechitCluster3Eta[i] = -999.;
        cscRechitCluster3Phi[i] = -999.;
        cscRechitCluster3JetVetoPt[i] = 0.0;
        cscRechitCluster3JetVetoEta[i] = 0.0;
        cscRechitCluster3JetVetoPhi[i] = 0.0;

        cscRechitCluster3JetVetoElectronEnergyFraction[i] = 0.0;
        cscRechitCluster3JetVetoPhotonEnergyFraction[i] = 0.0;
        cscRechitCluster3JetVetoChargedHadronEnergyFraction[i] = 0.0;
        cscRechitCluster3JetVetoNeutralHadronEnergyFraction[i] = 0.0;
        cscRechitCluster3JetVetoMuonEnergyFraction[i] = 0.0;
        cscRechitCluster3JetVetoE[i] = 0.0;
        cscRechitCluster3GenJetVetoPt[i] = 0.0;
        cscRechitCluster3GenJetVetoE[i] = 0.0;
        cscRechitCluster3MuonVetoPt[i] = 0.0;
        cscRechitCluster3MuonVetoE[i] = 0.0;
        cscRechitCluster3MuonVetoPhi[i] = 0.0;
        cscRechitCluster3MuonVetoEta[i] = 0.0;
        cscRechitCluster3MuonVetoLooseIso[i] = false;
        cscRechitCluster3MuonVetoTightIso[i] = false;
        cscRechitCluster3MuonVetoVTightIso[i] = false;
        cscRechitCluster3MuonVetoVVTightIso[i] = false;
        cscRechitCluster3MuonVetoTightId[i] = false;
        cscRechitCluster3MuonVetoLooseId[i] = false;
        cscRechitCluster3MuonVetoIso[i] = false;
        cscRechitCluster3IsoMuonVetoPt[i] = 0.0;
        cscRechitCluster3IsoMuonVetoE[i] = 0.0;
        cscRechitCluster3IsoMuonVetoPhi[i] = 0.0;
        cscRechitCluster3IsoMuonVetoEta[i] = 0.0;
        cscRechitCluster3GenMuonVetoE[i] = 0.0;
        cscRechitCluster3GenMuonVetoPt[i] = 0.0;
        cscRechitCluster3MuonVetoType[i] = 999;

        cscRechitCluster3GenMuonVetoProdX[i] = 0.0;
        cscRechitCluster3GenMuonVetoProdY[i] = 0.0;
        cscRechitCluster3GenMuonVetoProdZ[i] = 0.0;
        cscRechitCluster3GenMuonVetoLLPDist[i] = 999.;
        cscRechitCluster3GenMuonVetoLLPIndex[i] = 999;

        cscRechitCluster3JetVetoPt_0p6[i] = 0.0;
        cscRechitCluster3JetVetoPt_0p8[i] = 0.0;
        cscRechitCluster3JetVetoE_0p6[i] = 0.0;
        cscRechitCluster3JetVetoE_0p8[i] = 0.0;
        cscRechitCluster3MuonVetoPt_0p6[i] = 0.0;
        cscRechitCluster3MuonVetoPt_0p8[i] = 0.0;
        cscRechitCluster3MuonVetoE_0p6[i] = 0.0;
        cscRechitCluster3MuonVetoE_0p8[i] = 0.0;

        cscRechitCluster3ZLep1[i] = false;
        cscRechitCluster3ZLep2[i] = false;
        cscRechitCluster3ZLep2Tag[i] = false;
        cscRechitCluster3ZLep1Tag[i] = false;
        cscRechitCluster3ZLep1Id[i] = -999;
        cscRechitCluster3ZLep2Id[i] = -999;
        cscRechitCluster3ZLep1LooseIso[i] = false;
        cscRechitCluster3ZLep1TightIso[i] = false;
        cscRechitCluster3ZLep1VTightIso[i] = false;
        cscRechitCluster3ZLep1VVTightIso[i] = false;
        cscRechitCluster3ZLep1TightId[i] = false;
        cscRechitCluster3ZLep2LooseIso[i] = false;
        cscRechitCluster3ZLep2TightIso[i] = false;
        cscRechitCluster3ZLep2VTightIso[i] = false;
        cscRechitCluster3ZLep2VVTightIso[i] = false;
        cscRechitCluster3ZLep2TightId[i] = false;
        cscRechitCluster3NChamber[i] = -999;
        cscRechitCluster3MaxChamberRatio[i] = -999.;
        cscRechitCluster3MaxChamber[i] = -999;
        cscRechitCluster3NStation[i] = -999;
        cscRechitCluster3NStation5[i] = -999;
        cscRechitCluster3NStation10[i] = -999;
        cscRechitCluster3NStation10perc[i] = -999;
        cscRechitCluster3AvgStation[i] = -999.;
        cscRechitCluster3AvgStation5[i] = -999.;
        cscRechitCluster3AvgStation10[i] = -999.;
        cscRechitCluster3AvgStation10perc[i] = -999.;
        cscRechitCluster3MaxStationRatio[i] = -999.;
        cscRechitCluster3MaxStation[i] = -999;
        cscRechitCluster3Me11Ratio[i] = -999.;
        cscRechitCluster3Me12Ratio[i] = -999.;
        cscRechitCluster3NRechitChamberPlus11[i] = -999;
        cscRechitCluster3NRechitChamberPlus12[i] = -999;
        cscRechitCluster3NRechitChamberPlus13[i] = -999;
        cscRechitCluster3NRechitChamberPlus21[i] = -999;
        cscRechitCluster3NRechitChamberPlus22[i] = -999;
        cscRechitCluster3NRechitChamberPlus31[i] = -999;
        cscRechitCluster3NRechitChamberPlus32[i] = -999;
        cscRechitCluster3NRechitChamberPlus41[i] = -999;
        cscRechitCluster3NRechitChamberPlus42[i] = -999;
        cscRechitCluster3NRechitChamberMinus11[i] = -999;
        cscRechitCluster3NRechitChamberMinus12[i] = -999;
        cscRechitCluster3NRechitChamberMinus13[i] = -999;
        cscRechitCluster3NRechitChamberMinus21[i] = -999;
        cscRechitCluster3NRechitChamberMinus22[i] = -999;
        cscRechitCluster3NRechitChamberMinus31[i] = -999;
        cscRechitCluster3NRechitChamberMinus32[i] = -999;
        cscRechitCluster3NRechitChamberMinus41[i] = -999;
        cscRechitCluster3NRechitChamberMinus42[i] = -999;
        cscRechitCluster3Met_dPhi[i] = 999.;
        cscRechitCluster3MetXYCorr_dPhi[i] = 999.;

        cscRechitCluster3MetHEM_dPhi[i] = 999.;
        cscRechitCluster3MetHEMXYCorr_dPhi[i] = 999.;
        cscRechitCluster3MetEENoise_dPhi[i] = 999.;
        cscRechitCluster3MetEENoiseXYCorr_dPhi[i] = 999.;
        cscRechitCluster3MetJesUp_dPhi[i] = 999.;
        cscRechitCluster3MetJesDown_dPhi[i] = 999.;


        dtRechitCluster_match_gParticle_deltaR[i] = 999.;
        dtRechitCluster_match_gParticle_Id[i] = -999;
        dtRechitCluster_match_gParticle_Pt[i] = -999.;
        dtRechitCluster_match_gParticle_Eta[i] = -999.;
        dtRechitCluster_match_gParticle_Phi[i] = -999.;
        dtRechitCluster_match_gParticle_E[i] = -999.;
        dtRechitCluster_match_gParticle_Status[i] = -999;
        dtRechitCluster_match_gParticle_MotherId[i] = -999;
        dtRechitCluster_match_gParticle_deltaR[i] = -999.;

        dtRechitCluster_match_MB1hits_0p4[i] = 0;
        dtRechitCluster_match_MB1hits_0p5[i] = 0;
        dtRechitCluster_match_MB1hits_cosmics_plus[i] = 0;
        dtRechitCluster_match_MB1hits_cosmics_minus[i] = 0;
        dtRechitCluster_match_RPChits_dPhi0p5[i] = 0;
        dtRechitCluster_match_RPCBx_dPhi0p5[i] = 0;
        dtRechitCluster_match_RB1_0p4[i] = 0;
        dtRechitCluster_match_RB1_dPhi0p5[i] = 0;



        dtRechitCluster_match_gLLP[i] = false;
        dtRechitCluster_match_gLLP_minDeltaR[i] = 999;
        dtRechitCluster_match_gLLP_eta[i] = 999.;
        dtRechitCluster_match_gLLP_phi[i] = 999.;
        dtRechitCluster_match_gLLP_decay_x[i] = 999.;
        dtRechitCluster_match_gLLP_decay_y[i] = 999.;
        dtRechitCluster_match_gLLP_decay_z[i] = 999.;
        dtRechitCluster_match_gLLP_ctau[i] = 999.;
        dtRechitCluster_match_gLLP_beta[i] = 999.;
        dtRechitCluster_match_gLLP_csc[i] = false;
        dtRechitCluster_match_gLLP_dt[i] = false;



        dtRechitCluster_match_gLLP_e[i] = 999.;
        dtRechitCluster_match_gLLP_pt[i] = 999.;
        dtRechitClusterLep_dPhi[i] = 999.;



        dtRechitClusterSize[i] = -999;
        dtRechitClusterX[i] = -999.;
        dtRechitClusterY[i] = -999.;
        dtRechitClusterZ[i] = -999.;

        dtRechitClusterWheel[i] = -999;

        dtRechitClusterEta[i] = -999.;
        dtRechitClusterPhi[i] = -999.;

        dtRechitClusterJetVetoPt[i] = 0.0;
        dtRechitClusterJetVetoEta[i] = 0.0;
        dtRechitClusterJetVetoPhi[i] = 0.0;


        dtRechitClusterMuonVetoPt[i] = 0.0;
        dtRechitClusterMuonVetoE[i] = 0.0;
        dtRechitClusterMuonVetoPhi[i] = 0.0;
        dtRechitClusterMuonVetoEta[i] = 0.0;
        dtRechitClusterMuonVetoLooseId[i] = false;
        dtRechitClusterMuonVetoGlobal[i] = false;



        dtRechitClusterNChamber[i] = -999;
        dtRechitClusterMaxChamberRatio[i] = -999.;
        dtRechitClusterMaxChamber[i] = -999;
        dtRechitClusterNStation10[i] = -999;
        dtRechitClusterAvgStation10[i] = -999.;
        dtRechitClusterMaxStationRatio[i] = -999.;
        dtRechitClusterMaxStation[i] = -999;

        dtRechitClusterNSegmentStation1[i] = -999;
        dtRechitClusterNSegmentStation2[i] = -999;
        dtRechitClusterNSegmentStation3[i] = -999;
        dtRechitClusterNSegmentStation4[i] = -999;


        dtRechitClusterMetEENoise_dPhi[i] = 999.;




  }

  for(int i = 0;i<2;i++)
  {
    gLLP_multiplicity[i]= 0;
    gLLP_eta[i] = 0.0;
    gLLP_phi[i] = 0.0;
    gLLP_beta[i] = 0.0;
    gLLP_e[i] = 0.0;
    gLLP_pt[i] = 0.0;
    gLLP_lepdPhi[i] = 0.0;
    gLLP_csc[i] = 0.0;
    gLLP_dt[i] = 0.0;
    gLLP_ctau[i] = 0.0;
    gLLP_decay_vertex_r[i] = 0.0;
    gLLP_decay_vertex_x[i] = 0.0;
    gLLP_decay_vertex_y[i] = 0.0;
    gLLP_decay_vertex_z[i] = 0.0;
    gLLP_daughter_deltaR[i] = -999.0;

  }
  for(int i = 0;i<4;i++)
  {
    gLLP_daughter_pt[i] = -999.0;
    gLLP_daughter_eta[i] = -999.0;
    gLLP_daughter_phi[i] = -999.0;
    gLLP_daughter_e[i] = -999.0;
    gLLP_daughter_mass[i] = -999.0;
    gLLP_daughter_id[i] = 999;
  }


  genMetPtTrue = -999.;
  genMetPhiTrue = -999.;
  genMetPtCalo = -999.;
  genMetPhiCalo = -999.;
  // nGenParticle = 0;
  // nGenJets = 0;
  // for( int i = 0; i < N_MAX_GPARTICLES; i++ )
  // {
  //   gParticleId[i] = 0;
  //   gParticleStatus[i] = 999;
  //   gParticleMotherId[i] = 0;
  //   gParticlePt[i] = -999.;
  //   gParticleEta[i] = -999.;
  //   gParticlePhi[i] = -999.;
  //   gParticleE[i] = -999.;
  //   genJetE[i] = -999.;
  //   genJetPt[i] = -999.;
  //   genJetEta[i] = -999.;
  //   genJetPhi[i] = -999.;
  //   genJetMET[i] = -999.;
  //
  // }

  //leptons

  nMuons = 0;
  for( int i = 0; i < N_MAX_LEPTONS; i++ )
  {
    muonPt[i]     = -999.;
    muonEta[i]    = -999.;
    muonPhi[i]    = -999.;
  }
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
    lepPassVetoId[i] = false;
    lepPassId[i] = false;
    lepFromZ[i] = false;
    lepPassLooseIso[i] = false;
    lepPassTightIso[i] = false;
    lepPassVTightIso[i] = false;
    lepPassVTightIso[i] = false;
    lepEff[i] = 1.0;
    lepSF[i] = 1.0;
    lepTriggerSF[i] = 1.0;
    lepTightIdSF[i] = 1.0;
    lepLooseIdSF[i] = 1.0;
    lepTightIsoSF[i] = 1.0;
    lepLooseIsoSF[i] = 1.0;
    lepTriggerMCEfficiency[i] = 1.0;
    lepTightIdMCEfficiency[i] = 1.0;
    lepLooseIdMCEfficiency[i] = 1.0;
    lepTightIsoMCEfficiency[i] = 1.0;
    lepLooseIsoMCEfficiency[i] = 1.0;
    lepTag[i] = false;

  }
  //Z-candidate
  ZMass1 = -999.; ZMass = -999.; ZPt = -999.; ZEta = -999.; ZPhi = -999.;
  MT = -999.;
  ZleptonIndex1 = -999; ZleptonIndex2 = -999;
  //jets
  nJets = 0;
  for( int i = 0; i < N_MAX_JETS; i++ )
  {
    jetE[i]      = -999.;
    jetPt[i]     = -999.;
    jetEta[i]    = -999.;
    jetPhi[i]    = -999.;
    jetTime[i]   = -999.;
    // jetLoosePassId[i] = false;
    jetPassId[i] = false;
    jetPtJESUp[i] = -999.;
    jetPtJESDown[i] = -999.;
    jetEJESUp[i] = -999.;
    jetEJESDown[i] = -999.;
    JecUnc[i] = -999.;
    ecalNRechits[i] = -999.;
    ecalRechitE[i] = -999.;

    jetElectronEnergyFraction[i] = -999.;
    jetPhotonEnergyFraction[i] = -999.;
    jetChargedHadronEnergyFraction[i] = -999.;
    jetNeutralHadronEnergyFraction[i] = -999.;
    jetMuonEnergyFraction[i] = -999.;
    jet_match_genJet_minDeltaR[i] = -999.;
    jet_match_genJet_index[i] = -999;
    jet_match_genJet_pt[i] = -999.;
    jetTightPassId[i] = false;
  }

  for(int i = 0; i <NTriggersMAX; i++){
    HLTDecision[i] = false;
  }
  SingleMuonTrigger = false;
  SingleEleTrigger = false;
  SingleLepTrigger = false;
};

void HNLMuonSystemTree::InitTree()
{
  assert(tree_);
  InitVariables();

  tree_->SetBranchAddress("runNum",      &runNum);
  tree_->SetBranchAddress("lumiSec",     &lumiSec);
  tree_->SetBranchAddress("evtNum",      &evtNum);
  tree_->SetBranchAddress("category",    &category);
  tree_->SetBranchAddress("mX",      &mX);
  tree_->SetBranchAddress("mH",      &mH);
  tree_->SetBranchAddress("ctau",      &ctau);
  tree_->SetBranchAddress("ZCategory",      &ZCategory);

  tree_->SetBranchAddress("npv",         &npv);
  tree_->SetBranchAddress("npu",         &npu);
  tree_->SetBranchAddress("weight",      &weight);
  tree_->SetBranchAddress("scaleWeights",      &scaleWeights);


  tree_->SetBranchAddress("lepOverallSF",      &lepOverallSF);


  tree_->SetBranchAddress("sf_facScaleUp",      &sf_facScaleUp);
  tree_->SetBranchAddress("sf_facScaleDown",      &sf_facScaleDown);
  tree_->SetBranchAddress("sf_renScaleUp",      &sf_renScaleUp);
  tree_->SetBranchAddress("sf_renScaleDown",      &sf_renScaleDown);
  tree_->SetBranchAddress("sf_facRenScaleUp",      &sf_facRenScaleUp);
  tree_->SetBranchAddress("sf_facRenScaleDown",      &sf_facRenScaleDown);
  tree_->SetBranchAddress("metSF",      &metSF);

  tree_->SetBranchAddress("pileupWeight",      &pileupWeight);
  tree_->SetBranchAddress("pileupWeightUp",      &pileupWeightUp);
  tree_->SetBranchAddress("pileupWeightDown",      &pileupWeightDown);
  tree_->SetBranchAddress("Flag_HBHENoiseFilter",      &Flag_HBHENoiseFilter);
  tree_->SetBranchAddress("Flag_HBHEIsoNoiseFilter",      &Flag_HBHEIsoNoiseFilter);
  tree_->SetBranchAddress("Flag_BadPFMuonFilter",      &Flag_BadPFMuonFilter);
  tree_->SetBranchAddress("Flag_CSCTightHaloFilter",      &Flag_CSCTightHaloFilter);
  tree_->SetBranchAddress("Flag_BadChargedCandidateFilter",      &Flag_BadChargedCandidateFilter);
  tree_->SetBranchAddress("Flag_eeBadScFilter",      &Flag_eeBadScFilter);
  tree_->SetBranchAddress("Flag_globalSuperTightHalo2016Filter",      &Flag_globalSuperTightHalo2016Filter);
  tree_->SetBranchAddress("Flag_goodVertices",      &Flag_goodVertices);
  tree_->SetBranchAddress("Flag_ecalBadCalibFilter",      &Flag_ecalBadCalibFilter);
  // tree_->SetBranchAddress("Flag_all",      &Flag_all);

  tree_->SetBranchAddress("Flag2_HBHENoiseFilter",      &Flag2_HBHENoiseFilter);
  tree_->SetBranchAddress("Flag2_HBHEIsoNoiseFilter",      &Flag2_HBHEIsoNoiseFilter);
  tree_->SetBranchAddress("Flag2_BadPFMuonFilter",      &Flag2_BadPFMuonFilter);
  tree_->SetBranchAddress("Flag2_globalSuperTightHalo2016Filter",      &Flag2_globalSuperTightHalo2016Filter);
  tree_->SetBranchAddress("Flag2_globalTightHalo2016Filter",      &Flag2_globalTightHalo2016Filter);
  tree_->SetBranchAddress("Flag2_BadChargedCandidateFilter",      &Flag2_BadChargedCandidateFilter);
  tree_->SetBranchAddress("Flag2_EcalDeadCellTriggerPrimitiveFilter",      &Flag2_EcalDeadCellTriggerPrimitiveFilter);
  tree_->SetBranchAddress("Flag2_ecalBadCalibFilter",      &Flag2_ecalBadCalibFilter);
  tree_->SetBranchAddress("Flag2_eeBadScFilter",      &Flag2_eeBadScFilter);
  tree_->SetBranchAddress("Flag2_all",      &Flag2_all);
  tree_->SetBranchAddress("EE_prefiring",      &EE_prefiring);


  tree_->SetBranchAddress("rho",         &rho);
  tree_->SetBranchAddress("met",         &met);
  tree_->SetBranchAddress("HT",         &HT);
  tree_->SetBranchAddress("metNoMu",         &metNoMu);


  tree_->SetBranchAddress("metPhi",      &metPhi);
  tree_->SetBranchAddress("metXYCorr",      &metXYCorr);
  tree_->SetBranchAddress("metPhiXYCorr",      &metPhiXYCorr);

  tree_->SetBranchAddress("jetMet_dPhi",      &jetMet_dPhi);
  tree_->SetBranchAddress("jetMet_dPhiMin",      &jetMet_dPhiMin);
  tree_->SetBranchAddress("jetMet_dPhiMin4",      &jetMet_dPhiMin4);

  tree_->SetBranchAddress("metJESUp",      &metJESUp);
  tree_->SetBranchAddress("metJESDown",      &metJESDown);
  tree_->SetBranchAddress("metPhiJESUp",      &metPhiJESUp);
  tree_->SetBranchAddress("metPhiJESDown",      &metPhiJESDown);
  tree_->SetBranchAddress("metJESUpSF",      &metJESUpSF);
  tree_->SetBranchAddress("metJESDownSF",      &metJESDownSF);
  tree_->SetBranchAddress("metEENoise",      &metEENoise);
  tree_->SetBranchAddress("metPhiEENoise",      &metPhiEENoise);
  tree_->SetBranchAddress("metEENoiseXYCorr",      &metEENoiseXYCorr);
tree_->SetBranchAddress("metPhiEENoiseXYCorr",      &metPhiEENoiseXYCorr);

  tree_->SetBranchAddress("metHEM",      &metHEM);
  tree_->SetBranchAddress("metPhiHEM",      &metPhiHEM);

    tree_->SetBranchAddress("metHEMXYCorr",      &metHEMXYCorr);
    tree_->SetBranchAddress("metPhiHEMXYCorr",      &metPhiHEMXYCorr);


  tree_->SetBranchAddress("genMetPtTrue",         &genMetPtTrue);
  tree_->SetBranchAddress("genMetPhiTrue",      &genMetPhiTrue);
  tree_->SetBranchAddress("genMetPtCalo",      &genMetPtCalo);
  tree_->SetBranchAddress("genMetPhiCalo",      &genMetPhiCalo);

  // tree_->SetBranchAddress("nGenParticle",      &nGenParticle);
  // tree_->SetBranchAddress("gParticleId",      &gParticleId);
  // tree_->SetBranchAddress("gParticleStatus",      &gParticleStatus);
  // tree_->SetBranchAddress("gParticleMotherId",      &gParticleMotherId);
  // tree_->SetBranchAddress("gParticleE",      &gParticleE);
  // tree_->SetBranchAddress("gParticlePt",      &gParticlePt);
  // tree_->SetBranchAddress("gParticleEta",      &gParticleEta);
  // tree_->SetBranchAddress("gParticlePhi",      &gParticlePhi);
  //
  // tree_->SetBranchAddress("nGenJets",      &nGenJets);
  // tree_->SetBranchAddress("genJetE",      &genJetE);
  // tree_->SetBranchAddress("genJetPt",      &genJetPt);
  // tree_->SetBranchAddress("genJetEta",      &genJetEta);
  // tree_->SetBranchAddress("genJetPhi",      &genJetPhi);
  // tree_->SetBranchAddress("genJetMET",      &genJetMET);

  tree_->SetBranchAddress("gWPt",      &gWPt);

  tree_->SetBranchAddress("gLepId",      &gLepId);
  tree_->SetBranchAddress("gLepPt",      &gLepPt);
  tree_->SetBranchAddress("gLepPhi",      &gLepPhi);
  tree_->SetBranchAddress("gLepE",      &gLepE);
  tree_->SetBranchAddress("gLepEta",      &gLepEta);
  tree_->SetBranchAddress("gHiggsPt",      &gHiggsPt);
  tree_->SetBranchAddress("gHiggsPhi",      &gHiggsPhi);
  tree_->SetBranchAddress("gHiggsE",      &gHiggsE);
  tree_->SetBranchAddress("gHiggsEta",      &gHiggsEta);
  //CSC
  tree_->SetBranchAddress("nCscRechits",             &nCscRechits);
  tree_->SetBranchAddress("nCscPositiveYRechits",             &nCscPositiveYRechits);
  tree_->SetBranchAddress("nCscNegativeYRechits",             &nCscNegativeYRechits);
  tree_->SetBranchAddress("cscPosTpeak",             &cscPosTpeak);
  tree_->SetBranchAddress("cscNegTpeak",             &cscNegTpeak);


  tree_->SetBranchAddress("nEarlyCscRechits",             &nEarlyCscRechits);
  tree_->SetBranchAddress("nLateCscRechits",             &nLateCscRechits);
  tree_->SetBranchAddress("nEarly2CscRechits",             &nEarly2CscRechits);
  tree_->SetBranchAddress("nLate2CscRechits",             &nLate2CscRechits);
  tree_->SetBranchAddress("nLate2CscRechits",             &nLate2CscRechits);
  tree_->SetBranchAddress("nCscRings",             &nCscRings);
  tree_->SetBranchAddress("nCscRechitsChamberPlus11",           &nCscRechitsChamberPlus11);
  tree_->SetBranchAddress("nCscRechitsChamberPlus12",           &nCscRechitsChamberPlus12);
  tree_->SetBranchAddress("nCscRechitsChamberPlus13",           &nCscRechitsChamberPlus13);
  tree_->SetBranchAddress("nCscRechitsChamberPlus21",           &nCscRechitsChamberPlus21);
  tree_->SetBranchAddress("nCscRechitsChamberPlus22",           &nCscRechitsChamberPlus22);
  tree_->SetBranchAddress("nCscRechitsChamberPlus31",           &nCscRechitsChamberPlus31);
  tree_->SetBranchAddress("nCscRechitsChamberPlus32",           &nCscRechitsChamberPlus32);
  tree_->SetBranchAddress("nCscRechitsChamberPlus41",           &nCscRechitsChamberPlus41);
  tree_->SetBranchAddress("nCscRechitsChamberPlus42",           &nCscRechitsChamberPlus42);

  tree_->SetBranchAddress("nCscRechitsChamberMinus11",            &nCscRechitsChamberMinus11);
  tree_->SetBranchAddress("nCscRechitsChamberMinus12",            &nCscRechitsChamberMinus12);
  tree_->SetBranchAddress("nCscRechitsChamberMinus13",            &nCscRechitsChamberMinus13);
  tree_->SetBranchAddress("nCscRechitsChamberMinus21",            &nCscRechitsChamberMinus21);
  tree_->SetBranchAddress("nCscRechitsChamberMinus22",            &nCscRechitsChamberMinus22);
  tree_->SetBranchAddress("nCscRechitsChamberMinus31",            &nCscRechitsChamberMinus31);
  tree_->SetBranchAddress("nCscRechitsChamberMinus32",            &nCscRechitsChamberMinus32);
  tree_->SetBranchAddress("nCscRechitsChamberMinus41",            &nCscRechitsChamberMinus41);
  tree_->SetBranchAddress("nCscRechitsChamberMinus42",            &nCscRechitsChamberMinus42);


  tree_->SetBranchAddress("nDTRechits",            &nDTRechits);
  tree_->SetBranchAddress("nDTPositiveYRechits",            &nDTPositiveYRechits);
  tree_->SetBranchAddress("nDTNegativeYRechits",            &nDTNegativeYRechits);
  tree_->SetBranchAddress("nDtRings",             &nDtRings);
  tree_->SetBranchAddress("nDtWheels25",             &nDtWheels25);
  tree_->SetBranchAddress("nDtStations25",             &nDtStations25);

  tree_->SetBranchAddress("nDTRechitsWheelMinus2",             &nDTRechitsWheelMinus2);
  tree_->SetBranchAddress("nDTRechitsWheelMinus1",             &nDTRechitsWheelMinus1);
  tree_->SetBranchAddress("nDTRechitsWheel0",             &nDTRechitsWheel0);
  tree_->SetBranchAddress("nDTRechitsWheelPlus1",             &nDTRechitsWheelPlus1);
  tree_->SetBranchAddress("nDTRechitsWheelPlus2",             &nDTRechitsWheelPlus2);

  tree_->SetBranchAddress("nDTRechitsStation1",             &nDTRechitsStation1);
  tree_->SetBranchAddress("nDTRechitsStation2",             &nDTRechitsStation2);
  tree_->SetBranchAddress("nDTRechitsStation3",             &nDTRechitsStation3);
  tree_->SetBranchAddress("nDTRechitsStation4",             &nDTRechitsStation4);


  tree_->SetBranchAddress("nDTRechitsChamberMinus12",            &nDTRechitsChamberMinus12);
  tree_->SetBranchAddress("nDTRechitsChamberMinus11",            &nDTRechitsChamberMinus11);
  tree_->SetBranchAddress("nDTRechitsChamber10",            &nDTRechitsChamber10);
  tree_->SetBranchAddress("nDTRechitsChamberPlus11",            &nDTRechitsChamberPlus11);
  tree_->SetBranchAddress("nDTRechitsChamberPlus12",            &nDTRechitsChamberPlus12);
  tree_->SetBranchAddress("nDTRechitsChamberMinus22",            &nDTRechitsChamberMinus22);
  tree_->SetBranchAddress("nDTRechitsChamberMinus21",            &nDTRechitsChamberMinus21);
  tree_->SetBranchAddress("nDTRechitsChamber20",            &nDTRechitsChamber20);
  tree_->SetBranchAddress("nDTRechitsChamberPlus21",            &nDTRechitsChamberPlus21);
  tree_->SetBranchAddress("nDTRechitsChamberPlus22",            &nDTRechitsChamberPlus22);
  tree_->SetBranchAddress("nDTRechitsChamberMinus32",            &nDTRechitsChamberMinus32);
  tree_->SetBranchAddress("nDTRechitsChamberMinus31",            &nDTRechitsChamberMinus31);
  tree_->SetBranchAddress("nDTRechitsChamber30",            &nDTRechitsChamber30);

  tree_->SetBranchAddress("nDTRechitsChamberPlus31",            &nDTRechitsChamberPlus31);
  tree_->SetBranchAddress("nDTRechitsChamberPlus32",            &nDTRechitsChamberPlus32);
  tree_->SetBranchAddress("nDTRechitsChamberMinus42",            &nDTRechitsChamberMinus42);
  tree_->SetBranchAddress("nDTRechitsChamberMinus41",            &nDTRechitsChamberMinus41);
  tree_->SetBranchAddress("nDTRechitsChamber40",            &nDTRechitsChamber40);
  tree_->SetBranchAddress("nDTRechitsChamberPlus41",            &nDTRechitsChamberPlus41);
  tree_->SetBranchAddress("nDTRechitsChamberPlus42",            &nDTRechitsChamberPlus42);


  //DT CLUSTER

    tree_->SetBranchAddress("nDtRechitClusters",             &nDtRechitClusters);
    tree_->SetBranchAddress("dtRechitClusterX",             dtRechitClusterX);
    tree_->SetBranchAddress("dtRechitClusterY",             dtRechitClusterY);
    tree_->SetBranchAddress("dtRechitClusterZ",             dtRechitClusterZ);
    tree_->SetBranchAddress("dtRechitClusterWheel",             dtRechitClusterWheel);
    tree_->SetBranchAddress("dtRechitClusterEta",             dtRechitClusterEta);
    tree_->SetBranchAddress("dtRechitClusterPhi",             dtRechitClusterPhi);
    tree_->SetBranchAddress("dtRechitClusterSize",             dtRechitClusterSize);
    tree_->SetBranchAddress("dtRechitClusterMaxStation",             dtRechitClusterMaxStation);
    tree_->SetBranchAddress("dtRechitClusterMaxStationRatio",             dtRechitClusterMaxStationRatio);
    tree_->SetBranchAddress("dtRechitClusterNStation10",             dtRechitClusterNStation10);
    tree_->SetBranchAddress("dtRechitClusterAvgStation10",             dtRechitClusterAvgStation10);
    tree_->SetBranchAddress("dtRechitClusterMaxChamber",             dtRechitClusterMaxChamber);
    tree_->SetBranchAddress("dtRechitClusterMaxChamberRatio",             dtRechitClusterMaxChamberRatio);
    tree_->SetBranchAddress("dtRechitClusterNChamber",             dtRechitClusterNChamber);
    tree_->SetBranchAddress("dtRechitClusterJetVetoPt",             dtRechitClusterJetVetoPt);
    tree_->SetBranchAddress("dtRechitClusterJetVetoEta",             dtRechitClusterJetVetoEta);
    tree_->SetBranchAddress("dtRechitClusterJetVetoPhi",             dtRechitClusterJetVetoPhi);
    tree_->SetBranchAddress("dtRechitClusterJetVetoE",             dtRechitClusterJetVetoE);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoPt",             dtRechitClusterMuonVetoPt);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoE",             dtRechitClusterMuonVetoE);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoPhi",             dtRechitClusterMuonVetoPhi);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoEta",             dtRechitClusterMuonVetoEta);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoLooseId",             dtRechitClusterMuonVetoLooseId);
    tree_->SetBranchAddress("dtRechitClusterMuonVetoGlobal",             dtRechitClusterMuonVetoGlobal);
    tree_->SetBranchAddress("dtRechitClusterNSegmentStation1",             dtRechitClusterNSegmentStation1);
    tree_->SetBranchAddress("dtRechitClusterNSegmentStation2",             dtRechitClusterNSegmentStation2);
    tree_->SetBranchAddress("dtRechitClusterNSegmentStation3",             dtRechitClusterNSegmentStation3);
    tree_->SetBranchAddress("dtRechitClusterNSegmentStation4",             dtRechitClusterNSegmentStation4);
    tree_->SetBranchAddress("dtRechitCluster_match_gParticle_Id",             dtRechitCluster_match_gParticle_Id);
    tree_->SetBranchAddress("dtRechitCluster_match_gParticle_Pt",             dtRechitCluster_match_gParticle_Pt);
    tree_->SetBranchAddress("dtRechitCluster_match_gParticle_Eta",             dtRechitCluster_match_gParticle_Eta);
    tree_->SetBranchAddress("dtRechitCluster_match_gParticle_Phi",             dtRechitCluster_match_gParticle_Phi);
    tree_->SetBranchAddress("dtRechitCluster_match_gParticle_E",             dtRechitCluster_match_gParticle_E);
    tree_->SetBranchAddress("dtRechitCluster_match_gParticle_Status",             dtRechitCluster_match_gParticle_Status);
    tree_->SetBranchAddress("dtRechitCluster_match_gParticle_MotherId",             dtRechitCluster_match_gParticle_MotherId);
    tree_->SetBranchAddress("dtRechitCluster_match_gParticle_deltaR",             dtRechitCluster_match_gParticle_deltaR);
    tree_->SetBranchAddress("dtRechitCluster_match_gLLP",             dtRechitCluster_match_gLLP);
    tree_->SetBranchAddress("dtRechitCluster_match_gLLP_minDeltaR",             dtRechitCluster_match_gLLP_minDeltaR);
    tree_->SetBranchAddress("dtRechitCluster_match_gLLP_eta",             dtRechitCluster_match_gLLP_eta);
    tree_->SetBranchAddress("dtRechitCluster_match_gLLP_phi",             dtRechitCluster_match_gLLP_phi);
    tree_->SetBranchAddress("dtRechitCluster_match_gLLP_decay_x",             dtRechitCluster_match_gLLP_decay_x);
    tree_->SetBranchAddress("dtRechitCluster_match_gLLP_decay_y",             dtRechitCluster_match_gLLP_decay_y);
    tree_->SetBranchAddress("dtRechitCluster_match_gLLP_decay_z",             dtRechitCluster_match_gLLP_decay_z);
    tree_->SetBranchAddress("dtRechitCluster_match_gLLP_ctau",             dtRechitCluster_match_gLLP_ctau);
    tree_->SetBranchAddress("dtRechitCluster_match_gLLP_beta",             dtRechitCluster_match_gLLP_beta);
    tree_->SetBranchAddress("dtRechitCluster_match_gLLP_csc",             dtRechitCluster_match_gLLP_csc);
    tree_->SetBranchAddress("dtRechitCluster_match_gLLP_dt",             dtRechitCluster_match_gLLP_dt);
    tree_->SetBranchAddress("dtRechitCluster_match_gLLP_e",             dtRechitCluster_match_gLLP_e);
    tree_->SetBranchAddress("dtRechitCluster_match_gLLP_pt",             dtRechitCluster_match_gLLP_pt);
    tree_->SetBranchAddress("dtRechitClusterLep_dPhi",             dtRechitClusterLep_dPhi);
    tree_->SetBranchAddress("dtRechitCluster_match_MB1hits_0p4",             dtRechitCluster_match_MB1hits_0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_MB1hits_0p5",             dtRechitCluster_match_MB1hits_0p5);
    tree_->SetBranchAddress("dtRechitCluster_match_MB1hits_cosmics_plus",             dtRechitCluster_match_MB1hits_cosmics_plus);
    tree_->SetBranchAddress("dtRechitCluster_match_MB1hits_cosmics_minus",             dtRechitCluster_match_MB1hits_cosmics_minus);
    tree_->SetBranchAddress("dtRechitCluster_match_RPChits_dPhi0p5",             dtRechitCluster_match_RPChits_dPhi0p5);
    tree_->SetBranchAddress("dtRechitCluster_match_RPCBx_dPhi0p5",             dtRechitCluster_match_RPCBx_dPhi0p5);
    tree_->SetBranchAddress("dtRechitCluster_match_RB1_0p4",             dtRechitCluster_match_RB1_0p4);
    tree_->SetBranchAddress("dtRechitCluster_match_RB1_dPhi0p5",             dtRechitCluster_match_RB1_dPhi0p5);
    tree_->SetBranchAddress("dtRechitClusterMetEENoise_dPhi",             dtRechitClusterMetEENoise_dPhi);




  // CSC CLUSTER
  tree_->SetBranchAddress("nCscRechitClusters3",             &nCscRechitClusters3);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me1112_0p4",             &cscRechitCluster3_match_Me1112_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me1112_0p6",             &cscRechitCluster3_match_Me1112_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me1112_0p8",             &cscRechitCluster3_match_Me1112_0p8);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me11_0p4",             &cscRechitCluster3_match_Me11_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me11_0p6",             &cscRechitCluster3_match_Me11_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me11_0p8",             &cscRechitCluster3_match_Me11_0p8);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me12_0p4",             &cscRechitCluster3_match_Me12_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me12_0p6",             &cscRechitCluster3_match_Me12_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_Me12_0p8",             &cscRechitCluster3_match_Me12_0p8);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_id",             &cscRechitCluster3_match_gParticle_id);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle",             &cscRechitCluster3_match_gParticle);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_minDeltaR",             &cscRechitCluster3_match_gParticle_minDeltaR);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_index",             &cscRechitCluster3_match_gParticle_index);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_eta",             &cscRechitCluster3_match_gParticle_eta);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_phi",             &cscRechitCluster3_match_gParticle_phi);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_E",             &cscRechitCluster3_match_gParticle_E);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_pt",             &cscRechitCluster3_match_gParticle_pt);
  tree_->SetBranchAddress("cscRechitCluster3_match_gParticle_MotherId",             &cscRechitCluster3_match_gParticle_MotherId);

  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP",             &cscRechitCluster3_match_gLLP);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_index",             &cscRechitCluster3_match_gLLP_index);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_minDeltaR",             &cscRechitCluster3_match_gLLP_minDeltaR);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_eta",             &cscRechitCluster3_match_gLLP_eta);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_phi",             &cscRechitCluster3_match_gLLP_phi);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_decay_r",             &cscRechitCluster3_match_gLLP_decay_r);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_decay_x",             &cscRechitCluster3_match_gLLP_decay_x);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_decay_y",             &cscRechitCluster3_match_gLLP_decay_y);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_decay_z",             &cscRechitCluster3_match_gLLP_decay_z);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_ctau",             &cscRechitCluster3_match_gLLP_ctau);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_beta",             &cscRechitCluster3_match_gLLP_beta);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_csc",             &cscRechitCluster3_match_gLLP_csc);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_e",             &cscRechitCluster3_match_gLLP_e);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_pt",             &cscRechitCluster3_match_gLLP_pt);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_lepdPhi",             &cscRechitCluster3_match_gLLP_lepdPhi);

  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_daughter0_deltaR",             &cscRechitCluster3_match_gLLP_daughter0_deltaR);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_daughter1_deltaR",             &cscRechitCluster3_match_gLLP_daughter1_deltaR);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_daughter2_deltaR",             &cscRechitCluster3_match_gLLP_daughter2_deltaR);
  tree_->SetBranchAddress("cscRechitCluster3_match_gLLP_daughter3_deltaR",             &cscRechitCluster3_match_gLLP_daughter3_deltaR);


  tree_->SetBranchAddress("cscRechitCluster3Me11Ratio",             &cscRechitCluster3Me11Ratio);
  tree_->SetBranchAddress("cscRechitCluster3Me12Ratio",             &cscRechitCluster3Me12Ratio);
  tree_->SetBranchAddress("cscRechitCluster3MaxStation",             &cscRechitCluster3MaxStation);
  tree_->SetBranchAddress("cscRechitCluster3MaxStationRatio",             &cscRechitCluster3MaxStationRatio);
  tree_->SetBranchAddress("cscRechitCluster3NStation",             &cscRechitCluster3NStation);
  tree_->SetBranchAddress("cscRechitCluster3NStation5",             &cscRechitCluster3NStation5);
  tree_->SetBranchAddress("cscRechitCluster3NStation10",             &cscRechitCluster3NStation10);
  tree_->SetBranchAddress("cscRechitCluster3NStation10perc",             &cscRechitCluster3NStation10perc);
  tree_->SetBranchAddress("cscRechitCluster3AvgStation",             &cscRechitCluster3AvgStation);
  tree_->SetBranchAddress("cscRechitCluster3AvgStation5",             &cscRechitCluster3AvgStation5);
  tree_->SetBranchAddress("cscRechitCluster3AvgStation10",             &cscRechitCluster3AvgStation10);
  tree_->SetBranchAddress("cscRechitCluster3AvgStation10perc",             &cscRechitCluster3AvgStation10perc);
  tree_->SetBranchAddress("cscRechitCluster3MaxChamber",             &cscRechitCluster3MaxChamber);
  tree_->SetBranchAddress("cscRechitCluster3MaxChamberRatio",             &cscRechitCluster3MaxChamberRatio);
  tree_->SetBranchAddress("cscRechitCluster3NChamber",             &cscRechitCluster3NChamber);
  tree_->SetBranchAddress("cscRechitCluster3X",             cscRechitCluster3X);
  tree_->SetBranchAddress("cscRechitCluster3Y",             cscRechitCluster3Y);
  tree_->SetBranchAddress("cscRechitCluster3Z",             cscRechitCluster3Z);
  tree_->SetBranchAddress("cscRechitCluster3Time",             cscRechitCluster3Time);
  tree_->SetBranchAddress("cscRechitCluster3TimeWire",             cscRechitCluster3TimeWire);
  tree_->SetBranchAddress("cscRechitCluster3TimeWirePruned",             cscRechitCluster3TimeWirePruned);
  tree_->SetBranchAddress("cscRechitCluster3TimeTotal",             cscRechitCluster3TimeTotal);

  tree_->SetBranchAddress("cscRechitCluster3GenMuonDeltaR",             cscRechitCluster3GenMuonDeltaR);

  tree_->SetBranchAddress("cscRechitCluster3TimeSpread",             cscRechitCluster3TimeSpread);
  tree_->SetBranchAddress("cscRechitCluster3TimeTotalSpread",             cscRechitCluster3TimeTotalSpread);
  tree_->SetBranchAddress("cscRechitCluster3TimeTotalSpreadPruned",             cscRechitCluster3TimeTotalSpreadPruned);
  tree_->SetBranchAddress("cscRechitCluster3TimeWireSpread",             cscRechitCluster3TimeWireSpread);
  tree_->SetBranchAddress("cscRechitCluster3MajorAxis",             cscRechitCluster3MajorAxis);
  tree_->SetBranchAddress("cscRechitCluster3MinorAxis",             cscRechitCluster3MinorAxis);
  tree_->SetBranchAddress("cscRechitCluster3RSpread",             cscRechitCluster3RSpread);

  tree_->SetBranchAddress("cscRechitCluster3XSpread",             cscRechitCluster3XSpread);
  tree_->SetBranchAddress("cscRechitCluster3YSpread",             cscRechitCluster3YSpread);
  tree_->SetBranchAddress("cscRechitCluster3ZSpread",             cscRechitCluster3ZSpread);
  tree_->SetBranchAddress("cscRechitCluster3EtaPhiSpread",             cscRechitCluster3EtaPhiSpread);
  tree_->SetBranchAddress("cscRechitCluster3XYSpread",             cscRechitCluster3XYSpread);
  tree_->SetBranchAddress("cscRechitCluster3EtaSpread",             cscRechitCluster3EtaSpread);
  tree_->SetBranchAddress("cscRechitCluster3PhiSpread",             cscRechitCluster3PhiSpread);
  tree_->SetBranchAddress("cscRechitCluster3DeltaRSpread",             cscRechitCluster3DeltaRSpread);
  tree_->SetBranchAddress("cscRechitCluster3Eta",             cscRechitCluster3Eta);
  tree_->SetBranchAddress("cscRechitCluster3Phi",             cscRechitCluster3Phi);


  tree_->SetBranchAddress("cscRechitCluster3GenJetVetoPt",             cscRechitCluster3GenJetVetoPt);
  tree_->SetBranchAddress("cscRechitCluster3GenJetVetoE",             cscRechitCluster3GenJetVetoE);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoPt",             cscRechitCluster3JetVetoPt);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoEta",             cscRechitCluster3JetVetoEta);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoPhi",             cscRechitCluster3JetVetoPhi);

  tree_->SetBranchAddress("cscRechitCluster3JetVetoElectronEnergyFraction",             cscRechitCluster3JetVetoElectronEnergyFraction);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoPhotonEnergyFraction",             cscRechitCluster3JetVetoPhotonEnergyFraction);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoNeutralHadronEnergyFraction",             cscRechitCluster3JetVetoNeutralHadronEnergyFraction);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoChargedHadronEnergyFraction",             cscRechitCluster3JetVetoChargedHadronEnergyFraction);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoMuonEnergyFraction",             cscRechitCluster3JetVetoMuonEnergyFraction);


  tree_->SetBranchAddress("cscRechitCluster3JetVetoE",             cscRechitCluster3JetVetoE);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoPt",             cscRechitCluster3MuonVetoPt);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoE",             cscRechitCluster3MuonVetoE);

  tree_->SetBranchAddress("cscRechitCluster3JetVetoPt_0p6",             cscRechitCluster3JetVetoPt_0p6);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoPt_0p8",             cscRechitCluster3JetVetoPt_0p8);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoE_0p6",             cscRechitCluster3JetVetoE_0p6);
  tree_->SetBranchAddress("cscRechitCluster3JetVetoE_0p8",             cscRechitCluster3JetVetoE_0p8);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoPt_0p6",             cscRechitCluster3MuonVetoPt_0p6);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoPt_0p8",             cscRechitCluster3MuonVetoPt_0p8);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoE_0p6",             cscRechitCluster3MuonVetoE_0p6);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoE_0p8",             cscRechitCluster3MuonVetoE_0p8);

  tree_->SetBranchAddress("cscRechitCluster3ZLep1",             cscRechitCluster3ZLep1);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2",             cscRechitCluster3ZLep2);
  tree_->SetBranchAddress("cscRechitCluster3ZLep1Tag",             cscRechitCluster3ZLep1Tag);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2Tag",             cscRechitCluster3ZLep2Tag);
  tree_->SetBranchAddress("cscRechitCluster3ZLep1Id",             cscRechitCluster3ZLep1Id);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2Id",             cscRechitCluster3ZLep2Id);
  tree_->SetBranchAddress("cscRechitCluster3ZLep1LooseIso",             cscRechitCluster3ZLep1LooseIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep1TightIso",             cscRechitCluster3ZLep1TightIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep1VTightIso",             cscRechitCluster3ZLep1VTightIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep1VVTightIso",             cscRechitCluster3ZLep1VVTightIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep1TightId",             cscRechitCluster3ZLep1TightId);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2LooseIso",             cscRechitCluster3ZLep2LooseIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2TightIso",             cscRechitCluster3ZLep2TightIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2VTightIso",             cscRechitCluster3ZLep2VTightIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2VVTightIso",             cscRechitCluster3ZLep2VVTightIso);
  tree_->SetBranchAddress("cscRechitCluster3ZLep2TightId",             cscRechitCluster3ZLep2TightId);

  tree_->SetBranchAddress("cscRechitCluster3MuonVetoPhi",             cscRechitCluster3MuonVetoPhi);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoEta",             cscRechitCluster3MuonVetoEta);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoIso",             cscRechitCluster3MuonVetoIso);

    tree_->SetBranchAddress("cscRechitCluster3MuonVetoLooseIso",             cscRechitCluster3MuonVetoLooseIso);
    tree_->SetBranchAddress("cscRechitCluster3MuonVetoTightIso",             cscRechitCluster3MuonVetoTightIso);
    tree_->SetBranchAddress("cscRechitCluster3MuonVetoVTightIso",             cscRechitCluster3MuonVetoVTightIso);
    tree_->SetBranchAddress("cscRechitCluster3MuonVetoVVTightIso",             cscRechitCluster3MuonVetoVVTightIso);
    tree_->SetBranchAddress("cscRechitCluster3MuonVetoTightId",             cscRechitCluster3MuonVetoTightId);
    tree_->SetBranchAddress("cscRechitCluster3MuonVetoLooseId",             cscRechitCluster3MuonVetoLooseId);


  tree_->SetBranchAddress("cscRechitCluster3IsoMuonVetoPt",             cscRechitCluster3IsoMuonVetoPt);
  tree_->SetBranchAddress("cscRechitCluster3IsoMuonVetoE",             cscRechitCluster3IsoMuonVetoE);
  tree_->SetBranchAddress("cscRechitCluster3IsoMuonVetoPhi",             cscRechitCluster3IsoMuonVetoPhi);
  tree_->SetBranchAddress("cscRechitCluster3IsoMuonVetoEta",             cscRechitCluster3IsoMuonVetoEta);
  tree_->SetBranchAddress("cscRechitCluster3GenMuonVetoPt",             cscRechitCluster3GenMuonVetoPt);
  tree_->SetBranchAddress("cscRechitCluster3MuonVetoType",             cscRechitCluster3MuonVetoType);
  tree_->SetBranchAddress("cscRechitCluster3GenMuonVetoProdX",             cscRechitCluster3GenMuonVetoProdX);
  tree_->SetBranchAddress("cscRechitCluster3GenMuonVetoProdY",             cscRechitCluster3GenMuonVetoProdY);
  tree_->SetBranchAddress("cscRechitCluster3GenMuonVetoProdZ",             cscRechitCluster3GenMuonVetoProdZ);
  tree_->SetBranchAddress("cscRechitCluster3GenMuonVetoLLPDist",             cscRechitCluster3GenMuonVetoLLPDist);
  tree_->SetBranchAddress("cscRechitCluster3GenMuonVetoLLPIndex",             cscRechitCluster3GenMuonVetoLLPIndex);

  tree_->SetBranchAddress("cscRechitCluster3GenMuonVetoE",             cscRechitCluster3GenMuonVetoE);
  tree_->SetBranchAddress("cscRechitCluster3Size",             cscRechitCluster3Size);

  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus11",             cscRechitCluster3NRechitChamberPlus11);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus12",             cscRechitCluster3NRechitChamberPlus12);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus13",             cscRechitCluster3NRechitChamberPlus13);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus21",             cscRechitCluster3NRechitChamberPlus21);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus22",             cscRechitCluster3NRechitChamberPlus22);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus31",             cscRechitCluster3NRechitChamberPlus31);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus32",             cscRechitCluster3NRechitChamberPlus32);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus41",             cscRechitCluster3NRechitChamberPlus41);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberPlus42",             cscRechitCluster3NRechitChamberPlus42);

  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus11",             cscRechitCluster3NRechitChamberMinus11);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus12",             cscRechitCluster3NRechitChamberMinus12);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus13",             cscRechitCluster3NRechitChamberMinus13);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus21",             cscRechitCluster3NRechitChamberMinus21);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus22",             cscRechitCluster3NRechitChamberMinus22);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus31",             cscRechitCluster3NRechitChamberMinus31);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus32",             cscRechitCluster3NRechitChamberMinus32);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus41",             cscRechitCluster3NRechitChamberMinus41);
  tree_->SetBranchAddress("cscRechitCluster3NRechitChamberMinus42",             cscRechitCluster3NRechitChamberMinus42);

  tree_->SetBranchAddress("cscRechitCluster3_match_cscRechits_0p4",             cscRechitCluster3_match_cscRechits_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_cscSeg_0p4",             cscRechitCluster3_match_cscSeg_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_ME11Seg_0p4",             cscRechitCluster3_match_ME11Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_ME12Seg_0p4",             cscRechitCluster3_match_ME12Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_cscSeg_0p6",             cscRechitCluster3_match_cscSeg_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_ME11Seg_0p6",             cscRechitCluster3_match_ME11Seg_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_ME12Seg_0p6",             cscRechitCluster3_match_ME12Seg_0p6);

  tree_->SetBranchAddress("cscRechitCluster3_match_dtRechits_0p4",             cscRechitCluster3_match_dtRechits_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_MB1_0p4",             cscRechitCluster3_match_MB1_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_dtRechits_0p6",             cscRechitCluster3_match_dtRechits_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_dtRechits_phi0p2",             cscRechitCluster3_match_dtRechits_phi0p2);


  tree_->SetBranchAddress("cscRechitCluster3_match_MB1_0p6",             cscRechitCluster3_match_MB1_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_dtSeg_0p4",             cscRechitCluster3_match_dtSeg_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_MB1Seg_0p4",             cscRechitCluster3_match_MB1Seg_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_dtSeg_0p6",             cscRechitCluster3_match_dtSeg_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_MB1Seg_0p6",             cscRechitCluster3_match_MB1Seg_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_RB1_0p4",             cscRechitCluster3_match_RB1_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_RE12_0p4",             cscRechitCluster3_match_RE12_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_RB1_0p6",             cscRechitCluster3_match_RB1_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_RE12_0p6",             cscRechitCluster3_match_RE12_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_highEta_0p4",             cscRechitCluster3_match_highEta_0p4);
  tree_->SetBranchAddress("cscRechitCluster3_match_highEta_0p6",             cscRechitCluster3_match_highEta_0p6);
  tree_->SetBranchAddress("cscRechitCluster3_match_highEta_0p8",             cscRechitCluster3_match_highEta_0p8);
  tree_->SetBranchAddress("cscRechitCluster3_match_cluster_dR",             cscRechitCluster3_match_cluster_dR);
  tree_->SetBranchAddress("cscRechitCluster3_match_cluster_index",             cscRechitCluster3_match_cluster_index);

  tree_->SetBranchAddress("cscRechitCluster3Met_dPhi",             cscRechitCluster3Met_dPhi);
  tree_->SetBranchAddress("cscRechitCluster3MetXYCorr_dPhi",             cscRechitCluster3MetXYCorr_dPhi);


  tree_->SetBranchAddress("cscRechitCluster3MetHEM_dPhi",             cscRechitCluster3MetHEM_dPhi);
  tree_->SetBranchAddress("cscRechitCluster3MetHEMXYCorr_dPhi",             cscRechitCluster3MetHEMXYCorr_dPhi);
  tree_->SetBranchAddress("cscRechitCluster3MetEENoise_dPhi",             cscRechitCluster3MetEENoise_dPhi);
  tree_->SetBranchAddress("cscRechitCluster3MetEENoiseXYCorr_dPhi",             cscRechitCluster3MetEENoiseXYCorr_dPhi);
  tree_->SetBranchAddress("cscRechitCluster3MetJesUp_dPhi",             cscRechitCluster3MetJesUp_dPhi);
  tree_->SetBranchAddress("cscRechitCluster3MetJesDown_dPhi",             cscRechitCluster3MetJesDown_dPhi);

  tree_->SetBranchAddress("gLLP_multiplicity",    gLLP_multiplicity);

  tree_->SetBranchAddress("gLLP_eta",    gLLP_eta);
  tree_->SetBranchAddress("gLLP_phi",    gLLP_phi);
  tree_->SetBranchAddress("gLLP_csc",    gLLP_csc);
  tree_->SetBranchAddress("gLLP_dt",    gLLP_dt);
  tree_->SetBranchAddress("gLLP_ctau",    gLLP_ctau);
  tree_->SetBranchAddress("gLLP_beta",    gLLP_beta);
  tree_->SetBranchAddress("gLLP_e",    gLLP_e);
  tree_->SetBranchAddress("gLLP_pt",    gLLP_pt);
  tree_->SetBranchAddress("gLLP_lepdPhi",    gLLP_lepdPhi);

  tree_->SetBranchAddress("gLLP_decay_vertex_r",    gLLP_decay_vertex_r);
  tree_->SetBranchAddress("gLLP_decay_vertex_x",    gLLP_decay_vertex_x);
  tree_->SetBranchAddress("gLLP_decay_vertex_y",    gLLP_decay_vertex_y);
  tree_->SetBranchAddress("gLLP_decay_vertex_z",    gLLP_decay_vertex_z);

  tree_->SetBranchAddress("gLLP_daughter_deltaR",    gLLP_daughter_deltaR);



tree_->SetBranchAddress("gLLP_daughter_id",          gLLP_daughter_id);
tree_->SetBranchAddress("gLLP_daughter_pt",          gLLP_daughter_pt);
    tree_->SetBranchAddress("gLLP_daughter_eta",          gLLP_daughter_eta);
    tree_->SetBranchAddress("gLLP_daughter_phi",          gLLP_daughter_phi);
    tree_->SetBranchAddress("gLLP_daughter_e",          gLLP_daughter_e);
    tree_->SetBranchAddress("gLLP_daughter_mass",          gLLP_daughter_mass);

    tree_->SetBranchAddress("nMuons",    &nMuons);
    tree_->SetBranchAddress("muonPt",       muonPt);
    tree_->SetBranchAddress("muonEta",      muonEta);
    tree_->SetBranchAddress("muonPhi",      muonPhi);

  //Leptons
  tree_->SetBranchAddress("nLeptons",    &nLeptons);
  tree_->SetBranchAddress("lepE",        lepE);
  tree_->SetBranchAddress("lepPt",       lepPt);
  tree_->SetBranchAddress("lepEta",      lepEta);
  tree_->SetBranchAddress("lepPhi",      lepPhi);
  tree_->SetBranchAddress("lepPdgId",  lepPdgId);
  tree_->SetBranchAddress("lepDZ",     lepDZ);
  tree_->SetBranchAddress("lepEff", lepEff);
  tree_->SetBranchAddress("lepSF", lepSF);

  tree_->SetBranchAddress("lepTriggerSF", lepTriggerSF);
  tree_->SetBranchAddress("lepTightIdSF", lepTightIdSF);
  tree_->SetBranchAddress("lepLooseIdSF", lepLooseIdSF);
  tree_->SetBranchAddress("lepTightIsoSF", lepTightIsoSF);
  tree_->SetBranchAddress("lepLooseIsoSF", lepLooseIsoSF);
  tree_->SetBranchAddress("lepTriggerMCEfficiency", lepTriggerMCEfficiency);
  tree_->SetBranchAddress("lepTightIdMCEfficiency", lepTightIdMCEfficiency);
  tree_->SetBranchAddress("lepLooseIdMCEfficiency", lepLooseIdMCEfficiency);
  tree_->SetBranchAddress("lepTightIsoMCEfficiency", lepTightIsoMCEfficiency);
  tree_->SetBranchAddress("lepLooseIsoMCEfficiency", lepLooseIsoMCEfficiency);
  tree_->SetBranchAddress("lepTag", lepTag);



  // tree_->SetBranchAddress("lepLoosePassId", lepLoosePassId);
  // tree_->SetBranchAddress("lepMediumPassId", lepMediumPassId);
  // tree_->SetBranchAddress("lepTightPassId", lepTightPassId);
  tree_->SetBranchAddress("lepPassId", lepPassId);
  tree_->SetBranchAddress("lepFromZ", lepFromZ);

  tree_->SetBranchAddress("lepPassVetoId", lepPassVetoId);
  tree_->SetBranchAddress("lepPassLooseIso", lepPassLooseIso);
  tree_->SetBranchAddress("lepPassTightIso", lepPassTightIso);
  tree_->SetBranchAddress("lepPassVTightIso", lepPassVTightIso);
  tree_->SetBranchAddress("lepPassVTightIso", lepPassVTightIso);


  //Z-candidate

  tree_->SetBranchAddress("ZMass1",       &ZMass1);

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
  tree_->SetBranchAddress("jetPt",     jetPt);
  tree_->SetBranchAddress("jetEta",    jetEta);
  tree_->SetBranchAddress("jetPhi",    jetPhi);
  tree_->SetBranchAddress("jetTime",   jetTime);
  tree_->SetBranchAddress("jetPassId", jetPassId);

  tree_->SetBranchAddress("jetPtJESUp", jetPtJESUp);
  tree_->SetBranchAddress("jetPtJESDown", jetPtJESDown);
  tree_->SetBranchAddress("jetEJESUp", jetEJESUp);
  tree_->SetBranchAddress("jetEJESDown", jetEJESDown);
  tree_->SetBranchAddress("JecUnc", JecUnc);

  tree_->SetBranchAddress("jet_match_genJet_pt", jet_match_genJet_pt);
  tree_->SetBranchAddress("jet_match_genJet_index", jet_match_genJet_index);
  tree_->SetBranchAddress("jet_match_genJet_minDeltaR", jet_match_genJet_minDeltaR);

  // tree_->SetBranchAddress("ecalNRechits",   ecalNRechits);*/
  // tree_->SetBranchAddress("ecalRechitE", ecalRechitE);
  tree_->SetBranchAddress("jetElectronEnergyFraction", jetElectronEnergyFraction);
  tree_->SetBranchAddress("jetPhotonEnergyFraction", jetPhotonEnergyFraction);
  tree_->SetBranchAddress("jetChargedHadronEnergyFraction", jetChargedHadronEnergyFraction);
  tree_->SetBranchAddress("jetNeutralHadronEnergyFraction", jetNeutralHadronEnergyFraction);
  tree_->SetBranchAddress("jetMuonEnergyFraction", jetMuonEnergyFraction);

  // tree_->SetBranchAddress("jetLoosePassId", jetLoosePassId);
  tree_->SetBranchAddress("jetTightPassId", jetTightPassId);
  // triggers
  tree_->SetBranchAddress("HLTDecision",   HLTDecision);
  tree_->SetBranchAddress("SingleMuonTrigger",   &SingleMuonTrigger);
  tree_->SetBranchAddress("SingleEleTrigger",   &SingleEleTrigger);
  tree_->SetBranchAddress("SingleLepTrigger",   &SingleLepTrigger);
};

void HNLMuonSystemTree::LoadTree(const char* file)
{
  f_ = TFile::Open(file);
  assert(f_);
  tree_ = dynamic_cast<TTree*>(f_->Get("MuonSystem"));
  InitTree();
  assert(tree_);
};

void HNLMuonSystemTree::CreateTree()
{
  tree_ = new TTree("MuonSystem","MuonSystem");
  f_ = 0;

  tree_->Branch("runNum",      &runNum,     "runNum/i");      // event run number
  tree_->Branch("lumiSec",     &lumiSec,    "lumiSec/i");     // event lumi section
  tree_->Branch("evtNum",      &evtNum,     "evtNum/i");      // event number
  tree_->Branch("mH",      &mH,     "mH/I");      // event number
  tree_->Branch("mX",      &mX,     "mX/I");      // event number
  tree_->Branch("ctau",      &ctau,     "ctau/I");      // event number
  tree_->Branch("ZCategory",    &ZCategory,   "ZCategory/I");    // dilepton category

  tree_->Branch("category",    &category,   "category/i");    // dilepton category
  tree_->Branch("npv",         &npv,        "npv/i");         // number of primary vertices
  tree_->Branch("npu",         &npu,        "npu/i");         // number of in-time PU events (MC)
  tree_->Branch("weight",      &weight,     "weight/F");
  tree_->Branch("scaleWeights",      scaleWeights,     "scaleWeights[9]/F");
  tree_->Branch("lepOverallSF",      &lepOverallSF,     "lepOverallSF/F");


  tree_->Branch("sf_facScaleUp",      &sf_facScaleUp,     "sf_facScaleUp/F");
  tree_->Branch("sf_facScaleDown",      &sf_facScaleDown,     "sf_facScaleDown/F");
  tree_->Branch("sf_renScaleUp",      &sf_renScaleUp,     "sf_renScaleUp/F");
  tree_->Branch("sf_renScaleDown",      &sf_renScaleDown,     "sf_renScaleDown/F");
  tree_->Branch("sf_facRenScaleUp",      &sf_facRenScaleUp,     "sf_facRenScaleUp/F");
  tree_->Branch("sf_facRenScaleDown",      &sf_facRenScaleDown,     "sf_facRenScaleDown/F");
  tree_->Branch("metSF",      &metSF,     "metSF/F");

  tree_->Branch("pileupWeight",      &pileupWeight,     "pileupWeight/F");
  tree_->Branch("pileupWeightUp",      &pileupWeightUp,     "pileupWeightUp/F");
  tree_->Branch("pileupWeightDown",      &pileupWeightDown,     "pileupWeightDown/F");
  tree_->Branch("Flag_HBHENoiseFilter",      &Flag_HBHENoiseFilter,     "Flag_HBHENoiseFilter/O");
  tree_->Branch("Flag_BadPFMuonFilter",      &Flag_BadPFMuonFilter,     "Flag_BadPFMuonFilter/O");
  tree_->Branch("Flag_HBHEIsoNoiseFilter",      &Flag_HBHEIsoNoiseFilter,     "Flag_HBHEIsoNoiseFilter/O");
  tree_->Branch("Flag_CSCTightHaloFilter",      &Flag_CSCTightHaloFilter,     "Flag_CSCTightHaloFilter/O");
  tree_->Branch("Flag_globalSuperTightHalo2016Filter",      &Flag_globalSuperTightHalo2016Filter,     "Flag_globalSuperTightHalo2016Filter/O");
  tree_->Branch("Flag_goodVertices",      &Flag_goodVertices,     "Flag_goodVertices/O");
  tree_->Branch("Flag_ecalBadCalibFilter",      &Flag_ecalBadCalibFilter,     "Flag_ecalBadCalibFilter/O");
  tree_->Branch("Flag_BadChargedCandidateFilter",      &Flag_BadChargedCandidateFilter,     "Flag_BadChargedCandidateFilter/O");
  tree_->Branch("Flag_eeBadScFilter",      &Flag_eeBadScFilter,     "Flag_eeBadScFilter/O");
  tree_->Branch("Flag_all",      &Flag_all,     "Flag_all/O");

  tree_->Branch("Flag2_HBHENoiseFilter",      &Flag2_HBHENoiseFilter,     "Flag2_HBHENoiseFilter/O");
  tree_->Branch("Flag2_HBHEIsoNoiseFilter",      &Flag2_HBHEIsoNoiseFilter,     "Flag2_HBHEIsoNoiseFilter/O");
  tree_->Branch("Flag2_BadPFMuonFilter",      &Flag2_BadPFMuonFilter,     "Flag2_BadPFMuonFilter/O");
  tree_->Branch("Flag2_globalSuperTightHalo2016Filter",      &Flag2_globalSuperTightHalo2016Filter,     "Flag2_globalSuperTightHalo2016Filter/O");
  tree_->Branch("Flag2_globalTightHalo2016Filter",      &Flag2_globalTightHalo2016Filter,     "Flag2_globalTightHalo2016Filter/O");
  tree_->Branch("Flag2_BadChargedCandidateFilter",      &Flag2_BadChargedCandidateFilter,     "Flag2_BadChargedCandidateFilter/O");
  tree_->Branch("Flag2_EcalDeadCellTriggerPrimitiveFilter",      &Flag2_EcalDeadCellTriggerPrimitiveFilter,     "Flag2_EcalDeadCellTriggerPrimitiveFilter/O");
  tree_->Branch("Flag2_ecalBadCalibFilter",      &Flag2_ecalBadCalibFilter,     "Flag2_ecalBadCalibFilter/O");
  tree_->Branch("Flag2_eeBadScFilter",      &Flag2_eeBadScFilter,     "Flag2_eeBadScFilter/O");
  tree_->Branch("Flag2_all",      &Flag2_all,     "Flag2_all/O");
  tree_->Branch("EE_prefiring",      &EE_prefiring,     "EE_prefiring/O");



  tree_->Branch("rho",         &rho,        "rho/F");
  tree_->Branch("met",         &met,        "met/F");         // MET
  tree_->Branch("metNoMu",         &metNoMu,        "metNoMu/F");         // MET


  tree_->Branch("metPhi",      &metPhi,     "metPhi/F");      // phi(MET)
  tree_->Branch("metXYCorr",      &metXYCorr,     "metXYCorr/F");      // phi(MET)
  tree_->Branch("metPhiXYCorr",      &metPhiXYCorr,     "metPhiXYCorr/F");      // phi(MET)

  tree_->Branch("HT",      &HT,     "HT/F");      // phi(MET)

  tree_->Branch("jetMet_dPhi",      &jetMet_dPhi,     "jetMet_dPhi/F");      // phi(MET)
  tree_->Branch("jetMet_dPhiMin",      &jetMet_dPhiMin,     "jetMet_dPhiMin/F");      // phi(MET)
  tree_->Branch("jetMet_dPhiMin4",      &jetMet_dPhiMin4,     "jetMet_dPhiMin4/F");      // phi(MET)

  tree_->Branch("metJESUp",      &metJESUp,     "metJESUp/F");      // phi(MET)
  tree_->Branch("metJESDown",      &metJESDown,     "metJESDown/F");      // phi(MET)
  tree_->Branch("metPhiJESUp",      &metPhiJESUp,     "metPhiJESUp/F");      // phi(metPhi)
  tree_->Branch("metPhiJESDown",      &metPhiJESDown,     "metPhiJESDown/F");      // phi(metPhi)
  tree_->Branch("metJESUpSF",      &metJESUpSF,     "metJESUpSF/F");      // phi(MET)
  tree_->Branch("metJESDownSF",      &metJESDownSF,     "metJESDownSF/F");      // phi(MET)
  tree_->Branch("metEENoise",      &metEENoise,     "metEENoise/F");      // phi(MET)
  tree_->Branch("metPhiEENoise",      &metPhiEENoise,     "metPhiEENoise/F");      // phi(MET)
  tree_->Branch("metHEM",      &metHEM,     "metHEM/F");      // phi(MET)
  tree_->Branch("metPhiHEM",      &metPhiHEM,     "metPhiHEM/F");      // phi(MET)
  tree_->Branch("metEENoiseXYCorr",      &metEENoiseXYCorr,     "metEENoiseXYCorr/F");      // phi(MET)
  tree_->Branch("metPhiEENoiseXYCorr",      &metPhiEENoiseXYCorr,     "metPhiEENoiseXYCorr/F");      // phi(MET)
  tree_->Branch("metHEMXYCorr",      &metHEMXYCorr,     "metHEMXYCorr/F");      // phi(MET)
  tree_->Branch("metPhiHEMXYCorr",      &metPhiHEMXYCorr,     "metPhiHEMXYCorr/F");      // phi(MET)
  tree_->Branch("genMetPtTrue",         &genMetPtTrue,        "genMetPtTrue/F");         // MET
  tree_->Branch("genMetPhiTrue",      &genMetPhiTrue,     "genMetPhiTrue/F");      // phi(MET)
  tree_->Branch("genMetPtCalo",         &genMetPtCalo,        "genMetPtCalo/F");         // MET
  tree_->Branch("genMetPhiCalo",      &genMetPhiCalo,     "genMetPhiCalo/F");      // phi(MET)

  // tree_->Branch("nGenParticle",      &nGenParticle,   "nGenParticle/I");
  // tree_->Branch("gParticleId",      gParticleId,  "gParticleId[nGenParticle]/I");
  // tree_->Branch("gParticleStatus",      gParticleStatus,  "gParticleStatus[nGenParticle]/I");
  // tree_->Branch("gParticleMotherId",      gParticleMotherId,  "gParticleMotherId[nGenParticle]/I");
  // tree_->Branch("gParticleE",      gParticleE,  "gParticleE[nGenParticle]/F");
  // tree_->Branch("gParticlePt",      gParticlePt,  "gParticlePt[nGenParticle]/F");
  // tree_->Branch("gParticleEta",      gParticleEta,  "gParticleEta[nGenParticle]/F");
  // tree_->Branch("gParticlePhi",      gParticlePhi,  "gParticlePhi[nGenParticle]/F");

  // tree_->Branch("nGenJets",      &nGenJets,  "nGenJets/I");
  // tree_->Branch("genJetE",      genJetE,  "genJetE[nGenJets]/F");
  // tree_->Branch("genJetPt",      genJetPt,  "genJetPt[nGenJets]/F");
  // tree_->Branch("genJetEta",      genJetEta,  "genJetEta[nGenJets]/F");
  // tree_->Branch("genJetPhi",      genJetPhi,  "genJetPhi[nGenJets]/F");
  // tree_->Branch("genJetMET",      genJetMET,  "genJetMET[nGenJets]/F");

  tree_->Branch("gWPt",         &gWPt,        "gWPt/F");

  tree_->Branch("gLepId",      &gLepId,     "gLepId/I");      // phi(MET)
  tree_->Branch("gLepPt",      &gLepPt,     "gLepPt/F");      // phi(MET)
  tree_->Branch("gLepE",      &gLepE,     "gLepE/F");      // phi(MET)
  tree_->Branch("gLepEta",      &gLepEta,     "gLepEta/F");      // phi(MET)
  tree_->Branch("gLepPhi",      &gLepPhi,     "gLepPhi/F");      // phi(MET)
  tree_->Branch("gHiggsPt",      &gHiggsPt,     "gHiggsPt/F");      // phi(MET)
  tree_->Branch("gHiggsE",      &gHiggsE,     "gHiggsE/F");      // phi(MET)
  tree_->Branch("gHiggsEta",      &gHiggsEta,     "gHiggsEta/F");      // phi(MET)
  tree_->Branch("gHiggsPhi",      &gHiggsPhi,     "gHiggsPhi/F");      // phi(MET)
  //CSC

  tree_->Branch("nCscRechits",             &nCscRechits, "nCscRechits/I");
  tree_->Branch("nEarlyCscRechits",             &nEarlyCscRechits, "nEarlyCscRechits/I");
  tree_->Branch("nLateCscRechits",             &nLateCscRechits, "nLateCscRechits/I");
  tree_->Branch("nEarly2CscRechits",             &nEarly2CscRechits, "nEarly2CscRechits/I");
  tree_->Branch("nLate2CscRechits",             &nLate2CscRechits, "nLate2CscRechits/I");
  tree_->Branch("nCscRings",             &nCscRings, "nCscRings/I");
  tree_->Branch("nCscPositiveYRechits",             &nCscPositiveYRechits, "nCscPositiveYRechits/I");
  tree_->Branch("nCscNegativeYRechits",             &nCscNegativeYRechits, "nCscNegativeYRechits/I");
  tree_->Branch("cscPosTpeak",             &cscPosTpeak, "cscPosTpeak/F");
  tree_->Branch("cscNegTpeak",             &cscNegTpeak, "cscNegTpeak/F");

  tree_->Branch("nCscRechitsChamberPlus11",            &nCscRechitsChamberPlus11,             "nCscRechitsChamberPlus11/I");
  tree_->Branch("nCscRechitsChamberPlus12",            &nCscRechitsChamberPlus12,             "nCscRechitsChamberPlus12/I");
  tree_->Branch("nCscRechitsChamberPlus13",            &nCscRechitsChamberPlus13,             "nCscRechitsChamberPlus13/I");
  tree_->Branch("nCscRechitsChamberPlus21",            &nCscRechitsChamberPlus21,             "nCscRechitsChamberPlus21/I");
  tree_->Branch("nCscRechitsChamberPlus22",            &nCscRechitsChamberPlus22,             "nCscRechitsChamberPlus22/I");
  tree_->Branch("nCscRechitsChamberPlus31",            &nCscRechitsChamberPlus31,             "nCscRechitsChamberPlus31/I");
  tree_->Branch("nCscRechitsChamberPlus32",            &nCscRechitsChamberPlus32,             "nCscRechitsChamberPlus32/I");
  tree_->Branch("nCscRechitsChamberPlus41",            &nCscRechitsChamberPlus41,             "nCscRechitsChamberPlus41/I");
  tree_->Branch("nCscRechitsChamberPlus42",            &nCscRechitsChamberPlus42,             "nCscRechitsChamberPlus42/I");
  tree_->Branch("nCscRechitsChamberMinus11",            &nCscRechitsChamberMinus11,             "nCscRechitsChamberMinus11/I");
  tree_->Branch("nCscRechitsChamberMinus12",            &nCscRechitsChamberMinus12,             "nCscRechitsChamberMinus12/I");
  tree_->Branch("nCscRechitsChamberMinus13",            &nCscRechitsChamberMinus13,             "nCscRechitsChamberMinus13/I");
  tree_->Branch("nCscRechitsChamberMinus21",            &nCscRechitsChamberMinus21,             "nCscRechitsChamberMinus21/I");
  tree_->Branch("nCscRechitsChamberMinus22",            &nCscRechitsChamberMinus22,             "nCscRechitsChamberMinus22/I");
  tree_->Branch("nCscRechitsChamberMinus31",            &nCscRechitsChamberMinus31,             "nCscRechitsChamberMinus31/I");
  tree_->Branch("nCscRechitsChamberMinus32",            &nCscRechitsChamberMinus32,             "nCscRechitsChamberMinus32/I");
  tree_->Branch("nCscRechitsChamberMinus41",            &nCscRechitsChamberMinus41,             "nCscRechitsChamberMinus41/I");
  tree_->Branch("nCscRechitsChamberMinus42",            &nCscRechitsChamberMinus42,             "nCscRechitsChamberMinus42/I");

  tree_->Branch("nDTRechits",            &nDTRechits,             "nDTRechits/I");
  tree_->Branch("nDtRings",             &nDtRings, "nDtRings/I");
  tree_->Branch("nDTPositiveYRechits",             &nDTPositiveYRechits, "nDTPositiveYRechits/I");
  tree_->Branch("nDTNegativeYRechits",             &nDTNegativeYRechits, "nDTNegativeYRechits/I");

  tree_->Branch("nDTRechitsChamberMinus12",            &nDTRechitsChamberMinus12,             "nDTRechitsChamberMinus12/I");
  tree_->Branch("nDTRechitsChamberMinus11",            &nDTRechitsChamberMinus11,             "nDTRechitsChamberMinus11/I");
  tree_->Branch("nDTRechitsChamber10",            &nDTRechitsChamber10,             "nDTRechitsChamber10/I");
  tree_->Branch("nDTRechitsChamberPlus11",            &nDTRechitsChamberPlus11,             "nDTRechitsChamberPlus11/I");
  tree_->Branch("nDTRechitsChamberPlus12",            &nDTRechitsChamberPlus12,             "nDTRechitsChamberPlus12/I");
  tree_->Branch("nDTRechitsChamberMinus22",            &nDTRechitsChamberMinus22,             "nDTRechitsChamberMinus22/I");
  tree_->Branch("nDTRechitsChamberMinus21",            &nDTRechitsChamberMinus21,             "nDTRechitsChamberMinus21/I");
  tree_->Branch("nDTRechitsChamber20",            &nDTRechitsChamber20,             "nDTRechitsChamber20/I");
  tree_->Branch("nDTRechitsChamberPlus21",            &nDTRechitsChamberPlus21,             "nDTRechitsChamberPlus21/I");
  tree_->Branch("nDTRechitsChamberPlus22",            &nDTRechitsChamberPlus22,             "nDTRechitsChamberPlus22/I");
  tree_->Branch("nDTRechitsChamberMinus32",            &nDTRechitsChamberMinus32,             "nDTRechitsChamberMinus32/I");
  tree_->Branch("nDTRechitsChamberMinus31",            &nDTRechitsChamberMinus31,             "nDTRechitsChamberMinus31/I");
  tree_->Branch("nDTRechitsChamber30",            &nDTRechitsChamber30,             "nDTRechitsChamber30/I");

  tree_->Branch("nDTRechitsChamberPlus31",            &nDTRechitsChamberPlus31,             "nDTRechitsChamberPlus31/I");
  tree_->Branch("nDTRechitsChamberPlus32",            &nDTRechitsChamberPlus32,             "nDTRechitsChamberPlus32/I");
  tree_->Branch("nDTRechitsChamberMinus42",            &nDTRechitsChamberMinus42,             "nDTRechitsChamberMinus42/I");
  tree_->Branch("nDTRechitsChamberMinus41",            &nDTRechitsChamberMinus41,             "nDTRechitsChamberMinus41/I");
  tree_->Branch("nDTRechitsChamber40",            &nDTRechitsChamber40,             "nDTRechitsChamber40/I");
  tree_->Branch("nDTRechitsChamberPlus41",            &nDTRechitsChamberPlus41,             "nDTRechitsChamberPlus41/I");
  tree_->Branch("nDTRechitsChamberPlus42",            &nDTRechitsChamberPlus42,             "nDTRechitsChamberPlus42/I");

  tree_->Branch("nDtWheels25",             &nDtWheels25, "nDtWheels25/I");
  tree_->Branch("nDtStations25",             &nDtStations25, "nDtStations25/I");

  tree_->Branch("nDTRechitsWheelMinus2",             &nDTRechitsWheelMinus2, "nDTRechitsWheelMinus2/I");
  tree_->Branch("nDTRechitsWheelMinus1",             &nDTRechitsWheelMinus1, "nDTRechitsWheelMinus1/I");
  tree_->Branch("nDTRechitsWheel0",             &nDTRechitsWheel0, "nDTRechitsWheel0/I");
  tree_->Branch("nDTRechitsWheelPlus1",             &nDTRechitsWheelPlus1, "nDTRechitsWheelPlus1/I");
  tree_->Branch("nDTRechitsWheelPlus2",             &nDTRechitsWheelPlus2, "nDTRechitsWheelPlus2/I");

  tree_->Branch("nDTRechitsStation1",             &nDTRechitsStation1, "nDTRechitsStation1/I");
  tree_->Branch("nDTRechitsStation2",             &nDTRechitsStation2, "nDTRechitsStation2/I");
  tree_->Branch("nDTRechitsStation3",             &nDTRechitsStation3, "nDTRechitsStation3/I");
  tree_->Branch("nDTRechitsStation4",             &nDTRechitsStation4, "nDTRechitsStation4/I");





    // dt rechit cluster

    tree_->Branch("nDtRechitClusters",             &nDtRechitClusters, "nDtRechitClusters/I");
    tree_->Branch("dtRechitClusterX",             dtRechitClusterX,             "dtRechitClusterX[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterY",             dtRechitClusterY,             "dtRechitClusterY[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterZ",             dtRechitClusterZ,             "dtRechitClusterZ[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterWheel",             dtRechitClusterWheel,             "dtRechitClusterWheel[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterPhi",             dtRechitClusterPhi,             "dtRechitClusterPhi[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterEta",             dtRechitClusterEta,             "dtRechitClusterEta[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterJetVetoPt",             dtRechitClusterJetVetoPt,             "dtRechitClusterJetVetoPt[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterJetVetoEta",             dtRechitClusterJetVetoEta,             "dtRechitClusterJetVetoEta[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterJetVetoPhi",             dtRechitClusterJetVetoPhi,             "dtRechitClusterJetVetoPhi[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterJetVetoE",             dtRechitClusterJetVetoE,             "dtRechitClusterJetVetoE[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterMuonVetoPt",             dtRechitClusterMuonVetoPt,             "dtRechitClusterMuonVetoPt[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterMuonVetoE",             dtRechitClusterMuonVetoE,             "dtRechitClusterMuonVetoE[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterMuonVetoPhi",             dtRechitClusterMuonVetoPhi,             "dtRechitClusterMuonVetoPhi[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterMuonVetoEta",             dtRechitClusterMuonVetoEta,             "dtRechitClusterMuonVetoEta[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterMuonVetoLooseId",             dtRechitClusterMuonVetoLooseId,             "dtRechitClusterMuonVetoLooseId[nDtRechitClusters]/O");
    tree_->Branch("dtRechitClusterMuonVetoGlobal",             dtRechitClusterMuonVetoGlobal,             "dtRechitClusterMuonVetoGlobal[nDtRechitClusters]/O");
    tree_->Branch("dtRechitClusterSize",             dtRechitClusterSize,             "dtRechitClusterSize[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNStation10",             dtRechitClusterNStation10,             "dtRechitClusterNStation10[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterAvgStation10",             dtRechitClusterAvgStation10,             "dtRechitClusterAvgStation10[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterMaxStation",             dtRechitClusterMaxStation,             "dtRechitClusterMaxStation[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterMaxStationRatio",             dtRechitClusterMaxStationRatio,             "dtRechitClusterMaxStationRatio[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterNChamber",             dtRechitClusterNChamber,             "dtRechitClusterNChamber[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterMaxChamber",             dtRechitClusterMaxChamber,             "dtRechitClusterMaxChamber[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterMaxChamberRatio",             dtRechitClusterMaxChamberRatio,             "dtRechitClusterMaxChamberRatio[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterNSegmentStation1",             dtRechitClusterNSegmentStation1,             "dtRechitClusterNSegmentStation1[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNSegmentStation2",             dtRechitClusterNSegmentStation2,             "dtRechitClusterNSegmentStation2[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNSegmentStation3",             dtRechitClusterNSegmentStation3,             "dtRechitClusterNSegmentStation3[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterNSegmentStation4",             dtRechitClusterNSegmentStation4,             "dtRechitClusterNSegmentStation4[nDtRechitClusters]/I");
    tree_->Branch("dtRechitClusterMetEENoise_dPhi",             dtRechitClusterMetEENoise_dPhi,             "dtRechitClusterMetEENoise_dPhi[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gParticle_Id",             dtRechitCluster_match_gParticle_Id,             "dtRechitCluster_match_gParticle_Id[nDtRechitClusters]/I");
    tree_->Branch("dtRechitCluster_match_gParticle_Pt",             dtRechitCluster_match_gParticle_Pt,             "dtRechitCluster_match_gParticle_Pt[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gParticle_Eta",             dtRechitCluster_match_gParticle_Eta,             "dtRechitCluster_match_gParticle_Eta[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gParticle_Phi",             dtRechitCluster_match_gParticle_Phi,             "dtRechitCluster_match_gParticle_Phi[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gParticle_E",             dtRechitCluster_match_gParticle_E,             "dtRechitCluster_match_gParticle_E[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gParticle_Status",             dtRechitCluster_match_gParticle_Status,             "dtRechitCluster_match_gParticle_Status[nDtRechitClusters]/I");
    tree_->Branch("dtRechitCluster_match_gParticle_MotherId",             dtRechitCluster_match_gParticle_MotherId,             "dtRechitCluster_match_gParticle_MotherId[nDtRechitClusters]/I");
    tree_->Branch("dtRechitCluster_match_gParticle_deltaR",             dtRechitCluster_match_gParticle_deltaR,             "dtRechitCluster_match_gParticle_deltaR[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gLLP",             dtRechitCluster_match_gLLP,             "dtRechitCluster_match_gLLP[nDtRechitClusters]/O");
    tree_->Branch("dtRechitCluster_match_gLLP_minDeltaR",             dtRechitCluster_match_gLLP_minDeltaR,             "dtRechitCluster_match_gLLP_minDeltaR[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gLLP_eta",             dtRechitCluster_match_gLLP_eta, "dtRechitCluster_match_gLLP_eta[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gLLP_phi",             dtRechitCluster_match_gLLP_phi, "dtRechitCluster_match_gLLP_phi[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gLLP_decay_x",             dtRechitCluster_match_gLLP_decay_x, "dtRechitCluster_match_gLLP_decay_x[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gLLP_decay_y",             dtRechitCluster_match_gLLP_decay_y, "dtRechitCluster_match_gLLP_decay_y[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gLLP_decay_z",             dtRechitCluster_match_gLLP_decay_z, "dtRechitCluster_match_gLLP_decay_z[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gLLP_ctau",             dtRechitCluster_match_gLLP_ctau, "dtRechitCluster_match_gLLP_ctau[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gLLP_beta",             dtRechitCluster_match_gLLP_beta, "dtRechitCluster_match_gLLP_beta[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gLLP_csc",             dtRechitCluster_match_gLLP_csc, "dtRechitCluster_match_gLLP_csc[nDtRechitClusters]/O");
    tree_->Branch("dtRechitCluster_match_gLLP_dt",             dtRechitCluster_match_gLLP_dt, "dtRechitCluster_match_gLLP_dt[nDtRechitClusters]/O");
    tree_->Branch("dtRechitCluster_match_gLLP_e",             dtRechitCluster_match_gLLP_e, "dtRechitCluster_match_gLLP_e[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_gLLP_pt",             dtRechitCluster_match_gLLP_pt, "dtRechitCluster_match_gLLP_pt[nDtRechitClusters]/F");
    tree_->Branch("dtRechitClusterLep_dPhi",             dtRechitClusterLep_dPhi, "dtRechitClusterLep_dPhi[nDtRechitClusters]/F");
    tree_->Branch("dtRechitCluster_match_RPCBx_dPhi0p5",             dtRechitCluster_match_RPCBx_dPhi0p5,             "dtRechitCluster_match_RPCBx_dPhi0p5[nDtRechitClusters]/I");
    tree_->Branch("dtRechitCluster_match_RB1_0p4",             dtRechitCluster_match_RB1_0p4,             "dtRechitCluster_match_RB1_0p4[nDtRechitClusters]/I");
    tree_->Branch("dtRechitCluster_match_RB1_dPhi0p5",             dtRechitCluster_match_RB1_dPhi0p5,             "dtRechitCluster_match_RB1_dPhi0p5[nDtRechitClusters]/I");
    tree_->Branch("dtRechitCluster_match_RPChits_dPhi0p5",             dtRechitCluster_match_RPChits_dPhi0p5,             "dtRechitCluster_match_RPChits_dPhi0p5[nDtRechitClusters]/I");

    tree_->Branch("dtRechitCluster_match_MB1hits_0p4",             dtRechitCluster_match_MB1hits_0p4,             "dtRechitCluster_match_MB1hits_0p4[nDtRechitClusters]/I");
    tree_->Branch("dtRechitCluster_match_MB1hits_0p5",             dtRechitCluster_match_MB1hits_0p5,             "dtRechitCluster_match_MB1hits_0p5[nDtRechitClusters]/I");
    tree_->Branch("dtRechitCluster_match_MB1hits_cosmics_plus",             dtRechitCluster_match_MB1hits_cosmics_plus,             "dtRechitCluster_match_MB1hits_cosmics_plus[nDtRechitClusters]/I");
    tree_->Branch("dtRechitCluster_match_MB1hits_cosmics_minus",             dtRechitCluster_match_MB1hits_cosmics_minus,             "dtRechitCluster_match_MB1hits_cosmics_minus[nDtRechitClusters]/I");




    tree_->Branch("nCscRechitClusters3",             &nCscRechitClusters3, "nCscRechitClusters3/I");
    tree_->Branch("cscRechitCluster3_match_Me1112_0p4",             cscRechitCluster3_match_Me1112_0p4,             "cscRechitCluster3_match_Me1112_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me1112_0p6",             cscRechitCluster3_match_Me1112_0p6,             "cscRechitCluster3_match_Me1112_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me1112_0p8",             cscRechitCluster3_match_Me1112_0p8,             "cscRechitCluster3_match_Me1112_0p8[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me11_0p4",             cscRechitCluster3_match_Me11_0p4,             "cscRechitCluster3_match_Me11_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me11_0p6",             cscRechitCluster3_match_Me11_0p6,             "cscRechitCluster3_match_Me11_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me11_0p8",             cscRechitCluster3_match_Me11_0p8,             "cscRechitCluster3_match_Me11_0p8[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me12_0p4",             cscRechitCluster3_match_Me12_0p4,             "cscRechitCluster3_match_Me12_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me12_0p6",             cscRechitCluster3_match_Me12_0p6,             "cscRechitCluster3_match_Me12_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_Me12_0p8",             cscRechitCluster3_match_Me12_0p8,             "cscRechitCluster3_match_Me12_0p8[nCscRechitClusters3]/I");

    tree_->Branch("cscRechitCluster3_match_gParticle_id",             cscRechitCluster3_match_gParticle_id,             "cscRechitCluster3_match_gParticle_id[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_gParticle",             cscRechitCluster3_match_gParticle,             "cscRechitCluster3_match_gParticle[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3_match_gParticle_minDeltaR",             cscRechitCluster3_match_gParticle_minDeltaR,             "cscRechitCluster3_match_gParticle_minDeltaR[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gParticle_index",             cscRechitCluster3_match_gParticle_index,             "cscRechitCluster3_match_gParticle_index[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_gParticle_eta",             cscRechitCluster3_match_gParticle_eta,             "cscRechitCluster3_match_gParticle_eta[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gParticle_phi",             cscRechitCluster3_match_gParticle_phi,             "cscRechitCluster3_match_gParticle_phi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gParticle_E",             cscRechitCluster3_match_gParticle_E,             "cscRechitCluster3_match_gParticle_E[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gParticle_pt",             cscRechitCluster3_match_gParticle_pt,             "cscRechitCluster3_match_gParticle_pt[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gParticle_MotherId",             cscRechitCluster3_match_gParticle_MotherId,             "cscRechitCluster3_match_gParticle_MotherId[nCscRechitClusters3]/I");

    tree_->Branch("cscRechitCluster3_match_gLLP",             cscRechitCluster3_match_gLLP,             "cscRechitCluster3_match_gLLP[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3_match_gLLP_minDeltaR",             cscRechitCluster3_match_gLLP_minDeltaR,             "cscRechitCluster3_match_gLLP_minDeltaR[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_index",             cscRechitCluster3_match_gLLP_index,             "cscRechitCluster3_match_gLLP_index[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_gLLP_eta",             cscRechitCluster3_match_gLLP_eta, "cscRechitCluster3_match_gLLP_eta[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_phi",             cscRechitCluster3_match_gLLP_phi, "cscRechitCluster3_match_gLLP_phi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_decay_r",             cscRechitCluster3_match_gLLP_decay_r, "cscRechitCluster3_match_gLLP_decay_r[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_decay_x",             cscRechitCluster3_match_gLLP_decay_x, "cscRechitCluster3_match_gLLP_decay_x[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_decay_y",             cscRechitCluster3_match_gLLP_decay_y, "cscRechitCluster3_match_gLLP_decay_y[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_decay_z",             cscRechitCluster3_match_gLLP_decay_z, "cscRechitCluster3_match_gLLP_decay_z[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_ctau",             cscRechitCluster3_match_gLLP_ctau, "cscRechitCluster3_match_gLLP_ctau[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_beta",             cscRechitCluster3_match_gLLP_beta, "cscRechitCluster3_match_gLLP_beta[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_csc",             cscRechitCluster3_match_gLLP_csc, "cscRechitCluster3_match_gLLP_csc[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3_match_gLLP_e",             cscRechitCluster3_match_gLLP_e, "cscRechitCluster3_match_gLLP_e[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_pt",             cscRechitCluster3_match_gLLP_pt, "cscRechitCluster3_match_gLLP_pt[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_lepdPhi",             cscRechitCluster3_match_gLLP_lepdPhi, "cscRechitCluster3_match_gLLP_lepdPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_daughter0_deltaR",             cscRechitCluster3_match_gLLP_daughter0_deltaR, "cscRechitCluster3_match_gLLP_daughter0_deltaR[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_daughter1_deltaR",             cscRechitCluster3_match_gLLP_daughter1_deltaR, "cscRechitCluster3_match_gLLP_daughter1_deltaR[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_daughter2_deltaR",             cscRechitCluster3_match_gLLP_daughter2_deltaR, "cscRechitCluster3_match_gLLP_daughter2_deltaR[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_gLLP_daughter3_deltaR",             cscRechitCluster3_match_gLLP_daughter3_deltaR, "cscRechitCluster3_match_gLLP_daughter3_deltaR[nCscRechitClusters3]/F");

    tree_->Branch("cscRechitCluster3X",             cscRechitCluster3X,             "cscRechitCluster3X[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3Y",             cscRechitCluster3Y,             "cscRechitCluster3Y[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3Z",             cscRechitCluster3Z,             "cscRechitCluster3Z[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3Time",             cscRechitCluster3Time,             "cscRechitCluster3Time[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3TimeWire",             cscRechitCluster3TimeWire,             "cscRechitCluster3TimeWire[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3TimeWirePruned",             cscRechitCluster3TimeWirePruned,             "cscRechitCluster3TimeWirePruned[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3TimeTotal",             cscRechitCluster3TimeTotal,             "cscRechitCluster3TimeTotal[nCscRechitClusters3]/F");

    tree_->Branch("cscRechitCluster3TimeSpread",             cscRechitCluster3TimeSpread,             "cscRechitCluster3TimeSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3TimeTotalSpread",             cscRechitCluster3TimeTotalSpread,             "cscRechitCluster3TimeTotalSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3TimeTotalSpreadPruned",             cscRechitCluster3TimeTotalSpreadPruned,             "cscRechitCluster3TimeTotalSpreadPruned[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3TimeWireSpread",             cscRechitCluster3TimeWireSpread,             "cscRechitCluster3TimeWireSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3GenMuonDeltaR",             cscRechitCluster3GenMuonDeltaR,             "cscRechitCluster3GenMuonDeltaR[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3XYSpread",             cscRechitCluster3XYSpread,             "cscRechitCluster3XYSpread[nCscRechitClusters3]/F");

    tree_->Branch("cscRechitCluster3MajorAxis",             cscRechitCluster3MajorAxis,             "cscRechitCluster3MajorAxis[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MinorAxis",             cscRechitCluster3MinorAxis,             "cscRechitCluster3MinorAxis[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3EtaPhiSpread",             cscRechitCluster3EtaPhiSpread,             "cscRechitCluster3EtaPhiSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3PhiSpread",             cscRechitCluster3PhiSpread,             "cscRechitCluster3PhiSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3EtaSpread",             cscRechitCluster3EtaSpread,             "cscRechitCluster3EtaSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3DeltaRSpread",             cscRechitCluster3DeltaRSpread,             "cscRechitCluster3DeltaRSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3XSpread",             cscRechitCluster3XSpread,             "cscRechitCluster3XSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3RSpread",             cscRechitCluster3RSpread,             "cscRechitCluster3RSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3YSpread",             cscRechitCluster3YSpread,             "cscRechitCluster3YSpread[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3ZSpread",             cscRechitCluster3ZSpread,             "cscRechitCluster3ZSpread[nCscRechitClusters3]/F");


    tree_->Branch("cscRechitCluster3Phi",             cscRechitCluster3Phi,             "cscRechitCluster3Phi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3Eta",             cscRechitCluster3Eta,             "cscRechitCluster3Eta[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoPt",             cscRechitCluster3JetVetoPt,             "cscRechitCluster3JetVetoPt[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoEta",             cscRechitCluster3JetVetoEta,             "cscRechitCluster3JetVetoEta[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoPhi",             cscRechitCluster3JetVetoPhi,             "cscRechitCluster3JetVetoPhi[nCscRechitClusters3]/F");

    tree_->Branch("cscRechitCluster3JetVetoE",             cscRechitCluster3JetVetoE,             "cscRechitCluster3JetVetoE[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3GenJetVetoPt",             cscRechitCluster3GenJetVetoPt,             "cscRechitCluster3GenJetVetoPt[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3GenJetVetoE",             cscRechitCluster3GenJetVetoE,             "cscRechitCluster3GenJetVetoE[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoPt",             cscRechitCluster3MuonVetoPt,             "cscRechitCluster3MuonVetoPt[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoE",             cscRechitCluster3MuonVetoE,             "cscRechitCluster3MuonVetoE[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoPhi",             cscRechitCluster3MuonVetoPhi,             "cscRechitCluster3MuonVetoPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoEta",             cscRechitCluster3MuonVetoEta,             "cscRechitCluster3MuonVetoEta[nCscRechitClusters3]/F");

    tree_->Branch("cscRechitCluster3JetVetoElectronEnergyFraction",             cscRechitCluster3JetVetoElectronEnergyFraction,             "cscRechitCluster3JetVetoElectronEnergyFraction[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoPhotonEnergyFraction",             cscRechitCluster3JetVetoPhotonEnergyFraction,             "cscRechitCluster3JetVetoPhotonEnergyFraction[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoChargedHadronEnergyFraction",             cscRechitCluster3JetVetoChargedHadronEnergyFraction,             "cscRechitCluster3JetVetoChargedHadronEnergyFraction[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoNeutralHadronEnergyFraction",             cscRechitCluster3JetVetoNeutralHadronEnergyFraction,             "cscRechitCluster3JetVetoNeutralHadronEnergyFraction[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoMuonEnergyFraction",             cscRechitCluster3JetVetoMuonEnergyFraction,             "cscRechitCluster3JetVetoMuonEnergyFraction[nCscRechitClusters3]/F");


    tree_->Branch("cscRechitCluster3JetVetoPt_0p6",             cscRechitCluster3JetVetoPt_0p6,        "cscRechitCluster3JetVetoPt_0p6[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoPt_0p8",             cscRechitCluster3JetVetoPt_0p8,        "cscRechitCluster3JetVetoPt_0p8[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoE_0p6",             cscRechitCluster3JetVetoE_0p6,        "cscRechitCluster3JetVetoE_0p6[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3JetVetoE_0p8",             cscRechitCluster3JetVetoE_0p8,        "cscRechitCluster3JetVetoE_0p8[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoPt_0p6",             cscRechitCluster3MuonVetoPt_0p6,        "cscRechitCluster3MuonVetoPt_0p6[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoPt_0p8",             cscRechitCluster3MuonVetoPt_0p8,        "cscRechitCluster3MuonVetoPt_0p8[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoE_0p6",             cscRechitCluster3MuonVetoE_0p6,        "cscRechitCluster3MuonVetoE_0p6[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoE_0p8",             cscRechitCluster3MuonVetoE_0p8,        "cscRechitCluster3MuonVetoE_0p8[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoLooseIso",             cscRechitCluster3MuonVetoLooseIso,             "cscRechitCluster3MuonVetoLooseIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3MuonVetoTightIso",             cscRechitCluster3MuonVetoTightIso,             "cscRechitCluster3MuonVetoTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3MuonVetoVTightIso",             cscRechitCluster3MuonVetoVTightIso,             "cscRechitCluster3MuonVetoVTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3MuonVetoVVTightIso",             cscRechitCluster3MuonVetoVVTightIso,             "cscRechitCluster3MuonVetoVVTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3MuonVetoTightId",             cscRechitCluster3MuonVetoTightId,             "cscRechitCluster3MuonVetoTightId[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3MuonVetoLooseId",             cscRechitCluster3MuonVetoLooseId,             "cscRechitCluster3MuonVetoLooseId[nCscRechitClusters3]/O");


    tree_->Branch("cscRechitCluster3MuonVetoIso",             cscRechitCluster3MuonVetoIso,             "cscRechitCluster3MuonVetoIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3IsoMuonVetoPt",             cscRechitCluster3IsoMuonVetoPt,             "cscRechitCluster3IsoMuonVetoPt[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3IsoMuonVetoE",             cscRechitCluster3IsoMuonVetoE,             "cscRechitCluster3IsoMuonVetoE[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3IsoMuonVetoPhi",             cscRechitCluster3IsoMuonVetoPhi,             "cscRechitCluster3IsoMuonVetoPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3IsoMuonVetoEta",             cscRechitCluster3IsoMuonVetoEta,             "cscRechitCluster3IsoMuonVetoEta[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3GenMuonVetoPt",             cscRechitCluster3GenMuonVetoPt,             "cscRechitCluster3GenMuonVetoPt[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MuonVetoType",             cscRechitCluster3MuonVetoType,             "cscRechitCluster3MuonVetoType[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3GenMuonVetoE",             cscRechitCluster3GenMuonVetoE,             "cscRechitCluster3GenMuonVetoE[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3GenMuonVetoProdX",             cscRechitCluster3GenMuonVetoProdX,             "cscRechitCluster3GenMuonVetoProdX[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3GenMuonVetoProdY",             cscRechitCluster3GenMuonVetoProdY,             "cscRechitCluster3GenMuonVetoProdY[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3GenMuonVetoProdZ",             cscRechitCluster3GenMuonVetoProdZ,             "cscRechitCluster3GenMuonVetoProdZ[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3GenMuonVetoLLPDist",             cscRechitCluster3GenMuonVetoLLPDist,             "cscRechitCluster3GenMuonVetoLLPDist[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3GenMuonVetoLLPIndex",             cscRechitCluster3GenMuonVetoLLPIndex,             "cscRechitCluster3GenMuonVetoLLPIndex[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3ZLep1",             cscRechitCluster3ZLep1,             "cscRechitCluster3ZLep1[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep2",             cscRechitCluster3ZLep2,             "cscRechitCluster3ZLep2[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep1Tag",             cscRechitCluster3ZLep1Tag,             "cscRechitCluster3ZLep1Tag[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep2Tag",             cscRechitCluster3ZLep2Tag,             "cscRechitCluster3ZLep2Tag[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep1Id",             cscRechitCluster3ZLep1Id,             "cscRechitCluster3ZLep1Id[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3ZLep2Id",             cscRechitCluster3ZLep2Id,             "cscRechitCluster3ZLep2Id[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3ZLep1LooseIso",             cscRechitCluster3ZLep1LooseIso,             "cscRechitCluster3ZLep1LooseIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep1TightIso",             cscRechitCluster3ZLep1TightIso,             "cscRechitCluster3ZLep1TightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep1VTightIso",             cscRechitCluster3ZLep1VTightIso,             "cscRechitCluster3ZLep1VTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep1VVTightIso",             cscRechitCluster3ZLep1VVTightIso,             "cscRechitCluster3ZLep1VVTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep1TightId",             cscRechitCluster3ZLep1TightId,             "cscRechitCluster3ZLep1TightId[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep2LooseIso",             cscRechitCluster3ZLep2LooseIso,             "cscRechitCluster3ZLep2LooseIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep2TightIso",             cscRechitCluster3ZLep2TightIso,             "cscRechitCluster3ZLep2TightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep2VTightIso",             cscRechitCluster3ZLep2VTightIso,             "cscRechitCluster3ZLep2VTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep2VVTightIso",             cscRechitCluster3ZLep2VVTightIso,             "cscRechitCluster3ZLep2VVTightIso[nCscRechitClusters3]/O");
    tree_->Branch("cscRechitCluster3ZLep2TightId",             cscRechitCluster3ZLep2TightId,             "cscRechitCluster3ZLep2TightId[nCscRechitClusters3]/O");


    tree_->Branch("cscRechitCluster3Size",             cscRechitCluster3Size,             "cscRechitCluster3Size[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3Me11Ratio",             cscRechitCluster3Me11Ratio,             "cscRechitCluster3Me11Ratio[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3Me12Ratio",             cscRechitCluster3Me12Ratio,             "cscRechitCluster3Me12Ratio[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3NStation",             cscRechitCluster3NStation,             "cscRechitCluster3NStation[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NStation5",             cscRechitCluster3NStation5,             "cscRechitCluster3NStation5[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NStation10",             cscRechitCluster3NStation10,             "cscRechitCluster3NStation10[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NStation10perc",             cscRechitCluster3NStation10perc,             "cscRechitCluster3NStation10perc[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3AvgStation",             cscRechitCluster3AvgStation,             "cscRechitCluster3AvgStation[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3AvgStation5",             cscRechitCluster3AvgStation5,             "cscRechitCluster3AvgStation5[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3AvgStation10",             cscRechitCluster3AvgStation10,             "cscRechitCluster3AvgStation10[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3AvgStation10perc",             cscRechitCluster3AvgStation10perc,             "cscRechitCluster3AvgStation10perc[nCscRechitClusters3]/F");


    tree_->Branch("cscRechitCluster3MaxStation",             cscRechitCluster3MaxStation,             "cscRechitCluster3MaxStation[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3MaxStationRatio",             cscRechitCluster3MaxStationRatio,             "cscRechitCluster3MaxStationRatio[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3NChamber",             cscRechitCluster3NChamber,             "cscRechitCluster3NChamber[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3MaxChamber",             cscRechitCluster3MaxChamber,             "cscRechitCluster3MaxChamber[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3MaxChamberRatio",             cscRechitCluster3MaxChamberRatio,             "cscRechitCluster3MaxChamberRatio[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus11",             cscRechitCluster3NRechitChamberPlus11,             "cscRechitCluster3NRechitChamberPlus11[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus12",             cscRechitCluster3NRechitChamberPlus12,             "cscRechitCluster3NRechitChamberPlus12[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus13",             cscRechitCluster3NRechitChamberPlus13,             "cscRechitCluster3NRechitChamberPlus13[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus21",             cscRechitCluster3NRechitChamberPlus21,             "cscRechitCluster3NRechitChamberPlus21[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus22",             cscRechitCluster3NRechitChamberPlus22,             "cscRechitCluster3NRechitChamberPlus22[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus31",             cscRechitCluster3NRechitChamberPlus31,             "cscRechitCluster3NRechitChamberPlus31[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus32",             cscRechitCluster3NRechitChamberPlus32,             "cscRechitCluster3NRechitChamberPlus32[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus41",             cscRechitCluster3NRechitChamberPlus41,             "cscRechitCluster3NRechitChamberPlus41[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberPlus42",             cscRechitCluster3NRechitChamberPlus42,             "cscRechitCluster3NRechitChamberPlus42[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus11",             cscRechitCluster3NRechitChamberMinus11,             "cscRechitCluster3NRechitChamberMinus11[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus12",             cscRechitCluster3NRechitChamberMinus12,             "cscRechitCluster3NRechitChamberMinus12[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus13",             cscRechitCluster3NRechitChamberMinus13,             "cscRechitCluster3NRechitChamberMinus13[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus21",             cscRechitCluster3NRechitChamberMinus21,             "cscRechitCluster3NRechitChamberMinus21[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus22",             cscRechitCluster3NRechitChamberMinus22,             "cscRechitCluster3NRechitChamberMinus22[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus31",             cscRechitCluster3NRechitChamberMinus31,             "cscRechitCluster3NRechitChamberMinus31[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus32",             cscRechitCluster3NRechitChamberMinus32,             "cscRechitCluster3NRechitChamberMinus32[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus41",             cscRechitCluster3NRechitChamberMinus41,             "cscRechitCluster3NRechitChamberMinus41[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3NRechitChamberMinus42",             cscRechitCluster3NRechitChamberMinus42,             "cscRechitCluster3NRechitChamberMinus42[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3Met_dPhi",             cscRechitCluster3Met_dPhi,             "cscRechitCluster3Met_dPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MetXYCorr_dPhi",             cscRechitCluster3MetXYCorr_dPhi,             "cscRechitCluster3MetXYCorr_dPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_cscRechits_0p4",             cscRechitCluster3_match_cscRechits_0p4,             "cscRechitCluster3_match_cscRechits_0p4[nCscRechitClusters3]/I");


    tree_->Branch("cscRechitCluster3MetHEM_dPhi",             cscRechitCluster3MetHEM_dPhi,             "cscRechitCluster3MetHEM_dPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MetHEMXYCorr_dPhi",             cscRechitCluster3MetHEMXYCorr_dPhi,             "cscRechitCluster3MetHEMXYCorr_dPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MetEENoise_dPhi",             cscRechitCluster3MetEENoise_dPhi,             "cscRechitCluster3MetEENoise_dPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MetEENoiseXYCorr_dPhi",             cscRechitCluster3MetEENoiseXYCorr_dPhi,             "cscRechitCluster3MetEENoiseXYCorr_dPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MetJesUp_dPhi",             cscRechitCluster3MetJesUp_dPhi,             "cscRechitCluster3MetJesUp_dPhi[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3MetJesDown_dPhi",             cscRechitCluster3MetJesDown_dPhi,             "cscRechitCluster3MetJesDown_dPhi[nCscRechitClusters3]/F");


    tree_->Branch("cscRechitCluster3_match_cscSeg_0p4",             cscRechitCluster3_match_cscSeg_0p4,             "cscRechitCluster3_match_cscSeg_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_ME11Seg_0p4",             cscRechitCluster3_match_ME11Seg_0p4,             "cscRechitCluster3_match_ME11Seg_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_ME12Seg_0p4",             cscRechitCluster3_match_ME12Seg_0p4,             "cscRechitCluster3_match_ME12Seg_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_cscSeg_0p6",             cscRechitCluster3_match_cscSeg_0p6,             "cscRechitCluster3_match_cscSeg_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_ME11Seg_0p6",             cscRechitCluster3_match_ME11Seg_0p6,             "cscRechitCluster3_match_ME11Seg_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_ME12Seg_0p6",             cscRechitCluster3_match_ME12Seg_0p6,             "cscRechitCluster3_match_ME12Seg_0p6[nCscRechitClusters3]/I");

    tree_->Branch("cscRechitCluster3_match_dtRechits_phi0p2",             cscRechitCluster3_match_dtRechits_phi0p2,             "cscRechitCluster3_match_dtRechits_phi0p2[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_dtRechits_0p4",             cscRechitCluster3_match_dtRechits_0p4,             "cscRechitCluster3_match_dtRechits_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_MB1_0p4",             cscRechitCluster3_match_MB1_0p4,             "cscRechitCluster3_match_MB1_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_dtRechits_0p6",             cscRechitCluster3_match_dtRechits_0p6,             "cscRechitCluster3_match_dtRechits_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_MB1_0p6",             cscRechitCluster3_match_MB1_0p6,             "cscRechitCluster3_match_MB1_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_dtSeg_0p4",             cscRechitCluster3_match_dtSeg_0p4,             "cscRechitCluster3_match_dtSeg_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_MB1Seg_0p4",             cscRechitCluster3_match_MB1Seg_0p4,             "cscRechitCluster3_match_MB1Seg_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_dtSeg_0p6",             cscRechitCluster3_match_dtSeg_0p6,             "cscRechitCluster3_match_dtSeg_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_MB1Seg_0p6",             cscRechitCluster3_match_MB1Seg_0p6,             "cscRechitCluster3_match_MB1Seg_0p6[nCscRechitClusters3]/I");

    tree_->Branch("cscRechitCluster3_match_RB1_0p4",             cscRechitCluster3_match_RB1_0p4,             "cscRechitCluster3_match_RB1_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_RE12_0p4",             cscRechitCluster3_match_RE12_0p4,             "cscRechitCluster3_match_RE12_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_RB1_0p6",             cscRechitCluster3_match_RB1_0p6,             "cscRechitCluster3_match_RB1_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_RE12_0p6",             cscRechitCluster3_match_RE12_0p6,             "cscRechitCluster3_match_RE12_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_highEta_0p4",             cscRechitCluster3_match_highEta_0p4,             "cscRechitCluster3_match_highEta_0p4[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_highEta_0p6",             cscRechitCluster3_match_highEta_0p6,             "cscRechitCluster3_match_highEta_0p6[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_highEta_0p8",             cscRechitCluster3_match_highEta_0p8,             "cscRechitCluster3_match_highEta_0p8[nCscRechitClusters3]/I");
    tree_->Branch("cscRechitCluster3_match_cluster_dR",             cscRechitCluster3_match_cluster_dR,             "cscRechitCluster3_match_cluster_dR[nCscRechitClusters3]/F");
    tree_->Branch("cscRechitCluster3_match_cluster_index",             cscRechitCluster3_match_cluster_index,             "cscRechitCluster3_match_cluster_index[nCscRechitClusters3]/I");


  //gLLP branches
  tree_->Branch("gLLP_multiplicity",          gLLP_multiplicity,          "gLLP_multiplicity[1]/I");
  tree_->Branch("gLLP_eta",          gLLP_eta,          "gLLP_eta[1]/F");
  tree_->Branch("gLLP_phi",          gLLP_phi,          "gLLP_phi[1]/F");
  tree_->Branch("gLLP_csc",          gLLP_csc,          "gLLP_csc[1]/F");
  tree_->Branch("gLLP_dt",          gLLP_dt,          "gLLP_dt[1]/F");
  tree_->Branch("gLLP_beta",          gLLP_beta,          "gLLP_beta[1]/F");
  tree_->Branch("gLLP_e",          gLLP_e,          "gLLP_e[1]/F");
  tree_->Branch("gLLP_pt",          gLLP_pt,          "gLLP_pt[1]/F");
  tree_->Branch("gLLP_lepdPhi",          gLLP_lepdPhi,          "gLLP_lepdPhi[1]/F");

  tree_->Branch("gLLP_ctau",          gLLP_ctau,          "gLLP_ctau[1]/F");


  tree_->Branch("gLLP_decay_vertex_r",          gLLP_decay_vertex_r,          "gLLP_decay_vertex_r[1]/F");
  tree_->Branch("gLLP_decay_vertex_x",          gLLP_decay_vertex_x,          "gLLP_decay_vertex_x[1]/F");
  tree_->Branch("gLLP_decay_vertex_y",          gLLP_decay_vertex_y,          "gLLP_decay_vertex_y[1]/F");
  tree_->Branch("gLLP_decay_vertex_z",          gLLP_decay_vertex_z,          "gLLP_decay_vertex_z[1]/F");
  tree_->Branch("gLLP_daughter_deltaR",          gLLP_daughter_deltaR,          "gLLP_daughter_deltaR[1]/F");

  tree_->Branch("gLLP_daughter_pt",          gLLP_daughter_pt,          "gLLP_daughter_pt[3]/F");
  tree_->Branch("gLLP_daughter_id",          gLLP_daughter_id,          "gLLP_daughter_id[3]/I");
  tree_->Branch("gLLP_daughter_eta",          gLLP_daughter_eta,          "gLLP_daughter_eta[3]/F");
  tree_->Branch("gLLP_daughter_phi",          gLLP_daughter_phi,          "gLLP_daughter_phi[3]/F");
  tree_->Branch("gLLP_daughter_e",          gLLP_daughter_e,          "gLLP_daughter_e[3]/F");
  tree_->Branch("gLLP_daughter_mass",          gLLP_daughter_mass,          "gLLP_daughter_mass[3]/F");


  //leptons
  tree_->Branch("nMuons",  &nMuons, "nMuons/I");
  tree_->Branch("muonPt",     muonPt,     "muonPt[nMuons]/F");
  tree_->Branch("muonEta",    muonEta,    "muonEta[nMuons]/F");
  tree_->Branch("muonPhi",    muonPhi,    "muonPhi[nMuons]/F");


  //leptons
  tree_->Branch("nLeptons",  &nLeptons, "nLeptons/I");
  tree_->Branch("lepE",      lepE,      "lepE[nLeptons]/F");
  tree_->Branch("lepPt",     lepPt,     "lepPt[nLeptons]/F");
  tree_->Branch("lepEta",    lepEta,    "lepEta[nLeptons]/F");
  tree_->Branch("lepPhi",    lepPhi,    "lepPhi[nLeptons]/F");
  tree_->Branch("lepPdgId",  lepPdgId,  "lepPdgId[nLeptons]/I");
  tree_->Branch("lepDZ",     lepDZ,     "lepDZ[nLeptons]/F");
  tree_->Branch("lepPassId", lepPassId, "lepPassId[nLeptons]/O");
  tree_->Branch("lepFromZ", lepFromZ, "lepFromZ[nLeptons]/O");
  tree_->Branch("lepEff", lepEff, "lepEff[nLeptons]/F");
  tree_->Branch("lepSF", lepSF, "lepSF[nLeptons]/F");
  tree_->Branch("lepTriggerSF", lepTriggerSF, "lepTriggerSF[nLeptons]/F");
  tree_->Branch("lepTightIdSF", lepTightIdSF, "lepTightIdSF[nLeptons]/F");
  tree_->Branch("lepLooseIdSF", lepLooseIdSF, "lepLooseIdSF[nLeptons]/F");
  tree_->Branch("lepTightIsoSF", lepTightIsoSF, "lepTightIsoSF[nLeptons]/F");
  tree_->Branch("lepLooseIsoSF", lepLooseIsoSF, "lepLooseIsoSF[nLeptons]/F");
  tree_->Branch("lepTriggerMCEfficiency", lepTriggerMCEfficiency, "lepTriggerMCEfficiency[nLeptons]/F");
  tree_->Branch("lepTightIdMCEfficiency", lepTightIdMCEfficiency, "lepTightIdMCEfficiency[nLeptons]/F");
  tree_->Branch("lepLooseIdMCEfficiency", lepLooseIdMCEfficiency, "lepLooseIdMCEfficiency[nLeptons]/F");
  tree_->Branch("lepTightIsoMCEfficiency", lepTightIsoMCEfficiency, "lepTightIsoMCEfficiency[nLeptons]/F");
  tree_->Branch("lepLooseIsoMCEfficiency", lepLooseIsoMCEfficiency, "lepLooseIsoMCEfficiency[nLeptons]/F");
  tree_->Branch("lepTag", lepTag, "lepTag[nLeptons]/O");



  tree_->Branch("lepPassLooseIso", lepPassLooseIso, "lepPassLooseIso[nLeptons]/O");
  tree_->Branch("lepPassTightIso", lepPassTightIso, "lepPassTightIso[nLeptons]/O");
  tree_->Branch("lepPassVTightIso", lepPassVTightIso, "lepPassVTightIso[nLeptons]/O");
  tree_->Branch("lepPassVVTightIso", lepPassVVTightIso, "lepPassVVTightIso[nLeptons]/O");

  tree_->Branch("lepPassVetoId", lepPassVetoId, "lepPassVetoId[nLeptons]/O");

  // tree_->Branch("lepLoosePassId", lepLoosePassId, "lepLoosePassId[nLeptons]/O");
  // tree_->Branch("lepMediumPassId", lepMediumPassId, "lepMediumPassId[nLeptons]/O");
  // tree_->Branch("lepTightPassId", lepTightPassId, "lepTightPassId[nLeptons]/O");


  tree_->Branch("MT",      &MT,  "MT/F");
  //Z-candidate
  tree_->Branch("ZMass1",      &ZMass1,  "ZMass1/F");

  tree_->Branch("ZMass",      &ZMass,  "ZMass/F");
  tree_->Branch("ZPt",        &ZPt,    "ZPt/F");
  tree_->Branch("ZEta",       &ZEta,   "ZEta/F");
  tree_->Branch("ZPhi",       &ZPhi,   "ZPhi/F");
  tree_->Branch("ZleptonIndex1", &ZleptonIndex1, "ZleptonIndex1/I");
  tree_->Branch("ZleptonIndex2", &ZleptonIndex2, "ZleptonIndex2/I");

  //jets
  tree_->Branch("nJets",     &nJets,    "nJets/I");
  tree_->Branch("jetE",      jetE,      "jetE[nJets]/F");
  tree_->Branch("jetPt",     jetPt,     "jetPt[nJets]/F");
  tree_->Branch("jetEta",    jetEta,    "jetEta[nJets]/F");
  tree_->Branch("jetPhi",    jetPhi,    "jetPhi[nJets]/F");
  tree_->Branch("jetTime",   jetTime,   "jetTime[nJets]/F");
  tree_->Branch("jetPassId", jetPassId, "jetPassId[nJets]/O");

  tree_->Branch("jetPtJESUp",     jetPtJESUp,     "jetPtJESUp[nJets]/F");
  tree_->Branch("jetPtJESDown",     jetPtJESDown,     "jetPtJESDown[nJets]/F");
  tree_->Branch("jetEJESUp",      jetEJESUp,      "jetEJESUp[nJets]/F");
  tree_->Branch("jetEJESDown",      jetEJESDown,      "jetEJESDown[nJets]/F");
  tree_->Branch("JecUnc",      JecUnc,      "JecUnc[nJets]/F");


  tree_->Branch("jet_match_genJet_pt", jet_match_genJet_pt, "jet_match_genJet_pt[nJets]/F");
  tree_->Branch("jet_match_genJet_index", jet_match_genJet_index, "jet_match_genJet_index[nJets]/I");
  tree_->Branch("jet_match_genJet_minDeltaR", jet_match_genJet_minDeltaR, "jet_match_genJet_minDeltaR[nJets]/F");

  // tree_->Branch("ecalNRechits",   ecalNRechits,   "ecalNRechits[nJets]/F");
  // tree_->Branch("ecalRechitE", ecalRechitE, "ecalRechitE[nJets]/F");
  // tree_->Branch("jetLoosePassId", jetLoosePassId, "jetLoosePassId[nJets]/O");
  tree_->Branch("jetTightPassId", jetTightPassId, "jetTightPassId[nJets]/O");
  tree_->Branch("HLTDecision", HLTDecision, "HLTDecision[982]/O"); //hardcoded
  tree_->Branch("SingleMuonTrigger", SingleMuonTrigger, "SingleMuonTrigger/O");
  tree_->Branch("SingleEleTrigger", SingleEleTrigger, "SingleEleTrigger/O");
  tree_->Branch("SingleLepTrigger", SingleLepTrigger, "SingleLepTrigger/O");

  tree_->Branch("jetMuonEnergyFraction",   jetMuonEnergyFraction,   "jetMuonEnergyFraction[nJets]/F");
  tree_->Branch("jetElectronEnergyFraction",   jetElectronEnergyFraction,   "jetElectronEnergyFraction[nJets]/F");
  tree_->Branch("jetPhotonEnergyFraction",   jetPhotonEnergyFraction,   "jetPhotonEnergyFraction[nJets]/F");
  tree_->Branch("jetChargedHadronEnergyFraction",   jetChargedHadronEnergyFraction,   "jetChargedHadronEnergyFraction[nJets]/F");
  tree_->Branch("jetNeutralHadronEnergyFraction",   jetNeutralHadronEnergyFraction,   "jetNeutralHadronEnergyFraction[nJets]/F");
};
