// -*- C++ -*-
//
//
// Fill a tree for selected candidates.
//


// system include files
#include <cassert>
#include <memory>
#include <cmath>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/Framework/interface/LuminosityBlock.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <DataFormats/Common/interface/TriggerResults.h>
#include <FWCore/Common/interface/TriggerNames.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>

#include <DataFormats/Common/interface/View.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/MET.h>
#include <DataFormats/METReco/interface/PFMET.h>
#include <DataFormats/METReco/interface/PFMETCollection.h>
#include <ZZAnalysis/AnalysisStep/interface/METCorrectionHandler.h>
#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/JetReco/interface/PFJetCollection.h>
#include <DataFormats/Math/interface/LorentzVector.h>
#include <CommonTools/UtilAlgos/interface/TFileService.h>
#include <DataFormats/Common/interface/MergeableCounter.h>
#include <DataFormats/VertexReco/interface/Vertex.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include <SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h>
#include <SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h>

#include <ZZAnalysis/AnalysisStep/interface/DaughterDataHelpers.h>
#include <ZZAnalysis/AnalysisStep/interface/FinalStates.h>
#include <ZZAnalysis/AnalysisStep/interface/MCHistoryTools.h>
//#include <ZZAnalysis/AnalysisStep/interface/PileUpWeight.h>
#include <ZZAnalysis/AnalysisStep/interface/PileUpWeight.h>
#include "SimDataFormats/HTXS/interface/HiggsTemplateCrossSections.h"


#include "ZZAnalysis/AnalysisStep/interface/EwkCorrections.h"
#include "ZZAnalysis/AnalysisStep/src/kFactors.C"
#include <ZZAnalysis/AnalysisStep/interface/bitops.h>
#include <ZZAnalysis/AnalysisStep/interface/Fisher.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonIsoHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/PhotonIDHelper.h>
#include <ZZAnalysis/AnalysisStep/interface/JetCleaner.h>
#include <ZZAnalysis/AnalysisStep/interface/utils.h>
#include <ZZAnalysis/AnalysisStep/interface/miscenums.h>
#include <ZZAnalysis/AnalysisStep/interface/ggF_qcd_uncertainty_2017.h>
#include <ZZAnalysis/AnalysisStep/interface/LeptonSFHelper.h>

#include <MelaAnalytics/CandidateLOCaster/interface/MELACandidateRecaster.h>
#include <CommonLHETools/LHEHandler/interface/LHEHandler.h>
#include <ZZAnalysis/AnalysisStep/interface/utils.h>//clara

#include "ZZ4lConfigHelper.h"
//#include "./ZZ4lConfigHelper.h"
//#include <./ZZ4lConfigHelper.h>
//#include "/afs/cern.ch/user/m/mrkim/work/HZZ/201211/CMSSW_10_2_18/src/ZZAnalysis/AnalysisStep/test/Ntuplizers/ZZ4lConfigHelper.h"
//#include </afs/cern.ch/user/m/mrkim/work/HZZ/201211/CMSSW_10_2_18/src/ZZAnalysis/AnalysisStep/test/Ntuplizers/ZZ4lConfigHelper.h>
//#include <ZZAnalysis/AnalysisStep/test/Ntuplizers/ZZ4lConfigHelper.h>
#include "HZZ4lNtupleFactory.h"

#include <TRandom3.h>
#include <TH2D.h>
#include "TLorentzVector.h"
#include "TSpline.h"
#include "TGraphErrors.h"


using namespace zzanalysis;//clara

namespace {
  
  bool writeJets = true;     // Write jets in the tree. FIXME: make this configurable
  bool writePhotons = true; // Write photons in the tree. FIXME: make this configurable
  bool addKinRefit = true;
  bool addVtxFit = false;
  bool addFSRDetails = false;
  bool addQGLInputs = true;
  bool skipMuDataMCWeight = false; // skip computation of data/MC weight for mu 
  bool skipEleDataMCWeight = false; // skip computation of data/MC weight for ele
  bool skipFakeWeight = true;   // skip computation of fake rate weight for CRs
  bool skipHqTWeight = true;    // skip computation of hQT weight

  //List of variables with default values
  Int_t RunNumber  = 0;
  Long64_t EventNumber  = 0;
  Int_t LumiNumber  = 0;
  Short_t NRecoMu  = 0;
  Short_t NRecoEle  = 0;
  Short_t Nvtx  = 0;
  Short_t NObsInt  = 0;
  Float_t NTrueInt  = 0;
  Float_t PUWeight  = 0;
  Float_t PUWeight_Up  = 0;
  Float_t PUWeight_Dn  = 0;

  Float_t KFactor_QCD_ggZZ_Nominal = 0;
  Float_t KFactor_QCD_ggZZ_PDFScaleDn = 0;
  Float_t KFactor_QCD_ggZZ_PDFScaleUp = 0;
  Float_t KFactor_QCD_ggZZ_QCDScaleDn = 0;
  Float_t KFactor_QCD_ggZZ_QCDScaleUp = 0;
  Float_t KFactor_QCD_ggZZ_AsDn = 0;
  Float_t KFactor_QCD_ggZZ_AsUp = 0;
  Float_t KFactor_QCD_ggZZ_PDFReplicaDn = 0;
  Float_t KFactor_QCD_ggZZ_PDFReplicaUp = 0;
  Float_t KFactor_EW_qqZZ = 0;
  Float_t KFactor_EW_qqZZ_unc = 0;
  Float_t KFactor_QCD_qqZZ_dPhi = 0;
  Float_t KFactor_QCD_qqZZ_M = 0;
  Float_t KFactor_QCD_qqZZ_Pt = 0;
  // Generic MET object
  METObject metobj;
  METObject metobj_corrected;
  Float_t GenMET = -99;
  Float_t GenMETPhi = -99;
  // MET with no HF
  //Float_t PFMETNoHF  =  -99;
  //Float_t PFMETNoHFPhi  =  -99;
  Short_t nCleanedJets  =  0;
  Short_t nCleanedJetsPt30  = 0;
  Short_t nCleanedJetsPt30_jesUp  = 0;
  Short_t nCleanedJetsPt30_jesDn  = 0;
  Short_t nCleanedJetsPt30_jerUp  = 0;
  Short_t nCleanedJetsPt30_jerDn  = 0;
  Short_t nCleanedJetsPt30BTagged  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesUp  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jesDn  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jerUp  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSF_jerDn  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSFUp  = 0;
  Short_t nCleanedJetsPt30BTagged_bTagSFDn  = 0;
  Short_t trigWord  = 0;
  Float_t ZZMass  = 0;
  Float_t ZZMassErr  = 0;
  Float_t ZZMassErrCorr  = 0;
  Float_t ZZMassPreFSR  = 0;
  Short_t ZZsel  = 0;
  Float_t ZZPt  = 0;
  Float_t ZZEta  = 0;
  Float_t ZZPhi  = 0;
  Float_t ZZjjPt = 0;
  Int_t CRflag  = 0;
  Float_t Z1Mass  = 0;
  Float_t Z1Pt  = 0;
  Short_t Z1Flav  = 0;
  Float_t ZZMassRefit  = 0;
  Float_t ZZMassRefitErr  = 0;
  Float_t ZZMassUnrefitErr  = 0;
  Float_t ZZMassCFit  = 0;
  Float_t ZZChi2CFit  = 0;
  Float_t Z2Mass  = 0;
  Float_t Z2Pt  = 0;
  Short_t Z2Flav  = 0;
  Float_t costhetastar  = 0;
  Float_t helphi  = 0;
  Float_t helcosthetaZ1  = 0;
  Float_t helcosthetaZ2  = 0;
  Float_t phistarZ1  = 0;
  Float_t phistarZ2  = 0;
  Float_t xi  = 0;
  Float_t xistar  = 0;
  Float_t TLE_dR_Z = -1; // Delta-R between a TLE and the Z it does not belong to.
  Float_t TLE_min_dR_3l = 999; // Minimum DR between a TLE and any of the other leptons
  Short_t evtPassMETTrigger = 0;

  using namespace std;//Clara
  
 //DEFINITION FOR JETTINESS
  TLorentzVector v4_na(0.0,0.0,1.0,1.0);//xyzE, clara
  TLorentzVector v4_nb(0.0,0.0,-1.0,1.0);//xy-zE, clara
  Int_t j=0;//clara, Glover counter v4_PEPM; 

  std::vector<float> Jettiness;//clara
  std::vector<float> LepPt;
  std::vector<float> LepEta;
  std::vector<float> LepPhi;
  std::vector<float> LepPx;//clara
  std::vector<float> LepPy;//clara
  std::vector<float> LepPz;//clara
  std::vector<float> LepE;//clara
  std::vector<float> LepSCEta;
  std::vector<short> LepLepId;
  std::vector<float> LepSIP;
  std::vector<float> Lepdxy;
  std::vector<float> Lepdz;
  std::vector<float> LepTime;
  std::vector<bool> LepisID;
  std::vector<float> LepBDT;
  std::vector<bool> LepisCrack;
  std::vector<char> LepMissingHit;
  std::vector<float> LepChargedHadIso;
  std::vector<float> LepNeutralHadIso;
  std::vector<float> LepPhotonIso;
  std::vector<float> LepPUIsoComponent;
  std::vector<float> LepCombRelIsoPF;
  std::vector<short> LepisLoose;
  std::vector<float> LepSF;
  std::vector<float> LepSF_Unc;
  std::vector<float> LepScale_Total_Up;
  std::vector<float> LepScale_Total_Dn;
  std::vector<float> LepScale_Stat_Up;
  std::vector<float> LepScale_Stat_Dn;
  std::vector<float> LepScale_Syst_Up;
  std::vector<float> LepScale_Syst_Dn;
  std::vector<float> LepScale_Gain_Up;
  std::vector<float> LepScale_Gain_Dn;
  std::vector<float> LepSigma_Total_Up;
  std::vector<float> LepSigma_Total_Dn;
  std::vector<float> LepSigma_Rho_Up;
  std::vector<float> LepSigma_Rho_Dn;
  std::vector<float> LepSigma_Phi_Up;
  std::vector<float> LepSigma_Phi_Dn;
 
  std::vector<float> tau0;
  std::vector<float> tau1;
  std::vector<float> tau2;
  std::vector<float> tau3;
  std::vector<float> tau4;
  std::vector<float> mintau;
  std::vector<float> tausum;

  std::vector<float> fsrPt;
  std::vector<float> fsrEta;
  std::vector<float> fsrPhi;
  std::vector<float> fsrDR;
  std::vector<short> fsrLept;
  std::vector<short> fsrLeptID;
  std::vector<float> fsrGenPt;
  Bool_t passIsoPreFSR = 0;

  std::vector<float> JetPt ;
  std::vector<float> JetEta ;
  std::vector<float> JetPhi ;
  std::vector<float> JetMass ;
  std::vector<float> JetBTagger ;
  std::vector<float> JetIsBtagged;
  std::vector<float> JetIsBtaggedWithSF;
  std::vector<float> JetIsBtaggedWithSFUp;
  std::vector<float> JetIsBtaggedWithSFDn;
  std::vector<float> JetQGLikelihood;
  std::vector<float> JetAxis2;
  std::vector<float> JetMult;
  std::vector<float> JetPtD;
  std::vector<float> JetSigma ;
  std::vector<short> JetHadronFlavour;
  std::vector<short> JetPartonFlavour;
   
  std::vector<float> JetPtJEC_noJER;
  std::vector<float> JetRawPt;

  std::vector<float> JetPUValue;
  std::vector<short> JetPUID;
    
  std::vector<short> JetID;
   
  std::vector<float> JetJESUp ;
  std::vector<float> JetJESDown ;
   
  std::vector<float> JetJERUp ;
  std::vector<float> JetJERDown ;

  Float_t DiJetMass  = -99;
//   Float_t DiJetMassPlus  = -99;
//   Float_t DiJetMassMinus  = -99;
  Float_t DiJetDEta  = -99;
  Float_t DiJetFisher  = -99;
  Short_t nExtraLep  = 0;
  Short_t nExtraZ  = 0;
  std::vector<float> ExtraLepPt;
  std::vector<float> ExtraLepEta;
  std::vector<float> ExtraLepPhi ;
  std::vector<float> ExtraLepPx;//clara
  std::vector<float> ExtraLepPy;//clara
  std::vector<float> ExtraLepPz;//clara
  std::vector<float> ExtraLepE;//clara
  std::vector<short> ExtraLepLepId;
  
  // Photon info
  std::vector<float> PhotonPt ;
  std::vector<float> PhotonEta ;
  std::vector<float> PhotonPhi ;
  std::vector<bool> PhotonIsCutBasedLooseID;
   
   
  Short_t genFinalState  = 0;
  Int_t genProcessId  = 0;
  Float_t genHEPMCweight  = 0;
  Float_t genHEPMCweight_NNLO  = 0;
  Float_t genHEPMCweight_POWHEGonly = 0;
	

  std::vector<float> LHEMotherPz;
  std::vector<float> LHEMotherE;
  std::vector<short> LHEMotherId;
  std::vector<float> LHEDaughterPt;
  std::vector<float> LHEDaughterEta;
  std::vector<float> LHEDaughterPhi;
  std::vector<float> LHEDaughterMass;
  std::vector<short> LHEDaughterId;
  std::vector<float> LHEAssociatedParticlePt;
  std::vector<float> LHEAssociatedParticleEta;
  std::vector<float> LHEAssociatedParticlePhi;
  std::vector<float> LHEAssociatedParticleMass;
  std::vector<short> LHEAssociatedParticleId;

  Float_t LHEPDFScale = 0;
  Float_t LHEweight_QCDscale_muR1_muF1  = 0;
  Float_t LHEweight_QCDscale_muR1_muF2  = 0;
  Float_t LHEweight_QCDscale_muR1_muF0p5  = 0;
  Float_t LHEweight_QCDscale_muR2_muF1  = 0;
  Float_t LHEweight_QCDscale_muR2_muF2  = 0;
  Float_t LHEweight_QCDscale_muR2_muF0p5  = 0;
  Float_t LHEweight_QCDscale_muR0p5_muF1  = 0;
  Float_t LHEweight_QCDscale_muR0p5_muF2  = 0;
  Float_t LHEweight_QCDscale_muR0p5_muF0p5  = 0;
  Float_t LHEweight_PDFVariation_Up = 0;
  Float_t LHEweight_PDFVariation_Dn = 0;
  Float_t LHEweight_AsMZ_Up = 0;
  Float_t LHEweight_AsMZ_Dn = 0;

  Float_t PythiaWeight_isr_muRoneoversqrt2 = 0;
  Float_t PythiaWeight_fsr_muRoneoversqrt2 = 0;
  Float_t PythiaWeight_isr_muRsqrt2 = 0;
  Float_t PythiaWeight_fsr_muRsqrt2 = 0;
  Float_t PythiaWeight_isr_muR0p5 = 0;
  Float_t PythiaWeight_fsr_muR0p5 = 0;
  Float_t PythiaWeight_isr_muR2 = 0;
  Float_t PythiaWeight_fsr_muR2 = 0;
  Float_t PythiaWeight_isr_muR0p25 = 0;
  Float_t PythiaWeight_fsr_muR0p25 = 0;
  Float_t PythiaWeight_isr_muR4 = 0;
  Float_t PythiaWeight_fsr_muR4 = 0;


  Short_t genExtInfo  = 0;
  Float_t xsection  = 0;
  Float_t genxsection = 0;
  Float_t genbranchingratio = 0;
  Float_t dataMCWeight  = 0;
  Float_t trigEffWeight  = 0;
  Float_t HqTMCweight  = 0;
  Float_t ZXFakeweight  = 0;
  Float_t overallEventWeight  = 0;
  Float_t L1prefiringWeight = 0;
  Float_t L1prefiringWeightUp = 0;
  Float_t L1prefiringWeightDn = 0;
  Float_t GenHMass  = 0;
  Float_t GenHPt  = 0;
  Float_t GenHRapidity  = 0;
  Float_t GenZ1Mass  = 0;
  Float_t GenZ1Eta  = 0;
  Float_t GenZ1Pt  = 0;
  Float_t GenZ1Phi  = 0;
  Float_t GenZ1Flav  = 0;
  Float_t GenZ2Mass  = 0;
  Float_t GenZ2Eta  = 0;
  Float_t GenZ2Pt  = 0;
  Float_t GenZ2Phi  = 0;
  Float_t GenZ2Flav  = 0;
  Float_t GenLep1Pt  = 0;
  Float_t GenLep1Eta  = 0;
  Float_t GenLep1Phi  = 0;
  Short_t GenLep1Id  = 0;
  Float_t GenLep2Pt  = 0;
  Float_t GenLep2Eta  = 0;
  Float_t GenLep2Phi  = 0;
  Short_t GenLep2Id  = 0;
  Float_t GenLep3Pt  = 0;
  Float_t GenLep3Eta  = 0;
  Float_t GenLep3Phi  = 0;
  Short_t GenLep3Id  = 0;
  Float_t GenLep4Pt  = 0;
  Float_t GenLep4Eta  = 0;
  Float_t GenLep4Phi  = 0;
  Short_t GenLep4Id  = 0;
  Float_t GenAssocLep1Pt  = 0;
  Float_t GenAssocLep1Eta  = 0;
  Float_t GenAssocLep1Phi  = 0;
  Short_t GenAssocLep1Id  = 0;
  Float_t GenAssocLep2Pt  = 0;
  Float_t GenAssocLep2Eta  = 0;
  Float_t GenAssocLep2Phi  = 0;
  Short_t GenAssocLep2Id  = 0;
	
  Int_t   htxsNJets = -1;
  Float_t htxsHPt = 0;
  Int_t   htxs_errorCode=-1;
  Int_t   htxs_prodMode=-1;
  Int_t   htxs_stage0_cat = -1;
  Int_t   htxs_stage1p1_cat = -1;
  Int_t   htxs_stage1p0_cat = -1;
  Int_t   htxs_stage1p2_cat = -1;
  Float_t ggH_NNLOPS_weight = 0;
  Float_t ggH_NNLOPS_weight_unc = 0;
  std::vector<float> qcd_ggF_uncertSF;


//FIXME: temporary fix to the mismatch of charge() and sign(pdgId()) for muons with BTT=4
  int getPdgId(const reco::Candidate* p) {
    int id = p->pdgId();
    if (id!=22 && //for TLEs
	signbit(id) && p->charge()<0) id*=-1; // negative pdgId must be positive charge
    return id;
  }

}//End of namespace

using namespace std;
using namespace edm;
//
// class declaration
//
class HZZ4lNtupleMaker : public edm::EDAnalyzer {
public:
  explicit HZZ4lNtupleMaker(const edm::ParameterSet&);
  ~HZZ4lNtupleMaker();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  static float EvalSpline(TSpline3* const& sp, float xval);

  static void addweight(float &weight, float weighttoadd);

private:
  virtual void beginJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);

  void BookAllBranches();
  virtual void FillKFactors(edm::Handle<GenEventInfoProduct>& genInfo, std::vector<const reco::Candidate *>& genZLeps);
  virtual void FillLHECandidate();
  virtual void FillCandidate(const pat::CompositeCandidate& higgs, bool evtPass, const edm::Event&, const Int_t CRflag);
  //virtual void FillJettiness(const pat::Jet& jet,const pat::MuonCollection& muon,const pat::ElectronCollection& electron,const pat::CompositeCandidate& higgs,const edm::Event&);
  virtual void FillJet(const pat::Jet& jet);
  virtual void FillPhoton(int year, const pat::Photon& photon);
  virtual void endJob() ;

  void FillHGenInfo(const math::XYZTLorentzVector Hp, float w);
  void FillZGenInfo(Short_t Z1Id, Short_t Z2Id,
                    const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2);
  void FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id,
    const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2, const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4);
  void FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps);

  
  Float_t getAllWeight(const vector<const reco::Candidate*>& leptons);
  Float_t getHqTWeight(double mH, double genPt) const;
  Float_t getFakeWeight(Float_t LepPt, Float_t LepEta, Int_t LepID, Int_t LepZ1ID);
  Int_t FindBinValue(TGraphErrors *tgraph, double value);

  void getCheckedUserFloat(const pat::CompositeCandidate& cand, const std::string& strval, Float_t& setval, Float_t defaultval=0);

  void buildMELABranches();
  void computeMELABranches(MELACandidate* cand);
  void updateMELAClusters_Common(const string clustertype);
  void updateMELAClusters_NoInitialQ(const string clustertype);
  void updateMELAClusters_NoInitialG(const string clustertype);
  void updateMELAClusters_NoAssociatedG(const string clustertype);
  void updateMELAClusters_NoInitialGNoAssociatedG(const string clustertype);
  void updateMELAClusters_BestLOAssociatedZ(const string clustertype);
  void updateMELAClusters_BestLOAssociatedW(const string clustertype);
  void updateMELAClusters_BestLOAssociatedVBF(const string clustertype);
  void updateMELAClusters_BestNLOVHApproximation(const string clustertype);
  void updateMELAClusters_BestNLOVBFApproximation(const string clustertype);
  void pushRecoMELABranches(const pat::CompositeCandidate& cand);
  void pushLHEMELABranches();
  void clearMELABranches();
	
	

  // ----------member data ---------------------------
  ZZ4lConfigHelper myHelper;
  int theChannel;
  std::string theCandLabel;
  TString theFileName;

  HZZ4lNtupleFactory *myTree;
  TH1F *hCounter;

  bool isMC;
  bool is_loose_ele_selection; // Collection includes candidates with loose electrons/TLEs
  bool applySkim;       //   "     "      "         skim (if skipEmptyEvents=true)
  bool skipEmptyEvents; // Skip events whith no selected candidate (otherwise, gen info is preserved for all events; candidates not passing trigger&&skim are flagged with negative ZZsel)
  FailedTreeLevel failedTreeLevel;  //if/how events with no selected candidate are written to a separate tree (see miscenums.h for details)
  edm::InputTag metTag;
  bool applyTrigger;    // Keep only events passing trigger (overriden if skipEmptyEvents=False)
  bool applyTrigEffWeight;// apply trigger efficiency weight (concerns samples where trigger is not applied)
  Float_t xsec;
  Float_t genxsec;
  Float_t genbr;
  int year;
  double sqrts;
  double Hmass;

  Mela mela;
  std::vector<std::string> recoMElist;
  std::vector<MELAOptionParser*> recome_originalopts;
  std::vector<MELAOptionParser*> recome_copyopts;
  std::vector<std::string> lheMElist;
  //std::vector<MELAOptionParser*> lheme_originalopts;
  std::vector<MELAOptionParser*> lheme_copyopts;
  std::vector<MELAHypothesis*> lheme_units;
  std::vector<MELAHypothesis*> lheme_aliased_units;
  std::vector<MELAComputation*> lheme_computers;
  std::vector<MELACluster*> lheme_clusters;

  bool addLHEKinematics;
  LHEHandler* lheHandler;
  METCorrectionHandler* metCorrHandler;
  int apply_K_NNLOQCD_ZZGG; // 0: Do not; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
  bool apply_K_NNLOQCD_ZZQQB;
  bool apply_K_NLOEW_ZZQQB;
  bool apply_QCD_GGF_UNCERT;

  edm::EDGetTokenT<edm::View<reco::Candidate> > genParticleToken;
  edm::Handle<edm::View<reco::Candidate> > genParticles;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > candToken;
  edm::EDGetTokenT<edm::View<pat::CompositeCandidate> > lhecandToken;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultToken;
  edm::EDGetTokenT<vector<reco::Vertex> > vtxToken;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
  edm::EDGetTokenT<pat::PhotonCollection> photonToken;
  edm::EDGetTokenT<pat::METCollection> metToken;
  //edm::EDGetTokenT<pat::METCollection> metNoHFToken;
  edm::EDGetTokenT<pat::MuonCollection> muonToken;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken;
  edm::EDGetTokenT<HTXS::HiggsClassification> htxsToken;
  edm::EDGetTokenT<edm::MergeableCounter> preSkimToken;
  edm::EDGetTokenT<LHERunInfoProduct> lheRunInfoToken;
   
  edm::EDGetTokenT< double > prefweight_token;
  edm::EDGetTokenT< double > prefweightup_token;
  edm::EDGetTokenT< double > prefweightdown_token;

  PileUpWeight* pileUpReweight;

  //counters
  Float_t Nevt_Gen;
  Float_t Nevt_Gen_lumiBlock;

  Float_t gen_ZZ4mu;
  Float_t gen_ZZ4e;
  Float_t gen_ZZ2mu2e;
  Float_t gen_ZZ2l2tau;
  Float_t gen_ZZ2emu2tau;
  Float_t gen_ZZ4tau;
  Float_t gen_ZZ4mu_EtaAcceptance;
  Float_t gen_ZZ4mu_LeptonAcceptance;
  Float_t gen_ZZ4e_EtaAcceptance;
  Float_t gen_ZZ4e_LeptonAcceptance;
  Float_t gen_ZZ2mu2e_EtaAcceptance;
  Float_t gen_ZZ2mu2e_LeptonAcceptance;
  Float_t gen_BUGGY;
  Float_t gen_Unknown;

  Float_t gen_sumPUWeight;
  Float_t gen_sumGenMCWeight;
  Float_t gen_sumWeights;

  string sampleName;

  std::vector<const reco::Candidate *> genFSR;

  std::vector<std::vector<float> > ewkTable;
  TSpline3* spkfactor_ggzz_nnlo[9]; // Nominal, PDFScaleDn, PDFScaleUp, QCDScaleDn, QCDScaleUp, AsDn, AsUp, PDFReplicaDn, PDFReplicaUp
  TSpline3* spkfactor_ggzz_nlo[9]; // Nominal, PDFScaleDn, PDFScaleUp, QCDScaleDn, QCDScaleUp, AsDn, AsUp, PDFReplicaDn, PDFReplicaUp

  LeptonSFHelper *lepSFHelper;

  TH2D* h_weight; //HqT weights
  //TH2F *h_ZXWeightMuo;
  //TH2F *h_ZXWeightEle;
  TH2D* h_ZXWeight[4];
	
  TGraphErrors *gr_NNLOPSratio_pt_powheg_0jet;
  TGraphErrors *gr_NNLOPSratio_pt_powheg_1jet;
  TGraphErrors *gr_NNLOPSratio_pt_powheg_2jet;
  TGraphErrors *gr_NNLOPSratio_pt_powheg_3jet;

  bool printedLHEweightwarning;
};//End of class HZZ4lNtupleMaker

//
// constructors and destructor
//
HZZ4lNtupleMaker::HZZ4lNtupleMaker(const edm::ParameterSet& pset) :
  myHelper(pset),
  theChannel(myHelper.channel()), // Valid options: ZZ, ZLL, ZL
  theCandLabel(pset.getUntrackedParameter<string>("CandCollection")), // Name of input ZZ collection
  theFileName(pset.getUntrackedParameter<string>("fileName")),
  myTree(nullptr),
  skipEmptyEvents(pset.getParameter<bool>("skipEmptyEvents")), // Do not store events with no selected candidate (normally: true)
  failedTreeLevel(FailedTreeLevel(pset.getParameter<int>("failedTreeLevel"))),
  metTag(pset.getParameter<edm::InputTag>("metSrc")),
  applyTrigger(pset.getParameter<bool>("applyTrigger")), // Reject events that do not pass trigger (normally: true)
  applyTrigEffWeight(pset.getParameter<bool>("applyTrigEff")), //Apply an additional efficiency weights for MC samples where triggers are not present (normally: false)
  xsec(pset.getParameter<double>("xsec")),
  genxsec(pset.getParameter<double>("GenXSEC")),
  genbr(pset.getParameter<double>("GenBR")),
  year(pset.getParameter<int>("setup")),
  sqrts(SetupToSqrts(year)),
  Hmass(pset.getParameter<double>("superMelaMass")),
  mela(sqrts, Hmass, TVar::ERROR),
  recoMElist(pset.getParameter<std::vector<std::string>>("recoProbabilities")),

  lheMElist(pset.getParameter<std::vector<std::string>>("lheProbabilities")),
  addLHEKinematics(pset.getParameter<bool>("AddLHEKinematics")),
  lheHandler(nullptr),
  metCorrHandler(nullptr),
  apply_K_NNLOQCD_ZZGG(pset.getParameter<int>("Apply_K_NNLOQCD_ZZGG")),
  apply_K_NNLOQCD_ZZQQB(pset.getParameter<bool>("Apply_K_NNLOQCD_ZZQQB")),
  apply_K_NLOEW_ZZQQB(pset.getParameter<bool>("Apply_K_NLOEW_ZZQQB")),
  apply_QCD_GGF_UNCERT(pset.getParameter<bool>("Apply_QCD_GGF_UNCERT")),

  pileUpReweight(nullptr),
  sampleName(pset.getParameter<string>("sampleName")),
  h_weight(0),

  printedLHEweightwarning(false)
{
  //cout<< "Beginning Constructor\n\n\n" <<endl;
  consumesMany<std::vector< PileupSummaryInfo > >();
  genParticleToken = consumes<edm::View<reco::Candidate> >(edm::InputTag("prunedGenParticles"));
  genInfoToken = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  consumesMany<LHEEventProduct>();
  candToken = consumes<edm::View<pat::CompositeCandidate> >(edm::InputTag(theCandLabel));

  is_loose_ele_selection = false;
  if(pset.exists("is_loose_ele_selection")) { 
    is_loose_ele_selection = pset.getParameter<bool>("is_loose_ele_selection");
  }
  triggerResultToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults"));
  vtxToken = consumes<vector<reco::Vertex> >(edm::InputTag("goodPrimaryVertices"));
  jetToken = consumes<edm::View<pat::Jet> >(edm::InputTag("cleanJets"));
  photonToken = consumes<pat::PhotonCollection>(edm::InputTag("slimmedPhotons"));
  metToken = consumes<pat::METCollection>(metTag);
  //metNoHFToken = consumes<pat::METCollection>(edm::InputTag("slimmedMETsNoHF"));
  muonToken = consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons"));
  electronToken = consumes<pat::ElectronCollection>(edm::InputTag("slimmedElectrons"));
  preSkimToken = consumes<edm::MergeableCounter,edm::InLumi>(edm::InputTag("preSkimCounter"));
  lheRunInfoToken = consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer"));
   
  if (skipEmptyEvents) {
    applySkim=true;
  } else {
    applyTrigger=false; // This overrides the card applyTrigger
    applySkim=false;
    failedTreeLevel=noFailedTree; // This overrides the card failedTreeLevel
  }

  if (applyTrigEffWeight&&applyTrigger) {
    cout << "ERROR: cannot have applyTrigEffWeight == applyTrigger == true" << endl;
  }

  isMC = myHelper.isMC();
   
  if( isMC && (year == 2016 || year == 2017))
  {
     prefweight_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"));
     prefweightup_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
     prefweightdown_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown"));
  }
   
   
  addLHEKinematics = addLHEKinematics || !lheMElist.empty();
  if (isMC){
    lheHandler = new LHEHandler(
      ((MELAEvent::CandidateVVMode)(pset.getParameter<int>("VVMode")+1)), // FIXME: Need to pass strings and interpret them instead!
      pset.getParameter<int>("VVDecayMode"),
      (addLHEKinematics ? LHEHandler::doHiggsKinematics : LHEHandler::noKinematics),
      year, LHEHandler::tryNNPDF30, LHEHandler::tryNLO
    );
    metCorrHandler = new METCorrectionHandler(Form("%i", year));
    htxsToken = consumes<HTXS::HiggsClassification>(edm::InputTag("rivetProducerHTXS","HiggsClassification"));
    pileUpReweight = new PileUpWeight(myHelper.sampleType(), myHelper.setup());
  }

  Nevt_Gen = 0;
  Nevt_Gen_lumiBlock = 0;

  //For Efficiency studies
  gen_ZZ4mu = 0;
  gen_ZZ4e = 0;
  gen_ZZ2mu2e = 0;
  gen_ZZ2l2tau = 0;
  gen_ZZ2emu2tau = 0;
  gen_ZZ4tau = 0;
  gen_ZZ4mu_EtaAcceptance = 0;
  gen_ZZ4mu_LeptonAcceptance = 0;
  gen_ZZ4e_EtaAcceptance = 0;
  gen_ZZ4e_LeptonAcceptance = 0;
  gen_ZZ2mu2e_EtaAcceptance = 0;
  gen_ZZ2mu2e_LeptonAcceptance = 0;
  gen_BUGGY = 0;
  gen_Unknown = 0;

  gen_sumPUWeight = 0.f;
  gen_sumGenMCWeight = 0.f;
  gen_sumWeights =0.f;

  std::string fipPath;

  // Read EWK K-factor table from file
  edm::FileInPath ewkFIP("ZZAnalysis/AnalysisStep/data/kfactors/ZZ_EwkCorrections.dat");
  fipPath=ewkFIP.fullPath();
  ewkTable = EwkCorrections::readFile_and_loadEwkTable(fipPath.data());

  // Read the ggZZ k-factor shape from file
  TString strZZGGKFVar[9]={
    "Nominal", "PDFScaleDn", "PDFScaleUp", "QCDScaleDn", "QCDScaleUp", "AsDn", "AsUp", "PDFReplicaDn", "PDFReplicaUp"
  };
  edm::FileInPath ggzzFIP_NNLO("ZZAnalysis/AnalysisStep/data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NNLO_NNPDF_NarrowWidth_13TeV.root");
  fipPath=ggzzFIP_NNLO.fullPath();
  TFile* ggZZKFactorFile = TFile::Open(fipPath.data());
  for (unsigned int ikf=0; ikf<9; ikf++) spkfactor_ggzz_nnlo[ikf] = (TSpline3*)ggZZKFactorFile->Get(Form("sp_kfactor_%s", strZZGGKFVar[ikf].Data()))->Clone(Form("sp_kfactor_%s_NNLO", strZZGGKFVar[ikf].Data()));
  ggZZKFactorFile->Close();
  edm::FileInPath ggzzFIP_NLO("ZZAnalysis/AnalysisStep/data/kfactors/Kfactor_Collected_ggHZZ_2l2l_NLO_NNPDF_NarrowWidth_13TeV.root");
  fipPath=ggzzFIP_NLO.fullPath();
  ggZZKFactorFile = TFile::Open(fipPath.data());
  for (unsigned int ikf=0; ikf<9; ikf++) spkfactor_ggzz_nlo[ikf] = (TSpline3*)ggZZKFactorFile->Get(Form("sp_kfactor_%s", strZZGGKFVar[ikf].Data()))->Clone(Form("sp_kfactor_%s_NLO", strZZGGKFVar[ikf].Data()));
  ggZZKFactorFile->Close();
	
  edm::FileInPath NNLOPS_weight_path("ZZAnalysis/AnalysisStep/data/ggH_NNLOPS_Weights/NNLOPS_reweight.root");
  fipPath=NNLOPS_weight_path.fullPath();
  TFile* NNLOPS_weight_file = TFile::Open(fipPath.data());
  gr_NNLOPSratio_pt_powheg_0jet = (TGraphErrors*)NNLOPS_weight_file->Get("gr_NNLOPSratio_pt_powheg_0jet");
  gr_NNLOPSratio_pt_powheg_1jet = (TGraphErrors*)NNLOPS_weight_file->Get("gr_NNLOPSratio_pt_powheg_1jet");
  gr_NNLOPSratio_pt_powheg_2jet = (TGraphErrors*)NNLOPS_weight_file->Get("gr_NNLOPSratio_pt_powheg_2jet");
  gr_NNLOPSratio_pt_powheg_3jet = (TGraphErrors*)NNLOPS_weight_file->Get("gr_NNLOPSratio_pt_powheg_3jet");

  //Scale factors for data/MC efficiency
  if (!skipEleDataMCWeight && isMC) { lepSFHelper = new LeptonSFHelper(); }

  if (!skipHqTWeight) {
    //HqT weights
    edm::FileInPath HqTfip("ZZAnalysis/AnalysisStep/test/Macros/HqTWeights.root");
    std::string fipPath=HqTfip.fullPath();
    TFile *fHqt = TFile::Open(fipPath.data(),"READ");
    h_weight = (TH2D*)fHqt->Get("wH_400")->Clone();//FIXME: Ask simon to provide the 2D histo
    fHqt->Close();
  }

  if (!skipFakeWeight) {
    //CR fake rate weight
    TString filename;filename.Form("ZZAnalysis/AnalysisStep/test/Macros/FR2_2011_AA_electron.root");
    if(year==2015)filename.Form("ZZAnalysis/AnalysisStep/test/Macros/FR2_AA_ControlSample_ABCD.root");
    edm::FileInPath fipEleZX(filename.Data());
    std::string fipPath=fipEleZX.fullPath();
    TFile *FileZXWeightEle = TFile::Open(fipPath.data(),"READ");

    filename.Form("ZZAnalysis/AnalysisStep/test/Macros/FR2_2011_AA_muon.root");
    if(year==2015)filename.Form("ZZAnalysis/AnalysisStep/test/Macros/FR2_AA_muon.root");
    edm::FileInPath fipMuZX(filename.Data());
    fipPath=fipMuZX.fullPath();
    TFile *FileZXWeightMuo = TFile::Open(fipPath.data(),"READ");

    h_ZXWeight[0]=(TH2D*)FileZXWeightEle->Get("eff_Z1ee_plus_electron")->Clone();
    h_ZXWeight[1]=(TH2D*)FileZXWeightEle->Get("eff_Z1mumu_plus_electron")->Clone();
    h_ZXWeight[2]=(TH2D*)FileZXWeightMuo->Get("eff_Z1ee_plus_muon")->Clone();
    h_ZXWeight[3]=(TH2D*)FileZXWeightMuo->Get("eff_Z1mumu_plus_muon")->Clone();

    FileZXWeightEle->Close();
    FileZXWeightMuo->Close();
  }
}//End of No name parentheis

HZZ4lNtupleMaker::~HZZ4lNtupleMaker()
{
  clearMELABranches(); // Cleans LHE branches
  delete lheHandler;
  delete pileUpReweight;
  delete metCorrHandler;
}


// ------------ method called for each event  ------------
void HZZ4lNtupleMaker::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  myTree->InitializeVariables();

  //----------------------------------------------------------------------
  // Analyze MC truth; collect MC weights and update counters (this is done for all generated events,
  // including those that do not pass skim, trigger etc!)

  bool gen_ZZ4lInEtaAcceptance = false;   // All 4 gen leptons in eta acceptance
  bool gen_ZZ4lInEtaPtAcceptance = false; // All 4 gen leptons in eta,pT acceptance

  const reco::Candidate * genH = 0;
  std::vector<const reco::Candidate *> genZLeps;
  std::vector<const reco::Candidate *> genAssocLeps;

  edm::Handle<GenEventInfoProduct> genInfo;

  if (isMC) {
    // get PU weights
    vector<Handle<std::vector< PileupSummaryInfo > > >  PupInfos; //FIXME support for miniAOD v1/v2 where name changed; catch does not work...
    event.getManyByType(PupInfos);
    Handle<std::vector< PileupSummaryInfo > > PupInfo = PupInfos.front();
//     try {
//       cout << "TRY HZZ4lNtupleMaker" <<endl;
//       event.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
//     } catch (const cms::Exception& e){
//       cout << "FAIL HZZ4lNtupleMaker" <<endl;
//       event.getByLabel(edm::InputTag("slimmedAddPileupInfo"), PupInfo);
//     }

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(PVI->getBunchCrossing() == 0) {
        NObsInt  = PVI->getPU_NumInteractions();
        NTrueInt = PVI->getTrueNumInteractions();
        break;
      }
    }

    // get PU weight
    PUWeight = pileUpReweight->weight(NTrueInt);
    PUWeight_Up = pileUpReweight->weight(NTrueInt, PileUpWeight::PUvar::VARUP);
    PUWeight_Dn = pileUpReweight->weight(NTrueInt, PileUpWeight::PUvar::VARDOWN);
     
    // L1 prefiring weights
    if( year == 2016 || year == 2017 )
    {
       edm::Handle< double > theprefweight;
       event.getByToken(prefweight_token, theprefweight ) ;
       L1prefiringWeight =(*theprefweight);
        
       edm::Handle< double > theprefweightup;
       event.getByToken(prefweightup_token, theprefweightup ) ;
       L1prefiringWeightUp =(*theprefweightup);
        
       edm::Handle< double > theprefweightdown;
       event.getByToken(prefweightdown_token, theprefweightdown ) ;
       L1prefiringWeightDn =(*theprefweightdown);
    }
    else if ( year == 2018 )
    {
       L1prefiringWeight   = 1.;
       L1prefiringWeightUp = 1.;
       L1prefiringWeightDn = 1.;
    }

    event.getByToken(genParticleToken, genParticles);
    event.getByToken(genInfoToken, genInfo);

    edm::Handle<HTXS::HiggsClassification> htxs;
    event.getByToken(htxsToken,htxs);

    MCHistoryTools mch(event, sampleName, genParticles, genInfo);
    genFinalState = mch.genFinalState();
    genProcessId = mch.getProcessID();
    genHEPMCweight_NNLO = genHEPMCweight = mch.gethepMCweight(); // Overridden by LHEHandler if genHEPMCweight==1.
                                                                 // For 2017 MC, genHEPMCweight is reweighted later from NNLO to NLO
    const auto& genweights = genInfo->weights();
    if (genweights.size() > 1){
      if ((genweights.size() != 14 && genweights.size() != 46) || genweights[0] != genweights[1]){
        cms::Exception e("GenWeights");
        e << "Expected to find 1 gen weight, or 14 or 46 with the first two the same, found " << genweights.size() << ":\n";
        for (auto w : genweights) e << w << " ";
        throw e;
      }
      auto nominal = genweights[0];
      PythiaWeight_isr_muRoneoversqrt2 = genweights[2] / nominal;
      PythiaWeight_fsr_muRoneoversqrt2 = genweights[3] / nominal;
      PythiaWeight_isr_muRsqrt2 = genweights[4] / nominal;
      PythiaWeight_fsr_muRsqrt2 = genweights[5] / nominal;

      PythiaWeight_isr_muR0p5 = genweights[6] / nominal;
      PythiaWeight_fsr_muR0p5 = genweights[7] / nominal;
      PythiaWeight_isr_muR2 = genweights[8] / nominal;
      PythiaWeight_fsr_muR2 = genweights[9] / nominal;

      PythiaWeight_isr_muR0p25 = genweights[10] / nominal;
      PythiaWeight_fsr_muR0p25 = genweights[11] / nominal;
      PythiaWeight_isr_muR4 = genweights[12] / nominal;
      PythiaWeight_fsr_muR4 = genweights[13] / nominal;
    } else {
      PythiaWeight_isr_muRsqrt2 = PythiaWeight_isr_muRoneoversqrt2 = PythiaWeight_isr_muR2 =
      PythiaWeight_isr_muR0p5 = PythiaWeight_isr_muR4 = PythiaWeight_isr_muR0p25 =
      PythiaWeight_fsr_muRsqrt2 = PythiaWeight_fsr_muRoneoversqrt2 = PythiaWeight_fsr_muR2 =
      PythiaWeight_fsr_muR0p5 = PythiaWeight_fsr_muR4 = PythiaWeight_fsr_muR0p25 = 1;
    }

   htxsNJets = htxs->jets30.size();
   htxsHPt = htxs->higgs.Pt();
   htxs_stage0_cat = htxs->stage0_cat;
   htxs_stage1p0_cat = htxs->stage1_cat_pTjet30GeV;
   htxs_stage1p1_cat = htxs->stage1_1_cat_pTjet30GeV;
   htxs_stage1p2_cat = htxs->stage1_2_cat_pTjet30GeV;
   htxs_errorCode=htxs->errorCode;
   htxs_prodMode= htxs->prodMode;

   genExtInfo = mch.genAssociatedFS();

   //Information on generated candidates, will be used later
   genH = mch.genH();
   genZLeps     = mch.sortedGenZZLeps();
   genAssocLeps = mch.genAssociatedLeps();
   genFSR       = mch.genFSR();



    if(genH != 0){
      FillHGenInfo(genH->p4(),getHqTWeight(genH->p4().M(),genH->p4().Pt()));
    }
    else if(genZLeps.size()==4){ // for 4l events take the mass of the ZZ(4l) system
      FillHGenInfo((genZLeps.at(0)->p4()+genZLeps.at(1)->p4()+genZLeps.at(2)->p4()+genZLeps.at(3)->p4()),0);
    }

    if (genFinalState!=BUGGY) {

      if (genZLeps.size()==4) {

        // "generated Zs" defined with standard pairing applied on gen leptons (genZLeps is sorted by MCHistoryTools)
        FillZGenInfo(genZLeps.at(0)->pdgId()*genZLeps.at(1)->pdgId(),
                     genZLeps.at(2)->pdgId()*genZLeps.at(3)->pdgId(),
                     genZLeps.at(0)->p4()+genZLeps.at(1)->p4(),
                     genZLeps.at(2)->p4()+genZLeps.at(3)->p4());

        // Gen leptons
        FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), genZLeps.at(3)->pdgId(),
           genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), genZLeps.at(3)->p4());

      }

      if (genZLeps.size()==3) {
        FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), genZLeps.at(2)->pdgId(), 0,
                       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), genZLeps.at(2)->p4(), *(new math::XYZTLorentzVector));
      }
      if (genZLeps.size()==2) {
        FillLepGenInfo(genZLeps.at(0)->pdgId(), genZLeps.at(1)->pdgId(), 0, 0,
                       genZLeps.at(0)->p4(), genZLeps.at(1)->p4(), *(new math::XYZTLorentzVector), *(new math::XYZTLorentzVector));
      }

      if (genAssocLeps.size()==1 || genAssocLeps.size()==2) {
        FillAssocLepGenInfo(genAssocLeps);
      }

      // LHE information
      edm::Handle<LHEEventProduct> lhe_evt;
      vector<edm::Handle<LHEEventProduct> > lhe_handles;
      event.getManyByType(lhe_handles);
      if (!lhe_handles.empty()){
        lhe_evt = lhe_handles.front();
        lheHandler->setHandle(&lhe_evt);
        lheHandler->extract();
        FillLHECandidate(); // Also writes weights
        lheHandler->clear();
      }
      //else cerr << "lhe_handles.size()==0" << endl;

      // keep track of sum of weights
      addweight(gen_sumPUWeight, PUWeight);
      addweight(gen_sumGenMCWeight, genHEPMCweight);
      addweight(gen_sumWeights, PUWeight*genHEPMCweight);

      mch.genAcceptance(gen_ZZ4lInEtaAcceptance, gen_ZZ4lInEtaPtAcceptance);
    }//End of 

    addweight(Nevt_Gen_lumiBlock, 1); // Needs to be outside the if-block

    if (genFinalState == EEEE) {
      addweight(gen_ZZ4e, 1);
      if (gen_ZZ4lInEtaAcceptance) addweight(gen_ZZ4e_EtaAcceptance, 1);
      if (gen_ZZ4lInEtaPtAcceptance) addweight(gen_ZZ4e_LeptonAcceptance, 1);
    } else if (genFinalState == MMMM) {
      addweight(gen_ZZ4mu, 1);
      if (gen_ZZ4lInEtaAcceptance) addweight(gen_ZZ4mu_EtaAcceptance, 1);
      if (gen_ZZ4lInEtaPtAcceptance) addweight(gen_ZZ4mu_LeptonAcceptance, 1);
    } else if (genFinalState == EEMM) {
      addweight(gen_ZZ2mu2e, 1);
      if (gen_ZZ4lInEtaAcceptance) addweight(gen_ZZ2mu2e_EtaAcceptance, 1);
      if (gen_ZZ4lInEtaPtAcceptance) addweight(gen_ZZ2mu2e_LeptonAcceptance, 1);
    } else if (genFinalState == llTT){
      addweight(gen_ZZ2emu2tau, 1);
      addweight(gen_ZZ2l2tau, 1);
    } else if (genFinalState == TTTT){
      addweight(gen_ZZ4tau, 1);
      addweight(gen_ZZ2l2tau, 1);
    } else if (genFinalState == BUGGY){ // handle MCFM ZZ->4tau mZ<2mtau bug
      addweight(gen_BUGGY, 1);
      return; // BUGGY events are skipped
    } else {
      addweight(gen_Unknown, 1);
    }

// End of MC history analysis ------------------------------------------
  }//End of if(isMC) 
  else {
    ++Nevt_Gen_lumiBlock; // keep track of # events for data as well
  }



  // Get candidate collection
  edm::Handle<edm::View<pat::CompositeCandidate> > candHandle;
  event.getByToken(candToken, candHandle);
  if(candHandle.failedToGet()) {
    if(is_loose_ele_selection) return; // The collection can be missing in this case since we have a filter to skip the module when a regular candidate is present.
    else edm::LogError("") << "ZZ collection not found in non-loose electron flow. This should never happen";
  }
  const edm::View<pat::CompositeCandidate>* cands = candHandle.product();

  // For Z+L CRs, we want only events with exactly 1 Z+l candidate. FIXME: this has to be reviewed.
  if (theChannel==ZL && cands->size() != 1) return;


  // Retrieve trigger results
  Handle<edm::TriggerResults> triggerResults;
  event.getByToken(triggerResultToken, triggerResults);

  bool failed = false;

  // Apply MC filter (skip event)
  // Heshy note: I'm not turning return into failed = true because it looks like it's applied even if !skipEmptyEvents.
  //             It only does anything if the MCFILTER variable is set in the csv file, which is not currently the case.
  if (isMC && !(myHelper.passMCFilter(event,triggerResults))) return;

  // Apply skim
  bool evtPassSkim = myHelper.passSkim(event,triggerResults,trigWord);
  if (applySkim && !evtPassSkim) failed = true;       //but gen information will still be recorded if failedTreeLevel != 0

  // Apply trigger request (skip event)
  bool evtPassTrigger = myHelper.passTrigger(event,triggerResults,trigWord);
  if (applyTrigger && !evtPassTrigger) failed = true; //but gen information will still be recorded if failedTreeLevel != 0
	
  // Apply MET trigger request (skip event)
  evtPassMETTrigger = myHelper.passMETTrigger(event,triggerResults);

  if (skipEmptyEvents && !failedTreeLevel && (cands->size() == 0 || failed)) return; // Skip events with no candidate, unless skipEmptyEvents = false or failedTreeLevel != 0

  //Fill MC truth information
  if (isMC) FillKFactors(genInfo, genZLeps);

  // General event information
  RunNumber=event.id().run();
  LumiNumber=event.luminosityBlock();
  EventNumber=event.id().event();
  xsection=xsec;
  genxsection=genxsec;
  genbranchingratio=genbr;

  // Primary vertices
  Handle<vector<reco::Vertex> > vertices;
  event.getByToken(vtxToken,vertices);
  Nvtx=vertices->size();


  // Jets (cleaned wrt all tight isolated leptons)
  Handle<edm::View<pat::Jet> > CleanedJets;
  event.getByToken(jetToken, CleanedJets);
  vector<const pat::Jet*> cleanedJets;
  for(edm::View<pat::Jet>::const_iterator jet = CleanedJets->begin(); jet != CleanedJets->end(); ++jet){
    cleanedJets.push_back(&*jet);
  }
   
   // Photons
   Handle<pat::PhotonCollection> photonCands;
   event.getByToken(photonToken, photonCands);
   vector<const pat::Photon*> photons;
   
   for(unsigned int i = 0; i< photonCands->size(); ++i){
      const pat::Photon* photon = &((*photonCands)[i]);
      photons.push_back(&*photon);
   }

   
   if (writePhotons){
      for (unsigned i=0; i<photons.size(); ++i) {
            FillPhoton(year, *(photons.at(i)));
         }
   }
      
  // MET
  Handle<pat::METCollection> metHandle;
  event.getByToken(metToken, metHandle);

  GenMET=GenMETPhi=-99;
  if (metHandle.isValid()){
    const pat::MET &met = metHandle->front();

    metobj.extras.met = metobj.extras.met_original = metobj.extras.met_raw
      = metobj.extras.met_METup = metobj.extras.met_METdn
 



     = metobj.extras.met_JERup = metobj.extras.met_JERdn
      = metobj.extras.met_PUup = metobj.extras.met_PUdn

      = metobj_corrected.extras.met = metobj_corrected.extras.met_original = metobj_corrected.extras.met_raw
      = metobj_corrected.extras.met_METup = metobj_corrected.extras.met_METdn
      = metobj_corrected.extras.met_JERup = metobj_corrected.extras.met_JERdn
      = metobj_corrected.extras.met_PUup = metobj_corrected.extras.met_PUdn

      = met.pt();
    metobj.extras.phi = metobj.extras.phi_original = metobj.extras.phi_raw
      = metobj.extras.phi_METup = metobj.extras.phi_METdn
      = metobj.extras.phi_JECup = metobj.extras.phi_JECdn
      = metobj.extras.phi_JERup = metobj.extras.phi_JERdn
      = metobj.extras.phi_PUup = metobj.extras.phi_PUdn

      = metobj_corrected.extras.phi = metobj_corrected.extras.phi_original = metobj_corrected.extras.phi_raw
      = metobj_corrected.extras.phi_METup = metobj_corrected.extras.phi_METdn
      = metobj_corrected.extras.phi_JECup = metobj_corrected.extras.phi_JECdn
      = metobj_corrected.extras.phi_JERup = metobj_corrected.extras.phi_JERdn
      = metobj_corrected.extras.phi_PUup = metobj_corrected.extras.phi_PUdn

      = met.phi();

    metobj.extras.met_JECup = metobj_corrected.extras.met_JECup = met.shiftedPt(pat::MET::JetEnUp);
    metobj.extras.met_JECdn = metobj_corrected.extras.met_JECdn = met.shiftedPt(pat::MET::JetEnDown);
    metobj.extras.phi_JECup = metobj_corrected.extras.phi_JECup = met.shiftedPhi(pat::MET::JetEnUp);
    metobj.extras.phi_JECdn = metobj_corrected.extras.phi_JECdn = met.shiftedPhi(pat::MET::JetEnDown);

    if (isMC && metCorrHandler && met.genMET()){
      GenMET = met.genMET()->pt();
      GenMETPhi = met.genMET()->phi();
      metCorrHandler->correctMET(GenMET, GenMETPhi, &metobj_corrected, false); // FIXME: Last argument should be for isFastSim, but we don't have it yet
    }
    else if (isMC){
      cms::Exception e("METCorrectionHandler");
      e << "Either no met.genMET or metCorrHandler!";
      throw e;
    }
  }//End of if (metHandle.isValid())
  else{
    metobj.extras.met = metobj.extras.met_original = metobj.extras.met_raw
      = metobj.extras.met_METup = metobj.extras.met_METdn
      = metobj.extras.met_JECup = metobj.extras.met_JECdn
      = metobj.extras.met_JERup = metobj.extras.met_JERdn
      = metobj.extras.met_PUup = metobj.extras.met_PUdn

      = metobj_corrected.extras.met = metobj_corrected.extras.met_original = metobj_corrected.extras.met_raw
      = metobj_corrected.extras.met_METup = metobj_corrected.extras.met_METdn
      = metobj_corrected.extras.met_JECup = metobj_corrected.extras.met_JECdn
      = metobj_corrected.extras.met_JERup = metobj_corrected.extras.met_JERdn
      = metobj_corrected.extras.met_PUup = metobj_corrected.extras.met_PUdn

      = metobj.extras.phi = metobj.extras.phi_original = metobj.extras.phi_raw
      = metobj.extras.phi_METup = metobj.extras.phi_METdn
      = metobj.extras.phi_JECup = metobj.extras.phi_JECdn
      = metobj.extras.phi_JERup = metobj.extras.phi_JERdn
      = metobj.extras.phi_PUup = metobj.extras.phi_PUdn

      = metobj_corrected.extras.phi = metobj_corrected.extras.phi_original = metobj_corrected.extras.phi_raw
      = metobj_corrected.extras.phi_METup = metobj_corrected.extras.phi_METdn
      = metobj_corrected.extras.phi_JECup = metobj_corrected.extras.phi_JECdn
      = metobj_corrected.extras.phi_JERup = metobj_corrected.extras.phi_JERdn
      = metobj_corrected.extras.phi_PUup = metobj_corrected.extras.phi_PUdn

      = -99;
  }//End of else
  //Handle<pat::METCollection> metNoHFHandle;
  //event.getByToken(metNoHFToken, metNoHFHandle);
  //if(metNoHFHandle.isValid()){
  //  PFMETNoHF = metNoHFHandle->front().pt();
  //  PFMETNoHFPhi = metNoHFHandle->front().phi();
  //}


  // number of reconstructed leptons
  edm::Handle<pat::MuonCollection> muonHandle;
  event.getByToken(muonToken, muonHandle);
  for(unsigned int i = 0; i< muonHandle->size(); ++i){
    const pat::Muon* m = &((*muonHandle)[i]);
    if(m->pt()>5 && m->isPFMuon()) // these cuts are implicit in miniAOD
      NRecoMu++;
  }
  edm::Handle<pat::ElectronCollection> electronHandle;
  event.getByToken(electronToken, electronHandle);
  for(unsigned int i = 0; i< electronHandle->size(); ++i){
    const pat::Electron* e = &((*electronHandle)[i]);
    if(e->pt()>5) // this cut is implicit in miniAOD
      NRecoEle++;
  }

  if(isMC && apply_QCD_GGF_UNCERT)
  {

	  

    if (htxsNJets==0)
    {
      ggH_NNLOPS_weight = gr_NNLOPSratio_pt_powheg_0jet->Eval(min((double) htxsHPt, 125.0));
      ggH_NNLOPS_weight_unc=(gr_NNLOPSratio_pt_powheg_0jet->GetErrorY(FindBinValue(gr_NNLOPSratio_pt_powheg_0jet, min((double) htxsHPt, 125.0))))/ggH_NNLOPS_weight;

    }
	  else if (htxsNJets==1)
	  {
		  ggH_NNLOPS_weight = gr_NNLOPSratio_pt_powheg_1jet->Eval(min((double)htxsHPt,625.0));
		  ggH_NNLOPS_weight_unc=(gr_NNLOPSratio_pt_powheg_1jet->GetErrorY(FindBinValue(gr_NNLOPSratio_pt_powheg_1jet,min((double)htxsHPt,125.0))))/ggH_NNLOPS_weight;
	  }
	  else if (htxsNJets==2)
	  {
		  ggH_NNLOPS_weight = gr_NNLOPSratio_pt_powheg_2jet->Eval(min((double)htxsHPt,800.0));
		  ggH_NNLOPS_weight_unc=(gr_NNLOPSratio_pt_powheg_2jet->GetErrorY(FindBinValue(gr_NNLOPSratio_pt_powheg_2jet,min((double)htxsHPt,125.0))))/ggH_NNLOPS_weight;
	  }
	  else if (htxsNJets>=3)
	  {
		  ggH_NNLOPS_weight = gr_NNLOPSratio_pt_powheg_3jet->Eval(min((double)htxsHPt,925.0));
		  ggH_NNLOPS_weight_unc=(gr_NNLOPSratio_pt_powheg_3jet->GetErrorY(FindBinValue(gr_NNLOPSratio_pt_powheg_3jet,min((double)htxsHPt,125.0))))/ggH_NNLOPS_weight;
	  }
	  else
	  {
		  ggH_NNLOPS_weight = 1.0;
		  ggH_NNLOPS_weight_unc = 0.0;
	  }


	  std::vector<double> qcd_ggF_uncertSF_tmp;
	  qcd_ggF_uncertSF.clear();

    ////////////////////////////////////////////////////////////
    //////////////////     CHECK THIS!!!!!    //////////////////
    // Why is this done with the STXS 1.0 bins uncertainties? //
    //////////////////     CHECK THIS!!!!!    //////////////////
    ////////////////////////////////////////////////////////////

	  qcd_ggF_uncertSF_tmp = qcd_ggF_uncertSF_2017(htxsNJets, htxsHPt, htxs_stage1p0_cat);
	  qcd_ggF_uncertSF = std::vector<float>(qcd_ggF_uncertSF_tmp.begin(),qcd_ggF_uncertSF_tmp.end());


  }//End of if(isMC && apply_QCD_GGF_UNCERT) 
   
  //Loop on the candidates
  vector<Int_t> CRFLAG(cands->size());
  for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
    if (failed) break; //don't waste time on this
    size_t icand= cand-cands->begin();

    //    int candChannel = cand->userFloat("candChannel"); // This is currently the product of pdgId of leptons (eg 14641, 28561, 20449)

    if (theChannel==ZLL) {
      // Cross check region for Z + 1 loose electron + 1 loose TLE (defined only in loose_ele paths) 
      if (is_loose_ele_selection) {
        if (cand->userFloat("isBestCRZLL")&&cand->userFloat("CRZLL")) set_bit(CRFLAG[icand],ZLL);
      }

      // AA CRs
      if (cand->userFloat("isBestCRZLLss")&&cand->userFloat("CRZLLss")) set_bit(CRFLAG[icand],CRZLLss);

      // A CRs
      if (cand->userFloat("isBestCRZLLos_2P2F")&&cand->userFloat("CRZLLos_2P2F")) set_bit(CRFLAG[icand],CRZLLos_2P2F);
      if (cand->userFloat("isBestCRZLLos_3P1F")&&cand->userFloat("CRZLLos_3P1F")) set_bit(CRFLAG[icand],CRZLLos_3P1F);

      if (CRFLAG[icand]) { // This candidate belongs to one of the CRs: perform additional jet cleaning.
        // Note that this is (somewhat incorrectly) done per-event, so there could be some over-cleaning in events with >1 CR candidate.
        for (unsigned i=0; i<cleanedJets.size(); ++i) {
          if (cleanedJets[i]!=0  && (!jetCleaner::isGood(*cand, *(cleanedJets[i])))) {
            cleanedJets[i]=0;
          }
        }
      }
    }
  }//End of for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) 

  //v4_Jet.clear();//Clara
  // Count and store jets, after additional cleaning for CRs...
  for (unsigned i=0; i<cleanedJets.size(); ++i) {
    if (cleanedJets[i]==0) {
      continue; // Jet has been suppressed by additional cleaning
    }
    //v4_Jet.push_back(cleandJets->p4());//Clara
    ++nCleanedJets;

    // count jes up/down njets pt30
    float jes_unc = cleanedJets[i]->userFloat("jes_unc");

    float pt_nominal = cleanedJets[i]->pt();
    float pt_jes_up = pt_nominal * (1.0 + jes_unc);
    float pt_jes_dn = pt_nominal * (1.0 - jes_unc);

    if(pt_nominal>30){
      ++nCleanedJetsPt30;
      if(cleanedJets[i]->userFloat("isBtagged")) ++nCleanedJetsPt30BTagged;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF_Up")) ++nCleanedJetsPt30BTagged_bTagSFUp;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF_Dn")) ++nCleanedJetsPt30BTagged_bTagSFDn;
    }
    if(pt_jes_up>30){
      ++nCleanedJetsPt30_jesUp;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesUp;
    }
    if(pt_jes_dn>30){
      ++nCleanedJetsPt30_jesDn;
      if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jesDn;
    }

    // count jer up/down njets pt30
    float pt_jer_up = cleanedJets[i]->userFloat("pt_jerup");
    float pt_jer_dn = cleanedJets[i]->userFloat("pt_jerdn");
     
    if(pt_jer_up>30){
       ++nCleanedJetsPt30_jerUp;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jerUp;
    }
    if(pt_jer_dn>30){
       ++nCleanedJetsPt30_jerDn;
       if(cleanedJets[i]->userFloat("isBtaggedWithSF")) ++nCleanedJetsPt30BTagged_bTagSF_jerDn;
    }
     
    if (writeJets) FillJet(*(cleanedJets.at(i))); // No additional pT cut (for JEC studies)
  }//End of for (unsigned i=0; i<cleanedJets.size(); ++i)

  // Now we can write the variables for candidates
  int nFilled=0;
  for( edm::View<pat::CompositeCandidate>::const_iterator cand = cands->begin(); cand != cands->end(); ++cand) {
    if (failed) break; //don't waste time on this
    size_t icand= cand-cands->begin();

    if (!( theChannel==ZL || CRFLAG[icand] || (bool)(cand->userFloat("isBestCand")) )) continue; // Skip events other than the best cand (or CR candidates in the CR)

    //For the SR, also fold information about acceptance in CRflag.
    if (isMC && (theChannel==ZZ)) {
      if (gen_ZZ4lInEtaAcceptance)   set_bit(CRFLAG[icand],28);
      if (gen_ZZ4lInEtaPtAcceptance) set_bit(CRFLAG[icand],29);
    }
    FillCandidate(*cand, evtPassTrigger&&evtPassSkim, event, CRFLAG[icand]);

    // Fill the candidate as one entry in the tree. Do not reinitialize the event variables, as in CRs
    // there could be several candidates per event.
    myTree->FillCurrentTree(true);
    ++nFilled;
  }

  // If no candidate was filled but we still want to keep gen-level and weights, we need to fill one entry anyhow.
  if (nFilled==0) {
    if (skipEmptyEvents==false)
      myTree->FillCurrentTree(true);
    else
      myTree->FillCurrentTree(false); //puts it in the failed tree if there is one
  }
}//End of void HZZ4lNtupleMaker::analyze


void HZZ4lNtupleMaker::FillJet(const pat::Jet& jet)
{
   JetPt  .push_back( jet.pt());
   JetEta .push_back( jet.eta());
   JetPhi .push_back( jet.phi());
   JetMass .push_back( jet.p4().M());
   
   JetBTagger .push_back( jet.userFloat("bTagger"));
   JetIsBtagged .push_back( jet.userFloat("isBtagged"));
   JetIsBtaggedWithSF .push_back( jet.userFloat("isBtaggedWithSF"));
   JetIsBtaggedWithSFUp .push_back( jet.userFloat("isBtaggedWithSF_Up"));
   JetIsBtaggedWithSFDn .push_back( jet.userFloat("isBtaggedWithSF_Dn"));
   JetQGLikelihood .push_back( jet.userFloat("qgLikelihood"));
   if(addQGLInputs){
     JetAxis2 .push_back( jet.userFloat("axis2"));
     JetMult .push_back( jet.userFloat("mult"));
     JetPtD .push_back( jet.userFloat("ptD"));
   }
   JetSigma .push_back(jet.userFloat("jes_unc"));
   
   JetRawPt  .push_back( jet.userFloat("RawPt"));
   JetPtJEC_noJER .push_back( jet.userFloat("pt_JEC_noJER"));
   
   JetJESUp .push_back(jet.userFloat("pt_jesup"));
   JetJESDown .push_back(jet.userFloat("pt_jesdn"));

   JetJERUp .push_back(jet.userFloat("pt_jerup"));
   JetJERDown .push_back(jet.userFloat("pt_jerdn"));
    
   JetID.push_back(jet.userFloat("JetID"));
   JetPUID.push_back(jet.userFloat("PUjetID"));

   if (jet.hasUserFloat("pileupJetIdUpdated:fullDiscriminant")) { // if JEC is reapplied, we set this
     JetPUValue.push_back(jet.userFloat("pileupJetIdUpdated:fullDiscriminant"));
   } else {
     JetPUValue.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));
   }
   

   JetHadronFlavour .push_back(jet.hadronFlavour());
   JetPartonFlavour .push_back(jet.partonFlavour());
}//End of void HZZ4lNtupleMaker::FillJet(const pat::Jet& jet)

void HZZ4lNtupleMaker::FillPhoton(int year, const pat::Photon& photon)
{
   PhotonPt  .push_back( photon.pt());
   PhotonEta .push_back( photon.eta());
   PhotonPhi .push_back( photon.phi());
   
   PhotonIsCutBasedLooseID .push_back( PhotonIDHelper::isCutBasedID_Loose(year, photon) );
}

float HZZ4lNtupleMaker::EvalSpline(TSpline3* const& sp, float xval){
  double xmin = sp->GetXmin();
  double xmax = sp->GetXmax();
  double res=0;
  if (xval<xmin){
    res=sp->Eval(xmin);
    double deriv=sp->Derivative(xmin);
    res += deriv*(xval-xmin);
  }
  else if (xval>xmax){
    res=sp->Eval(xmax);
    double deriv=sp->Derivative(xmax);
    res += deriv*(xval-xmax);
  }
  else res=sp->Eval(xval);
  return res;
}//End of float HZZ4lNtupleMaker::EvalSpline()

void HZZ4lNtupleMaker::FillKFactors(edm::Handle<GenEventInfoProduct>& genInfo, std::vector<const reco::Candidate *>& genZLeps){
  KFactor_QCD_ggZZ_Nominal=1;
  KFactor_QCD_ggZZ_PDFScaleDn=1;
  KFactor_QCD_ggZZ_PDFScaleUp=1;
  KFactor_QCD_ggZZ_QCDScaleDn=1;
  KFactor_QCD_ggZZ_QCDScaleUp=1;
  KFactor_QCD_ggZZ_AsDn=1;
  KFactor_QCD_ggZZ_AsUp=1;
  KFactor_QCD_ggZZ_PDFReplicaDn=1;
  KFactor_QCD_ggZZ_PDFReplicaUp=1;
  KFactor_QCD_qqZZ_dPhi=1;
  KFactor_QCD_qqZZ_M=1;
  KFactor_QCD_qqZZ_Pt=1;
  KFactor_EW_qqZZ=1;
  KFactor_EW_qqZZ_unc=0;

  if (isMC){
    GenEventInfoProduct  genInfoP = *(genInfo.product());
    if (apply_K_NNLOQCD_ZZGG>0 && apply_K_NNLOQCD_ZZGG!=3){
      if (spkfactor_ggzz_nnlo[0]!=0) KFactor_QCD_ggZZ_Nominal = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[0], GenHMass);
      if (spkfactor_ggzz_nnlo[1]!=0) KFactor_QCD_ggZZ_PDFScaleDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[1], GenHMass);
      if (spkfactor_ggzz_nnlo[2]!=0) KFactor_QCD_ggZZ_PDFScaleUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[2], GenHMass);
      if (spkfactor_ggzz_nnlo[3]!=0) KFactor_QCD_ggZZ_QCDScaleDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[3], GenHMass);
      if (spkfactor_ggzz_nnlo[4]!=0) KFactor_QCD_ggZZ_QCDScaleUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[4], GenHMass);
      if (spkfactor_ggzz_nnlo[5]!=0) KFactor_QCD_ggZZ_AsDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[5], GenHMass);
      if (spkfactor_ggzz_nnlo[6]!=0) KFactor_QCD_ggZZ_AsUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[6], GenHMass);
      if (spkfactor_ggzz_nnlo[7]!=0) KFactor_QCD_ggZZ_PDFReplicaDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[7], GenHMass);
      if (spkfactor_ggzz_nnlo[8]!=0) KFactor_QCD_ggZZ_PDFReplicaUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nnlo[8], GenHMass);
      if (apply_K_NNLOQCD_ZZGG==2){
        if (spkfactor_ggzz_nlo[0]!=0){
          float divisor = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[0], GenHMass);
          KFactor_QCD_ggZZ_Nominal /= divisor;
          KFactor_QCD_ggZZ_PDFScaleDn /= divisor;
          KFactor_QCD_ggZZ_PDFScaleUp /= divisor;
          KFactor_QCD_ggZZ_QCDScaleDn /= divisor;
          KFactor_QCD_ggZZ_QCDScaleUp /= divisor;
          KFactor_QCD_ggZZ_AsDn /= divisor;
          KFactor_QCD_ggZZ_AsUp /= divisor;
          KFactor_QCD_ggZZ_PDFReplicaDn /= divisor;
          KFactor_QCD_ggZZ_PDFReplicaUp /= divisor;
        }
        else{
          KFactor_QCD_ggZZ_Nominal=0;
          KFactor_QCD_ggZZ_PDFScaleDn=0;
          KFactor_QCD_ggZZ_PDFScaleUp=0;
          KFactor_QCD_ggZZ_QCDScaleDn=0;
          KFactor_QCD_ggZZ_QCDScaleUp=0;
          KFactor_QCD_ggZZ_AsDn=0;
          KFactor_QCD_ggZZ_AsUp=0;
          KFactor_QCD_ggZZ_PDFReplicaDn=0;
          KFactor_QCD_ggZZ_PDFReplicaUp=0;
        }
      }
    }
    else if (apply_K_NNLOQCD_ZZGG==3){
      if (spkfactor_ggzz_nlo[0]!=0) KFactor_QCD_ggZZ_Nominal = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[0], GenHMass);
      if (spkfactor_ggzz_nlo[1]!=0) KFactor_QCD_ggZZ_PDFScaleDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[1], GenHMass);
      if (spkfactor_ggzz_nlo[2]!=0) KFactor_QCD_ggZZ_PDFScaleUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[2], GenHMass);
      if (spkfactor_ggzz_nlo[3]!=0) KFactor_QCD_ggZZ_QCDScaleDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[3], GenHMass);
      if (spkfactor_ggzz_nlo[4]!=0) KFactor_QCD_ggZZ_QCDScaleUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[4], GenHMass);
      if (spkfactor_ggzz_nlo[5]!=0) KFactor_QCD_ggZZ_AsDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[5], GenHMass);
      if (spkfactor_ggzz_nlo[6]!=0) KFactor_QCD_ggZZ_AsUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[6], GenHMass);
      if (spkfactor_ggzz_nlo[7]!=0) KFactor_QCD_ggZZ_PDFReplicaDn = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[7], GenHMass);
      if (spkfactor_ggzz_nlo[8]!=0) KFactor_QCD_ggZZ_PDFReplicaUp = HZZ4lNtupleMaker::EvalSpline(spkfactor_ggzz_nlo[8], GenHMass);
    }

    if (genFinalState!=BUGGY){
      if (genZLeps.size()==4) {
        // Calculate NNLO/NLO QCD K factors for qqZZ
        if (apply_K_NNLOQCD_ZZQQB){
          bool sameflavor=(genZLeps.at(0)->pdgId()*genZLeps.at(1)->pdgId() == genZLeps.at(2)->pdgId()*genZLeps.at(3)->pdgId());
          KFactor_QCD_qqZZ_dPhi = kfactor_qqZZ_qcd_dPhi(fabs(GenZ1Phi-GenZ2Phi), (sameflavor) ? 1 : 2);
          KFactor_QCD_qqZZ_M    = kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 2) / kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 1);
          KFactor_QCD_qqZZ_Pt   = kfactor_qqZZ_qcd_Pt(GenHPt, (sameflavor) ? 1 : 2);
        }
        // Calculate NLO EWK K factors for qqZZ
        if (apply_K_NLOEW_ZZQQB){
          TLorentzVector GENZ1Vec, GENZ2Vec, GENZZVec;
          GENZ1Vec.SetPtEtaPhiM(GenZ1Pt, GenZ1Eta, GenZ1Phi, GenZ1Mass);
          GENZ2Vec.SetPtEtaPhiM(GenZ2Pt, GenZ2Eta, GenZ2Phi, GenZ2Mass);
          GENZZVec = GENZ1Vec + GENZ2Vec;

          KFactor_EW_qqZZ = EwkCorrections::getEwkCorrections(genParticles, ewkTable, genInfoP, GENZ1Vec, GENZ2Vec);

          bool sameflavor=(genZLeps.at(0)->pdgId()*genZLeps.at(1)->pdgId() == genZLeps.at(2)->pdgId()*genZLeps.at(3)->pdgId());
          float K_NNLO_LO = kfactor_qqZZ_qcd_M(GenHMass, (sameflavor) ? 1 : 2, 2);
          float rho = GENZZVec.Pt()/(GenLep1Pt+GenLep2Pt+GenLep3Pt+GenLep4Pt);
          if (rho<0.3) KFactor_EW_qqZZ_unc = fabs((K_NNLO_LO-1.)*(1.-KFactor_EW_qqZZ));
          else KFactor_EW_qqZZ_unc = fabs(1.-KFactor_EW_qqZZ);
        }
      }
    }
  }

}//End of void HZZ4lNtupleMaker::FillKFactors(edm::Handle<GenEventInfoProduct>& genInfo, std::vector<const reco::Candidate *>& genZLeps){


void HZZ4lNtupleMaker::FillLHECandidate(){
  LHEMotherPz.clear();
  LHEMotherE.clear();
  LHEMotherId.clear();
  LHEDaughterPt.clear();
  LHEDaughterEta.clear();
  LHEDaughterPhi.clear();
  LHEDaughterMass.clear();
  LHEDaughterId.clear();
  LHEAssociatedParticlePt.clear();
  LHEAssociatedParticleEta.clear();
  LHEAssociatedParticlePhi.clear();
  LHEAssociatedParticleMass.clear();
  LHEAssociatedParticleId.clear();

  LHEPDFScale = 0;
  genHEPMCweight_POWHEGonly = 0;
  LHEweight_QCDscale_muR1_muF1=0;
  LHEweight_QCDscale_muR1_muF2=0;
  LHEweight_QCDscale_muR1_muF0p5=0;
  LHEweight_QCDscale_muR2_muF1=0;
  LHEweight_QCDscale_muR2_muF2=0;
  LHEweight_QCDscale_muR2_muF0p5=0;
  LHEweight_QCDscale_muR0p5_muF1=0;
  LHEweight_QCDscale_muR0p5_muF2=0;
  LHEweight_QCDscale_muR0p5_muF0p5=0;

  MELACandidate* cand = lheHandler->getBestCandidate();
  if (cand!=0 && addLHEKinematics){
    for (int imot=0; imot<cand->getNMothers(); imot++){
      MELAParticle* apart = cand->getMother(imot);
      if (apart==0){ LHEMotherPz.clear(); LHEMotherE.clear(); LHEMotherId.clear(); break; } // Something went wrong
      LHEMotherPz.push_back(apart->z());
      LHEMotherE.push_back(apart->t());
      LHEMotherId.push_back((short)apart->id);
    }

    for (int iV=0; iV<min(2, cand->getNSortedVs()); iV++){
      MELAParticle* Vi = cand->getSortedV(iV);
      if (Vi!=0){
        for (int iVj=0; iVj<Vi->getNDaughters(); iVj++){
          MELAParticle* Vij = Vi->getDaughter(iVj);
          if (Vij!=0){
            LHEDaughterPt.push_back(Vij->pt());
            if (abs(Vij->pt() / Vij->z()) > 2e-8) {
              LHEDaughterEta.push_back(Vij->eta());
            } else {
              edm::LogWarning("ZeroPt") << "pt = 0!  Using eta = +/-1e10\n" << Vij->id << " " << Vij->x() << " " << Vij->y() << " " << Vij->z() << " " << Vij->t();
              LHEDaughterEta.push_back(copysign(1e10, Vij->z()));
            }
            LHEDaughterPhi.push_back(Vij->phi());
            LHEDaughterMass.push_back(Vij->m());
            LHEDaughterId.push_back((short)Vij->id);
          }
        }
      }
    }

    vector<MELAParticle*> AssociatedParticle;
    vector<MELAParticle*> tmpAssociatedParticle;
    for (int aa=0; aa<cand->getNAssociatedJets(); aa++){
      MELAParticle* apart = cand->getAssociatedJet(aa);
      tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<cand->getNAssociatedLeptons(); aa++){
      MELAParticle* apart = cand->getAssociatedLepton(aa);
      if (!PDGHelpers::isANeutrino(apart->id)) tmpAssociatedParticle.push_back(apart);
    }
    for (int aa=0; aa<cand->getNAssociatedNeutrinos(); aa++){
      MELAParticle* apart = cand->getAssociatedNeutrino(aa);
      tmpAssociatedParticle.push_back(apart);
    }
    while (!tmpAssociatedParticle.empty()){ // Re-sort all associated particles by leading pT (categories are individually sorted, but mixing categories loses this sorting)
      MELAParticle* tmpPart=0;
      int pos=0;
      for (unsigned int el=0; el<tmpAssociatedParticle.size(); el++){
        if (tmpPart==0){ tmpPart = tmpAssociatedParticle.at(el); pos=el; }
        else if (tmpPart->pt()<tmpAssociatedParticle.at(el)->pt()){ tmpPart = tmpAssociatedParticle.at(el); pos=el; } // Safer to do in two steps
      }
      AssociatedParticle.push_back(tmpPart);
      tmpAssociatedParticle.erase(tmpAssociatedParticle.begin()+pos);
    }
    for (unsigned int aa=0; aa<AssociatedParticle.size(); aa++){
      MELAParticle* apart = AssociatedParticle.at(aa);
      if (apart!=0){
        LHEAssociatedParticlePt.push_back(apart->pt());
        if (abs(apart->pt() / apart->z()) > 2e-8) {
          LHEAssociatedParticleEta.push_back(apart->eta());
        } else {
          edm::LogWarning("ZeroPt") << "pt = 0!  Using eta = +/-1e10\n" << apart->id << " " << apart->x() << " " << apart->y() << " " << apart->z() << " " << apart->t();
          LHEAssociatedParticleEta.push_back(copysign(1e10, apart->z()));
        }
        LHEAssociatedParticlePhi.push_back(apart->phi());
        LHEAssociatedParticleMass.push_back(apart->m());
        LHEAssociatedParticleId.push_back((short)apart->id);
      }
    }

    /*
    cout << "NEW EVENT:" << endl;
    cout << "Mothers:" << endl;
    for (unsigned int ipart=0; ipart<LHEMotherId.size(); ipart++) cout << "\t Mot" << ipart << " (pz, E, id) = " << LHEMotherPz.at(ipart) << " " << LHEMotherE.at(ipart) << " " << LHEMotherId.at(ipart) << endl;
    cout << "Daughters:" << endl;
    for (unsigned int ipart=0; ipart<LHEDaughterId.size(); ipart++) cout << "\t Dau" << ipart << " (pt, eta, phi, m, id) = " << LHEDaughterPt.at(ipart) << " " << LHEDaughterEta.at(ipart) << " " << LHEDaughterPhi.at(ipart) << " " << LHEDaughterMass.at(ipart) << " " << LHEDaughterId.at(ipart) << endl;
    cout << "Associated:" << endl;
    for (unsigned int ipart=0; ipart<LHEAssociatedParticleId.size(); ipart++) cout << "\t APart" << ipart << " (pt, eta, phi, m, id) = " << LHEAssociatedParticlePt.at(ipart) << " " << LHEAssociatedParticleEta.at(ipart) << " " << LHEAssociatedParticlePhi.at(ipart) << " " << LHEAssociatedParticleMass.at(ipart) << " " << LHEAssociatedParticleId.at(ipart) << endl;
    cout << endl;
    */

    computeMELABranches(cand);
    pushLHEMELABranches();
  }

  LHEPDFScale = lheHandler->getPDFScale();
  if (genHEPMCweight==1.){
    genHEPMCweight_NNLO = genHEPMCweight = lheHandler->getLHEOriginalWeight();
    if (!printedLHEweightwarning && genHEPMCweight!=1.) {
      printedLHEweightwarning = true;
      edm::LogWarning("InconsistentWeights") << "Gen weight is 1, LHE weight is " << genHEPMCweight;
    }
  }
  genHEPMCweight *= lheHandler->getWeightRescale();

  genHEPMCweight_POWHEGonly = lheHandler->getMemberZeroWeight();
  LHEweight_QCDscale_muR1_muF1 = lheHandler->getLHEWeight(0, 1.);
  LHEweight_QCDscale_muR1_muF2 = lheHandler->getLHEWeight(1, 1.);
  LHEweight_QCDscale_muR1_muF0p5 = lheHandler->getLHEWeight(2, 1.);
  LHEweight_QCDscale_muR2_muF1 = lheHandler->getLHEWeight(3, 1.);
  LHEweight_QCDscale_muR2_muF2 = lheHandler->getLHEWeight(4, 1.);
  LHEweight_QCDscale_muR2_muF0p5 = lheHandler->getLHEWeight(5, 1.);
  LHEweight_QCDscale_muR0p5_muF1 = lheHandler->getLHEWeight(6, 1.);
  LHEweight_QCDscale_muR0p5_muF2 = lheHandler->getLHEWeight(7, 1.);
  LHEweight_QCDscale_muR0p5_muF0p5 = lheHandler->getLHEWeight(8, 1.);
  LHEweight_PDFVariation_Up = lheHandler->getLHEWeight_PDFVariationUpDn(1, 1.);
  LHEweight_PDFVariation_Dn = lheHandler->getLHEWeight_PDFVariationUpDn(-1, 1.);
  LHEweight_AsMZ_Up = lheHandler->getLHEWeigh_AsMZUpDn(1, 1.);
  LHEweight_AsMZ_Dn = lheHandler->getLHEWeigh_AsMZUpDn(-1, 1.);
}//End of void HZZ4lNtupleMaker::FillLHECandidate()


void HZZ4lNtupleMaker::FillCandidate(const pat::CompositeCandidate& cand, bool evtPass, const edm::Event& event, Int_t CRFLAG)
{
  //Initialize a new candidate into the tree
  //myTree->createNewCandidate(); // this doesn't do anything anymore

  //Reinitialize the per-candidate vectors (necessary because in CRs we can store more than 1 candidate per event)
  Jettiness.clear();//clara
  
  LepPt.clear();
  LepEta.clear();
  LepPhi.clear();

  LepPx.clear();//clara
  LepPy.clear();//clara
  LepPz.clear();//clara
  LepE.clear();//clara
  
  LepSCEta.clear();
  LepLepId.clear();
  LepSIP.clear();
  Lepdxy.clear();
  Lepdz.clear();
  LepTime.clear();
  LepisID.clear();
  LepBDT.clear();
  LepisCrack.clear();
  LepMissingHit.clear();
  LepChargedHadIso.clear();
  LepNeutralHadIso.clear();
  LepPhotonIso.clear();
  LepPUIsoComponent.clear();
  LepCombRelIsoPF.clear();
  LepSF.clear();
  LepSF_Unc.clear();
	
  LepScale_Total_Up.clear();
  LepScale_Total_Dn.clear();
  LepScale_Stat_Up.clear();
  LepScale_Stat_Dn.clear();
  LepScale_Syst_Up.clear();
  LepScale_Syst_Dn.clear();
  LepScale_Gain_Up.clear();
  LepScale_Gain_Dn.clear();
  LepSigma_Total_Up.clear();
  LepSigma_Total_Dn.clear();
  LepSigma_Rho_Up.clear();
  LepSigma_Rho_Dn.clear();
  LepSigma_Phi_Up.clear();
  LepSigma_Phi_Dn.clear();

  TLE_dR_Z = -1;
  fsrPt.clear();
  fsrEta.clear();
  fsrPhi.clear();
  fsrLept.clear();
  fsrLeptID.clear();
  fsrDR.clear();
  fsrGenPt.clear();
  ExtraLepPt.clear();
  ExtraLepEta.clear();
  ExtraLepPhi.clear();
  ExtraLepPx.clear();//clara
  ExtraLepPy.clear();//clara
  ExtraLepPz.clear();//clara
  ExtraLepE.clear();//clara
  ExtraLepLepId.clear();

  CRflag = CRFLAG;

  if(theChannel!=ZL){
    //Fill the info on the Higgs candidate
    ZZMass = cand.p4().mass();
    ZZMassErr = cand.userFloat("massError");
    ZZMassErrCorr = cand.userFloat("massErrorCorr");
    ZZMassPreFSR = cand.userFloat("m4l");

    ZZPt  = cand.p4().pt();
    ZZEta = cand.p4().eta();
    ZZPhi = cand.p4().phi();
    ZZjjPt = cand.userFloat("ZZjjPt");
     
    if(addKinRefit){
      if (cand.hasUserFloat("ZZMassRefit")) {
  ZZMassRefit = cand.userFloat("ZZMassRefit");
  ZZMassRefitErr = cand.userFloat("ZZMassRefitErr");
  ZZMassUnrefitErr = cand.userFloat("ZZMassUnrefitErr");
      }
    }
    if(addVtxFit){
      ZZMassCFit = cand.userFloat("CFitM");
      ZZChi2CFit = cand.userFloat("CFitChi2");
    }

    DiJetMass  = cand.userFloat("DiJetMass");
    DiJetDEta  = cand.userFloat("DiJetDEta");
    DiJetFisher  = cand.userFloat("DiJetFisher");
    //    DiJetMassPlus  = cand.userFloat("DiJetMassPlus");
    //    DiJetMassMinus  = cand.userFloat("DiJetMassMinus");

    //Fill the angular variables
    helcosthetaZ1 = cand.userFloat("costheta1");
    helcosthetaZ2 = cand.userFloat("costheta2");
    helphi       = cand.userFloat("phi");
    costhetastar = cand.userFloat("costhetastar");
    phistarZ1      = cand.userFloat("phistar1");
    //phistarZ2      = cand.userFloat("phistar2");
    xi            = cand.userFloat("xi");
    xistar        = cand.userFloat("xistar");

    // Get MELA probabilities
    pushRecoMELABranches(cand);

  }//End of if(theChannel!=ZL)

  //Z1 and Z2 variables
  const reco::Candidate* Z1;
  const reco::Candidate* Z2;
  vector<const reco::Candidate*> jet;//Clara
  vector<const reco::Candidate*> leptons;
  vector<const reco::Candidate*> fsrPhot;
  vector<short> fsrIndex;
  vector<string> labels;

  if (theChannel!=ZL) { // Regular 4l candidates
    Z1   = cand.daughter("Z1");
    Z2   = cand.daughter("Z2");
    userdatahelpers::getSortedLeptons(cand, leptons, labels, fsrPhot, fsrIndex);
  } else {              // Special handling of Z+l candidates
    Z1   = cand.daughter(0); // the Z
    Z2   = cand.daughter(1); // This is actually the additional lepton!
    userdatahelpers::getSortedLeptons(cand, leptons, labels, fsrPhot, fsrIndex, false); // note: we get just 3 leptons in this case.
  }

  Z1Mass = Z1->mass();
  Z1Pt =   Z1->pt();
  Z1Flav =  getPdgId(Z1->daughter(0)) * getPdgId(Z1->daughter(1));
 
  Z2Mass = Z2->mass();
  Z2Pt =   Z2->pt();
  Z2Flav = theChannel==ZL ? getPdgId(Z2) : getPdgId(Z2->daughter(0)) * getPdgId(Z2->daughter(1));

  const reco::Candidate* non_TLE_Z = nullptr;
  size_t TLE_index = 999;
  if(abs(Z1Flav) == 11*11 || abs(Z1Flav) == 13*13) non_TLE_Z = Z1;
  if(abs(Z2Flav) == 11*11 || abs(Z2Flav) == 13*13) non_TLE_Z = Z2;
  for (size_t i=0; i<leptons.size(); ++i){
    if(abs(leptons[i]->pdgId()) == 22) TLE_index = i;
  }
  if(TLE_index < 999 && non_TLE_Z != nullptr) {
    TLE_dR_Z = reco::deltaR(non_TLE_Z->p4(), leptons[TLE_index]->p4()); 
  }

  Int_t sel = 0;
  if(theChannel==ZZ){

    // Precomputed selections
    bool candPass70Z2Loose = cand.userFloat("Z2Mass") &&
                             cand.userFloat("MAllComb") &&
                             cand.userFloat("pt1")>20 && cand.userFloat("pt2")>10. &&
                             ZZMass>70.;
    bool candPassFullSel70 = cand.userFloat("SR");
    bool candPassFullSel   = cand.userFloat("FullSel");
    bool candIsBest = cand.userFloat("isBestCand");
    bool passMz_zz = (Z1Mass>60. && Z1Mass<120. && Z2Mass>60. && Z2Mass<120.);   //FIXME hardcoded cut

    if (candIsBest) {
      //    sel = 10; //FIXME see above
      if (candPass70Z2Loose) sel=70;
      if (candPassFullSel70){ // includes MZ2 > 12
        sel = 90;
        if (candPassFullSel){
          sel=100;
          if (passMz_zz) sel = 120;
        }
      }
    }
  } else if(theChannel==ZLL) sel = 20; 
  else if(theChannel==ZL) sel = 10;
 

  if (!(evtPass)) {sel = -sel;} // avoid confusion when we write events which do not pass trigger/skim

  // Retrieve the userFloat of the leptons in vectors ordered in the same way.
  vector<float> SIP(4);
  vector<float> combRelIsoPF(4);
  passIsoPreFSR = true;
  

  for (unsigned int i=0; i<leptons.size(); ++i){
    float curr_dR = 999;
    if(i != TLE_index && TLE_index < 999)
      curr_dR = reco::deltaR(leptons[i]->p4(), leptons[TLE_index]->p4());
    if(curr_dR < TLE_min_dR_3l) TLE_min_dR_3l = curr_dR;

    short lepFlav = std::abs(leptons[i]->pdgId());

    SIP[i]             = userdatahelpers::getUserFloat(leptons[i],"SIP");
    passIsoPreFSR      = passIsoPreFSR&&(userdatahelpers::getUserFloat(leptons[i],"combRelIsoPF")<LeptonIsoHelper::isoCut(leptons[i]));

    //in the Legacy approach,  FSR-corrected iso is attached to the Z, not to the lepton!
    if (theChannel!=ZL) {
      combRelIsoPF[i]    = cand.userFloat(labels[i]+"combRelIsoPFFSRCorr"); // Note: the
      assert(SIP[i] == cand.userFloat(labels[i]+"SIP")); // Check that I don't mess up with labels[] and leptons[]
    } else {
      //FIXME: at the moment,  FSR-corrected iso is not computed for Z+L CRs
      combRelIsoPF[i]    = userdatahelpers::getUserFloat(leptons[i],"combRelIsoPF");
    }

    //Fill the info on the lepton candidates
    
    LepPt .push_back( leptons[i]->pt() );
    LepEta.push_back( leptons[i]->eta() );
    LepPhi.push_back( leptons[i]->phi() );
    LepPx.push_back( leptons[i]->px() );//clara
    LepPy.push_back( leptons[i]->py() );//clara
    LepPz.push_back( leptons[i]->pz() );//clara
    LepE.push_back( leptons[i]->p4().E());//clara
    
    LepSCEta.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"SCeta") : -99. );
    int id =  leptons[i]->pdgId();
    if(id == 22 && (i == 1 || i == 3)) id=-22; //FIXME this assumes a standard ordering of leptons.
    LepLepId.push_back( id );
    LepSIP  .push_back( SIP[i] );
    Lepdxy  .push_back( userdatahelpers::getUserFloat(leptons[i],"dxy") );
    Lepdz   .push_back( userdatahelpers::getUserFloat(leptons[i],"dz") );
    LepTime .push_back( lepFlav==13 ? userdatahelpers::getUserFloat(leptons[i],"time") : 0. );
    LepisID .push_back( userdatahelpers::getUserFloat(leptons[i],"ID") );
    LepBDT  .push_back( userdatahelpers::getUserFloat(leptons[i],"BDT") );
    LepisCrack.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"isCrack") : 0 );
    LepMissingHit.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"missingHit") : 0 );
    LepScale_Total_Up.push_back( userdatahelpers::getUserFloat(leptons[i],"scale_total_up") );
    LepScale_Total_Dn.push_back( userdatahelpers::getUserFloat(leptons[i],"scale_total_dn") );
    LepScale_Stat_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_stat_up") : -99. );
    LepScale_Stat_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_stat_dn") : -99. );
    LepScale_Syst_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_syst_up") : -99. );
    LepScale_Syst_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_syst_dn") : -99. );
    LepScale_Gain_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_gain_up") : -99. );
    LepScale_Gain_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"scale_gain_dn") : -99. );
    LepSigma_Total_Up.push_back( userdatahelpers::getUserFloat(leptons[i],"sigma_total_up") );
    LepSigma_Total_Dn.push_back( userdatahelpers::getUserFloat(leptons[i],"sigma_total_dn") );
    LepSigma_Rho_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"sigma_rho_up") : -99. );
    LepSigma_Rho_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"sigma_rho_dn") : -99. );
    LepSigma_Phi_Up.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"sigma_phi_up") : -99. );
    LepSigma_Phi_Dn.push_back( lepFlav==11 ? userdatahelpers::getUserFloat(leptons[i],"sigma_phi_dn") : -99. );
    LepChargedHadIso.push_back( userdatahelpers::getUserFloat(leptons[i],"PFChargedHadIso") );
    LepNeutralHadIso.push_back( userdatahelpers::getUserFloat(leptons[i],"PFNeutralHadIso") );
    LepPhotonIso.push_back( userdatahelpers::getUserFloat(leptons[i],"PFPhotonIso") );
	  LepPUIsoComponent.push_back( lepFlav==13 ? userdatahelpers::getUserFloat(leptons[i],"PFPUChargedHadIso") : 0. );
    LepCombRelIsoPF.push_back( combRelIsoPF[i] );
    LepisLoose.push_back(userdatahelpers::hasUserFloat(leptons[i],"isLoose") == 1 ? userdatahelpers::getUserFloat(leptons[i],"isLoose") : -2);

  }//End of for (unsigned int i=0; i<leptons.size(); ++i) 

  // FSR
  for (unsigned i=0; i<fsrPhot.size(); ++i) {
    math::XYZTLorentzVector fsr = fsrPhot[i]->p4();
    fsrPt.push_back(fsr.pt());
    fsrEta.push_back(fsr.eta());
    fsrPhi.push_back(fsr.phi());
    fsrLept.push_back(fsrIndex[i]+1);
    fsrLeptID.push_back(leptons[fsrIndex[i]]->pdgId());
    fsrDR.push_back(ROOT::Math::VectorUtil::DeltaR(leptons[fsrIndex[i]]->momentum(), fsrPhot[i]->momentum()));
    int igen = MCHistoryTools::fsrMatch(fsrPhot[i], genFSR);
    double dRGenVsReco = -1.;
    double genpT = -1.;

    if (igen>=0) {
      dRGenVsReco = ROOT::Math::VectorUtil::DeltaR(genFSR[igen]->momentum(), fsrPhot[i]->momentum());
//       pTGen = genFSR[igen]->pt();
//       etaGen = genFSR[igen]->eta();
//       phiGen = genFSR[igen]->phi();
      if (dRGenVsReco<0.3) {// matching cut - FIXME
        genpT=genFSR[igen]->pt();
      }
    }
    fsrGenPt.push_back(genpT);


  }

  //convention: 0 -> 4mu   1 -> 4e   2 -> 2mu2e
  if(CRFLAG){
    ZXFakeweight = 1.;
    for(int izx=0;izx<2;izx++)
      ZXFakeweight *= getFakeWeight(Z2->daughter(izx)->pt(),Z2->daughter(izx)->eta(),Z2->daughter(izx)->pdgId(),Z1->daughter(0)->pdgId());
  }
  ZZsel = sel;

  if (theChannel!=ZL) {

    //Fill the info on categorization
    nExtraLep = cand.userFloat("nExtraLep"); // Why is this still a float at this point?
    nExtraZ=cand.userFloat("nExtraZ");

    //Fill the info on the extra leptons
    for (int iExtraLep=1; iExtraLep<=(int)nExtraLep; iExtraLep++){
      TString extraString;extraString.Form("ExtraLep%d",iExtraLep);
      if (cand.hasUserCand(extraString.Data())){
        //for(int iextra=0;iextra<4;iextra++)varExtra[iextra].Prepend(extraString.Data());
        reco::CandidatePtr candPtr=cand.userCand(extraString.Data());

        ExtraLepPt.push_back(candPtr->pt());
        ExtraLepEta.push_back(candPtr->eta());
        ExtraLepPhi.push_back(candPtr->phi());
        ExtraLepPx.push_back(candPtr->px());//clara
        ExtraLepPy.push_back(candPtr->py());//clara
        ExtraLepPz.push_back(candPtr->pz());//clara
        ExtraLepE.push_back(candPtr->p4().E());//clara
        }
    }
  }

   //Compute the data/MC weight
   dataMCWeight = 1.;
 //When the trigger is not applied in the MC, apply a trigger efficiency factor instead (FIXME: here hardcoding the efficiencies computed for ICHEP2016)
  trigEffWeight = 1.;
  if(isMC) {
    
    dataMCWeight = getAllWeight(leptons);
    
    if (applyTrigEffWeight){
      Int_t ZZFlav = abs(Z1Flav*Z2Flav);
      if(ZZFlav==121*121 || ZZFlav==121*242) //4e
	trigEffWeight = 0.992;
      if(ZZFlav==169*169) //4mu
	trigEffWeight = 0.996;
      if(ZZFlav==169*121 || ZZFlav==169*242) //2e2mu
	trigEffWeight = 0.995;
    }
  }
  //Store an overall event weight (which is normalized by gen_sumWeights)
  overallEventWeight = PUWeight * genHEPMCweight * dataMCWeight * trigEffWeight;

  /* // debug printout for weights
  cout<<"Event "<<event.id().run()<<":"<<event.luminosityBlock()<<":"<<event.id().event()<<endl;
 cout<<" pileup weight =         "<<PUWeight<<endl;
  cout<<" sign of gen. weight =   "<<genHEPMCweight/fabs(genHEPMCweight)<<endl;
  cout<<" lepton data/MC weight = "<<dataMCWeight<<endl;
  for(unsigned int i=0; i<leptons.size(); ++i)
    cout<<"   lepton ID="<<leptons[i]->pdgId()<<", pT="<<leptons[i]->pt()<<", weight="<<getAllWeight(leptons[i])<<endl;
  cout<<" trigger eff. weight =   "<<trigEffWeight<<endl;
  cout<<"product of all =         "<<overallEventWeight/fabs(genHEPMCweight)<<endl;
  cout<<endl;
 */ 
  
  
// ******************************************
// *       <Jettiness Algorithem>           *
// *  Author: KIM Mi Ran                    *
// *  Sunkyunkwan University,South Korea    *
// *  26th,March,2021                       *
// ******************************************
  Int_t nJet=JetPt.size();
  Int_t nLep=LepPx.size();
  Int_t nExtraLep=ExtraLepPx.size();
  Float_t e_mass=0.00051;//GeV/c^2
  Float_t mu_mass=0.10566;//GeV/c^2
 
  TLorentzVector v4ExtraLep[nExtraLep];
  TLorentzVector v4Lep[nLep];
  TLorentzVector v4Jet[nJet];
  TLorentzVector v4_allJetSum[nJet+1];    
  
  TLorentzVector  v4_ExtraLepAdd;
  TLorentzVector  v4_LepAddTot;
  TLorentzVector  v4_allLepSum;
  TLorentzVector  v4_allLepJetSum[nJet+1];    
 
  //Initialization
  v4_allLepSum.SetPxPyPzE(0,0,0,0);    
  for(Int_t i=0;i<nJet+1;i++){
  v4_allJetSum[i].SetPxPyPzE(0,0,0,0);     
  v4_allLepJetSum[i].SetPxPyPzE(0,0,0,0);     
  }
  
  //----Beam defintion----//
  Float_t Tau[nJet+1][nJet+2];//Tau elements array
  Float_t TauSum[nJet+1];//min valus Sum of Tau elements
  TLorentzVector v4_qa[nJet+1];//Beam=Scalar products v4 sum of jet+lepton and the direction v4_na:Z
  TLorentzVector v4_qb[nJet+1];//Beam=Scalra products v4 sum of jet+lepton and the direction v4_nb:-Z
  Float_t QQ[nJet+1];//Energy Normalization Factor

/*
 //--v4ExtraLepton---//
  for(Int_t i=0;i<nExtraLep;i++){
      v4ExtraLep[i].SetPxPyPzE(ExtraLepPx.at(i),ExtraLepPy.at(i),ExtraLepPz.at(i),ExtraLepE.at(i));
   v4_ExtraLepAdd +=v4ExtraLep[i];//ExtraLep Sum
  }//End of for(Int_t i=0;i<nExtraLep;i++)
 */

  //---Original Method for ExtraLepton---//
  for(Int_t i=0;i<nExtraLep;i++){
    if(abs(ExtraLepLepId.at(i))==11){
      TLorentzVector v4_PEPM;//v4 Container of PtEtaPhiM values
     v4_PEPM.Clear();
      v4_PEPM.SetPtEtaPhiM(ExtraLepPt.at(i),ExtraLepEta.at(i),ExtraLepPhi.at(i),e_mass);
      v4ExtraLep[i].SetPxPyPzE(v4_PEPM.Px(),v4_PEPM.Py(),v4_PEPM.Pz(),v4_PEPM.E());  
  }//End of (abs(ExtraLepLepId->at(i))==11)

    if(abs(ExtraLepLepId.at(i))==13){
      TLorentzVector v4_PEPM;//v4 Container of PtEtaPhiM values
      v4_PEPM.Clear();
      v4_PEPM.SetPtEtaPhiM(ExtraLepPt.at(i),ExtraLepEta.at(i),ExtraLepPhi.at(i),mu_mass);
      v4ExtraLep[i].SetPxPyPzE(v4_PEPM.Px(),v4_PEPM.Py(),v4_PEPM.Pz(),v4_PEPM.E());
    }//End of if(abs(ExtraLepLepId.at(i))==13)
  
   v4_ExtraLepAdd +=v4ExtraLep[i];//ExtraLep Sum
 }//End of for(Int_t i=0;i<nExtraLep;i++)


/*
//--4LeptionsSum--// 
  for(Int_t i=0;i<nLep;i++){
      v4Lep[i].SetPxPyPzE(LepPx.at(i),LepPy.at(i),LepPz.at(i),LepE.at(i));
      v4_allLepSum +=v4Lep[i];      
  }//End of for(Int_t i=0;i<nExtraLep;i++)
*/


  //---Original Method for Lepton---//
  for(Int_t i=0;i<nLep;i++){
      if(abs(LepLepId.at(i))==11){
      TLorentzVector v4_PEPM;//v4 Container of PtEtaPhiM values
      v4_PEPM.Clear();
      v4_PEPM.SetPtEtaPhiM(LepPt.at(i),LepEta.at(i),LepPhi.at(i),e_mass);
      v4Lep[i].SetPxPyPzE(v4_PEPM.Px(),v4_PEPM.Py(),v4_PEPM.Pz(),v4_PEPM.E());
    }//End of (abs(ExtraLepLepId->at(i))==11)

    if(abs(LepLepId.at(i))==13){
      TLorentzVector v4_PEPM;//v4 Container of PtEtaPhiM values
      v4_PEPM.Clear();
      v4_PEPM.SetPtEtaPhiM(LepPt.at(i),LepEta.at(i),LepPhi.at(i),mu_mass);
      v4Lep[i].SetPxPyPzE(v4_PEPM.Px(),v4_PEPM.Py(),v4_PEPM.Pz(),v4_PEPM.E());
    }//End of if(abs(ExtraLepLepId->at(i))==13)
  v4_allLepSum +=v4Lep[i];//Lepton Sum
  }//End of for (Int_t i=0;i<nLep;i++)

//---v4_Jet original Method-----//
  for(Int_t i=0;i<nJet;i++){
      TLorentzVector v4_PEPM;//v4 Containner of PtEtaPhiM values
      v4_PEPM.Clear();
      v4_PEPM.SetPtEtaPhiM(JetPt.at(i),JetEta.at(i),JetPhi.at(i),JetMass.at(i));
      v4Jet[i].SetPxPyPzE(v4_PEPM.Px(),v4_PEPM.Py(),v4_PEPM.Pz(),v4_PEPM.E());
      
        if(i==0){v4_allJetSum[0]+=v4Jet[0];}
        if(i>0){v4_allJetSum[i]+=v4_allJetSum[i-1]+v4Jet[i];}
  }//End of for (Int_t i=0;i<nJet;i++)
 
 //---4Leptons+Jet---//
  v4_LepAddTot=v4_ExtraLepAdd+v4_allLepSum;//USE IT ONLY WHEN INCLUDE ADDITIONAL LEPTONS FOR QQ
  v4_allLepJetSum[0]=v4_allLepSum;//nJet=0
 
  for(Int_t i=1;i<nJet+1;i++){//nJet!=0
    v4_allLepJetSum[i]=v4_allLepSum +v4_allJetSum[i-1];//All signal Lepton+Njets
  }

  //---qa Beam vector---//
  v4_qa[0]=(v4_allLepJetSum[0]*v4_nb)*v4_na;//Beam
  for(Int_t i=1;i<nJet+1;i++){
    v4_qa[i]=(v4_allLepJetSum[i]*v4_nb)*v4_na;
  }
  //---qb Beam vector---//
  v4_qb[0]=(v4_allLepJetSum[0]*v4_na)*v4_nb;//Beam
  for(Int_t i=1;i<nJet+1;i++){
    v4_qb[i]=(v4_allLepJetSum[i]*v4_na)*v4_nb;
  }
 
  //---QQ:Normalization of Energy Factor---//
  QQ[0]=(v4_allLepJetSum[0]*v4_nb)*(v4_na*v4_allLepJetSum[0]);
  for(Int_t i=1;i<nJet+1;i++){
     QQ[i]=(v4_allLepJetSum[i]*v4_nb)*(v4_na*v4_allLepJetSum[i]);
  }

  //---Initialization of TauSum[nJet+1]---//
  for(Int_t i=0; i<nJet+1;i++){
    TauSum[i]=0.0;
  }//End of for(Int_t i=0; i<nJet+1;i++)

//---------(nExtraLe=0) calculation of Jettiness.at(0) to Jettiness.at(2)-----
  if(nExtraLep==0){
    for(Int_t jn=0;jn<nJet+1;jn++){
      for(Int_t j=jn;j<nJet;j++){
        Tau[jn][0]=v4_qa[jn]*v4Jet[j];
        Tau[jn][1]=v4_qb[jn]*v4Jet[j];
        for(Int_t k=0;k<jn;k++){
          Tau[jn][k+2]=v4Jet[k]*v4Jet[j];
        }//End of for (Int_t k=0;k<jn;k++)  
	Float_t min3;
        min3=Tau[jn][0];
        //-----Find min elemets of Tau------
        for (Int_t m=0;m<jn+2;m++){
          if(Tau[jn][m]<min3) min3=Tau[jn][m];
        }//End of for (Int_t m=0;m<jn+2;m++)
        TauSum[jn] +=min3;
      }//End of for (Int_t j=jn;j<nJet;j++)
      if(jn==nJet) TauSum[jn]=0.;
      
      Jettiness.push_back(TauSum[jn]/QQ[jn]);
    }//End of for (Int_t jn=0;jn<nJet+1;jn++) 
  }//End of if(nExtraLep==0)
 
//---------(nExtraLep!=0) calculation of Jettiness.at(0) to Jettiness.at(2)-----

  if(nExtraLep!=0){ 
    for(Int_t jn=0;jn<nJet+1;jn++){
      for(Int_t j=jn;j<nJet;j++){
        Tau[jn][0]=v4_qa[jn]*v4Jet[j];
        Tau[jn][1]=v4_qb[jn]*v4Jet[j];
        for(Int_t k=0;k<jn;k++){
          Tau[jn][k+2]=v4Jet[k]*v4Jet[j];
        }//End of for (Int_t k=0;k<jn;k++)  
        Float_t min3;
        min3=Tau[jn][0];
        //-----Find min elemets of Tau------
        for (Int_t m=0;m<jn+2;m++){
          if(Tau[jn][m]<min3) min3=Tau[jn][m];
        }//End of for (Int_t m=0;m<jn+2;m++)
        TauSum[jn] +=min3;
      }//End of for (Int_t j=jn;j<nJet;j++)
      for (Int_t j=0;j<nExtraLep;j++){
        Tau[jn][0]=v4_qa[jn]*v4ExtraLep[j];
        Tau[jn][1]=v4_qb[jn]*v4ExtraLep[j];
        for(Int_t k=0;k<jn;k++){
          Tau[jn][k+2]=v4Jet[k]*v4ExtraLep[j];
        }//End of for(Int_t k=0;k<jn;k++)
      Float_t min3;
      min3=Tau[jn][0];
        for(Int_t m=0;m<jn+2;m++){
          if(Tau[jn][m]<min3) min3=Tau[jn][m];
          //mintau[jn].push_back(min3);
        }//End of for (Int_t m=0;m<jn+2;m++)
      TauSum[jn]+=min3; 
      }//End of for(Int_t j=0;j<ExtraLep_N;j++)
      
      Jettiness.push_back(TauSum[jn]/QQ[jn]);
    }//End of for (Int_t jn=0;jn<nJet+1;jn++) 
  }//End of if(nExtraLep!=0)
 //}//End of if(EventNumber==68417 && LumiNumber==12)
 

}//End of void HZZ4lNtupleMaker::FillCandidat

void HZZ4lNtupleMaker::getCheckedUserFloat(const pat::CompositeCandidate& cand, const std::string& strval, Float_t& setval, Float_t defaultval){
  if (cand.hasUserFloat(strval)) setval = cand.userFloat(strval);
  else setval = defaultval;
}



// ------------ method called once each job just before starting event loop  ------------
void HZZ4lNtupleMaker::beginJob()
{
  edm::Service<TFileService> fs;
  TTree *candTree = fs->make<TTree>(theFileName,"Event Summary");
  TTree *candTree_failed = 0;
  if (failedTreeLevel)
    candTree_failed = fs->make<TTree>(theFileName+"_failed","Event Summary");
  myTree = new HZZ4lNtupleFactory(candTree, candTree_failed);
  const int nbins = 45;
  hCounter = fs->make<TH1F>("Counters", "Counters", nbins, 0., nbins);
  BookAllBranches();
  buildMELABranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void HZZ4lNtupleMaker::endJob()
{
  hCounter->SetBinContent(0 ,gen_sumWeights); // also stored in bin 40
  hCounter->SetBinContent(1 ,Nevt_Gen-gen_BUGGY);
  hCounter->SetBinContent(2 ,gen_ZZ4mu);
  hCounter->SetBinContent(3 ,gen_ZZ4e);
  hCounter->SetBinContent(4 ,gen_ZZ2mu2e);
  hCounter->SetBinContent(5 ,gen_ZZ2l2tau);
  hCounter->SetBinContent(6 ,gen_ZZ4mu_EtaAcceptance);
  hCounter->SetBinContent(7 ,gen_ZZ4mu_LeptonAcceptance);
  hCounter->SetBinContent(8 ,gen_ZZ2emu2tau);
  hCounter->SetBinContent(9 ,gen_ZZ4tau);
  hCounter->SetBinContent(10,gen_ZZ4e_EtaAcceptance);
  hCounter->SetBinContent(11,gen_ZZ4e_LeptonAcceptance);
  hCounter->SetBinContent(14,gen_ZZ2mu2e_EtaAcceptance);
  hCounter->SetBinContent(15,gen_ZZ2mu2e_LeptonAcceptance);
  hCounter->SetBinContent(19,gen_BUGGY);
  hCounter->SetBinContent(20,gen_Unknown);

  hCounter->SetBinContent(40,gen_sumWeights); // Also stored in underflow bin; added here for convenience
  hCounter->SetBinContent(41,gen_sumGenMCWeight);
  hCounter->SetBinContent(42,gen_sumPUWeight);

  TH1 *h[1] ={ hCounter };
  for (int i = 0; i < 1; i++) {
    h[i]->GetXaxis()->SetBinLabel(1 ,"Nevt_Gen");
    h[i]->GetXaxis()->SetBinLabel(2 ,"gen_ZZ4mu");
    h[i]->GetXaxis()->SetBinLabel(3 ,"gen_ZZ4e");
    h[i]->GetXaxis()->SetBinLabel(4 ,"gen_ZZ2mu2e");
    h[i]->GetXaxis()->SetBinLabel(5 ,"gen_ZZ2l2tau");
    h[i]->GetXaxis()->SetBinLabel(6 ,"gen_ZZ4mu_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(7 ,"gen_ZZ4mu_LeptonAcceptance");
    h[i]->GetXaxis()->SetBinLabel(8 ,"gen_ZZ2emu2tau");
    h[i]->GetXaxis()->SetBinLabel(9 ,"gen_ZZ4tau");
    h[i]->GetXaxis()->SetBinLabel(10,"gen_ZZ4e_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(11,"gen_ZZ4e_LeptonAcceptance");
    h[i]->GetXaxis()->SetBinLabel(14,"gen_ZZ2mu2e_EtaAcceptance");
    h[i]->GetXaxis()->SetBinLabel(15,"gen_ZZ2mu2e_LeptonAcceptance");
    h[i]->GetXaxis()->SetBinLabel(19,"gen_BUGGY");
    h[i]->GetXaxis()->SetBinLabel(20,"gen_Unknown");

    h[i]->GetXaxis()->SetBinLabel(40,"gen_sumWeights");
    h[i]->GetXaxis()->SetBinLabel(41,"gen_sumGenMCWeight");
    h[i]->GetXaxis()->SetBinLabel(42,"gen_sumPUWeight");
  }

  delete myTree;

  return;
}

// ------------ method called when starting to processes a run  ------------
void HZZ4lNtupleMaker::beginRun(edm::Run const& iRun, edm::EventSetup const&)
{
  static bool firstRun=true;
  if (firstRun){
    if (lheHandler){
      edm::Handle<LHERunInfoProduct> lhe_runinfo;
      iRun.getByLabel(edm::InputTag("externalLHEProducer"), lhe_runinfo);
      lheHandler->setHeaderFromRunInfo(&lhe_runinfo);
    }
    firstRun=false;
  }
}

// ------------ method called when ending the processing of a run  ------------
void HZZ4lNtupleMaker::endRun(edm::Run const& iRun, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void HZZ4lNtupleMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
  Nevt_Gen_lumiBlock = 0;
}

// ------------ method called when ending the processing of a luminosity block  ------------
void HZZ4lNtupleMaker::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
  Float_t Nevt_preskim = -1.;
  edm::Handle<edm::MergeableCounter> preSkimCounter;
  if (iLumi.getByToken(preSkimToken, preSkimCounter)) { // Counter before skim. Does not exist for non-skimmed samples.
    Nevt_preskim = preSkimCounter->value;
    // We do not use a filtering skim for the time being; so this is just left as a check in case we need it again in the future.
    if (!std::uncaught_exception() && Nevt_preskim>=0.) assert(Nevt_preskim == Nevt_Gen_lumiBlock);
  }

  Nevt_Gen += Nevt_Gen_lumiBlock;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HZZ4lNtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


Float_t HZZ4lNtupleMaker::getAllWeight(const vector<const reco::Candidate*>& leptons) 
{
  Float_t totWeight = 1.;
	
  for(unsigned int i=0; i<leptons.size(); ++i){ 
    Int_t   myLepID = abs(leptons[i]->pdgId());
    if (skipMuDataMCWeight&& myLepID==13) return 1.;
    if (skipEleDataMCWeight&& myLepID==11) return 1.;
    
    float SF = 1.0;
    float SF_Unc = 0.0;
  
    Float_t myLepPt = leptons[i]->pt();
    Float_t myLepEta = leptons[i]->eta();
     
    Float_t SCeta;
    if (myLepID == 11) SCeta = userdatahelpers::getUserFloat(leptons[i],"SCeta");
    else SCeta = myLepEta;
    
    Float_t mySCeta;
     
    // Deal with very rare cases when SCeta is out of 2.5 bonds
    if ( myLepEta <= 2.5 && SCeta >= 2.5) mySCeta = 2.49;
    else if ( myLepEta >= -2.5 && SCeta <= -2.5) mySCeta = -2.49;
    else mySCeta = SCeta;
     
    bool isCrack;
    if (myLepID == 11) isCrack = userdatahelpers::getUserFloat(leptons[i],"isCrack");
    else isCrack = false;
     
 
    SF = lepSFHelper->getSF(year,myLepID,myLepPt,myLepEta, mySCeta, isCrack);
    SF_Unc = lepSFHelper->getSFError(year,myLepID,myLepPt,myLepEta, mySCeta, isCrack);

    LepSF.push_back(SF);
    LepSF_Unc.push_back(SF_Unc);

    totWeight *= SF;
  } 

  return totWeight;
}


Float_t HZZ4lNtupleMaker::getHqTWeight(double mH, double genPt) const
{
  if (skipHqTWeight) return 1.;

  //cout<<"mH = "<<mH<<", genPt = "<<genPt<<endl;
  if (mH<400 || genPt>250) return 1.;

  double weight = 1.;

  const int masses[4] = {400,600,800,1000};
  double massDiff = 1000;
  int iMass = -1;
  for (int i=0; i<4; ++i){
    double massDiffTmp = std::fabs(mH-masses[i]);
    if (massDiffTmp<massDiff){
      massDiff = massDiffTmp;
      iMass = i;
    }
  }

  if (iMass>=0) {
    weight = h_weight->GetBinContent(h_weight->FindBin(genPt));
  }
  return weight;
}


// Added by CO
Float_t HZZ4lNtupleMaker::getFakeWeight(Float_t LepPt, Float_t LepEta, Int_t LepID, Int_t LepZ1ID)
{
  if (skipFakeWeight) return 1.;

  // year 0 = 2011
  // year 1 = 2012

  Float_t weight  = 1.;

  Int_t nHisto=0;

  Float_t myLepPt   = LepPt;
  Float_t myLepEta  = fabs(LepEta);
  Int_t   myLepID   = abs(LepID);
  Int_t   myZ1LepID = abs(LepZ1ID);

  //cout << " pt = " << myLepPt << " eta = " << myLepEta << " ZID = " << myZ1LepID << " LepID = " << myLepID << endl;

  //avoid to go out of the TH boundary
  if(myLepPt > 79.) myLepPt = 79.;
  if(myLepID==13)nHisto+=2;
  if(myZ1LepID==13)nHisto+=1;


  TString Z1flavor = "Z1ee";     if(myZ1LepID==13) Z1flavor = "Z1mumu";
  TString Z2flavor = "electron"; if(myLepID==13)   Z2flavor = "muon";
  TString histo_name = "eff_"+Z1flavor+"_plus_"+Z2flavor;

  //cout << " histo = " << histo_name << endl;

  weight = h_ZXWeight[nHisto]->GetBinContent(h_ZXWeight[nHisto]->GetXaxis()->FindBin(myLepPt), h_ZXWeight[nHisto]->GetYaxis()->FindBin(myLepEta));
  /*
  h_ZXWeight[0]=FileZXWeightEle->Get("eff_Z1ee_plus_electron")->Clone();
  h_ZXWeight[2]=FileZXWeightEle->Get("eff_Z1mumu_plus_electron")->Clone();

  h_ZXWeight[5]=FileZXWeightMuo->Get("eff_Z1ee_plus_muon")->Clone();
  h_ZXWeight[7]=FileZXWeightMuo->Get("eff_Z1mumu_plus_muon")->Clone();
*/
  return weight;

} // end of getFakeWeight
void HZZ4lNtupleMaker::FillZGenInfo(Short_t Z1Id, Short_t Z2Id,
                                    const math::XYZTLorentzVector pZ1, const math::XYZTLorentzVector pZ2)
{
  GenZ1Mass= pZ1.M();
  GenZ1Pt= pZ1.Pt();
  GenZ1Eta= pZ1.Eta();
  GenZ1Phi= pZ1.Phi();
  GenZ1Flav= Z1Id;

  GenZ2Mass= pZ2.M();
  GenZ2Pt= pZ2.Pt();
  GenZ2Eta= pZ2.Eta();
  GenZ2Phi= pZ2.Phi();
  GenZ2Flav= Z2Id;

  return;
}

void HZZ4lNtupleMaker::FillLepGenInfo(Short_t Lep1Id, Short_t Lep2Id, Short_t Lep3Id, Short_t Lep4Id,
                                      const math::XYZTLorentzVector Lep1, const math::XYZTLorentzVector Lep2,
                                      const math::XYZTLorentzVector Lep3, const math::XYZTLorentzVector Lep4)
{
  GenLep1Pt=Lep1.Pt();
  GenLep1Eta=Lep1.Eta();
  GenLep1Phi=Lep1.Phi();
  GenLep1Id=Lep1Id;

  GenLep2Pt=Lep2.Pt();
  GenLep2Eta=Lep2.Eta();
  GenLep2Phi=Lep2.Phi();
  GenLep2Id=Lep2Id;

  GenLep3Pt=Lep3.Pt();
  GenLep3Eta=Lep3.Eta();
  GenLep3Phi=Lep3.Phi();
  GenLep3Id=Lep3Id;

  GenLep4Pt=Lep4.Pt();
  GenLep4Eta=Lep4.Eta();
  GenLep4Phi=Lep4.Phi();
  GenLep4Id=Lep4Id;

  //can comment this back in if Gen angles are needed for any reason...
  //TUtil::computeAngles(zzanalysis::tlv(Lep1), Lep1Id, zzanalysis::tlv(Lep2), Lep2Id, zzanalysis::tlv(Lep3), Lep3Id, zzanalysis::tlv(Lep4), Lep4Id, Gencosthetastar, GenhelcosthetaZ1, GenhelcosthetaZ2, Genhelphi, GenphistarZ1);

  return;
}

void HZZ4lNtupleMaker::FillAssocLepGenInfo(std::vector<const reco::Candidate *>& AssocLeps)
{
  if (AssocLeps.size() >= 1) {
    GenAssocLep1Pt =AssocLeps.at(0)->p4().Pt();
    GenAssocLep1Eta=AssocLeps.at(0)->p4().Eta();
    GenAssocLep1Phi=AssocLeps.at(0)->p4().Phi();
    GenAssocLep1Id =AssocLeps.at(0)->pdgId();
  }
  if (AssocLeps.size() >= 2) {
    GenAssocLep2Pt =AssocLeps.at(1)->p4().Pt();
    GenAssocLep2Eta=AssocLeps.at(1)->p4().Eta();
    GenAssocLep2Phi=AssocLeps.at(1)->p4().Phi();
    GenAssocLep2Id =AssocLeps.at(1)->pdgId();
  }

  return;
}


void HZZ4lNtupleMaker::FillHGenInfo(const math::XYZTLorentzVector pH, float w)
{
  GenHMass=pH.M();
  GenHPt=pH.Pt();
  GenHRapidity=pH.Rapidity();

  HqTMCweight=w;

  return;
}


void HZZ4lNtupleMaker::BookAllBranches(){
   //Event variables
  myTree->Book("RunNumber",RunNumber, failedTreeLevel >= minimalFailedTree);
  myTree->Book("EventNumber",EventNumber, failedTreeLevel >= minimalFailedTree);
  myTree->Book("LumiNumber",LumiNumber, failedTreeLevel >= minimalFailedTree);
  myTree->Book("NRecoMu",NRecoMu, failedTreeLevel >= fullFailedTree);
  myTree->Book("NRecoEle",NRecoEle, failedTreeLevel >= fullFailedTree);
  myTree->Book("Nvtx",Nvtx, failedTreeLevel >= fullFailedTree);
  myTree->Book("NObsInt",NObsInt, failedTreeLevel >= fullFailedTree);
  myTree->Book("NTrueInt",NTrueInt, failedTreeLevel >= fullFailedTree);

  myTree->Book("GenMET", GenMET, failedTreeLevel >= minimalFailedTree);
  myTree->Book("GenMETPhi", GenMETPhi, failedTreeLevel >= minimalFailedTree);
  myTree->Book("PFMET", metobj.extras.met, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMET_jesUp", metobj.extras.met_JECup, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMET_jesDn", metobj.extras.met_JECdn, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi", metobj.extras.phi, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi_jesUp", metobj.extras.phi_JECup, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi_jesDn", metobj.extras.phi_JECdn, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMET_corrected", metobj_corrected.extras.met, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMET_corrected_jesUp", metobj_corrected.extras.met_JECup, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMET_corrected_jesDn", metobj_corrected.extras.met_JECdn, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMET_corrected_jerUp", metobj_corrected.extras.met_JERup, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMET_corrected_jerDn", metobj_corrected.extras.met_JERdn, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMET_corrected_puUp", metobj_corrected.extras.met_PUup, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMET_corrected_puDn", metobj_corrected.extras.met_PUdn, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMET_corrected_metUp", metobj_corrected.extras.met_METup, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMET_corrected_metDn", metobj_corrected.extras.met_METdn, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi_corrected", metobj_corrected.extras.phi, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi_corrected_jesUp", metobj_corrected.extras.phi_JECup, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi_corrected_jesDn", metobj_corrected.extras.phi_JECdn, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi_corrected_jerUp", metobj_corrected.extras.phi_JERup, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi_corrected_jerDn", metobj_corrected.extras.phi_JERdn, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi_corrected_puUp", metobj_corrected.extras.phi_PUup, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi_corrected_puDn", metobj_corrected.extras.phi_PUdn, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi_corrected_metUp", metobj_corrected.extras.phi_METup, failedTreeLevel >= fullFailedTree);
  myTree->Book("PFMETPhi_corrected_metDn", metobj_corrected.extras.phi_METdn, failedTreeLevel >= fullFailedTree);
  //myTree->Book("PFMETNoHF",PFMETNoHF, failedTreeLevel >= fullFailedTree);
  //myTree->Book("PFMETNoHFPhi",PFMETNoHFPhi, failedTreeLevel >= fullFailedTree);

  myTree->Book("nCleanedJets",nCleanedJets, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30",nCleanedJetsPt30, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesUp",nCleanedJetsPt30_jesUp, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jesDn",nCleanedJetsPt30_jesDn, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jerUp",nCleanedJetsPt30_jerUp, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30_jerDn",nCleanedJetsPt30_jerDn, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged",nCleanedJetsPt30BTagged, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF",nCleanedJetsPt30BTagged_bTagSF, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesUp",nCleanedJetsPt30BTagged_bTagSF_jesUp, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jesDn",nCleanedJetsPt30BTagged_bTagSF_jesDn, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jerUp",nCleanedJetsPt30BTagged_bTagSF_jerUp, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSF_jerDn",nCleanedJetsPt30BTagged_bTagSF_jerDn, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSFUp",nCleanedJetsPt30BTagged_bTagSFUp, failedTreeLevel >= fullFailedTree);
  myTree->Book("nCleanedJetsPt30BTagged_bTagSFDn",nCleanedJetsPt30BTagged_bTagSFDn, failedTreeLevel >= fullFailedTree);
  myTree->Book("trigWord",trigWord, failedTreeLevel >= minimalFailedTree);
  myTree->Book("evtPassMETFilter",evtPassMETTrigger, failedTreeLevel >= minimalFailedTree);
  myTree->Book("ZZMass",ZZMass, false);
  myTree->Book("ZZMassErr",ZZMassErr, false);
  myTree->Book("ZZMassErrCorr",ZZMassErrCorr, false);
  myTree->Book("ZZMassPreFSR",ZZMassPreFSR, false);
  myTree->Book("ZZsel",ZZsel, false);
  myTree->Book("ZZPt",ZZPt, false);
  myTree->Book("ZZEta",ZZEta, false);
  myTree->Book("ZZPhi",ZZPhi, false);
  myTree->Book("ZZjjPt",ZZjjPt, false);
  myTree->Book("CRflag",CRflag, false);
  myTree->Book("Z1Mass",Z1Mass, false);
  myTree->Book("Z1Pt",Z1Pt, false);
  myTree->Book("Z1Flav",Z1Flav, false);

  //Kin refitted info
  if (addKinRefit) {
    myTree->Book("ZZMassRefit",ZZMassRefit, false);
    myTree->Book("ZZMassRefitErr",ZZMassRefitErr, false);
    myTree->Book("ZZMassUnrefitErr",ZZMassUnrefitErr, false);
  }
  if (addVtxFit){
    myTree->Book("ZZMassCFit",ZZMassCFit, false);
    myTree->Book("ZZChi2CFit",ZZChi2CFit, false);
  }

  //Z2 variables
  myTree->Book("Z2Mass",Z2Mass, false);
  myTree->Book("Z2Pt",Z2Pt, false);
  myTree->Book("Z2Flav",Z2Flav, false);
  myTree->Book("costhetastar",costhetastar, false);
  myTree->Book("helphi",helphi, false);
  myTree->Book("helcosthetaZ1",helcosthetaZ1, false);
  myTree->Book("helcosthetaZ2",helcosthetaZ2, false);
  myTree->Book("phistarZ1",phistarZ1, false);
  myTree->Book("phistarZ2",phistarZ2, false);
  myTree->Book("xi",xi, false);
  myTree->Book("xistar",xistar, false);

  if (is_loose_ele_selection) {
    myTree->Book("TLE_dR_Z",TLE_dR_Z, false);
    myTree->Book("TLE_min_dR_3l",TLE_min_dR_3l, false);
  }

  myTree->Book("Jettiness",Jettiness, false);//clara
  myTree->Book("LepPx",LepPx, false);//clara
  myTree->Book("LepPy",LepPy, false);//clara
  myTree->Book("LepPz",LepPz, false);//clara
  myTree->Book("LepE",LepE, false);//clara
  myTree->Book("LepPt",LepPt, false);
  myTree->Book("LepEta",LepEta, false);
  myTree->Book("LepPhi",LepPhi, false);
  myTree->Book("LepSCEta",LepSCEta, false);
  myTree->Book("LepLepId",LepLepId, false);
  myTree->Book("LepSIP",LepSIP, false);
  myTree->Book("Lepdxy",Lepdxy, false);
  myTree->Book("Lepdz",Lepdz, false);
  myTree->Book("LepTime",LepTime, false);
  myTree->Book("LepisID",LepisID, false);
  myTree->Book("LepisLoose",LepisLoose, false);
  myTree->Book("LepBDT",LepBDT, false);
  myTree->Book("LepisCrack",LepisCrack, false);
  myTree->Book("LepMissingHit",LepMissingHit, false);
  myTree->Book("LepChargedHadIso",LepChargedHadIso, false);
  myTree->Book("LepNeutralHadIso",LepNeutralHadIso, false);
  myTree->Book("LepPhotonIso",LepPhotonIso, false);
  myTree->Book("LepPUIsoComponent",LepPUIsoComponent, false);
  myTree->Book("LepCombRelIsoPF",LepCombRelIsoPF, false);
  myTree->Book("LepSF",LepSF, false);
  myTree->Book("LepSF_Unc",LepSF_Unc, false);
  myTree->Book("LepScale_Total_Up",LepScale_Total_Up, false);
  myTree->Book("LepScale_Total_Dn",LepScale_Total_Dn, false);
  myTree->Book("LepScale_Stat_Up",LepScale_Stat_Up, false);
  myTree->Book("LepScale_Stat_Dn",LepScale_Stat_Dn, false);
  myTree->Book("LepScale_Syst_Up",LepScale_Syst_Up, false);
  myTree->Book("LepScale_Syst_Dn",LepScale_Syst_Dn, false);
  myTree->Book("LepScale_Gain_Up",LepScale_Gain_Up, false);
  myTree->Book("LepScale_Gain_Dn",LepScale_Gain_Dn, false);
  myTree->Book("LepSigma_Total_Up",LepSigma_Total_Up, false);
  myTree->Book("LepSigma_Total_Dn",LepSigma_Total_Dn, false);
  myTree->Book("LepSigma_Rho_Up",LepSigma_Rho_Up, false);
  myTree->Book("LepSigma_Rho_Dn",LepSigma_Rho_Dn, false);
  myTree->Book("LepSigma_Phi_Up",LepSigma_Phi_Up, false);
  myTree->Book("LepSigma_Phi_Dn",LepSigma_Phi_Up, false);

  myTree->Book("fsrPt",fsrPt, false);
  myTree->Book("fsrEta",fsrEta, false);
  myTree->Book("fsrPhi",fsrPhi, false);
  myTree->Book("fsrLept",fsrLept, false);
  myTree->Book("passIsoPreFSR",passIsoPreFSR, false);
  if (addFSRDetails) {
    myTree->Book("fsrDR",fsrDR, false);
    myTree->Book("fsrLeptId",fsrLeptID, false);
    myTree->Book("fsrGenPt",fsrGenPt, false);
  }


  myTree->Book("JetPt",JetPt, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetEta",JetEta, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPhi",JetPhi, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetMass",JetMass, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetBTagger",JetBTagger, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetIsBtagged",JetIsBtagged, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetIsBtaggedWithSF",JetIsBtaggedWithSF, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetIsBtaggedWithSFUp",JetIsBtaggedWithSFUp, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetIsBtaggedWithSFDn",JetIsBtaggedWithSFDn, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetQGLikelihood",JetQGLikelihood, failedTreeLevel >= fullFailedTree);
  if(addQGLInputs){
    myTree->Book("JetAxis2",JetAxis2, failedTreeLevel >= fullFailedTree);
    myTree->Book("JetMult",JetMult, failedTreeLevel >= fullFailedTree);
    myTree->Book("JetPtD",JetPtD, failedTreeLevel >= fullFailedTree);
  }
  myTree->Book("JetSigma",JetSigma, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetHadronFlavour",JetHadronFlavour, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPartonFlavour",JetPartonFlavour, failedTreeLevel >= fullFailedTree);

  myTree->Book("JetRawPt",JetRawPt, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPtJEC_noJER",JetPtJEC_noJER, failedTreeLevel >= fullFailedTree);
  
  myTree->Book("JetPt_JESUp",JetJESUp, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JESDown",JetJESDown, failedTreeLevel >= fullFailedTree);
   
  myTree->Book("JetPt_JERUp",JetJERUp, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPt_JERDown",JetJERDown, failedTreeLevel >= fullFailedTree);

  myTree->Book("JetID", JetID, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPUID", JetPUID, failedTreeLevel >= fullFailedTree);
  myTree->Book("JetPUValue", JetPUValue, failedTreeLevel >= fullFailedTree);

  myTree->Book("DiJetMass",DiJetMass, false);
//   myTree->Book("DiJetMassPlus",DiJetMassPlus, false); // FIXME: add back once filled again
//   myTree->Book("DiJetMassMinus",DiJetMassMinus, false);
  myTree->Book("DiJetDEta",DiJetDEta, false);
  myTree->Book("DiJetFisher",DiJetFisher, false);
  
  //Photon variables
  myTree->Book("PhotonPt",PhotonPt, failedTreeLevel >= fullFailedTree);
  myTree->Book("PhotonEta",PhotonEta, failedTreeLevel >= fullFailedTree);
  myTree->Book("PhotonPhi",PhotonPhi, failedTreeLevel >= fullFailedTree);
  myTree->Book("PhotonIsCutBasedLooseID",PhotonIsCutBasedLooseID, failedTreeLevel >= fullFailedTree);
   
  myTree->Book("nExtraLep",nExtraLep, false);
  myTree->Book("nExtraZ",nExtraZ, false);
  myTree->Book("ExtraLepPt",ExtraLepPt, false);
  myTree->Book("ExtraLepEta",ExtraLepEta, false);
  myTree->Book("ExtraLepPhi",ExtraLepPhi, false);
  myTree->Book("ExtraLepLepId",ExtraLepLepId, false);
  myTree->Book("ExtraLepPx",ExtraLepPx, false);//clara
  myTree->Book("ExtraLepPy",ExtraLepPy, false);//clara
  myTree->Book("ExtraLepPz",ExtraLepPz, false);//clara
  myTree->Book("ExtraLepE",ExtraLepE, false);//clara
  myTree->Book("ZXFakeweight", ZXFakeweight, false);

  if (isMC){
    if (apply_K_NNLOQCD_ZZGG>0){
      myTree->Book("KFactor_QCD_ggZZ_Nominal", KFactor_QCD_ggZZ_Nominal, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_PDFScaleDn", KFactor_QCD_ggZZ_PDFScaleDn, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_PDFScaleUp", KFactor_QCD_ggZZ_PDFScaleUp, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_QCDScaleDn", KFactor_QCD_ggZZ_QCDScaleDn, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_QCDScaleUp", KFactor_QCD_ggZZ_QCDScaleUp, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_AsDn", KFactor_QCD_ggZZ_AsDn, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_AsUp", KFactor_QCD_ggZZ_AsUp, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_PDFReplicaDn", KFactor_QCD_ggZZ_PDFReplicaDn, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_ggZZ_PDFReplicaUp", KFactor_QCD_ggZZ_PDFReplicaUp, failedTreeLevel >= minimalFailedTree);
    }
    if (apply_K_NLOEW_ZZQQB){
      myTree->Book("KFactor_EW_qqZZ", KFactor_EW_qqZZ, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_EW_qqZZ_unc", KFactor_EW_qqZZ_unc, failedTreeLevel >= minimalFailedTree);
    }
    if (apply_K_NNLOQCD_ZZQQB){
      myTree->Book("KFactor_QCD_qqZZ_dPhi", KFactor_QCD_qqZZ_dPhi, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_qqZZ_M", KFactor_QCD_qqZZ_M, failedTreeLevel >= minimalFailedTree);
      myTree->Book("KFactor_QCD_qqZZ_Pt", KFactor_QCD_qqZZ_Pt, failedTreeLevel >= minimalFailedTree);
    }

    myTree->Book("genFinalState", genFinalState, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genProcessId", genProcessId, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genHEPMCweight", genHEPMCweight, failedTreeLevel >= minimalFailedTree);
    if (year == 2017 || year == 2018) myTree->Book("genHEPMCweight_NNLO", genHEPMCweight_NNLO, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genHEPMCweight_POWHEGonly", genHEPMCweight_POWHEGonly, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PUWeight", PUWeight, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PUWeight_Dn", PUWeight_Dn, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PUWeight_Up", PUWeight_Up, failedTreeLevel >= minimalFailedTree);
    myTree->Book("dataMCWeight", dataMCWeight, false);
    myTree->Book("trigEffWeight", trigEffWeight, false);
    myTree->Book("overallEventWeight", overallEventWeight, false);
    myTree->Book("L1prefiringWeight", L1prefiringWeight, false);
    myTree->Book("L1prefiringWeightUp", L1prefiringWeightUp, false);
    myTree->Book("L1prefiringWeightDn", L1prefiringWeightDn, false);
    myTree->Book("HqTMCweight", HqTMCweight, failedTreeLevel >= minimalFailedTree);
    myTree->Book("xsec", xsection, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genxsec", genxsection, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genBR", genbranchingratio, failedTreeLevel >= minimalFailedTree);
    myTree->Book("genExtInfo", genExtInfo, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenHMass", GenHMass, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenHPt", GenHPt, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenHRapidity", GenHRapidity, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenZ1Mass", GenZ1Mass, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ1Pt", GenZ1Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ1Phi", GenZ1Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ1Flav", GenZ1Flav, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenZ2Mass", GenZ2Mass, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ2Pt", GenZ2Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ2Phi", GenZ2Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenZ2Flav", GenZ2Flav, failedTreeLevel >= minimalFailedTree);
    myTree->Book("GenLep1Pt", GenLep1Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep1Eta", GenLep1Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep1Phi", GenLep1Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep1Id", GenLep1Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep2Pt", GenLep2Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep2Eta", GenLep2Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep2Phi", GenLep2Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep2Id", GenLep2Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep3Pt", GenLep3Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep3Eta", GenLep3Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep3Phi", GenLep3Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep3Id", GenLep3Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep4Pt", GenLep4Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep4Eta", GenLep4Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep4Phi", GenLep4Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenLep4Id", GenLep4Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep1Pt", GenAssocLep1Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep1Eta", GenAssocLep1Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep1Phi", GenAssocLep1Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep1Id", GenAssocLep1Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep2Pt", GenAssocLep2Pt, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep2Eta", GenAssocLep2Eta, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep2Phi", GenAssocLep2Phi, failedTreeLevel >= fullFailedTree);
    myTree->Book("GenAssocLep2Id", GenAssocLep2Id, failedTreeLevel >= fullFailedTree);
    myTree->Book("htxs_errorCode", htxs_errorCode, failedTreeLevel >= minimalFailedTree);
    myTree->Book("htxs_prodMode", htxs_prodMode, failedTreeLevel >= minimalFailedTree);
    myTree->Book("htxsNJets", htxsNJets, failedTreeLevel >= minimalFailedTree);
    myTree->Book("htxsHPt", htxsHPt, failedTreeLevel >= minimalFailedTree);
    myTree->Book("htxs_stage0_cat", htxs_stage0_cat, failedTreeLevel >= minimalFailedTree);
    myTree->Book("htxs_stage1p1_cat", htxs_stage1p1_cat, failedTreeLevel >= minimalFailedTree);
    myTree->Book("htxs_stage1p2_cat", htxs_stage1p2_cat, failedTreeLevel >= minimalFailedTree);
    if(apply_QCD_GGF_UNCERT)
      {
	myTree->Book("ggH_NNLOPS_weight", ggH_NNLOPS_weight, failedTreeLevel >= minimalFailedTree);
	myTree->Book("ggH_NNLOPS_weight_unc", ggH_NNLOPS_weight_unc, failedTreeLevel >= minimalFailedTree);
	myTree->Book("qcd_ggF_uncertSF", qcd_ggF_uncertSF, failedTreeLevel >= minimalFailedTree);
      }
	  

    if (addLHEKinematics){
      myTree->Book("LHEMotherPz", LHEMotherPz, failedTreeLevel >= LHEFailedTree);
      myTree->Book("LHEMotherE", LHEMotherE, failedTreeLevel >= LHEFailedTree);
      myTree->Book("LHEMotherId", LHEMotherId, failedTreeLevel >= LHEFailedTree);
      myTree->Book("LHEDaughterPt", LHEDaughterPt, failedTreeLevel >= LHEFailedTree);
      myTree->Book("LHEDaughterEta", LHEDaughterEta, failedTreeLevel >= LHEFailedTree);
      myTree->Book("LHEDaughterPhi", LHEDaughterPhi, failedTreeLevel >= LHEFailedTree);
      myTree->Book("LHEDaughterMass", LHEDaughterMass, failedTreeLevel >= LHEFailedTree);
      myTree->Book("LHEDaughterId", LHEDaughterId, failedTreeLevel >= LHEFailedTree);
      myTree->Book("LHEAssociatedParticlePt", LHEAssociatedParticlePt, failedTreeLevel >= LHEFailedTree);
      myTree->Book("LHEAssociatedParticleEta", LHEAssociatedParticleEta, failedTreeLevel >= LHEFailedTree);
      myTree->Book("LHEAssociatedParticlePhi", LHEAssociatedParticlePhi, failedTreeLevel >= LHEFailedTree);
      myTree->Book("LHEAssociatedParticleMass", LHEAssociatedParticleMass, failedTreeLevel >= LHEFailedTree);
      myTree->Book("LHEAssociatedParticleId", LHEAssociatedParticleId, failedTreeLevel >= LHEFailedTree);
    }

    myTree->Book("LHEPDFScale", LHEPDFScale, failedTreeLevel >= minimalFailedTree);
    myTree->Book("LHEweight_QCDscale_muR1_muF1", LHEweight_QCDscale_muR1_muF1, failedTreeLevel >= minimalFailedTree);
    myTree->Book("LHEweight_QCDscale_muR1_muF2", LHEweight_QCDscale_muR1_muF2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("LHEweight_QCDscale_muR1_muF0p5", LHEweight_QCDscale_muR1_muF0p5, failedTreeLevel >= minimalFailedTree);
    myTree->Book("LHEweight_QCDscale_muR2_muF1", LHEweight_QCDscale_muR2_muF1, failedTreeLevel >= minimalFailedTree);
    myTree->Book("LHEweight_QCDscale_muR2_muF2", LHEweight_QCDscale_muR2_muF2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("LHEweight_QCDscale_muR2_muF0p5", LHEweight_QCDscale_muR2_muF0p5, failedTreeLevel >= minimalFailedTree);
    myTree->Book("LHEweight_QCDscale_muR0p5_muF1", LHEweight_QCDscale_muR0p5_muF1, failedTreeLevel >= minimalFailedTree);
    myTree->Book("LHEweight_QCDscale_muR0p5_muF2", LHEweight_QCDscale_muR0p5_muF2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("LHEweight_QCDscale_muR0p5_muF0p5", LHEweight_QCDscale_muR0p5_muF0p5, failedTreeLevel >= minimalFailedTree);
    myTree->Book("LHEweight_PDFVariation_Up", LHEweight_PDFVariation_Up, failedTreeLevel >= minimalFailedTree);
    myTree->Book("LHEweight_PDFVariation_Dn", LHEweight_PDFVariation_Dn, failedTreeLevel >= minimalFailedTree);
    myTree->Book("LHEweight_AsMZ_Up", LHEweight_AsMZ_Up, failedTreeLevel >= minimalFailedTree);
    myTree->Book("LHEweight_AsMZ_Dn", LHEweight_AsMZ_Dn, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_isr_muR4", PythiaWeight_isr_muR4, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_isr_muR2", PythiaWeight_isr_muR2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_isr_muRsqrt2", PythiaWeight_isr_muRsqrt2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_isr_muRoneoversqrt2", PythiaWeight_isr_muRoneoversqrt2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_isr_muR0p5", PythiaWeight_isr_muR0p5, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_isr_muR0p25", PythiaWeight_isr_muR0p25, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_fsr_muR4", PythiaWeight_fsr_muR4, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_fsr_muR2", PythiaWeight_fsr_muR2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_fsr_muRsqrt2", PythiaWeight_fsr_muRsqrt2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_fsr_muRoneoversqrt2", PythiaWeight_fsr_muRoneoversqrt2, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_fsr_muR0p5", PythiaWeight_fsr_muR0p5, failedTreeLevel >= minimalFailedTree);
    myTree->Book("PythiaWeight_fsr_muR0p25", PythiaWeight_fsr_muR0p25, failedTreeLevel >= minimalFailedTree);
  }

  // MELA branches are booked under buildMELA
}

void HZZ4lNtupleMaker::addweight(float &weight, float weighttoadd) {
  weight += weighttoadd;
}


void HZZ4lNtupleMaker::buildMELABranches(){
  /***********************/
  /***********************/
  /**   Reco branches   **/
  /***********************/
  /***********************/
  for (unsigned int it=0; it<recoMElist.size(); it++){
    MELAOptionParser* me_opt = new MELAOptionParser(recoMElist.at(it));
    if (recoMElist.at(it).find("Copy")!=string::npos) recome_copyopts.push_back(me_opt);
    else recome_originalopts.push_back(me_opt);
  }
  // Resolve original options
  for (unsigned int it=0; it<recome_originalopts.size(); it++){
    MELAOptionParser* me_opt = recome_originalopts.at(it);
    myTree->BookMELABranches(me_opt, false, 0);
  }
  // Resolve copy options
  for (unsigned int it=0; it<recome_copyopts.size(); it++){
    MELAOptionParser* me_opt = recome_copyopts.at(it);
    MELAOptionParser* original_opt=0;
    // Find the original options
    for (unsigned int ih=0; ih<recome_originalopts.size(); ih++){
      if (me_opt->testCopyAlias(recome_originalopts.at(ih)->getAlias())){
        original_opt = recome_originalopts.at(ih);
        break;
      }
    }
    if (original_opt==0) continue;
    else me_opt->pickOriginalOptions(original_opt);
    myTree->BookMELABranches(me_opt, false, 0);
  }

  /**********************/
  /**********************/
  /**   LHE branches   **/
  /**********************/
  /**********************/
  for (unsigned int it=0; it<lheMElist.size(); it++){
    MELAOptionParser* lheme_opt;
    // First find out if the option has a copy specification
    // These copy options will be evaulated in a separate loop
    if (lheMElist.at(it).find("Copy")!=string::npos){
      lheme_opt = new MELAOptionParser(lheMElist.at(it));
      lheme_copyopts.push_back(lheme_opt);
      continue;
    }

    // Create a hypothesis for each option
    MELAHypothesis* lheme_hypo = new MELAHypothesis(&mela, lheMElist.at(it));
    lheme_units.push_back(lheme_hypo);

    lheme_opt = lheme_hypo->getOption();
    if (lheme_opt->isAliased()) lheme_aliased_units.push_back(lheme_hypo);

    // Create a computation for each hypothesis
    MELAComputation* lheme_computer = new MELAComputation(lheme_hypo);
    lheme_computers.push_back(lheme_computer);

    // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
    GMECHelperFunctions::addToMELACluster(lheme_computer, lheme_clusters);

    // Create the necessary branches for each computation
    myTree->BookMELABranches(lheme_opt, true, lheme_computer);
  }
  // Resolve copy options
  for (unsigned int it=0; it<lheme_copyopts.size(); it++){
    MELAOptionParser* lheme_opt = lheme_copyopts.at(it);
    MELAHypothesis* original_hypo=0;
    MELAOptionParser* original_opt=0;
    // Find the original options
    for (unsigned int ih=0; ih<lheme_aliased_units.size(); ih++){
      if (lheme_opt->testCopyAlias(lheme_aliased_units.at(ih)->getOption()->getAlias())){
        original_hypo = lheme_aliased_units.at(ih);
        original_opt = original_hypo->getOption();
        break;
      }
    }
    if (original_opt==0) continue;
    else lheme_opt->pickOriginalOptions(original_opt);
    // Create a new computation for the copy options
    MELAComputation* lheme_computer = new MELAComputation(original_hypo);
    lheme_computer->setOption(lheme_opt);
    lheme_computers.push_back(lheme_computer);

    // The rest is the same story...
    // Add the computation to a named cluster to keep track of JECUp/JECDn, or for best-pWH_SM Lep_WH computations
    GMECHelperFunctions::addToMELACluster(lheme_computer, lheme_clusters);

    // Create the necessary branches for each computation
    myTree->BookMELABranches(lheme_opt, true, lheme_computer);
  }
  // Loop over the computations to add any contingencies to aliased hypotheses
  for (unsigned int it=0; it<lheme_computers.size(); it++) lheme_computers.at(it)->addContingencies(lheme_aliased_units);

  if (DEBUG_MB){
    std::vector<MELABranch*>* lheme_branches = myTree->getLHEMELABranches();
    for (unsigned int ib=0; ib<lheme_branches->size(); ib++) lheme_branches->at(ib)->Print();
    for (unsigned int icl=0; icl<lheme_clusters.size(); icl++) cout << "LHE ME cluster " << lheme_clusters.at(icl)->getName() << " is present in " << lheme_clusters.size() << " clusters with #Computations = " << lheme_clusters.at(icl)->getComputations()->size() << endl;
  }
}

void HZZ4lNtupleMaker::computeMELABranches(MELACandidate* cand){
  mela.setCurrentCandidate(cand);
  // Sequantial computation
  updateMELAClusters_Common("Common");
  updateMELAClusters_NoInitialQ("NoInitialQ");
  updateMELAClusters_NoInitialG("NoInitialG");
  updateMELAClusters_BestLOAssociatedZ("BestLOAssociatedZ");
  updateMELAClusters_BestLOAssociatedW("BestLOAssociatedW");
  updateMELAClusters_BestLOAssociatedVBF("BestLOAssociatedVBF");
  updateMELAClusters_BestNLOVHApproximation("BestNLOZHApproximation");
  updateMELAClusters_BestNLOVHApproximation("BestNLOWHApproximation");
  updateMELAClusters_BestNLOVBFApproximation("BestNLOVBFApproximation");
  updateMELAClusters_NoAssociatedG("NoAssociatedG");
  updateMELAClusters_NoInitialGNoAssociatedG("NoInitialGNoAssociatedG");
  // Reverse sequence
  updateMELAClusters_NoInitialGNoAssociatedG("NoInitialGNoAssociatedGLast");
  updateMELAClusters_NoAssociatedG("NoAssociatedGLast");
  updateMELAClusters_BestNLOVBFApproximation("BestNLOVBFApproximationLast");
  updateMELAClusters_BestNLOVHApproximation("BestNLOWHApproximationLast");
  updateMELAClusters_BestNLOVHApproximation("BestNLOZHApproximationLast");
  updateMELAClusters_BestLOAssociatedVBF("BestLOAssociatedVBFLast");
  updateMELAClusters_BestLOAssociatedW("BestLOAssociatedWLast");
  updateMELAClusters_BestLOAssociatedZ("BestLOAssociatedZLast");
  updateMELAClusters_NoInitialG("NoInitialGLast");
  updateMELAClusters_NoInitialQ("NoInitialQLast");
  updateMELAClusters_Common("CommonLast");
  // Reset mela
  mela.resetInputEvent();
}
// Common ME computations that do not manipulate the LHE candidate
void HZZ4lNtupleMaker::updateMELAClusters_Common(const string clustertype){
  MELACandidate* melaCand = mela.getCurrentCandidate();
  if (melaCand==0) return;

  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
    MELACluster* theCluster = lheme_clusters.at(ic);
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }
}
// ME computations that require no quark initial state
void HZZ4lNtupleMaker::updateMELAClusters_NoInitialQ(const string clustertype){
  MELACandidate* melaCand = mela.getCurrentCandidate();
  if (melaCand==0) return;

  // Manipulate the candidate
  // Assign 0 to the id of quark mothers
  std::vector<int> motherIds;
  for (int imot=0; imot<melaCand->getNMothers(); imot++){
    motherIds.push_back(melaCand->getMother(imot)->id);
    if (PDGHelpers::isAQuark(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
  }

  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
    MELACluster* theCluster = lheme_clusters.at(ic);
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
}
// ME computations that require no gluon initial state
void HZZ4lNtupleMaker::updateMELAClusters_NoInitialG(const string clustertype){
  MELACandidate* melaCand = mela.getCurrentCandidate();
  if (melaCand==0) return;

  // Manipulate the candidate
  // Assign 0 to the id of gluon mothers
  std::vector<int> motherIds;
  for (int imot=0; imot<melaCand->getNMothers(); imot++){
    motherIds.push_back(melaCand->getMother(imot)->id);
    if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
  }

  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
    MELACluster* theCluster = lheme_clusters.at(ic);
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
}
// ME computations that require no gluons as associated particles
void HZZ4lNtupleMaker::updateMELAClusters_NoAssociatedG(const string clustertype){
  MELACandidate* melaCand = mela.getCurrentCandidate();
  if (melaCand==0) return;

  // Manipulate the candidate
  // Assign 0 to the id of gluon mothers
  std::vector<int> ajetIds;
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
    ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
    if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
  }

  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
    MELACluster* theCluster = lheme_clusters.at(ic);
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
}
// ME computations that require no gluon initial state and no gluons as associated particles
void HZZ4lNtupleMaker::updateMELAClusters_NoInitialGNoAssociatedG(const string clustertype){
  MELACandidate* melaCand = mela.getCurrentCandidate();
  if (melaCand==0) return;

  // Manipulate the candidate
  // Assign 0 to the id of gluon mothers
  std::vector<int> motherIds;
  std::vector<int> ajetIds;
  for (int imot=0; imot<melaCand->getNMothers(); imot++){
    motherIds.push_back(melaCand->getMother(imot)->id);
    if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
  }
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
    ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
    if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
  }

  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
    MELACluster* theCluster = lheme_clusters.at(ic);
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
}
// ME computations that require best Z, W or VBF topology at LO (no gluons)
void HZZ4lNtupleMaker::updateMELAClusters_BestLOAssociatedZ(const string clustertype){
  MELACandidate* melaCand = mela.getCurrentCandidate();
  if (melaCand==0) return;

  // Manipulate the candidate
  // Assign 0 to the id of gluon mothers
  std::vector<int> motherIds;
  std::vector<int> ajetIds;
  for (int imot=0; imot<melaCand->getNMothers(); imot++){
    motherIds.push_back(melaCand->getMother(imot)->id);
    if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
  }
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
    ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
    if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
  }
  // Give precedence to leptonic V decays
  bool hasALepV=false;
  for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
    MELAParticle* Vtmp = melaCand->getSortedV(iv);
    if (Vtmp!=0 && PDGHelpers::isAZBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
      if (
        PDGHelpers::isALepton(Vtmp->getDaughter(0)->id)
        ||
        PDGHelpers::isANeutrino(Vtmp->getDaughter(0)->id)
        ){
        hasALepV=true;
      }
    }
  }
  int bestVbyMass=-1;
  float bestVMassDiff=1e5;
  for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
    MELAParticle* Vtmp = melaCand->getSortedV(iv);
    if (Vtmp!=0 && PDGHelpers::isAZBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
      if (
        PDGHelpers::isAJet(Vtmp->getDaughter(0)->id)
        && hasALepV
        ){
        for (int idau=0; idau<Vtmp->getNDaughters(); idau++) Vtmp->getDaughter(idau)->setSelected(false);
      }
      else if (fabs(Vtmp->m()-PDGHelpers::Zmass)<bestVMassDiff){
        bestVMassDiff=fabs(Vtmp->m()-PDGHelpers::Zmass);
        bestVbyMass = iv;
      }
    }
  }
  for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
    MELAParticle* Vtmp = melaCand->getSortedV(iv);
    if (Vtmp!=0 && PDGHelpers::isAZBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
      for (int idau=0; idau<Vtmp->getNDaughters(); idau++) Vtmp->getDaughter(idau)->setSelected((iv==bestVbyMass));
    }
  }

  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
    MELACluster* theCluster = lheme_clusters.at(ic);
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
  for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
    MELAParticle* Vtmp = melaCand->getSortedV(iv);
    if (Vtmp!=0 && PDGHelpers::isAZBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
      for (int idau=0; idau<Vtmp->getNDaughters(); idau++) Vtmp->getDaughter(idau)->setSelected(true);
    }
  }
}
void HZZ4lNtupleMaker::updateMELAClusters_BestLOAssociatedW(const string clustertype){
  MELACandidate* melaCand = mela.getCurrentCandidate();
  if (melaCand==0) return;

  // Manipulate the candidate
  // Assign 0 to the id of gluon mothers
  std::vector<int> motherIds;
  std::vector<int> ajetIds;
  for (int imot=0; imot<melaCand->getNMothers(); imot++){
    motherIds.push_back(melaCand->getMother(imot)->id);
    if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
  }
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
    ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
    if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
  }
  // Give precedence to leptonic V decays
  bool hasALepV=false;
  for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
    MELAParticle* Vtmp = melaCand->getSortedV(iv);
    if (Vtmp!=0 && PDGHelpers::isAWBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
      if (
        PDGHelpers::isALepton(Vtmp->getDaughter(0)->id)
        ||
        PDGHelpers::isANeutrino(Vtmp->getDaughter(0)->id)
        ){
        hasALepV=true;
      }
    }
  }
  int bestVbyMass=-1;
  float bestVMassDiff=1e5;
  for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
    MELAParticle* Vtmp = melaCand->getSortedV(iv);
    if (Vtmp!=0 && PDGHelpers::isAWBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
      if (
        PDGHelpers::isAJet(Vtmp->getDaughter(0)->id)
        && hasALepV
        ){
        for (int idau=0; idau<Vtmp->getNDaughters(); idau++) Vtmp->getDaughter(idau)->setSelected(false);
      }
      else if (fabs(Vtmp->m()-PDGHelpers::Wmass)<bestVMassDiff){
        bestVMassDiff=fabs(Vtmp->m()-PDGHelpers::Wmass);
        bestVbyMass = iv;
      }
    }
  }
  for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
    MELAParticle* Vtmp = melaCand->getSortedV(iv);
    if (Vtmp!=0 && PDGHelpers::isAWBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
      for (int idau=0; idau<Vtmp->getNDaughters(); idau++) Vtmp->getDaughter(idau)->setSelected((iv==bestVbyMass));
    }
  }

  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
    MELACluster* theCluster = lheme_clusters.at(ic);
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
  for (int iv=2; iv<melaCand->getNSortedVs(); iv++){
    MELAParticle* Vtmp = melaCand->getSortedV(iv);
    if (Vtmp!=0 && PDGHelpers::isAWBoson(Vtmp->id) && Vtmp->getNDaughters()>=1){
      for (int idau=0; idau<Vtmp->getNDaughters(); idau++) Vtmp->getDaughter(idau)->setSelected(true);
    }
  }
}
void HZZ4lNtupleMaker::updateMELAClusters_BestLOAssociatedVBF(const string clustertype){
  // Same as updateMELAClusters_NoInitialGNoAssociatedG, but keep a separate function for future studies
  MELACandidate* melaCand = mela.getCurrentCandidate();
  if (melaCand==0) return;

  // Manipulate the candidate
  // Assign 0 to the id of gluon mothers
  std::vector<int> motherIds;
  std::vector<int> ajetIds;
  for (int imot=0; imot<melaCand->getNMothers(); imot++){
    motherIds.push_back(melaCand->getMother(imot)->id);
    if (PDGHelpers::isAGluon(melaCand->getMother(imot)->id)) melaCand->getMother(imot)->id = 0;
  }
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++){
    ajetIds.push_back(melaCand->getAssociatedJet(ijet)->id);
    if (PDGHelpers::isAGluon(melaCand->getAssociatedJet(ijet)->id)) melaCand->getAssociatedJet(ijet)->id = 0;
  }

  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
    MELACluster* theCluster = lheme_clusters.at(ic);
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  // Restore the candidate properties
  for (int imot=0; imot<melaCand->getNMothers(); imot++) melaCand->getMother(imot)->id = motherIds.at(imot); // Restore all mother ids
  for (int ijet=0; ijet<melaCand->getNAssociatedJets(); ijet++) melaCand->getAssociatedJet(ijet)->id = ajetIds.at(ijet); // Restore all jets
}
// ME computations that can approximate the NLO QCD (-/+ MiNLO extra jet) phase space to LO QCD in signal VBF or VH
// Use these for POWHEG samples
// MELACandidateRecaster has very specific use cases, so do not use these functions for other cases.
void HZZ4lNtupleMaker::updateMELAClusters_BestNLOVHApproximation(const string clustertype){
  MELACandidate* melaCand = mela.getCurrentCandidate();
  if (melaCand==0) return;

  // Check if any clusters request this computation
  bool clustersRequest=false;
  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
    MELACluster* theCluster = lheme_clusters.at(ic);
    if (theCluster->getName()==clustertype){
      clustersRequest=true;
      break;
    }
  }
  if (!clustersRequest) return;

  // Need one recaster for each of ZH and WH, so distinguish by the cluster name
  TVar::Production candScheme;
  if (clustertype.find("BestNLOZHApproximation")!=string::npos) candScheme = TVar::Had_ZH;
  else if (clustertype.find("BestNLOWHApproximation")!=string::npos) candScheme = TVar::Had_WH;
  else return;

  MELACandidateRecaster recaster(candScheme);
  MELACandidate* candModified=nullptr;
  MELAParticle* bestAV = MELACandidateRecaster::getBestAssociatedV(melaCand, candScheme);
  if (bestAV){
    recaster.copyCandidate(melaCand, candModified);
    recaster.deduceLOVHTopology(candModified);
    mela.setCurrentCandidate(candModified);
  }
  else return; // No associated Vs found. The algorithm won't work.

  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
    MELACluster* theCluster = lheme_clusters.at(ic);
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  delete candModified;
  mela.setCurrentCandidate(melaCand); // Go back to the original candidate
}
void HZZ4lNtupleMaker::updateMELAClusters_BestNLOVBFApproximation(const string clustertype){
  MELACandidate* melaCand = mela.getCurrentCandidate();
  if (melaCand==0) return;

  // Check if any clusters request this computation
  bool clustersRequest=false;
  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
    MELACluster* theCluster = lheme_clusters.at(ic);
    if (theCluster->getName()==clustertype){
      clustersRequest=true;
      break;
    }
  }
  if (!clustersRequest) return;

  // Need one recaster for VBF
  TVar::Production candScheme;
  if (clustertype.find("BestNLOVBFApproximation")!=string::npos) candScheme = TVar::JJVBF;
  else return;

  MELACandidateRecaster recaster(candScheme);
  MELACandidate* candModified=nullptr;
  recaster.copyCandidate(melaCand, candModified);
  recaster.reduceJJtoQuarks(candModified);
  mela.setCurrentCandidate(candModified);

  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++){
    MELACluster* theCluster = lheme_clusters.at(ic);
    if (theCluster->getName()==clustertype){
      // Re-compute all related hypotheses first...
      theCluster->computeAll();
      // ...then update the cluster
      theCluster->update();
    }
  }

  delete candModified;
  mela.setCurrentCandidate(melaCand); // Go back to the original candidate
}


void HZZ4lNtupleMaker::pushRecoMELABranches(const pat::CompositeCandidate& cand){
  std::vector<MELABranch*>* recome_branches = myTree->getRecoMELABranches();
 //Pull+push...
  for (unsigned int ib=0; ib<recome_branches->size(); ib++){
    std::string branchname = recome_branches->at(ib)->bname.Data();
    if (cand.hasUserFloat(branchname)) recome_branches->at(ib)->setValue((Float_t)cand.userFloat(branchname));
    else cerr << "HZZ4lNtupleMaker::pushRecoMELABranches: Candidate does not contain the reco ME " << branchname << " it should have calculated!" << endl;
  }
}
void HZZ4lNtupleMaker::pushLHEMELABranches(){
  std::vector<MELABranch*>* lheme_branches = myTree->getLHEMELABranches();
  // Pull + push...
  for (unsigned int ib=0; ib<lheme_branches->size(); ib++) lheme_branches->at(ib)->setVal();
   // ...then reset
  for (unsigned int ic=0; ic<lheme_clusters.size(); ic++) lheme_clusters.at(ic)->reset();
}
void HZZ4lNtupleMaker::clearMELABranches(){
  for (unsigned int it=0; it<lheme_clusters.size(); it++) delete lheme_clusters.at(it);
  for (unsigned int it=0; it<lheme_computers.size(); it++) delete lheme_computers.at(it);
  for (unsigned int it=0; it<lheme_copyopts.size(); it++) delete lheme_copyopts.at(it);
  //for (unsigned int it=0; it<lheme_aliased_units.size(); it++) delete lheme_aliased_units.at(it); // DO NOT DELETE THIS, WILL BE DELETED WITH lheme_units!
  for (unsigned int it=0; it<lheme_units.size(); it++) delete lheme_units.at(it);
  
  for (unsigned int it=0; it<recome_copyopts.size(); it++) delete recome_copyopts.at(it);
  for (unsigned int it=0; it<recome_originalopts.size(); it++) delete recome_originalopts.at(it);
}

Int_t HZZ4lNtupleMaker::FindBinValue(TGraphErrors *tgraph, double value)
{
   Double_t x_prev,x,y;
   Int_t bin = 0;
   x_prev = 0.;
   for(int i=0;i<tgraph->GetN();i++){
      tgraph->GetPoint(i,x,y);
      if(value > x_prev && value < x){
         bin = i;
         break;
      }
      else x_prev = x;
   }
   if (bin == 0) bin = 1;
   return bin-1;
}

//define this as a plug-in
DEFINE_FWK_MODULE(HZZ4lNtupleMaker);
