#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <math.h>
#include <fstream>

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>

#include "StyleFunc.h"

#include "RooFit.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooCategory.h"
#include "RooPlot.h"

using namespace std;
using namespace RooFit;

class TreePFCandEventData {
public:
  // ===== Class Methods =====
  void Init();

  unsigned int runNb;
  unsigned int eventNb;
  unsigned int LS;

  // -- centrality variables --
  Int_t           CentBin;
  Int_t           Npix, NpixelTracks, Ntracks, NtracksPtCut, NtracksEtaCut, NtracksEtaPtCut;
  Float_t         SumET_HF, SumET_HFplus, SumET_HFminus, SumET_HFplusEta4, SumET_HFminusEta4;
  Float_t         SumET_HFhit, SumET_HFhitPlus, SumET_HFhitMinus;
  Float_t         SumET_ZDC, SumET_ZDCplus, SumET_ZDCminus;
  Float_t         SumET_EEplus, SumET_EEminus;
  Float_t         SumET_EE, SumET_EB, SumET_ET;

  // -- Primary Vetex --
  Float_t         nPV;
  Float_t         RefVtx_x, RefVtx_y, RefVtx_z;
  Float_t         RefVtx_xError, RefVtx_yError, RefVtx_zError;

  // -- particle info & GEN info --
  Int_t                        nPFpart, nGENpart, njets;
  std::vector<Int_t>           *pfId, *genPDGId;
  std::vector<Float_t>         *pfEnergy, *jetEnergy;
  std::vector<Float_t>         *pfPt, *genPt,  *jetPt;
  std::vector<Float_t>         *pfEta, *genEta,  *jetEta;
  std::vector<Float_t>         *pfPhi, *genPhi,  *jetPhi;
  std::vector<Float_t>         *jetMass, *jetPU;
  std::vector<Int_t>           *jetMatchIndex, *pfCharge;
  std::vector<Float_t>         *pfTheta, *pfEt;
  std::vector<Float_t>         *pfVsPt, *pfVsPtInitial, *pfArea;
  Float_t                      vn[200];
  Float_t                      psin[200];
  Float_t                      sumpt[20];
  Float_t                      ueraw[1200];
  // (particle flow charged hadrons and muons)
  std::vector<Float_t>         *pfMuonPx, *pfMuonPy, *pfMuonPz;
  std::vector<bool>            *pfTrackerMuon;
  std::vector<Float_t>         *pfTrackerMuonPt;
  std::vector<Int_t>           *pfTrackHits;
  std::vector<Float_t>         *pfDxy, *pfDz, *pfChi2;
  std::vector<Float_t>         *pfGlobalMuonPt;
  std::vector<Float_t>         *pfChargedPx, *pfChargedPy, *pfChargedPz;
  std::vector<Float_t>         *pfChargedTrackRefPt;

  // -- generalTracks info --
  Int_t                        nTRACKpart;
  std::vector<Int_t>           *traQual, *traCharge;
  std::vector<Float_t>         *traPt,   *traEta,  *traPhi;
  std::vector<Int_t>           *traAlgo, *traHits;
  // track algorithm enum:
  // https://cmssdt.cern.ch/SDT/doxygen/CMSSW_8_0_24/doc/html/da/d0c/TrackBase_8h_source.html#l00099

  // -- MET info --
  Float_t                      recoPFMET, recoPFMETPhi, recoPFMETsumEt;
  Float_t                      recoPFMETmEtSig, recoPFMETSig;

  // -- Muon info (pat::muons) --
  Int_t                        nMUpart;
  std::vector<Float_t>         *muPx, *muPy, *muPz;
  std::vector<Float_t>         *muMt, *muPt, *muEta, *muPhi;
  std::vector<Int_t>           *muCharge, *mu_SelectionType;
  std::vector<Float_t>         *muTrackIso, *muCaloIso, *muEcalIso, *muHcalIso; //R0.3 default pp isolation setup
  std::vector<Float_t>         *muSumChargedHadronPt, *muSumNeutralHadronEt, *muSumPhotonEt, *muSumPUPt, *muPFBasedDBetaIso; //R0.4 (*mu POG default PF-based isolation)
  std::vector<bool>            *muHighPurity, *muIsTightMuon, *muIsGoodMuon, *muTrkMuArb, *muTMOneStaTight;
  std::vector<Int_t>           *muSelectionType, *muNTrkHits, *muNPixValHits, *muNPixWMea, *muNTrkWMea, *muStationsMatched, *muNMuValHits;
  std::vector<Float_t>         *muDxy, *muDxyErr, *muDz, *muDzErr;
  std::vector<Float_t>         *muPtInner, *muPtErrInner, *muPtGlobal, *muPtErrGlobal;
  std::vector<Float_t>         *muNormChi2Inner, *muNormChi2Global;
  std::vector<Float_t>         *muIso03_sumPt, *muIso04_sumPt, *muIso05_sumPt, *muIso03_emEt, *muIso04_emEt, *muIso05_emEt, *muIso03_hadEt, *muIso04_hadEt, *muIso05_hadEt;
  std::vector<Int_t>           *muIso03_nTracks, *muIso04_nTracks, *muIso05_nTracks;
  std::vector<bool>            *muNotPFMuon;

  // -- Trigger info --
  std::vector<ULong64_t>       *muTrig;
  ULong64_t                    HLTriggers;
  std::vector<Int_t>           *trigPrescale;
};

void TreePFCandEventData::Init()
{
  pfId = 0;
  pfPt = 0;
  pfEnergy = 0;
  pfVsPt = 0;
  pfVsPtInitial = 0;
  pfArea = 0;
  pfEta = 0;
  pfPhi = 0;
  pfCharge = 0;
  pfTheta = 0;
  pfEt = 0;
  pfMuonPx = 0;
  pfMuonPy = 0;
  pfMuonPz = 0;
  pfTrackerMuon = 0;
  pfTrackerMuonPt = 0;
  pfTrackHits = 0;
  pfDxy = 0;
  pfDz = 0;
  pfChi2 = 0;
  pfGlobalMuonPt = 0;
  pfChargedPx = 0;
  pfChargedPy = 0;
  pfChargedPz = 0;
  pfChargedTrackRefPt = 0;
  genPDGId = 0;
  genPt = 0;
  genEta = 0;
  genPhi = 0;
  traQual = 0;
  traCharge = 0;
  traPt = 0;
  traEta = 0;
  traPhi = 0;
  traAlgo = 0;
  traHits = 0;
  muPx = 0;
  muPy = 0;
  muPz = 0;
  muMt = 0;
  muPt = 0;
  muEta = 0;
  muPhi = 0;
  muCharge = 0;
  muSelectionType = 0;
  muTrackIso = 0;
  muCaloIso = 0;
  muEcalIso = 0;
  muHcalIso = 0;
  muSumChargedHadronPt = 0;
  muSumNeutralHadronEt = 0;
  muSumPhotonEt = 0;
  muSumPUPt = 0;
  muPFBasedDBetaIso = 0;
  muHighPurity = 0;
  muIsTightMuon = 0;
  muIsGoodMuon = 0;
  muTrkMuArb = 0;
  muTMOneStaTight = 0;
  muNTrkHits = 0;
  muNPixValHits = 0;
  muNPixWMea = 0;
  muNTrkWMea = 0;
  muStationsMatched = 0;
  muNMuValHits = 0;
  muDxy = 0;
  muDxyErr = 0;
  muDz = 0;
  muDzErr = 0;
  muPtInner = 0;
  muPtErrInner = 0;
  muPtGlobal = 0;
  muPtErrGlobal = 0;
  muNormChi2Inner = 0;
  muNormChi2Global = 0;
  muIso03_sumPt = 0;
  muIso04_sumPt = 0;
  muIso05_sumPt = 0;
  muIso03_emEt = 0;
  muIso04_emEt = 0;
  muIso05_emEt = 0;
  muIso03_hadEt = 0;
  muIso04_hadEt = 0;
  muIso05_hadEt = 0;
  muIso03_nTracks = 0;
  muIso04_nTracks = 0;
  muIso05_nTracks = 0;
  muNotPFMuon = 0;
  muTrig = 0;
  trigPrescale = 0;
}


class TreeToDataset : public exception {
public :
  //!pointer to the analyzed TTree or TChain
  TChain          *fChain;

  // List of branches
  TBranch        *b_runNb;   //!
  TBranch        *b_eventNb;   //!
  TBranch        *b_LS;   //!
  TBranch        *b_CentBin;   //!
  TBranch        *b_Npix;   //!
  TBranch        *b_NpixelTracks;   //!
  TBranch        *b_Ntracks;   //!
  TBranch        *b_NtracksPtCut;   //!
  TBranch        *b_NtracksEtaCut;   //!
  TBranch        *b_NtracksEtaPtCut;   //!
  TBranch        *b_SumET_HF;   //!
  TBranch        *b_SumET_HFplus;   //!
  TBranch        *b_SumET_HFminus;   //!
  TBranch        *b_SumET_HFplusEta4;   //!
  TBranch        *b_SumET_HFminusEta4;   //!
  TBranch        *b_SumET_HFhit;   //!
  TBranch        *b_SumET_HFhitPlus;   //!
  TBranch        *b_SumET_HFhitMinus;   //!
  TBranch        *b_SumET_ZDC;   //!
  TBranch        *b_SumET_ZDCplus;   //!
  TBranch        *b_SumET_ZDCminus;   //!
  TBranch        *b_SumET_EEplus;   //!
  TBranch        *b_SumET_EEminus;   //!
  TBranch        *b_SumET_EE;   //!
  TBranch        *b_SumET_EB;   //!
  TBranch        *b_SumET_ET;   //!
  TBranch        *b_nPV;   //!
  TBranch        *b_RefVtx_x;   //!
  TBranch        *b_RefVtx_y;   //!
  TBranch        *b_RefVtx_z;   //!
  TBranch        *b_RefVtx_xError;   //!
  TBranch        *b_RefVtx_yError;   //!
  TBranch        *b_RefVtx_zError;   //!
  TBranch        *b_nPFpart;   //!
  TBranch        *b_pfId;   //!
  TBranch        *b_pfPt;   //!
  TBranch        *b_pfEnergy;   //!
  TBranch        *b_pfVsPt;   //!
  TBranch        *b_pfVsPtInitial;   //!
  TBranch        *b_pfArea;   //!
  TBranch        *b_pfEta;   //!
  TBranch        *b_pfPhi;   //!
  TBranch        *b_pfCharge;   //!
  TBranch        *b_pfTheta;   //!
  TBranch        *b_pfEt;   //!
  TBranch        *b_vn;   //!
  TBranch        *b_vpsi;   //!
  TBranch        *b_sumpt;   //!
  TBranch        *b_pfMuonPx;   //!
  TBranch        *b_pfMuonPy;   //!
  TBranch        *b_pfMuonPz;   //!
  TBranch        *b_pfTrackerMuon;   //!
  TBranch        *b_pfTrackerMuonPt;   //!
  TBranch        *b_pfTrackHits;   //!
  TBranch        *b_pfDxy;   //!
  TBranch        *b_pfDz;   //!
  TBranch        *b_pfChi2;   //!
  TBranch        *b_pfGlobalMuonPt;   //!
  TBranch        *b_pfChargedPx;   //!
  TBranch        *b_pfChargedPy;   //!
  TBranch        *b_pfChargedPz;   //!
  TBranch        *b_pfChargedTrackRefPt;   //!
  TBranch        *b_nGENpart;   //!
  TBranch        *b_genPDGId;   //!
  TBranch        *b_genPt;   //!
  TBranch        *b_genEta;   //!
  TBranch        *b_genPhi;   //!
  TBranch        *b_nTRACKpart;   //!
  TBranch        *b_traQual;   //!
  TBranch        *b_traCharge;   //!
  TBranch        *b_traPt;   //!
  TBranch        *b_traEta;   //!
  TBranch        *b_traPhi;   //!
  TBranch        *b_traAlgo;   //!
  TBranch        *b_traHits;   //!
  TBranch        *b_recoPFMET;   //!
  TBranch        *b_recoPFMETPhi;   //!
  TBranch        *b_recoPFMETsumEt;   //!
  TBranch        *b_recoPFMETmEtSig;   //!
  TBranch        *b_recoPFMETSig;   //!
  TBranch        *b_nMUpart;   //!
  TBranch        *b_muPx;   //!
  TBranch        *b_muPy;   //!
  TBranch        *b_muPz;   //!
  TBranch        *b_muMt;   //!
  TBranch        *b_muPt;   //!
  TBranch        *b_muEta;   //!
  TBranch        *b_muPhi;   //!
  TBranch        *b_muCharge;   //!
  TBranch        *b_muSelectionType;   //!
  TBranch        *b_muTrackIso;   //!
  TBranch        *b_muCaloIso;   //!
  TBranch        *b_muEcalIso;   //!
  TBranch        *b_muHcalIso;   //!
  TBranch        *b_muSumChargedHadronPt;
  TBranch        *b_muSumNeutralHadronEt;
  TBranch        *b_muSumPhotonEt;
  TBranch        *b_muSumPUPt;
  TBranch        *b_muPFBasedDBetaIso;
  TBranch        *b_muHighPurity;   //!
  TBranch        *b_muIsTightMuon;   //!
  TBranch        *b_muIsGoodMuon;   //!
  TBranch        *b_muTrkMuArb;   //!
  TBranch        *b_muTMOneStaTight;   //!
  TBranch        *b_muNTrkHits;   //!
  TBranch        *b_muNPixValHits;   //!
  TBranch        *b_muNPixWMea;   //!
  TBranch        *b_muNTrkWMea;   //!
  TBranch        *b_muStationsMatched;   //!
  TBranch        *b_muNMuValHits;   //!
  TBranch        *b_muDxy;   //!
  TBranch        *b_muDxyErr;   //!
  TBranch        *b_muDz;   //!
  TBranch        *b_muDzErr;   //!
  TBranch        *b_muPtInner;   //!
  TBranch        *b_muPtErrInner;   //!
  TBranch        *b_muPtGlobal;   //!
  TBranch        *b_muPtErrGlobal;   //!
  TBranch        *b_muNormChi2Inner;   //!
  TBranch        *b_muNormChi2Global;   //!
  TBranch        *b_muIso03_sumPt;
  TBranch        *b_muIso04_sumPt;
  TBranch        *b_muIso05_sumPt;
  TBranch        *b_muIso03_emEt;
  TBranch        *b_muIso04_emEt;
  TBranch        *b_muIso05_emEt;
  TBranch        *b_muIso03_hadEt;
  TBranch        *b_muIso04_hadEt;
  TBranch        *b_muIso05_hadEt;
  TBranch        *b_muIso03_nTracks;
  TBranch        *b_muIso04_nTracks;
  TBranch        *b_muIso05_nTracks;
  TBranch        *b_muNotPFMuon;
  TBranch        *b_muTrig;   //!
  TBranch        *b_HLTriggers;   //!
  TBranch        *b_trigPrescale;   //!

  vector<string> filename;  // input file names

  bool doMC;
  int trigIdx;
  int isoCut;
  float cutValue;
  
  TreePFCandEventData pfEvt_;
  
  RooRealVar *TMass;
  RooRealVar *MET;
  RooRealVar *Pt;
  RooRealVar *Eta;
  
  RooDataSet *dataset;

  TreeToDataset(vector<string> _filelist, bool _doMC, int _trigIdx, int _isoCut, float _cutValue);
  virtual ~TreeToDataset();
  virtual string   OpenInputs();
  virtual void     SetBranches();
  virtual void     MakeRooDataset();
  virtual int      Loop();
  bool CheckIsolation(int i_mu);
};


TreeToDataset::TreeToDataset(vector<string> _filelist, bool _doMC, int _trigIdx=5, int _isoCut=13, float _cutValue=0.1)
{
  doMC = _doMC;
  trigIdx = _trigIdx;
  isoCut = _isoCut;
  cutValue = _cutValue;

  // Copy filenames from a list
  for (vector<string>::size_type idx=0; idx!=_filelist.size(); idx++) {
    filename.push_back(_filelist[idx]);
  }
}

string TreeToDataset::OpenInputs() {
  // Check if input files are valid
  for (vector<string>::size_type idx=0; idx!=filename.size(); idx++) {
    if (!TFile::Open(filename[idx].c_str()))
      return string("Cannot open input file: ")+ filename[idx];
  }

  fChain = new TChain("pfcandAnalyzer/pfTree");
  // Load files
  for (vector<string>::size_type idx=0; idx!=filename.size(); idx++) {
    fChain->AddFile(filename[idx].c_str());
    cout << "Loading : " << filename[idx] << endl;
  }

  // Initalize tree/chain
  pfEvt_.Init();
  SetBranches();

  return "";
}


TreeToDataset::~TreeToDataset()
{
  delete TMass;
  delete MET;
  delete Pt;
  delete Eta;
  delete dataset;

  if (!fChain) return;
  delete fChain;
}


void TreeToDataset::SetBranches() {
  // Set branch addresses and branch pointers
  if (!fChain) return;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("runNb", &pfEvt_.runNb, &b_runNb);
  fChain->SetBranchAddress("eventNb", &pfEvt_.eventNb, &b_eventNb);
  fChain->SetBranchAddress("LS", &pfEvt_.LS, &b_LS);
  fChain->SetBranchAddress("CentBin", &pfEvt_.CentBin, &b_CentBin);
  fChain->SetBranchAddress("Npix", &pfEvt_.Npix, &b_Npix);
  fChain->SetBranchAddress("NpixelTracks", &pfEvt_.NpixelTracks, &b_NpixelTracks);
  fChain->SetBranchAddress("Ntracks", &pfEvt_.Ntracks, &b_Ntracks);
  fChain->SetBranchAddress("NtracksPtCut", &pfEvt_.NtracksPtCut, &b_NtracksPtCut);
  fChain->SetBranchAddress("NtracksEtaCut", &pfEvt_.NtracksEtaCut, &b_NtracksEtaCut);
  fChain->SetBranchAddress("NtracksEtaPtCut", &pfEvt_.NtracksEtaPtCut, &b_NtracksEtaPtCut);
  fChain->SetBranchAddress("SumET_HF", &pfEvt_.SumET_HF, &b_SumET_HF);
  fChain->SetBranchAddress("SumET_HFplus", &pfEvt_.SumET_HFplus, &b_SumET_HFplus);
  fChain->SetBranchAddress("SumET_HFminus", &pfEvt_.SumET_HFminus, &b_SumET_HFminus);
  fChain->SetBranchAddress("SumET_HFplusEta4", &pfEvt_.SumET_HFplusEta4, &b_SumET_HFplusEta4);
  fChain->SetBranchAddress("SumET_HFminusEta4", &pfEvt_.SumET_HFminusEta4, &b_SumET_HFminusEta4);
  fChain->SetBranchAddress("SumET_HFhit", &pfEvt_.SumET_HFhit, &b_SumET_HFhit);
  fChain->SetBranchAddress("SumET_HFhitPlus", &pfEvt_.SumET_HFhitPlus, &b_SumET_HFhitPlus);
  fChain->SetBranchAddress("SumET_HFhitMinus", &pfEvt_.SumET_HFhitMinus, &b_SumET_HFhitMinus);
  fChain->SetBranchAddress("SumET_ZDC", &pfEvt_.SumET_ZDC, &b_SumET_ZDC);
  fChain->SetBranchAddress("SumET_ZDCplus", &pfEvt_.SumET_ZDCplus, &b_SumET_ZDCplus);
  fChain->SetBranchAddress("SumET_ZDCminus", &pfEvt_.SumET_ZDCminus, &b_SumET_ZDCminus);
  fChain->SetBranchAddress("SumET_EEplus", &pfEvt_.SumET_EEplus, &b_SumET_EEplus);
  fChain->SetBranchAddress("SumET_EEminus", &pfEvt_.SumET_EEminus, &b_SumET_EEminus);
  fChain->SetBranchAddress("SumET_EE", &pfEvt_.SumET_EE, &b_SumET_EE);
  fChain->SetBranchAddress("SumET_EB", &pfEvt_.SumET_EB, &b_SumET_EB);
  fChain->SetBranchAddress("SumET_ET", &pfEvt_.SumET_ET, &b_SumET_ET);
  fChain->SetBranchAddress("nPV", &pfEvt_.nPV, &b_nPV);
  fChain->SetBranchAddress("RefVtx_x", &pfEvt_.RefVtx_x, &b_RefVtx_x);
  fChain->SetBranchAddress("RefVtx_y", &pfEvt_.RefVtx_y, &b_RefVtx_y);
  fChain->SetBranchAddress("RefVtx_z", &pfEvt_.RefVtx_z, &b_RefVtx_z);
  fChain->SetBranchAddress("RefVtx_xError", &pfEvt_.RefVtx_xError, &b_RefVtx_xError);
  fChain->SetBranchAddress("RefVtx_yError", &pfEvt_.RefVtx_yError, &b_RefVtx_yError);
  fChain->SetBranchAddress("RefVtx_zError", &pfEvt_.RefVtx_zError, &b_RefVtx_zError);
  fChain->SetBranchAddress("nPFpart", &pfEvt_.nPFpart, &b_nPFpart);
  fChain->SetBranchAddress("pfId", &pfEvt_.pfId, &b_pfId);
  fChain->SetBranchAddress("pfPt", &pfEvt_.pfPt, &b_pfPt);
  fChain->SetBranchAddress("pfEnergy", &pfEvt_.pfEnergy, &b_pfEnergy);
  fChain->SetBranchAddress("pfVsPt", &pfEvt_.pfVsPt, &b_pfVsPt);
  fChain->SetBranchAddress("pfVsPtInitial", &pfEvt_.pfVsPtInitial, &b_pfVsPtInitial);
  fChain->SetBranchAddress("pfArea", &pfEvt_.pfArea, &b_pfArea);
  fChain->SetBranchAddress("pfEta", &pfEvt_.pfEta, &b_pfEta);
  fChain->SetBranchAddress("pfPhi", &pfEvt_.pfPhi, &b_pfPhi);
  fChain->SetBranchAddress("pfCharge", &pfEvt_.pfCharge, &b_pfCharge);
  fChain->SetBranchAddress("pfTheta", &pfEvt_.pfTheta, &b_pfTheta);
  fChain->SetBranchAddress("pfEt", &pfEvt_.pfEt, &b_pfEt);
  fChain->SetBranchAddress("vn", &pfEvt_.vn, &b_vn);
  fChain->SetBranchAddress("psin", &pfEvt_.psin, &b_vpsi);
  fChain->SetBranchAddress("sumpt", &pfEvt_.sumpt, &b_sumpt);
  fChain->SetBranchAddress("pfMuonPx", &pfEvt_.pfMuonPx, &b_pfMuonPx);
  fChain->SetBranchAddress("pfMuonPy", &pfEvt_.pfMuonPy, &b_pfMuonPy);
  fChain->SetBranchAddress("pfMuonPz", &pfEvt_.pfMuonPz, &b_pfMuonPz);
  fChain->SetBranchAddress("pfTrackerMuon", &pfEvt_.pfTrackerMuon, &b_pfTrackerMuon);
  fChain->SetBranchAddress("pfTrackerMuonPt", &pfEvt_.pfTrackerMuonPt, &b_pfTrackerMuonPt);
  fChain->SetBranchAddress("pfTrackHits", &pfEvt_.pfTrackHits, &b_pfTrackHits);
  fChain->SetBranchAddress("pfDxy", &pfEvt_.pfDxy, &b_pfDxy);
  fChain->SetBranchAddress("pfDz", &pfEvt_.pfDz, &b_pfDz);
  fChain->SetBranchAddress("pfChi2", &pfEvt_.pfChi2, &b_pfChi2);
  fChain->SetBranchAddress("pfGlobalMuonPt", &pfEvt_.pfGlobalMuonPt, &b_pfGlobalMuonPt);
  fChain->SetBranchAddress("pfChargedPx", &pfEvt_.pfChargedPx, &b_pfChargedPx);
  fChain->SetBranchAddress("pfChargedPy", &pfEvt_.pfChargedPy, &b_pfChargedPy);
  fChain->SetBranchAddress("pfChargedPz", &pfEvt_.pfChargedPz, &b_pfChargedPz);
  fChain->SetBranchAddress("pfChargedTrackRefPt", &pfEvt_.pfChargedTrackRefPt, &b_pfChargedTrackRefPt);
  fChain->SetBranchAddress("nGENpart", &pfEvt_.nGENpart, &b_nGENpart);
  fChain->SetBranchAddress("genPDGId", &pfEvt_.genPDGId, &b_genPDGId);
  fChain->SetBranchAddress("genPt", &pfEvt_.genPt, &b_genPt);
  fChain->SetBranchAddress("genEta", &pfEvt_.genEta, &b_genEta);
  fChain->SetBranchAddress("genPhi", &pfEvt_.genPhi, &b_genPhi);
  fChain->SetBranchAddress("nTRACKpart", &pfEvt_.nTRACKpart, &b_nTRACKpart);
  fChain->SetBranchAddress("traQual", &pfEvt_.traQual, &b_traQual);
  fChain->SetBranchAddress("traCharge", &pfEvt_.traCharge, &b_traCharge);
  fChain->SetBranchAddress("traPt", &pfEvt_.traPt, &b_traPt);
  fChain->SetBranchAddress("traEta", &pfEvt_.traEta, &b_traEta);
  fChain->SetBranchAddress("traPhi", &pfEvt_.traPhi, &b_traPhi);
  fChain->SetBranchAddress("traAlgo", &pfEvt_.traAlgo, &b_traAlgo);
  fChain->SetBranchAddress("traHits", &pfEvt_.traHits, &b_traHits);
  fChain->SetBranchAddress("recoPFMET", &pfEvt_.recoPFMET, &b_recoPFMET);
  fChain->SetBranchAddress("recoPFMETPhi", &pfEvt_.recoPFMETPhi, &b_recoPFMETPhi);
  fChain->SetBranchAddress("recoPFMETsumEt", &pfEvt_.recoPFMETsumEt, &b_recoPFMETsumEt);
  fChain->SetBranchAddress("recoPFMETmEtSig", &pfEvt_.recoPFMETmEtSig, &b_recoPFMETmEtSig);
  fChain->SetBranchAddress("recoPFMETSig", &pfEvt_.recoPFMETSig, &b_recoPFMETSig);
  fChain->SetBranchAddress("nMUpart", &pfEvt_.nMUpart, &b_nMUpart);
  fChain->SetBranchAddress("muPx", &pfEvt_.muPx, &b_muPx);
  fChain->SetBranchAddress("muPy", &pfEvt_.muPy, &b_muPy);
  fChain->SetBranchAddress("muPz", &pfEvt_.muPz, &b_muPz);
  fChain->SetBranchAddress("muMt", &pfEvt_.muMt, &b_muMt);
  fChain->SetBranchAddress("muPt", &pfEvt_.muPt, &b_muPt);
  fChain->SetBranchAddress("muEta", &pfEvt_.muEta, &b_muEta);
  fChain->SetBranchAddress("muPhi", &pfEvt_.muPhi, &b_muPhi);
  fChain->SetBranchAddress("muCharge", &pfEvt_.muCharge, &b_muCharge);
  fChain->SetBranchAddress("muSelectionType", &pfEvt_.muSelectionType, &b_muSelectionType);
  fChain->SetBranchAddress("muTrackIso", &pfEvt_.muTrackIso, &b_muTrackIso);
  fChain->SetBranchAddress("muCaloIso", &pfEvt_.muCaloIso, &b_muCaloIso);
  fChain->SetBranchAddress("muEcalIso", &pfEvt_.muEcalIso, &b_muEcalIso);
  fChain->SetBranchAddress("muHcalIso", &pfEvt_.muHcalIso, &b_muHcalIso);
  fChain->SetBranchAddress("muSumChargedHadronPt",&pfEvt_.muSumChargedHadronPt, &b_muSumChargedHadronPt);
  fChain->SetBranchAddress("muSumNeutralHadronEt",&pfEvt_.muSumNeutralHadronEt, &b_muSumNeutralHadronEt);
  fChain->SetBranchAddress("muSumPhotonEt",&pfEvt_.muSumPhotonEt, &b_muSumPhotonEt);
  fChain->SetBranchAddress("muSumPUPt",&pfEvt_.muSumPUPt, &b_muSumPUPt);
  fChain->SetBranchAddress("muPFBasedDBetaIso",&pfEvt_.muPFBasedDBetaIso, &b_muPFBasedDBetaIso);
  fChain->SetBranchAddress("muHighPurity", &pfEvt_.muHighPurity, &b_muHighPurity);
  fChain->SetBranchAddress("muIsTightMuon", &pfEvt_.muIsTightMuon, &b_muIsTightMuon);
  fChain->SetBranchAddress("muIsGoodMuon", &pfEvt_.muIsGoodMuon, &b_muIsGoodMuon);
  fChain->SetBranchAddress("muTrkMuArb", &pfEvt_.muTrkMuArb, &b_muTrkMuArb);
  fChain->SetBranchAddress("muTMOneStaTight", &pfEvt_.muTMOneStaTight, &b_muTMOneStaTight);
  fChain->SetBranchAddress("muNTrkHits", &pfEvt_.muNTrkHits, &b_muNTrkHits);
  fChain->SetBranchAddress("muNPixValHits", &pfEvt_.muNPixValHits, &b_muNPixValHits);
  fChain->SetBranchAddress("muNPixWMea", &pfEvt_.muNPixWMea, &b_muNPixWMea);
  fChain->SetBranchAddress("muNTrkWMea", &pfEvt_.muNTrkWMea, &b_muNTrkWMea);
  fChain->SetBranchAddress("muStationsMatched", &pfEvt_.muStationsMatched, &b_muStationsMatched);
  fChain->SetBranchAddress("muNMuValHits", &pfEvt_.muNMuValHits, &b_muNMuValHits);
  fChain->SetBranchAddress("muDxy", &pfEvt_.muDxy, &b_muDxy);
  fChain->SetBranchAddress("muDxyErr", &pfEvt_.muDxyErr, &b_muDxyErr);
  fChain->SetBranchAddress("muDz", &pfEvt_.muDz, &b_muDz);
  fChain->SetBranchAddress("muDzErr", &pfEvt_.muDzErr, &b_muDzErr);
  fChain->SetBranchAddress("muPtInner", &pfEvt_.muPtInner, &b_muPtInner);
  fChain->SetBranchAddress("muPtErrInner", &pfEvt_.muPtErrInner, &b_muPtErrInner);
  fChain->SetBranchAddress("muPtGlobal", &pfEvt_.muPtGlobal, &b_muPtGlobal);
  fChain->SetBranchAddress("muPtErrGlobal", &pfEvt_.muPtErrGlobal, &b_muPtErrGlobal);
  fChain->SetBranchAddress("muNormChi2Inner", &pfEvt_.muNormChi2Inner, &b_muNormChi2Inner);
  fChain->SetBranchAddress("muNormChi2Global", &pfEvt_.muNormChi2Global, &b_muNormChi2Global);
  fChain->SetBranchAddress("muIso03_sumPt",&pfEvt_.muIso03_sumPt, &b_muIso03_sumPt);
  fChain->SetBranchAddress("muIso04_sumPt",&pfEvt_.muIso04_sumPt, &b_muIso04_sumPt);
  fChain->SetBranchAddress("muIso05_sumPt",&pfEvt_.muIso05_sumPt, &b_muIso05_sumPt);
  fChain->SetBranchAddress("muIso03_emEt",&pfEvt_.muIso03_emEt, &b_muIso03_emEt);
  fChain->SetBranchAddress("muIso04_emEt",&pfEvt_.muIso04_emEt, &b_muIso04_emEt);
  fChain->SetBranchAddress("muIso05_emEt",&pfEvt_.muIso05_emEt, &b_muIso05_emEt);
  fChain->SetBranchAddress("muIso03_hadEt",&pfEvt_.muIso03_hadEt, &b_muIso03_hadEt);
  fChain->SetBranchAddress("muIso04_hadEt",&pfEvt_.muIso04_hadEt, &b_muIso04_hadEt);
  fChain->SetBranchAddress("muIso05_hadEt",&pfEvt_.muIso05_hadEt, &b_muIso05_hadEt);
  fChain->SetBranchAddress("muIso03_nTracks",&pfEvt_.muIso03_nTracks, &b_muIso03_nTracks);
  fChain->SetBranchAddress("muIso04_nTracks",&pfEvt_.muIso04_nTracks, &b_muIso04_nTracks);
  fChain->SetBranchAddress("muIso05_nTracks",&pfEvt_.muIso05_nTracks, &b_muIso05_nTracks);
  fChain->SetBranchAddress("muNotPFMuon",&pfEvt_.muNotPFMuon, &b_muNotPFMuon);
  fChain->SetBranchAddress("muTrig", &pfEvt_.muTrig, &b_muTrig);
  fChain->SetBranchAddress("HLTriggers", &pfEvt_.HLTriggers, &b_HLTriggers);
  fChain->SetBranchAddress("trigPrescale", &pfEvt_.trigPrescale, &b_trigPrescale);
}



void TreeToDataset::MakeRooDataset() {
  TMass = new RooRealVar("TMass","Transverse Mass",0,160,"GeV/c^{2}");
  MET = new RooRealVar("MET","Missing E_{T}",0,100,"GeV");
  Pt = new RooRealVar("Pt","p_{T}",0,200,"GeV/c");
  Eta = new RooRealVar("Eta","#eta",-2.4,2.4,"");

  RooArgList varlist(*TMass,*MET,*Pt,*Eta);

  dataset = new RooDataSet("dataset","WDataSet",varlist);
}

 
