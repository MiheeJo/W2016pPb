#include "Inputs.h"
#include "TreeToDataset.h"

bool TreeToDataset::CheckIsolation(int i_mu) {
  bool isolation = false;

  if (isoCut==0) return true;

  if (isoCut==13) { // PF-based combined relative isolation without this muon, without duplication (cone=0.3)
    float sumEtInCone = ( pfEvt_.muIso03_sumPt->at(i_mu) + pfEvt_.muIso03_emEt->at(i_mu) + pfEvt_.muIso03_hadEt->at(i_mu) ) / pfEvt_.muPt->at(i_mu);
    if ( sumEtInCone < cutValue )
      isolation = true;
  }
  else if (isoCut==14) { // PF-based combined relative isolation without this muon, without duplication (cone=0.3)
    float sumEtInCone = ( pfEvt_.muIso04_sumPt->at(i_mu) + pfEvt_.muIso04_emEt->at(i_mu) + pfEvt_.muIso04_hadEt->at(i_mu) ) / pfEvt_.muPt->at(i_mu);
    if ( sumEtInCone < cutValue )
      isolation = true;
  }
  else if (isoCut==15) { // PF-based combined relative isolation without this muon, without duplication (cone=0.3)
    float sumEtInCone = ( pfEvt_.muIso05_sumPt->at(i_mu) + pfEvt_.muIso05_emEt->at(i_mu) + pfEvt_.muIso05_hadEt->at(i_mu) ) / pfEvt_.muPt->at(i_mu);
    if ( sumEtInCone < cutValue )
      isolation = true;
  }
  else if (isoCut==2) { // PF-based combined relative isolation with delta beta correction (cone=0.4)
    if ( pfEvt_.muPFBasedDBetaIso->at(i_mu) < cutValue )
      isolation = true;
  }
  else if (isoCut==21) { // PF-based combined relative isolation without delta beta correction (cone=0.4)
    float sumEtInCone = pfEvt_.muSumChargedHadronPt->at(i_mu);
    float Et = pfEvt_.muSumNeutralHadronEt->at(i_mu) + pfEvt_.muSumPhotonEt->at(i_mu);
    if (0.<Et) sumEtInCone += Et;
    sumEtInCone /= pfEvt_.muPt->at(i_mu);
    if ( sumEtInCone < cutValue )
      isolation = true;
  }
  else if (isoCut==3) { // Tracker-based relative isolation
    if ( pfEvt_.muTrackIso->at(i_mu)/pfEvt_.muPt->at(i_mu) < cutValue )
      isolation = true;
  }

  return isolation;
}


int TreeToDataset::Loop()
{
  if (fChain == 0) return -1;

  Long64_t nentries = fChain->GetEntries();
  for (Long64_t evt=0; evt<nentries; evt++) {
    if ( evt%100000 == 0 ) cout << "Event: " << evt  << " / " << nentries << endl;
    fChain->GetEntry(evt);

    if ( pfEvt_.nMUpart != pfEvt_.muPt->size() ) {
      cout << "pfEvt_.nMUpart != muPt->size() AT " << evt << endl;
      return -1;
    }

    for ( int i_mu=0; i_mu < pfEvt_.nMUpart; i_mu++) {
      bool tightSelection = false, triggerSelection = false, isolation = false;

      if ( pfEvt_.muIsTightMuon->at(i_mu) ) tightSelection = true;

      if ( ( pfEvt_.muTrig->at(i_mu)&((ULong64_t)pow(2,trigIdx)) )==( (ULong64_t)pow(2,trigIdx) ) &&
           ( pfEvt_.HLTriggers&((ULong64_t)pow(2,trigIdx)) )==( (ULong64_t)pow(2,trigIdx) ) )
        triggerSelection = true;

      isolation = CheckIsolation(i_mu);

      if (tightSelection && triggerSelection && isolation) {
        MET->setVal(pfEvt_.recoPFMET);
        TMass->setVal(pfEvt_.muMt->at(i_mu));
        Pt->setVal(pfEvt_.muPt->at(i_mu));
        Eta->setVal(pfEvt_.muEta->at(i_mu));

        RooArgList varlist(*TMass,*MET,*Pt,*Eta);

        dataset->add(varlist);
      }
    } // end of i_mu loop
   
  } // end of evt loop

  return 0;
}

int main(int argc, char* argv[]) {

  /// *** Parse input options
  Inputs Opt(argc, argv);
  if (Opt.ParseOptions()) return -1; // When wrong inputs received
  Opt.ShowOptions();

  /// *** Load Files/Trees
  TreeToDataset *ITrees = new TreeToDataset(Opt.sources, Opt.doMC, 5, Opt.isoCut, 0.1); //trigIdx=5
  string out = ITrees->OpenInputs(); 
  if (out!="") {
    cout << out << endl;
    delete ITrees;
    return -1;
  }

  /// *** RooDataSet to be written
  ITrees->MakeRooDataset();
  if (ITrees->Loop()) {
    cout << "Problem while reading events\n";
    return -1;
  }

  /// *** Output TFile with RooDataSet
  TFile* Out = new TFile(Opt.outputname.c_str(),"RECREATE");
  Out->cd();
  ITrees->dataset->Write();
  Out->Close();


  return 0;
}
