// -*- C++ -*-
//
// Package:    Run3ScoutingAnalysisTools/ScoutingPackedCandPhotonTreeMakerRun3Monitor
// Class:      ScoutingPackedCandPhotonTreeMakerRun3Monitor
//
/**\class ScoutingPackedCandPhotonTreeMakerRun3Monitor ScoutingPackedCandPhotonTreeMakerRun3Monitor.cc Run3ScoutingAnalysisTools/ScoutingTreeMakerRun3/plugins/ScoutingPackedCandPhotonTreeMakerRun3Monitor.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  David Sperka
//         Created:  Sat, 11 Feb 2023 14:15:08 GMT
//
//

// system include files
#include <memory>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
//
// class declaration
//

class ScoutingPackedCandPhotonTreeMakerRun3Monitor : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit ScoutingPackedCandPhotonTreeMakerRun3Monitor(const edm::ParameterSet&);
  ~ScoutingPackedCandPhotonTreeMakerRun3Monitor() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  const edm::InputTag triggerResultsTag;
  const edm::EDGetTokenT<edm::TriggerResults>                 triggerResultsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon> >      muonsToken;
  const edm::EDGetTokenT<std::vector<pat::Muon> >             offlineMuonsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron> >  electronsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >    primaryVerticesToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> >    verticesToken;
  const edm::EDGetTokenT<double>                              rhoToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPhoton> >    photonsToken;
  const edm::EDGetTokenT<std::vector<pat::Photon> >           offlinePhotonsToken;
  //const edm::EDGetTokenT<std::vector<pat::PackedCandidateCollection> > pfCandidatesToken;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pfLabel_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle> >  pfcandsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPFJet> >     pfjetsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingTrack> >     tracksToken;
  
  std::vector<std::string> triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;
  
  bool doL1;
  triggerExpression::Data triggerCache_;
  
  edm::InputTag                algInputTag_;
  edm::InputTag                extInputTag_;
  edm::EDGetToken              algToken_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string>     l1Seeds_;
  std::vector<std::string>     l1MonitorSeeds_;
  std::vector<bool>            l1Result_;
  std::vector<bool>            l1Result_mon_;
  
  TTree* tree;

  //Defining photons variables
  int   nScoutingPhotons;
  std::vector<float> pho_scout_pt;
  std::vector<float> pho_scout_eta;
  std::vector<float> pho_scout_phi;
  int   nOfflinePhotons;
  std::vector<float> pho_pt;
  std::vector<float> pho_eta;
  std::vector<float> pho_phi;
  std::vector<float> pho_r9;
  std::vector<float> pho_full5x5_r9;
  std::vector<float> pho_hoe;
  std::vector<float> pho_sieie;
  int   nPfPhotons;
  std::vector<float> pho_pf_pt;
  std::vector<float> pho_pf_eta;
  std::vector<float> pho_pf_phi;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ScoutingPackedCandPhotonTreeMakerRun3Monitor::ScoutingPackedCandPhotonTreeMakerRun3Monitor(const edm::ParameterSet& iConfig):
  triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerresults")),
  triggerResultsToken      (consumes<edm::TriggerResults>                    (triggerResultsTag)),
  muonsToken               (consumes<std::vector<Run3ScoutingMuon> >         (iConfig.getParameter<edm::InputTag>("muons"))),
  offlineMuonsToken        (consumes<std::vector<pat::Muon> >                (iConfig.getUntrackedParameter<edm::InputTag>("offlineMuons"))),
  electronsToken           (consumes<std::vector<Run3ScoutingElectron> >     (iConfig.getParameter<edm::InputTag>("electrons"))),
  primaryVerticesToken     (consumes<std::vector<Run3ScoutingVertex> >       (iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  verticesToken            (consumes<std::vector<Run3ScoutingVertex> >       (iConfig.getParameter<edm::InputTag>("displacedVertices"))),
  rhoToken                 (consumes<double>                                 (iConfig.getParameter<edm::InputTag>("rho"))), 
  photonsToken             (consumes<std::vector<Run3ScoutingPhoton> >       (iConfig.getParameter<edm::InputTag>("photons"))),
  offlinePhotonsToken      (consumes<std::vector<pat::Photon> >              (iConfig.getUntrackedParameter<edm::InputTag>("offlinePhotons"))),
  pfLabel_                 (consumes<pat::PackedCandidateCollection>         (iConfig.getUntrackedParameter<edm::InputTag>("pfLabel"))),
  pfcandsToken             (consumes<std::vector<Run3ScoutingParticle> >     (iConfig.getParameter<edm::InputTag>("pfcands"))),
  pfjetsToken              (consumes<std::vector<Run3ScoutingPFJet> >        (iConfig.getParameter<edm::InputTag>("pfjets"))),
  tracksToken              (consumes<std::vector<Run3ScoutingTrack> >        (iConfig.getParameter<edm::InputTag>("tracks"))),
  doL1                     (iConfig.existsAs<bool>("doL1")               ?    iConfig.getParameter<bool>  ("doL1") : false)
{
  usesResource("TFileService");
  if (doL1) {
    algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
    extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
    algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
    l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
    l1MonitorSeeds_ = iConfig.getParameter<std::vector<std::string> >("l1MonitorSeeds");
    l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
  }
  else {
    l1Seeds_ = std::vector<std::string>();
    l1MonitorSeeds_ = std::vector<std::string>();
    l1GtUtils_ = 0;
  }
}

ScoutingPackedCandPhotonTreeMakerRun3Monitor::~ScoutingPackedCandPhotonTreeMakerRun3Monitor() {
}

//
// member functions
//

// ------------ method called for each event  ------------
void ScoutingPackedCandPhotonTreeMakerRun3Monitor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;
  
  edm::Handle<edm::TriggerResults> triggerResultsH;
  iEvent.getByToken(triggerResultsToken, triggerResultsH);

  bool passDST=false;
  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
      if (triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) passDST=true;
  }

  // Scouting path true => scouting reconstruction enabled
  // -----------------------------------------------------
  if (!passDST) {/*cout<<"failed DST"<<endl;*/ return;}
  
  l1Result_.clear();

  // NOTE: No selection on the Monitoring seeds
  // -------------------------------------------
  //bool passMonitor=false; 
  //bool passScouting=false; 
  if (doL1) {
      l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);

      for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
          bool l1htbit = 0;
          l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
          l1Result_.push_back( l1htbit );
          //if (l1htbit) passScouting=true; // No selection on the Scouting seeds, selection stores
      }
      for( unsigned int iseed = 0; iseed < l1MonitorSeeds_.size(); iseed++ ) {
          bool l1htbit = 0;
          l1GtUtils_->getFinalDecisionByName(string(l1MonitorSeeds_[iseed]), l1htbit);
	  l1Result_mon_.push_back( l1htbit );
          //if (l1htbit) passMonitor=true; // No selection on the Monitoring seeds, selection stored
      }
  }
  //if (!passMonitor) {/*cout<<"failed L1 seed"<<endl;*/ return;} // No selection on the Monitoring seeds
  //if (!passScouting) {/*cout<<"failed L1 seed"<<endl;*/ return;} // No selection on the Scouting seeds
    
 


  // --------------------- //
  // PHOTONS PF CANDIDATES //
  // --------------------- //
  edm::Handle<pat::PackedCandidateCollection> pfPhotonsH;
  iEvent.getByToken(pfLabel_,pfPhotonsH);
  //Handle<vector<pat::PackedCandidateCollection> > pfPhotonsH;
  //iEvent.getByToken(pfCandidatesToken, pfPhotonsH);  

  nPfPhotons=0;
  vector<int> pfpho_idx;
  int pf=0;

  pho_pf_pt.clear();
  pho_pf_eta.clear();
  pho_pf_phi.clear();

  // Minimal selection on the offline photons
  // ----------------------------------------
  for (auto pfpho_iter = pfPhotonsH->begin(); pfpho_iter != pfPhotonsH->end(); ++pfpho_iter) {
    //for (auto pat::PackedCandidateCollection::const_iterator pfpho_iter = pfPhotonsH->begin(); pfpho_iter!=pfPhotonsH->end(); ++pfpho_iter){
    if (pfpho_iter -> pdgId()== 22 && pfpho_iter->et()>1 && abs(pfpho_iter->eta())<2.5) { 
      //cout<<"offline photon pt: "<<photons_iter->pt()<<" eta: "<<photons_iter->eta()<<" phi: "<<photons_iter->phi()<<endl;                          
      nPfPhotons+=1;
      pho_pf_pt.push_back(pfpho_iter->et());
      pho_pf_eta.push_back(pfpho_iter->eta());
      pho_pf_phi.push_back(pfpho_iter->phi());
    }
    pf+=1;
  }

  cout << "ParticleFlow PHOTONS = " << nPfPhotons << endl;



  // --------------- //
  // OFFLINE PHOTONS //
  // --------------- //
  Handle<vector<pat::Photon> > offlinePhotonsH;
  iEvent.getByToken(offlinePhotonsToken, offlinePhotonsH);  
  //if (offlinePhotonsH->size()<1) return;

  nOfflinePhotons=0;
  vector<int> pho_idx;
  int p=0;

  pho_pt.clear();
  pho_eta.clear();
  pho_phi.clear();
  pho_r9.clear();
  pho_full5x5_r9.clear();
  pho_hoe.clear();
  pho_sieie.clear();

  // Minimal selection on the offline photons
  // ----------------------------------------
  for (auto photons_iter = offlinePhotonsH->begin(); photons_iter != offlinePhotonsH->end(); ++photons_iter) {
    if (photons_iter->pt()>1 &&  abs(photons_iter->eta())<2.5) { 
      //cout<<"offline photon pt: "<<photons_iter->pt()<<" eta: "<<photons_iter->eta()<<" phi: "<<photons_iter->phi()<<endl;                          
      nOfflinePhotons+=1;
      pho_pt.push_back(photons_iter->pt());
      pho_eta.push_back(photons_iter->eta());
      pho_phi.push_back(photons_iter->phi());
      pho_r9.push_back(photons_iter->r9());
      pho_full5x5_r9.push_back(photons_iter->full5x5_r9());
      pho_hoe.push_back(photons_iter->hadronicOverEm());
      pho_sieie.push_back(photons_iter->sigmaIetaIeta());
    }
    p+=1;
  }
  //cout << "OFFLINE PHOTONS = " << nOfflinePhotons << endl;

  // ---------------- //
  // SCOUTING PHOTONS //
  // ---------------- //
  Handle<vector<Run3ScoutingPhoton> > photonsH;
  iEvent.getByToken(photonsToken, photonsH);
  nScoutingPhotons=0;
  vector<int> pho_idx_scout;
  int p_s=0;

  pho_scout_pt.clear();
  pho_scout_eta.clear();
  pho_scout_phi.clear();

 // Minimal selection on the scouting photons
  // ---------------------------------------
  for (auto photons_scout_iter = photonsH->begin(); photons_scout_iter != photonsH->end(); ++photons_scout_iter) {
    if (photons_scout_iter->pt()>1 && abs(photons_scout_iter->eta())<2.5) {
      //cout<<"scouting photon pt: "<<photons_scout_iter->pt()<<" eta: "<<photons_scout_iter->eta()<<" phi: "<<photons_scout_iter->phi()<<endl;                
      nScoutingPhotons+=1;
      pho_idx_scout.push_back(p_s);
      pho_scout_pt.push_back(photons_scout_iter->pt());
      pho_scout_eta.push_back(photons_scout_iter->eta());
      pho_scout_phi.push_back(photons_scout_iter->phi());
    }
    p_s+=1;
  }

  /* //DEBUGGING
  std::cout << "##############################################################################" << endl;
  std::cout << "nOfflinePhotons = " << nOfflinePhotons << "   and pho_pt.size() = " << pho_pt.size() << endl;
  if (pho_pt.size() > 0){
    for (long unsigned int n=0; n < pho_pt.size(); n++)
      {
	std::cout << "---------------------" << endl;
	std::cout << "* Pho pt = " << pho_pt.at(n) << endl;
	std::cout << "* Pho eta = " << pho_eta.at(n) << endl;
	std::cout << "* Pho phi = " << pho_phi.at(n) << endl;
      }
  }
  std::cout << "##############################################################################" << endl;
  std::cout << "nScoutingPhotons = " << nScoutingPhotons << "   and pho_scout_pt.size() = " << pho_scout_pt.size() << endl;
  if (pho_scout_pt.size() > 0){
    for (long unsigned int f=0; f < pho_scout_pt.size(); f++)
      {
	std::cout << "---------------------" << endl;
	std::cout << "* Pho scout pt = " << pho_scout_pt.at(f) << endl;
	std::cout << "* Pho scout eta = " << pho_scout_eta.at(f) << endl;
	std::cout << "* Pho scout phi = " << pho_scout_phi.at(f) << endl;
      }
  }
  */
  /*if (pho_idx_scout.size()==1){
    cout << "SCOUTING PHOTONS = " << nScoutingPhotons << endl;
    pho_scout_pt=photonsH->at(pho_idx_scout[0]).pt();
    pho_scout_eta=photonsH->at(pho_idx_scout[0]).eta();
    pho_scout_phi=photonsH->at(pho_idx_scout[0]).phi();
    pho_pt=offlinePhotonsH->at(pho_idx[0]).pt();
    pho_eta=offlinePhotonsH->at(pho_idx[0]).eta();
    pho_phi=offlinePhotonsH->at(pho_idx[0]).phi();
    }*/
  
  tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void ScoutingPackedCandPhotonTreeMakerRun3Monitor::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"      , "tree");
    tree->Branch("nOfflinePhotons"     , &nOfflinePhotons              , "nOfflinePhotons/I");
    tree->Branch("nScoutingPhotons"    , &nScoutingPhotons             , "nScoutingPhotons/I");
    tree->Branch("nPfPhotons"          , &nPfPhotons                   , "nPfPhotons/I");
    tree->Branch("l1Result"            , "std::vector<bool>"           , &l1Result_, 32000, 0);
    tree->Branch("l1Result_mon"        , "std::vector<bool>"           , &l1Result_mon_, 32000, 0);
    tree->Branch("pho_scout_pt"        , "std::vector<float>"          , &pho_scout_pt, 32000, 0);                                    
    tree->Branch("pho_scout_eta"       , "std::vector<float>"          , &pho_scout_eta, 32000, 0);                                    
    tree->Branch("pho_scout_phi"       , "std::vector<float>"          , &pho_scout_phi, 32000, 0);                                    
    tree->Branch("pho_pt"              , "std::vector<float>"          , &pho_pt, 32000, 0);                                    
    tree->Branch("pho_eta"             , "std::vector<float>"          , &pho_eta, 32000, 0);                                    
    tree->Branch("pho_phi"             , "std::vector<float>"          , &pho_phi, 32000, 0);                                    
    tree->Branch("pho_r9"              , "std::vector<float>"          , &pho_r9, 32000, 0);                                    
    tree->Branch("pho_full5x5_r9"      , "std::vector<float>"          , &pho_full5x5_r9, 32000, 0);                                      
    tree->Branch("pho_hoe"           , "std::vector<float>"          , &pho_hoe, 32000, 0);                                    
    tree->Branch("pho_sieie"           , "std::vector<float>"          , &pho_sieie, 32000, 0);                                    
    tree->Branch("pho_pf_pt"           , "std::vector<float>"          , &pho_pf_pt, 32000, 0);                                    
    tree->Branch("pho_pf_eta"          , "std::vector<float>"          , &pho_pf_eta, 32000, 0);                                    
    tree->Branch("pho_pf_phi"          , "std::vector<float>"          , &pho_pf_phi, 32000, 0);                                    
}

// ------------ method called once each job just after ending the event loop  ------------
void ScoutingPackedCandPhotonTreeMakerRun3Monitor::endJob() {
  // please remove this method if not needed
}

void ScoutingPackedCandPhotonTreeMakerRun3Monitor::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
    // HLT paths
    triggerPathsVector.push_back("DST_Run3_PFScoutingPixelTracking_v*");

    HLTConfigProvider hltConfig;
    bool changedConfig = false;
    hltConfig.init(iRun, iSetup, triggerResultsTag.process(), changedConfig);

    for (size_t i = 0; i < triggerPathsVector.size(); i++) {
        triggerPathsMap[triggerPathsVector[i]] = -1;
    }

    for(size_t i = 0; i < triggerPathsVector.size(); i++){
        TPRegexp pattern(triggerPathsVector[i]);
        for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
            std::string pathName = hltConfig.triggerNames()[j];
            if(TString(pathName).Contains(pattern)){
                triggerPathsMap[triggerPathsVector[i]] = j;
            }
        }
    }
}


void ScoutingPackedCandPhotonTreeMakerRun3Monitor::endRun(edm::Run const&, edm::EventSetup const&) {
}

void ScoutingPackedCandPhotonTreeMakerRun3Monitor::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void ScoutingPackedCandPhotonTreeMakerRun3Monitor::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ScoutingPackedCandPhotonTreeMakerRun3Monitor::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ScoutingPackedCandPhotonTreeMakerRun3Monitor);
