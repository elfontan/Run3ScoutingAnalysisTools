// -*- C++ -*-
//
// Package:    Run3ScoutingAnalysisTools/gen_standardRecoTreeMakerRun3
// Class:      gen_standardRecoTreeMakerRun3
//
/**\class gen_standardRecoTreeMakerRun3 gen_standardRecoTreeMakerRun3.cc Run3ScoutingAnalysisTools/ScoutingTreeMakerRun3/plugins/gen_standardRecoTreeMakerRun3.cc

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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
//
// class declaration
//

class gen_standardRecoTreeMakerRun3 : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit gen_standardRecoTreeMakerRun3(const edm::ParameterSet&);
  ~gen_standardRecoTreeMakerRun3() override;

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
  const edm::EDGetTokenT<edm::TriggerResults>                  triggerResultsToken;
  const edm::EDGetTokenT<std::vector<reco::GenParticle> >      genParticleToken;
  const edm::EDGetTokenT<std::vector<pat::PackedGenParticle> > packedGenParticleToken;
  const edm::EDGetTokenT<std::vector<pat::Muon> >              offlineMuonsToken;
  const edm::EDGetTokenT<std::vector<pat::Photon> >            offlinePhotonsToken;
  const edm::EDGetTokenT<pat::PackedCandidateCollection>       pfLabel_;
  
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
  int nEta_mumugamma_gen = 0;
  int nEta_mumugamma_pho10_gen = 0;
  int nEta_mumugamma_gen_twoRecoOsMu = 0;
  int nEta_mumugamma_gen_twoRecoOsMu_recoPho = 0;
  int nEta_mumugamma_gen_twoRecoOsMu_recoPhoMatched = 0;
  int nEta_mumugamma_gen_twoRecoOsMu_recoPfPho = 0;

  int nEta;
  int nMuPlus;
  int nMuMinus; 
  int nGamma;

  // Defining gen variables
  float etagenmup_pt;
  float etagenmum_pt;
  float etagenmup_eta;
  float etagenmum_eta;
  float etagenmup_phi;
  float etagenmum_phi;
  float etagenpho_pt; 
  float etagenpho_eta;
  float etagenpho_phi;

  float gen_dimu_mass;
  float gen_mmg_mass;
  float gen_mmg_dr;
  float gen_mmg_pt;
  float gen_mmg_pho10_mass;

  //Defining reco variables
  float recoMu_mm_mass;
  float recoMu_recoPho_mmg_mass;
  float recoMu_recoPho_mmg_dr;
  float recoMu_recoPhoMatched_mmg_mass;
  float recoMu_recoPhoMatched_mmg_dr;
  float recoMu_recoPfPhoMatched_mmg_mass;
  float recoMu_recoPfPhoMatched_mmg_dr;

  float sel_recomup_pt;
  float sel_recomup_eta;
  float sel_recomup_phi;
  float sel_recomum_pt;
  float sel_recomum_eta;
  float sel_recomum_phi;
  float sel_recopfPho_pt;
  float sel_recopfPho_eta;
  float sel_recopfPho_phi;
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
gen_standardRecoTreeMakerRun3::gen_standardRecoTreeMakerRun3(const edm::ParameterSet& iConfig):
  triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerresults")),
  triggerResultsToken      (consumes<edm::TriggerResults>                    (triggerResultsTag)),
  genParticleToken         (consumes<std::vector<reco::GenParticle> >        (iConfig.getUntrackedParameter<edm::InputTag>("genP"))),
  packedGenParticleToken   (consumes<std::vector<pat::PackedGenParticle> >   (iConfig.getUntrackedParameter<edm::InputTag>("packedGenP"))),
  offlineMuonsToken        (consumes<std::vector<pat::Muon> >                (iConfig.getUntrackedParameter<edm::InputTag>("offlineMuons"))),
  offlinePhotonsToken      (consumes<std::vector<pat::Photon> >              (iConfig.getUntrackedParameter<edm::InputTag>("offlinePhotons"))),
  pfLabel_                 (consumes<pat::PackedCandidateCollection>         (iConfig.getUntrackedParameter<edm::InputTag>("pfLabel"))),
  doL1                     (iConfig.existsAs<bool>("doL1")              ?    iConfig.getParameter<bool>  ("doL1") : false)
{
  usesResource("TFileService");
  if (doL1) {
    algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
    extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
    algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
    l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
    l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
  }
  else {
    l1Seeds_ = std::vector<std::string>();
    l1MonitorSeeds_ = std::vector<std::string>();
    l1GtUtils_ = 0;
  }
}

gen_standardRecoTreeMakerRun3::~gen_standardRecoTreeMakerRun3() {
}

// ---------------- //
// Member functions //
// ---------------- //

// ------------ Method called for each event  ------------
int ev = 0;

void gen_standardRecoTreeMakerRun3::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  using namespace reco;

  std::cout << "------------------------------- "  << std::endl;
  std::cout << "------------EVENT " << ev << std::endl;
  std::cout << "------------------------------- "  << std::endl;
  ev+=1;

  // ------------- //
  // GEN PARTICLES //
  // ------------- //
  Handle<vector<reco::GenParticle> > genP;
  iEvent.getByToken(genParticleToken, genP);
  genP.product();
  Handle<vector<pat::PackedGenParticle> > packedGenP;
  iEvent.getByToken(packedGenParticleToken, packedGenP);
 
  vector<float> pt_genmup;
  vector<float> pt_genmum; 
  vector<float> pt_genpho;
  vector<float> mup_motherPt;
  vector<float> mum_motherPt;
  vector<float> pho_motherPt;
  vector<float> eta_genmup;
  vector<float> eta_genmum;
  vector<float> eta_genpho;
  vector<float> phi_genmup;
  vector<float> phi_genmum;
  vector<float> phi_genpho;
  vector<float> pfpho_recopt;
  vector<float> pfpho_recoeta;
  vector<float> pfpho_recophi;

  nEta = 0;
  nMuPlus = 0;
  nMuMinus = 0; 
  nGamma = 0;
  //gen_dimu_mass = 0;

  for (auto gen_iter = genP->begin(); gen_iter != genP->end(); ++gen_iter) {
    if (gen_iter->pdgId()==221 && gen_iter->numberOfDaughters()==3 && gen_iter->status()==2)
      {
    	nEta+=1;
      }
  }
  /*
  for (auto gen_iter = genP->begin(); gen_iter != genP->end(); ++gen_iter) {
    if (gen_iter->pdgId()==221 && gen_iter->numberOfDaughters()==3 && gen_iter->status()==2)
      {
    	nEta+=1;
      }
    if (gen_iter->pdgId()==13 && gen_iter->status()==1 && gen_iter->motherRef()->pdgId() == 221)
      //if (gen_iter->pdgId()==13 && gen_iter->status()==1)
      {
	nMuPlus+=1;
	pt_genmup.push_back(gen_iter->pt());
	eta_genmup.push_back(gen_iter->eta());
	phi_genmup.push_back(gen_iter->phi());
      }
    if (gen_iter->pdgId()==-13 && gen_iter->status()==1 && gen_iter->motherRef()->pdgId() == 221)
      //if (gen_iter->pdgId()==-13 && gen_iter->status()==1)
      {
	nMuMinus+=1;
	pt_genmum.push_back(gen_iter->pt());
	eta_genmum.push_back(gen_iter->eta());
	phi_genmum.push_back(gen_iter->phi());
      }
    if (gen_iter->pdgId()==22 && gen_iter->status()==1 && gen_iter->motherRef()->pdgId() == 221)
      //if (abs(gen_iter->pdgId())==22 && gen_iter->status()==1)
      {
	nGamma+=1;
	pt_genpho.push_back(gen_iter->pt());
	eta_genpho.push_back(gen_iter->eta());
	phi_genpho.push_back(gen_iter->phi());
      }
  }
  */
  for (auto packedgen_iter = packedGenP->begin(); packedgen_iter != packedGenP->end(); ++packedgen_iter) {
    //if (packedgen_iter->pdgId()==221 && packedgen_iter->numberOfDaughters()==3 && packedgen_iter->status()==2)
    //  {
    //	nEta+=1;
    //  }
    if (packedgen_iter->pdgId()==13 && packedgen_iter->status()==1 && packedgen_iter->motherRef()->pdgId() == 221 && packedgen_iter->pt() > 3)
      //if (packedgen_iter->pdgId()==13 && packedgen_iter->status()==1)
      {
	nMuPlus+=1;
	pt_genmup.push_back(packedgen_iter->pt());
	eta_genmup.push_back(packedgen_iter->eta());
	phi_genmup.push_back(packedgen_iter->phi());
	mup_motherPt.push_back(packedgen_iter->motherRef()->pt());
      }
    if (packedgen_iter->pdgId()==-13 && packedgen_iter->status()==1 && packedgen_iter->motherRef()->pdgId() == 221 && packedgen_iter->pt() > 3)
      //if (packedgen_iter->pdgId()==-13 && packedgen_iter->status()==1)
      {
	nMuMinus+=1;
	pt_genmum.push_back(packedgen_iter->pt());
	eta_genmum.push_back(packedgen_iter->eta());
	phi_genmum.push_back(packedgen_iter->phi());
	mum_motherPt.push_back(packedgen_iter->motherRef()->pt());
      }
    if (packedgen_iter->pdgId()==22 && packedgen_iter->status()==1 && packedgen_iter->motherRef()->pdgId() == 221 && packedgen_iter->pt() > 1)
      //if (abs(packedgen_iter->pdgId())==22 && packedgen_iter->status()==1)
      {
	nGamma+=1;
	//std::cout << "PHOTON from " << packedgen_iter->pdgId() << " with pT = " << packedgen_iter->pt() << std::endl;
	pt_genpho.push_back(packedgen_iter->pt());
	eta_genpho.push_back(packedgen_iter->eta());
	phi_genpho.push_back(packedgen_iter->phi());
	pho_motherPt.push_back(packedgen_iter->motherRef()->pt());
      }
  }

  std::cout << "-------------------------------" << nEta << std::endl;
  std::cout << "nEta = " << nEta << std::endl;
  std::cout << "nMuPlus = " << nMuPlus << std::endl;
  std::cout << "nMuMinus = " << nMuMinus << std::endl;
  std::cout << "nGamma = " << nGamma << std::endl;
  std::cout << "-------------------------------" << nEta << std::endl;

  if (!(nEta == 1 && nMuPlus > 0  && nMuMinus > 0 && nGamma > 0)) {return;}
 
  TLorentzVector reco_mu1;
  TLorentzVector reco_mu2;
  TLorentzVector reco_dimu;
  TLorentzVector reco_pho;
  TLorentzVector reco_phoMatched;
  TLorentzVector reco_mmg;
  TLorentzVector reco_mmg_matched;
  TLorentzVector reco_pfPho;
  TLorentzVector reco_mmg_pf;

  for (long unsigned int mp = 0; mp < pt_genmup.size(); ++mp)
    {
      for (long unsigned int mm = 0; mm < pt_genmum.size(); ++mm)
	{
	  for (long unsigned int g = 0; g < pt_genpho.size(); ++g)
	    {
	      TLorentzVector m1, m2; 
	      TLorentzVector gen_dimu;
	      TLorentzVector pho;
	      TLorentzVector gen_mmg;
	      std::cout << "MOTHER PT = " << mup_motherPt.at(mp) << "  , " << mum_motherPt.at(mm) << "  , " << pho_motherPt.at(g) << endl;
	      if (!(mup_motherPt.at(mp) == mum_motherPt.at(mm) && mup_motherPt.at(mp) == pho_motherPt.at(g))){continue;}
	      m1.SetPtEtaPhiM(pt_genmup.at(mp), eta_genmup.at(mp), phi_genmup.at(mp), 0.106);	
	      m2.SetPtEtaPhiM(pt_genmum.at(mm), eta_genmum.at(mm), phi_genmum.at(mm), 0.106);
	      pho.SetPtEtaPhiM(pt_genpho.at(g), eta_genpho.at(g), phi_genpho.at(g), 0.);
	      std::cout << "pt_genmup.at(mp) = " << pt_genmup.at(mp) << " eta_genmup.at(mp) = " << eta_genmup.at(mp) << std::endl;
	      std::cout << "pt_genmum.at(mm) = " << pt_genmum.at(mm) << " eta_genmum.at(mm) = " << eta_genmum.at(mm) << std::endl;
	      std::cout << "pt_genpho.at(mm) = " << pt_genpho.at(g) << " eta_genpho.at(g) = " << eta_genpho.at(g) << std::endl;
	      gen_dimu = m1 + m2;
	      //if (gen_dimu.M() > 0.548) continue;
	      gen_mmg = gen_dimu + pho;
	      //std::cout << "pt_genpho.at(g) = " << pt_genpho.at(g) << " eta_genpho.at(g) = " << eta_genpho.at(g) << std::endl;	  
	      //std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> mmg mass = " << gen_mmg.M() << std::endl; 

 	      if (abs(gen_mmg.M() - 0.548) < abs(gen_dimu.M() - 0.548)) 
		{
		  etagenmup_pt = pt_genmup.at(mp); 
		  etagenmup_eta = eta_genmup.at(mp); 
		  etagenmup_phi = phi_genmup.at(mp);
		  etagenmum_pt = pt_genmum.at(mm); 
		  etagenmum_eta = eta_genmum.at(mm); 
		  etagenmum_phi = phi_genmum.at(mm);
		  etagenpho_pt = pt_genpho.at(g); 
		  etagenpho_eta = eta_genpho.at(g); 
		  etagenpho_phi = phi_genpho.at(g);
		  nEta_mumugamma_gen+=1;
		  gen_dimu_mass = gen_dimu.M();
		  gen_mmg_mass = gen_mmg.M();
		  gen_mmg_pt = gen_mmg.Pt();
		  gen_mmg_dr = gen_dimu.DeltaR(pho);

		  std::cout << "########## GEN LEVEL ##########" << endl;		 
		  std::cout << "mm mass = " << gen_dimu.M() << std::endl; 
		  std::cout << "gen_mmg_mass = " << gen_mmg_mass << std::endl;	 
		  std::cout << "etagenmup_pt = " << etagenmup_pt << std::endl;	 
		  std::cout << "etagenmum_pt = " << etagenmum_pt << std::endl;	 
		  std::cout << "etagenpho_pt = " << etagenpho_pt << std::endl;	 
		  if (etagenpho_pt > 10) 
		    {
		      nEta_mumugamma_pho10_gen+=1;
		      gen_mmg_pho10_mass = gen_mmg.M();
		      std::cout << ">10 GeV! etagenpho_pt = " << etagenpho_pt << std::endl;	 
		      std::cout << "gen_mmg_pho10_mass = " << gen_mmg_pho10_mass << std::endl;	 
		    }
		  //std::cout << "etagenmup_pt = " << etagenmup_pt << std::endl;	 
		  //std::cout << "etagenmum_pt = " << etagenmum_pt << std::endl;	 
		  //std::cout << "etagenpho_pt = " << etagenpho_pt << std::endl;	 
		}
	    }
	}
    }

  // ------------- //
  // OFFLINE MUONS //
  // ------------- //
  Handle<vector<pat::Muon> > offlineMuonsH;
  iEvent.getByToken(offlineMuonsToken, offlineMuonsH);	   
  if (offlineMuonsH->size() < 2) return;
  
  std::cout << "########## MUONS ##########" << endl;		 
  std::cout << "offlineMuonsH->size() = " << offlineMuonsH->size() << endl;
  vector<int> idx;
  int j=0;
  float dr_genmup = 999.;
  float dr_genmum = 999.;
  
  // Minimal selection on the offline muons
  // --------------------------------------
  for (auto muons_iter = offlineMuonsH->begin(); muons_iter != offlineMuonsH->end(); ++muons_iter) {
    if (muons_iter->pt() > 3 && abs(muons_iter->eta()) < 2.4) { //removing the MediumMuon ID requirement, saved later
      //std::cout << " --- muons_iter->pt() = " << muons_iter->pt() << std::endl;
      if (muons_iter->pdgId() == 13)
	{
	  dr_genmup = sqrt( (muons_iter->eta() - etagenmup_eta)*(muons_iter->eta() - etagenmup_eta) + (muons_iter->phi() - etagenmup_phi)*(muons_iter->phi() - etagenmup_phi) ); 
	  //std::cout << " --- dr_genmup " << dr_genmup << std::endl;
	  if (dr_genmup < 0.01) 
	    {
	      idx.push_back(j);
	      //std::cout << " --- etagenmup_pt = " << etagenmup_pt << std::endl;
	      reco_mu1.SetPtEtaPhiM(muons_iter->pt(), muons_iter->eta(), muons_iter->phi(), 0.106);
	    }
	}
      if (muons_iter->pdgId() == -13)
	{
	  dr_genmum = sqrt( (muons_iter->eta() - etagenmum_eta)*(muons_iter->eta() - etagenmum_eta) + (muons_iter->phi() - etagenmum_phi)*(muons_iter->phi() - etagenmum_phi) ); 
	  //std::cout << " --- dr_genmum " << dr_genmum << std::endl;
	  if (dr_genmum < 0.01) 
	    {
	      idx.push_back(j);
	      //std::cout << " --- etagenmum_pt = " << etagenmum_pt << std::endl;
	      reco_mu2.SetPtEtaPhiM(muons_iter->pt(), muons_iter->eta(), muons_iter->phi(), 0.106);
	    }
	} 
    }
    j+=1;
  }
  
  std::cout << "idx.size() = " << idx.size() << std::endl;
  if (!(idx.size()==2)) {/*cout<<"Two reco muons matched to gen muons"<<endl;*/ return;}  
  // Checking opposite charge                                                                                     
  // ------------------------
  int checkOSCharge = offlineMuonsH->at(idx[0]).charge() * offlineMuonsH->at(idx[1]).charge();
  cout<<"Check OS requirement for offline muons = " << checkOSCharge << endl; 
  cout<<"My muons are = " << offlineMuonsH->at(idx[0]).pt() << " and " << offlineMuonsH->at(idx[1]).pt()<< endl; 
  if (checkOSCharge > 0){/*cout<<"not OS pair"<<endl;*/return;}
  //cout<<"Check OS requirement for offline muons = PASSED!" << endl; 
  nEta_mumugamma_gen_twoRecoOsMu+=1;		  
  reco_dimu = reco_mu1 + reco_mu2;
  recoMu_mm_mass = reco_dimu.M();
  cout<<"** RECO DIMUON MASS = " << recoMu_mm_mass << endl; 
  sel_recomup_pt = offlineMuonsH->at(idx[0]).pt();
  sel_recomum_pt = offlineMuonsH->at(idx[1]).pt();
  sel_recomup_eta = offlineMuonsH->at(idx[0]).eta();
  sel_recomum_eta = offlineMuonsH->at(idx[1]).eta();
  sel_recomup_phi = offlineMuonsH->at(idx[0]).phi();
  sel_recomum_phi = offlineMuonsH->at(idx[1]).phi();
  

  // -------------- //
  // OFFLINE PHOTON //
  // -------------- //
  Handle<vector<pat::Photon> > offlinePhotonsH;
  iEvent.getByToken(offlinePhotonsToken, offlinePhotonsH);	   
  std::cout << "########## PHOTONS ##########" << endl;		 
  std::cout << "offlinePhotonsH->size() = " << offlinePhotonsH->size() << endl;
  int p=0;
  vector<int> pho_idx;
  vector<int> matchedPho_idx;
  float dr_genpho = 999.;
  
  // Minimal selection on the offline photon
  // ----------------------------------------
  for (auto pho_iter = offlinePhotonsH->begin(); pho_iter != offlinePhotonsH->end(); ++pho_iter) {
    std::cout << "pho_iter->pt() = " << pho_iter->pt() << ", pho_iter->eta() = " << pho_iter->eta() << " , pho_iter->phi() = " << pho_iter->phi() << std::endl;
    if (pho_iter->pt()>10 && abs(pho_iter->eta())<2.5) {
      dr_genpho = sqrt( (pho_iter->eta() - etagenpho_eta)*(pho_iter->eta() - etagenpho_eta) + (pho_iter->phi() - etagenpho_phi)*(pho_iter->phi() - etagenpho_phi) );
      pho_idx.push_back(p);
      //reco_pho.SetPtEtaPhiM(pho_iter->pt(), pho_iter->eta(), pho_iter->phi(), 0.);
      if (dr_genpho < 0.2){
	std::cout << "--- dR matching = " << dr_genpho << std::endl;
	std::cout << "--- MATCHED PHOTON pt = " << pho_iter->pt() << "  with gen PHOTON pt = " << etagenpho_pt << std::endl;
	matchedPho_idx.push_back(p);
	reco_phoMatched.SetPtEtaPhiM(pho_iter->pt(), pho_iter->eta(), pho_iter->phi(), 0.);
      }
    }
    p+=1;
  }
  
  std::cout << "--- pho_idx.size() = " << pho_idx.size() << std::endl;
  std::cout << "--- matchedPho_idx.size() = " << matchedPho_idx.size() << std::endl;
  
  if (pho_idx.size()>0){
    float dr_reco_dimu_pho = 999.;
    float dr_min = 999.;
    for (long unsigned int pidx = 0; pidx < pho_idx.size(); ++pidx){          		    
      dr_min = sqrt( (reco_dimu.Eta() - offlinePhotonsH->at(pho_idx[pidx]).eta())*(reco_dimu.Eta() - offlinePhotonsH->at(pho_idx[pidx]).eta()) + (reco_dimu.Phi() - offlinePhotonsH->at(pho_idx[pidx]).phi())*(reco_dimu.Phi() - offlinePhotonsH->at(pho_idx[pidx]).phi()) );
      //std::cout << "@@ dr_min = " << dr_min << std::endl;	 
      if (dr_min < dr_reco_dimu_pho)
	{
	  //std::cout << "@@ dr_reco_dimu_pho = " << dr_reco_dimu_pho << std::endl;	 
	  reco_pho.SetPtEtaPhiM(offlinePhotonsH->at(pho_idx[pidx]).pt(), offlinePhotonsH->at(pho_idx[pidx]).eta(), offlinePhotonsH->at(pho_idx[pidx]).phi(), 0.);
	  dr_reco_dimu_pho = dr_min;
	}
    }
    
    reco_mmg = reco_dimu + reco_pho;
    recoMu_recoPho_mmg_dr = reco_mmg.DeltaR(reco_pho);
    std::cout << "@@ recoMu_recoPho_mmg_dr = " << recoMu_recoPho_mmg_dr << std::endl;  
    recoMu_recoPho_mmg_mass = reco_mmg.M();
    std::cout << "recoMu_recoPho_mmg_mass = " << recoMu_recoPho_mmg_mass << std::endl;	 
    nEta_mumugamma_gen_twoRecoOsMu_recoPho+=1;
  }
  
  if (matchedPho_idx.size()>0) 
    {
      nEta_mumugamma_gen_twoRecoOsMu_recoPhoMatched+=1;
      reco_mmg_matched = reco_dimu + reco_phoMatched;
      recoMu_recoPhoMatched_mmg_mass = reco_mmg_matched.M();		      
      std::cout << "recoMu_recoPhoMatched_mmg_mass = " << recoMu_recoPhoMatched_mmg_mass << std::endl;	 
      recoMu_recoPhoMatched_mmg_dr = reco_mmg.DeltaR(reco_phoMatched);		      
    }		  		  
  
  // --------------------- //                                                            
  // PHOTONS PF CANDIDATES //                                                                           
  // --------------------- //                                                                               
  edm::Handle<pat::PackedCandidateCollection> pfPhotonsH;
  iEvent.getByToken(pfLabel_,pfPhotonsH);
  
  int pf = 0;
  vector<int> pfpho_idx;
  float dr_genpho_pf = 999;
  for (auto pfpho_iter = pfPhotonsH->begin(); pfpho_iter != pfPhotonsH->end(); ++pfpho_iter) {
    if (pfpho_iter -> pdgId() == 22 && pfpho_iter->pt() > 1 && abs(pfpho_iter->eta()) < 2.5) {
      dr_genpho_pf = sqrt( (pfpho_iter->eta() - etagenpho_eta)*(pfpho_iter->eta() - etagenpho_eta) + (pfpho_iter->phi() - etagenpho_phi)*(pfpho_iter->phi() - etagenpho_phi) );
      //std::cout << "dr_genpho_pf = " << dr_genpho_pf << endl;                                                             
      if (dr_genpho_pf < 0.2){
	pfpho_idx.push_back(pf);
	pfpho_recopt.push_back(pfpho_iter->pt());
	pfpho_recoeta.push_back(pfpho_iter->eta());
	pfpho_recophi.push_back(pfpho_iter->phi());
      }
    }
    pf+=1;
  }
  if (pfpho_idx.size() > 0){
    std::cout << "pfpho_idx.size() = " << pfpho_idx.size() << endl;
    float dr_reco_dimu_pfpho = 999.;
    float dr_min_pf = 999.;
    for (long unsigned int gpf = 0; gpf < pfpho_idx.size(); ++gpf)
      {
	dr_min_pf = sqrt( (reco_dimu.Eta() - pfpho_recoeta.at(gpf))*(reco_dimu.Eta() - pfpho_recoeta.at(gpf)) + (reco_dimu.Phi() - pfpho_recophi.at(gpf))*(reco_dimu.Phi() - pfpho_recophi.at(gpf)));
	//std::cout << "@@ dr_min = " << dr_min << std::endl;	 
	if (dr_min_pf < dr_reco_dimu_pfpho)
	  {
	    //std::cout << "@@ dr_reco_dimu_pho = " << dr_reco_dimu_pho << std::endl;	 
	    reco_pfPho.SetPtEtaPhiM(pfpho_recopt.at(gpf), pfpho_recoeta.at(gpf), pfpho_recophi.at(gpf),0.);
	    dr_reco_dimu_pfpho = dr_min_pf;
	  }
      }
    nEta_mumugamma_gen_twoRecoOsMu_recoPfPho+=1;
    reco_mmg_pf = reco_dimu + reco_pfPho;
    recoMu_recoPfPhoMatched_mmg_mass = reco_mmg_pf.M();		      
    recoMu_recoPfPhoMatched_mmg_dr = reco_mmg_pf.DeltaR(reco_pfPho);		      
    std::cout << "%% recoMu_recoPfPhoMatched_mmg_mass = " << recoMu_recoPfPhoMatched_mmg_mass << std::endl; 	       
    std::cout << "%% with photon  = " << reco_pfPho.Pt() << std::endl; 	       
    sel_recopfPho_pt  = reco_pfPho.Pt();
    sel_recopfPho_eta = reco_pfPho.Eta();
    sel_recopfPho_phi = reco_pfPho.Phi();
  }
  
  
  std::cout << "-------------------------------------------" << std::endl;	 
  std::cout << "nEta_mumugamma_gen = " << nEta_mumugamma_gen << std::endl;	 
  std::cout << "nEta_mumugamma_pho10_gen = " << nEta_mumugamma_pho10_gen << std::endl;	 
  std::cout << "nEta_mumugamma_gen_twoRecoOsMu = " << nEta_mumugamma_gen_twoRecoOsMu << std::endl;	 
  std::cout << "nEta_mumugamma_gen_twoRecoOsMu_recoPho = " << nEta_mumugamma_gen_twoRecoOsMu_recoPho << std::endl;	 
  std::cout << "nEta_mumugamma_gen_twoRecoOsMu_recoPhoMatched = " << nEta_mumugamma_gen_twoRecoOsMu_recoPhoMatched << std::endl;	 
  std::cout << "nEta_mumugamma_gen_twoRecoOsMu_recoPfPho = " << nEta_mumugamma_gen_twoRecoOsMu_recoPfPho << std::endl;	 
  std::cout << "-------------------------------------------" << std::endl;	 
  
  tree->Fill();
  std::cout << ">>>>>>>>>>>>>>>> FILL!!!!!!!!!!!!!!!!!!!!" << std::endl;  
  std::cout << "Note: recoMu_mm_mass is " << recoMu_mm_mass << std::endl;  
  pt_genmup.clear();
  pt_genmum.clear();
  pt_genpho.clear();
  eta_genmup.clear();
  eta_genmum.clear();
  eta_genpho.clear();
  phi_genmup.clear();
  phi_genmum.clear();
  phi_genpho.clear();
}

// ------------ method called once each job just before starting event loop  ------------
void gen_standardRecoTreeMakerRun3::beginJob() {
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree"      , "tree");
    tree->Branch("nEta"                , &nEta                , "nEta/I");
    tree->Branch("nMuPlus"             , &nMuPlus             , "nMuPlus/I");
    tree->Branch("nMuMinus"            , &nMuMinus            , "nMuMinus/I");
    tree->Branch("nGamma"              , &nGamma              , "nGamma/I");
    tree->Branch("gen_dimu_mass"       , &gen_dimu_mass       , "gen_dimu_mass/F");
    tree->Branch("gen_mmg_mass"        , &gen_mmg_mass        , "gen_mmg_mass/F");
    tree->Branch("gen_mmg_pt"          , &gen_mmg_pt          , "gen_mmg_pt/F");
    tree->Branch("gen_mmg_dr"          , &gen_mmg_dr          , "gen_mmg_dr/F");
    tree->Branch("gen_mmg_pho10_mass"  , &gen_mmg_pho10_mass  , "gen_mmg_pho10_mass/F");

    tree->Branch("etagenmup_pt"        , &etagenmup_pt          , "etagenmup_pt/F");
    tree->Branch("etagenmup_eta"       , &etagenmup_eta         , "etagenmup_eta/F");
    tree->Branch("etagenmup_phi"       , &etagenmup_phi         , "etagenmup_phi/F");
    tree->Branch("etagenmum_pt"        , &etagenmum_pt          , "etagenmum_pt/F");
    tree->Branch("etagenmum_eta"       , &etagenmum_eta         , "etagenmum_eta/F");
    tree->Branch("etagenmum_phi"       , &etagenmum_phi         , "etagenmum_phi/F");
    tree->Branch("etagenpho_pt"        , &etagenpho_pt          , "etagenpho_pt/F");
    tree->Branch("etagenpho_eta"       , &etagenpho_eta         , "etagenpho_eta/F");
    tree->Branch("etagenpho_phi"       , &etagenpho_phi         , "etagenpho_phi/F");

    tree->Branch("sel_recomup_pt"          , &sel_recomup_pt          , "sel_recomup_pt/F");
    tree->Branch("sel_recomup_eta"         , &sel_recomup_eta         , "sel_recomup_eta/F");
    tree->Branch("sel_recomup_phi"         , &sel_recomup_phi         , "sel_recomup_phi/F");
    tree->Branch("sel_recomum_pt"          , &sel_recomum_pt          , "sel_recomum_pt/F");
    tree->Branch("sel_recomum_eta"         , &sel_recomum_eta         , "sel_recomum_eta/F");
    tree->Branch("sel_recomum_phi"         , &sel_recomum_phi         , "sel_recomum_phi/F");
    tree->Branch("sel_recopfPho_pt"        , &sel_recopfPho_pt        , "sel_recopfPho_pt/F");
    tree->Branch("sel_recopfPho_eta"       , &sel_recopfPho_eta       , "sel_recopfPho_eta/F");
    tree->Branch("sel_recopfPho_phi"       , &sel_recopfPho_phi       , "sel_recopfPho_phi/F");

    tree->Branch("recoMu_mm_mass"                   , &recoMu_mm_mass                   , "recoMu_mm_mass/F");
    tree->Branch("recoMu_recoPho_mmg_mass"          , &recoMu_recoPho_mmg_mass          , "recoMu_recoPho_mmg_mass/F");
    tree->Branch("recoMu_recoPhoMatched_mmg_mass"   , &recoMu_recoPhoMatched_mmg_mass   , "recoMu_recoPhoMatched_mmg_mass/F");
    tree->Branch("recoMu_recoPho_mmg_dr"            , &recoMu_recoPho_mmg_dr            , "recoMu_recoPho_mmg_dr/F");
    tree->Branch("recoMu_recoPhoMatched_mmg_dr"     , &recoMu_recoPhoMatched_mmg_dr     , "recoMu_recoPhoMatched_mmg_dr/F");
    tree->Branch("recoMu_recoPfPhoMatched_mmg_mass" , &recoMu_recoPfPhoMatched_mmg_mass , "recoMu_recoPfPhoMatched_mmg_mass/F");
    tree->Branch("recoMu_recoPfPhoMatched_mmg_dr"   , &recoMu_recoPfPhoMatched_mmg_dr   , "recoMu_recoPfPhoMatched_mmg_dr/F");    
}

// ------------ method called once each job just after ending the event loop  ------------
void gen_standardRecoTreeMakerRun3::endJob() {
  // please remove this method if not needed
}

void gen_standardRecoTreeMakerRun3::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
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


void gen_standardRecoTreeMakerRun3::endRun(edm::Run const&, edm::EventSetup const&) {
}

void gen_standardRecoTreeMakerRun3::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void gen_standardRecoTreeMakerRun3::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void gen_standardRecoTreeMakerRun3::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(gen_standardRecoTreeMakerRun3);
