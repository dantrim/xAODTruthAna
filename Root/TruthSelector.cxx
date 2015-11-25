#include "xAODTruthAna/TruthSelector.h"

//std/stl
#include <iostream>
#include <cstdio> // printf
#include <map>
using namespace std;

//ROOT
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TNamed.h"

#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODEventInfo/EventInfo.h"

//EDM
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

ClassImp(TruthSelector)


const double mev2gev = 1.0e-3;
typedef map<int, int> JftMultiplicity; // first=flavor, second=multiplicity

#undef CHECK
#define CHECK( ARG )    \
    do {    \
        const bool result = ARG; \
        if( ! result ) {    \
            ::Error("TruthSelector", "Failed to execute: \"%s\"", \
                #ARG ); \
            exit(1);    \
        }   \
    } while (false)


// ---------------------------------------------------------
TruthSelector::TruthSelector() :
    m_inputSampleName(""),
    m_outputFileName("superTruth.root"),
    m_outTreeFile(NULL),
    m_outputTree(NULL),
    m_tree(NULL),
    m_event(new xAOD::TEvent(xAOD::TEvent::kBranchAccess) ),
    m_dbg(false)
{
    cout << "TruthSelector" << endl;
}
TruthSelector::~TruthSelector()
{
    delete m_event;
}

// ---------------------------------------------------------
void TruthSelector::Init(TTree *tree)
{
    cout << "TruthSelector::Init" << endl;
    cout << "Input sample : " << m_inputSampleName << endl;
    if(tree) {
        cout << " >> Found tree with " << tree->GetEntries() << " entries" << endl;
    }
    else {
        cout << " >> Init received a null pointer" << endl;
        exit(1);
    }
    m_tree = tree;
    return;
}
// ---------------------------------------------------------
Bool_t TruthSelector::Notify()
{
    if(!m_tree) {
        cout << "TruthSelector::Notify    Tree pointer is null!" << endl;
        return kTRUE;
    }
    cout << "TruthSelector::Notify   Opening a new file" << endl;
    TFile* ifile = 0;
    TChain* chain = dynamic_cast<TChain*>(m_tree);
    if(chain) {
        ifile = chain->GetFile();
        cout << "chain : " << ifile->GetName() << endl;
    }
    else {
        ifile = m_tree->GetCurrentFile();
        cout << "not chain : " << ifile << endl;
        //cout << "not chain : " << ifile->GetName() << endl;
    }

    RETURN_CHECK(GetName(), m_event->readFrom(ifile, kFALSE) );
    return kTRUE;
}
// ---------------------------------------------------------
void TruthSelector::SlaveBegin(TTree*)
{
    cout << "TruthSelector::SlaveBegin" << endl;
    //xAOD::Init("TruthSelector").ignore();
    //m_event.readFrom(tree);
    initializeOutputTree();
    initializeTreeBranches();

}
void TruthSelector::Begin(TTree*)
{
    cout << "TruthSelector::Begin" << endl;
    return;
}
// ---------------------------------------------------------
void TruthSelector::initializeOutputTree()
{
    if(m_dbg > 10) cout << "TruthSelector::initializeOutputTree" << endl;

    m_outTreeFile = new TFile(m_outputFileName.c_str(), "recreate");
    m_outputTree = new TTree("superTruth", "superTruth");
    m_outputTree->SetAutoSave(10000000);
    m_outputTree->SetMaxTreeSize(3000000000u);
}
// ---------------------------------------------------------
void TruthSelector::initializeTreeBranches()
{
    if(m_dbg > 10) cout << "TruthSelector::initializeTreeBranches" << endl;

    // event vars
    m_outputTree->Branch("runNumber",
                              &m_br_runNumber             );
    m_outputTree->Branch("eventNumber",
                              &m_br_eventNumber           );
    m_outputTree->Branch("eventweight",
                              &m_br_eventWeight           );

    // lepton vars
    m_outputTree->Branch("nLeptons",
                              &m_br_lepton_n              );
    m_outputTree->Branch("l_pt",
                              &m_br_lepton_pt             );
    m_outputTree->Branch("l_eta",
                              &m_br_lepton_eta            );
    m_outputTree->Branch("l_phi",
                              &m_br_lepton_phi            );
    m_outputTree->Branch("l_pdgId",
                              &m_br_lepton_pdgId          );
    m_outputTree->Branch("l_q",
                              &m_br_lepton_charge         );
    m_outputTree->Branch("l_type",
                              &m_br_lepton_type           );
    m_outputTree->Branch("l_origin",
                              &m_br_lepton_origin         );
    m_outputTree->Branch("isSF",
                              &m_br_isSF                  );
    m_outputTree->Branch("isDF",
                              &m_br_isDF                  );
    m_outputTree->Branch("isEE",
                              &m_br_isEE                  );
    m_outputTree->Branch("isMM",
                              &m_br_isMM                  );
    m_outputTree->Branch("isEM",
                              &m_br_isEM                  );

    // jet vars
    m_outputTree->Branch("j_n",     &m_br_jet_n                 );
    m_outputTree->Branch("j_pt",    &m_br_jet_pt                );
    m_outputTree->Branch("j_eta",   &m_br_jet_eta               );
    m_outputTree->Branch("j_phi",   &m_br_jet_phi               );
    m_outputTree->Branch("j_m",     &m_br_jet_m                 );
    m_outputTree->Branch("j_flav",  &m_br_jet_flav              );

    // ljet vars
    m_outputTree->Branch("lj_n",     &m_br_ljet_n                 );
    m_outputTree->Branch("lj_pt",    &m_br_ljet_pt                );
    m_outputTree->Branch("lj_eta",   &m_br_ljet_eta               );
    m_outputTree->Branch("lj_phi",   &m_br_ljet_phi               );
    m_outputTree->Branch("lj_m",     &m_br_ljet_m                 );
    m_outputTree->Branch("lj_flav",  &m_br_ljet_flav              );

    // cjet vars
    m_outputTree->Branch("cj_n",     &m_br_cjet_n                 );
    m_outputTree->Branch("cj_pt",    &m_br_cjet_pt                );
    m_outputTree->Branch("cj_eta",   &m_br_cjet_eta               );
    m_outputTree->Branch("cj_phi",   &m_br_cjet_phi               );
    m_outputTree->Branch("cj_m",     &m_br_cjet_m                 );
    m_outputTree->Branch("cj_flav",  &m_br_cjet_flav              );

    // bjet vars
    m_outputTree->Branch("bj_n",     &m_br_bjet_n                 );
    m_outputTree->Branch("bj_pt",    &m_br_bjet_pt                );
    m_outputTree->Branch("bj_eta",   &m_br_bjet_eta               );
    m_outputTree->Branch("bj_phi",   &m_br_bjet_phi               );
    m_outputTree->Branch("bj_m",     &m_br_bjet_m                 );
    m_outputTree->Branch("bj_flav",  &m_br_bjet_flav              );

    // hjet vars
    m_outputTree->Branch("hj_n",     &m_br_hjet_n                 );
    m_outputTree->Branch("hj_pt",    &m_br_hjet_pt                );
    m_outputTree->Branch("hj_eta",   &m_br_hjet_eta               );
    m_outputTree->Branch("hj_phi",   &m_br_hjet_phi               );
    m_outputTree->Branch("hj_m",     &m_br_hjet_m                 );
    m_outputTree->Branch("hj_flav",  &m_br_hjet_flav              );

    // gjet vars
    m_outputTree->Branch("gj_n",     &m_br_gjet_n                 );
    m_outputTree->Branch("gj_pt",    &m_br_gjet_pt                );
    m_outputTree->Branch("gj_eta",   &m_br_gjet_eta               );
    m_outputTree->Branch("gj_phi",   &m_br_gjet_phi               );
    m_outputTree->Branch("gj_m",     &m_br_gjet_m                 );
    m_outputTree->Branch("gj_flav",  &m_br_gjet_flav              );




}
void TruthSelector::clearContainers()
{

    m_br_runNumber = 0;
    m_br_eventNumber = 0;
    m_br_eventWeight = 0;
    m_br_eventWeights.clear();

    m_br_lepton_pt.clear();
    m_br_lepton_eta.clear();
    m_br_lepton_phi.clear();
    m_br_lepton_pdgId.clear();
    m_br_lepton_charge.clear();
    m_br_lepton_type.clear();
    m_br_lepton_origin.clear();
    m_br_isSF = 0;
    m_br_isDF = 0;
    m_br_isEE = 0;
    m_br_isMM = 0;
    m_br_isEM = 0;

    m_br_jet_n = 0;
    m_br_jet_pt.clear();
    m_br_jet_eta.clear();
    m_br_jet_phi.clear();
    m_br_jet_m.clear();
    m_br_jet_flav.clear();

    m_br_ljet_n = 0;
    m_br_ljet_pt.clear();
    m_br_ljet_eta.clear();
    m_br_ljet_phi.clear();
    m_br_ljet_m.clear();
    m_br_ljet_flav.clear();

    m_br_bjet_n = 0;
    m_br_bjet_pt.clear();
    m_br_bjet_eta.clear();
    m_br_bjet_phi.clear();
    m_br_bjet_m.clear();
    m_br_bjet_flav.clear();

    m_br_cjet_n = 0;
    m_br_cjet_pt.clear();
    m_br_cjet_eta.clear();
    m_br_cjet_phi.clear();
    m_br_cjet_m.clear();
    m_br_cjet_flav.clear();

    m_br_hjet_n = 0;
    m_br_hjet_pt.clear();
    m_br_hjet_eta.clear();
    m_br_hjet_phi.clear();
    m_br_hjet_m.clear();
    m_br_hjet_flav.clear();

    m_br_gjet_n = 0;
    m_br_gjet_pt.clear();
    m_br_gjet_eta.clear();
    m_br_gjet_phi.clear();
    m_br_gjet_m.clear();
    m_br_gjet_flav.clear();

}
// ---------------------------------------------------------
void TruthSelector::saveOutputTree()
{
    cout << "TruthSelector::saveOutputTree" <<endl;

    m_outTreeFile = m_outputTree->GetCurrentFile();
    if(m_outTreeFile) {
        TDirectory* current_directory = gROOT->CurrentDirectory();
        m_outTreeFile->cd();
        TNamed inputSampleName("inputSample", m_inputSampleName.c_str());
        inputSampleName.Write();
        current_directory->cd();
    }
    m_outputTree->Write(0, TObject::kOverwrite);
    cout << "output tree saved to " << m_outTreeFile->GetName() << endl;
    m_outTreeFile->Close();
}

// ------------------------------------------------------------
// Define some convenience methods

template <class Container>
double sumPt(Container &container)
{
    double totPt = 0.0;
    for(const auto *e : container)
        totPt += (e->pt() * mev2gev);
    return totPt;
}

struct IsPdg {
    int pdg;
    IsPdg(int p) : pdg(p) {}
    bool operator()(const xAOD::TruthParticle* l) { return (l && l->absPdgId()==pdg); }
};
IsPdg isElectron(11), isMuon(13), isTau(15);

struct PtGreater {
    bool operator()(const xAOD::TruthParticle* a, const xAOD::TruthParticle* b) { return a->pt()>b->pt(); }
} byPt;

struct PtGreaterJet {
    bool operator()(const xAOD::Jet* a, const xAOD::Jet* b) { return a->pt()>b->pt(); }
} byPtJet;

bool is_ee(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1) { return isElectron(l0) and isElectron(l1); }
bool is_em(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1) { return ((isElectron(l0) and isMuon(l1)) or
                                                                                   (isMuon(l0) and isElectron(l1))); }
bool is_mm(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1) { return isMuon(l0) and isMuon(l1); }

// ---------------------------------------------------------
Bool_t TruthSelector::Process(Long64_t entry)
{
    if(m_dbg > 10) cout << "TruthSelector::Process" << endl;
    static Long64_t chainEntry = -1;
    chainEntry++;

    m_event->getEntry(chainEntry);

    // clear the variables that we store
    clearContainers();

    // ---------------------------------
    //  Event information
    // ---------------------------------
    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), m_event->retrieve(ei, "EventInfo") );

    m_br_runNumber   = ei->runNumber();
    m_br_eventNumber = ei->eventNumber();
    m_br_eventWeight = ei->mcEventWeight();

    if(m_dbg || chainEntry%5000==0){
        cout << "*************************************************************" << endl;
        cout << "**** Processing entry " << setw(6) << chainEntry
             << "  run " << setw(6) << m_br_runNumber
             << "  event " << setw(7) << m_br_eventNumber<< " ****" << endl;
        cout << "*************************************************************" << endl;
    }

    // ---------------------------------
    //  Jets
    // ---------------------------------
    const xAOD::JetContainer *jets = 0;
    RETURN_CHECK( GetName(), m_event->retrieve( jets, "AntiKt4TruthJets" ));

    std::vector<xAOD::Jet*> v_jet;

    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();
    for( ; jet_itr != jet_end; ++jet_itr) {
        if( ((*jet_itr)->pt() * mev2gev) < 20 ) continue;
        xAOD::Jet* jet = new xAOD::Jet();
        jet->makePrivateStore( **jet_itr );
        v_jet.push_back(jet);
    }
    sort(v_jet.begin(), v_jet.end(), byPtJet);

    // ---------------------------------
    //  Electrons
    // ---------------------------------
    const xAOD::TruthParticleContainer* electrons = 0;
    RETURN_CHECK( GetName(), m_event->retrieve( electrons, "TruthElectrons" ));

    std::vector<xAOD::TruthParticle*> v_electron;

    xAOD::TruthParticleContainer::const_iterator ele_itr = electrons->begin();
    xAOD::TruthParticleContainer::const_iterator ele_end = electrons->end();
    for( ; ele_itr != ele_end ; ++ele_itr ) {
        if( ((*ele_itr)->pt() * mev2gev) < 10 ) continue;
        if( fabs((*ele_itr)->eta()) > 2.47 ) continue;
        bool is_prompt = ((*ele_itr)->auxdata<unsigned int>( "classifierParticleType" ) == 2);
        if(not is_prompt) continue;
        xAOD::TruthParticle* electron = new xAOD::TruthParticle();
        electron->makePrivateStore(**ele_itr);
        v_electron.push_back(electron);
    }
    sort(v_electron.begin(), v_electron.end(), byPt);

    // ---------------------------------
    //  Muons
    // ---------------------------------
    const xAOD::TruthParticleContainer* muons = 0;
    RETURN_CHECK( GetName(), m_event->retrieve( muons, "TruthMuons" ));

    std::vector<xAOD::TruthParticle*> v_muon;

    xAOD::TruthParticleContainer::const_iterator mu_itr = muons->begin();
    xAOD::TruthParticleContainer::const_iterator mu_end = muons->end();
    for( ; mu_itr != mu_end ; ++mu_itr) {
        if( ((*mu_itr)->pt() * mev2gev) < 10 ) continue;
        if( fabs((*mu_itr)->eta()) > 2.4 ) continue;
        bool is_prompt = ((*mu_itr)->auxdata<unsigned int>( "classifierParticleType" ) == 6);
        if(not is_prompt) continue;
        xAOD::TruthParticle* muon = new xAOD::TruthParticle();
        muon->makePrivateStore(**mu_itr);
        v_muon.push_back(muon);
    } 
    sort(v_muon.begin(), v_muon.end(), byPt);

    // ---------------------------------
    //  MET
    // ---------------------------------
    const xAOD::MissingETContainer* met = 0;
    RETURN_CHECK( GetName(), m_event->retrieve( met, "MET_Truth" ));

    xAOD::MissingETContainer::const_iterator met_it = met->begin();


    // ----------------------------
    // Overlap Removal
    // ----------------------------
  
    // Jet Removal
    for(int i_jet = 0 ; i_jet < (int)v_jet.size() ; i_jet++ ){
      for(int i_el = 0 ; i_el < (int)v_electron.size() ; i_el++ ){
        if( v_jet.at(i_jet)->p4().DeltaR(v_electron.at(i_el)->p4()) < 0.2) { // Remove Jet element and move backward
          v_jet.erase( v_jet.begin()+i_jet);
          i_jet-- ;
          break;
        }
      }
    }
  
    // Electron Removal
    for(int i_el = 0 ; i_el < (int)v_electron.size() ; i_el++ ){
      for(int i_jet = 0 ; i_jet < (int)v_jet.size() ; i_jet++ ){
        if( v_electron.at(i_el)->p4().DeltaR(v_jet.at(i_jet)->p4()) < 0.4) { // Remove Electron element and move backward
          v_electron.erase( v_electron.begin()+i_el);
          i_el-- ;
          break;
        }
      }
    }
  
    // Muon Removal
    for(int i_mu = 0 ; i_mu < (int)v_muon.size() ; i_mu++ ){
      for(int i_jet = 0 ; i_jet < (int)v_jet.size() ; i_jet++ ){
        if( v_muon.at(i_mu)->p4().DeltaR(v_jet.at(i_jet)->p4()) < 0.4) { // Remove Muon element and move backward
          v_muon.erase( v_muon.begin()+i_mu);
          i_mu-- ;
          break;
        }
      }
    }

    // ----------------------------
    // Leptons
    // ----------------------------
    std::vector<xAOD::TruthParticle*> v_lepton;

    for(int i = 0; i < (int)v_electron.size(); i++)
        v_lepton.push_back(v_electron.at(i));
    for(int i = 0; i < (int)v_muon.size(); i++)
        v_lepton.push_back(v_muon.at(i));
    sort(v_lepton.begin(), v_lepton.end(), byPt);

    if(not (v_lepton.size() == 2) ) return kTRUE; // only store dilepton events

    // ---------------------------
    //  B-Jets
    // ---------------------------
    JftMultiplicity truthJetFlavorMultiplicity;
    std::vector<xAOD::Jet*> v_bjet, v_cjet, v_ljet, v_hjet, v_gjet;

    for(const auto jet : v_jet) {
        int flavor = abs( jet->auxdata<int>("PartonTruthLabelID"));
        if(flavor <4)
            v_ljet.push_back(jet);
        if(flavor==4)
            v_cjet.push_back(jet);
        if(flavor==5)
            v_bjet.push_back(jet);
        if(flavor==4 || flavor==5)
            v_hjet.push_back(jet);
        if(flavor==21)
            v_gjet.push_back(jet);

        if(truthJetFlavorMultiplicity.find(flavor)!=truthJetFlavorMultiplicity.end())
            truthJetFlavorMultiplicity[flavor] += 1;
        else
            truthJetFlavorMultiplicity[flavor] = 1;
    } // jet

    JftMultiplicity::iterator fj_it = truthJetFlavorMultiplicity.begin();
    JftMultiplicity::iterator fj_end = truthJetFlavorMultiplicity.end();
    for( ; fj_it != fj_end; ++fj_it ){
        cout << "flav: " << fj_it->first << "  mult: " << fj_it->second << endl;
    } 

    // ------------------------------
    //  Calculate lepton variables
    // ------------------------------
    m_br_lepton_n = (int)v_lepton.size();
    for(const auto l : v_lepton) {
        m_br_lepton_pt.push_back        (l->pt()*mev2gev);
        m_br_lepton_eta.push_back       (l->eta());
        m_br_lepton_phi.push_back       (l->phi());
        m_br_lepton_pdgId.push_back     (l->pdgId());
        m_br_lepton_charge.push_back    (l->charge());
        m_br_lepton_type.push_back      (l->auxdata<unsigned int>("classifierParticleType"));
        m_br_lepton_origin.push_back    (l->auxdata<unsigned int>("classifierParticleOrigin"));
    }

    m_br_isSF = (is_ee(v_lepton.at(0), v_lepton.at(1)) || is_mm(v_lepton.at(0), v_lepton.at(1)));
    m_br_isDF = !m_br_isSF;
    m_br_isEE = is_ee(v_lepton.at(0), v_lepton.at(1)) ? 1 : 0;
    m_br_isMM = is_mm(v_lepton.at(0), v_lepton.at(1)) ? 1 : 0;
    m_br_isEM = is_em(v_lepton.at(0), v_lepton.at(1)) ? 1 : 0;

    // ------------------------------
    //  Calculate jet variables
    // ------------------------------
    m_br_jet_n  = (int)v_jet.size();
    m_br_ljet_n = (int)v_ljet.size();
    m_br_cjet_n = (int)v_cjet.size();
    m_br_bjet_n = (int)v_bjet.size();
    m_br_hjet_n = (int)v_hjet.size();
    m_br_gjet_n = (int)v_gjet.size();

    for(const auto j : v_jet) {
        m_br_jet_pt.push_back           (j->pt()*mev2gev);
        m_br_jet_eta.push_back          (j->eta());
        m_br_jet_phi.push_back          (j->phi());
        m_br_jet_m.push_back            (j->m());
        m_br_jet_flav.push_back         (abs(j->auxdata<int>("PartonTruthLabelID")));
    }
    for(const auto j : v_ljet) {
        m_br_ljet_pt.push_back           (j->pt()*mev2gev);
        m_br_ljet_eta.push_back          (j->eta());
        m_br_ljet_phi.push_back          (j->phi());
        m_br_ljet_m.push_back            (j->m()*mev2gev);
        m_br_ljet_flav.push_back         (abs(j->auxdata<int>("PartonTruthLabelID")));
    }
    for(const auto j : v_cjet) {
        m_br_cjet_pt.push_back           (j->pt()*mev2gev);
        m_br_cjet_eta.push_back          (j->eta());
        m_br_cjet_phi.push_back          (j->phi());
        m_br_cjet_m.push_back            (j->m()*mev2gev);
        m_br_cjet_flav.push_back         (abs(j->auxdata<int>("PartonTruthLabelID")));
    }
    for(const auto j : v_bjet) {
        m_br_bjet_pt.push_back           (j->pt()*mev2gev);
        m_br_bjet_eta.push_back          (j->eta());
        m_br_bjet_phi.push_back          (j->phi());
        m_br_bjet_m.push_back            (j->m()*mev2gev);
        m_br_bjet_flav.push_back         (abs(j->auxdata<int>("PartonTruthLabelID")));
    }
    for(const auto j : v_gjet) {
        m_br_gjet_pt.push_back           (j->pt()*mev2gev);
        m_br_gjet_eta.push_back          (j->eta());
        m_br_gjet_phi.push_back          (j->phi());
        m_br_gjet_m.push_back            (j->m()*mev2gev);
        m_br_gjet_flav.push_back         (abs(j->auxdata<int>("PartonTruthLabelID")));
    }
    for(const auto j : v_hjet) {
        m_br_hjet_pt.push_back           (j->pt()*mev2gev);
        m_br_hjet_eta.push_back          (j->eta());
        m_br_hjet_phi.push_back          (j->phi());
        m_br_hjet_m.push_back            (j->m()*mev2gev);
        m_br_hjet_flav.push_back         (abs(j->auxdata<int>("PartonTruthLabelID")));
    }

    // -----------------------------
    //  Fill the output ntuple
    // -----------------------------

    m_outputTree->Fill();

    return kTRUE;
    
}
// ---------------------------------------------------------
void TruthSelector::Terminate()
{
    cout << "TruthSelector::Terminate" << endl;
    saveOutputTree();
}



        
