#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <xAODTruthAna/TruthSelector.h>

// STD include(s):
#include <stdio.h>      // sscanf

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"

// EDM include(s):
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"


#include <algorithm>
#include <iostream>
#include <cstdio> // printf
using namespace std;

#include "TFile.h"

const double mev2gev = 1.0e-3;

#define EL_RETURN_CHECK( CONTEXT, EXP )                  \
    do {                                                 \
        if( !EXP.isSuccess() ) {                         \
            Error( CONTEXT,                              \
                XAOD_MESSAGE( "failed to execute: %s" ), \
                #EXP );                                  \
            return EL::StatusCode::FAILURE;              \
        }                                                \
    } while ( false )                                    \

// this is needed to distribute the algorithm to the workers
ClassImp(TruthSelector)

TruthSelector :: TruthSelector() :
    verbose(false)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}

EL::StatusCode TruthSelector :: setupJob (EL::Job& job)
{
    // Here you put code that sets up the job on the submission object
    // so that it is ready to work with your algorithm, e.g. you can
    // request the D3PDReader service or add output files.  Any code you
    // put here could instead also go into the submission script.  The
    // sole advantage of putting it here is that it gets automatically
    // activated/deactivated when you add/remove the algorithm from your
    // job, which may or may not be of value to you.

    // let's initialize the algorithm to use the xAODRootAccess package
    job.useXAOD();
    EL_RETURN_CHECK( "setupJob()", xAOD::Init() );
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode TruthSelector :: histInitialize ()
{
    // Here you do everything that needs to be done at the very
    // beginning on each worker node, e.g. create histograms and output
    // trees.  This method gets called before any input files are
    // connected.

//    TFile *outputFile = wk()->getOutputFile (*m_outFileName);
//    m_outTree = new TTree("superTruth","superTruth");
//    m_outTree->SetDirectory(outputFile);
//
//    m_outTree->Branch("runNumber",
//                            &m_br_runNumber                         );
//    m_outTree->Branch("eventNumber",
//                            &m_br_eventNumber                       );
//    m_outTree->Branch("eventweight",
//                            &m_br_eventWeight                       );
//    m_outTree->Branch("l_pt", /* lepton pTs */
//                            &m_br_lpt                               );

    return EL::StatusCode::SUCCESS;

}

EL::StatusCode TruthSelector :: fileExecute ()
{
    // Here you do everything that needs to be done exactly once for every
    // single file, e.g. collect a list of all lumi-blocks processed
    return EL::StatusCode::SUCCESS;
}


EL::StatusCode TruthSelector :: initialize ()
{
    // Here you do everything that you need to do after the first input
    // file has been connected and before the first event is processed,
    // e.g. create additional histograms based on which variables are
    // available in the input files.  You can also create all of your
    // histograms and trees in here, but be aware that this method
    // doesn't get called if no events are processed.  So any objects
    // you create here won't be available in the output if you have no
    // input events.

//    xAOD::TEvent* event = wk()->xaodEvent();
//
//    Info("initialize()", "Number of events = %lli", event->getEntries() ); // long long int

    return EL::StatusCode::SUCCESS;
}
/*
template <class Container>
double sumPt(Container &container)
{
    double toPt = 0.0;
    for(const auto *e : container)
        toPt += (e->pt() * mev2gev);
    return toPt;
}

typedef map<int, int> JftMultiplicity_t; // first = flavor, second = multiplicity

struct isPdg {
    int pdg;
    isPdg(int p) : pdg(p) {}
    bool operator()(const xAOD::TruthParticle* l) { return (l && l->absPdgId()==pdg);}
};
isPdg isElectron(11), isMuon(13), isTau(15);

struct PtGreater {
    bool operator()(const xAOD::TruthParticle* a, const xAOD::TruthParticle* b) { return a->pt()>b->pt(); } 
} byPt;

bool is_ee(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1) { return isElectron(l0) and isElectron(l1); }
bool is_em(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1) { return ((isElectron(l0) and isMuon(l1)) or (isMuon(l0) and isElectron(l1))); }
bool is_mm(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1) { return isMuon(l0) and isMuon(l1); }
*/
EL::StatusCode TruthSelector :: execute()
{
    // Here you do everything that needs to be done on every single
    // events, e.g. read input variables, apply cuts, and fill
    // histograms and trees.  This is where most of your actual analysis
    // code will go.
/*
    xAOD::TEvent* event = wk()->xaodEvent();
    double eventWeight = 1.0;


    // --------------------------------------
    //  Event Information
    // --------------------------------------
    const xAOD::EventInfo* eventinfo = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( eventinfo, "EventInfo" ));

    double runNumber = eventinfo->runNumber();
    double eventNumber = eventinfo->eventNumber();
    eventWeight *= eventinfo->mcEventWeight();
    

    const xAOD::TruthEventContainer* truthEvents = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( truthEvents, "TruthEvents" ));

    // --------------------------------------
    //  Electrons
    // --------------------------------------
    const xAOD::TruthParticleContainer* electrons = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( electrons, "TruthElectrons" ));

    std::vector<xAOD::TruthParticle*> v_electron;

    xAOD::TruthParticleContainer::const_iterator electron_itr = electrons->begin();
    xAOD::TruthParticleContainer::const_iterator electron_end = electrons->end();
    for( ; electron_itr != electron_end; ++electron_itr ) {
        if( ((*electron_itr)->pt() * mev2gev) < 10 ) continue;
        if( fabs((*electron_itr)->eta()) > 2.47 ) continue;
        bool is_prompt = ((*electron_itr)->auxdata<unsigned int>("classiferParticleType") == 2);
        if(not is_prompt) continue;
        xAOD::TruthParticle* electron = new xAOD::TruthParticle();
        electron->makePrivateStore( **electron_itr );
        v_electron.push_back(electron);
    }
    sort(v_electron.begin(), v_electron.end(), byPt);

    // --------------------------------------
    //  Muons
    // --------------------------------------
    const xAOD::TruthParticleContainer* muons = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( muons, "TruthMuons" ));

    std::vector<xAOD::TruthParticle*> v_muon;

    xAOD::TruthParticleContainer::const_iterator muon_itr = muons->begin();
    xAOD::TruthParticleContainer::const_iterator muon_end = muons->end();
    for( ; muon_itr != muon_end; ++muon_itr ) {
        if( ((*muon_itr)->pt() * mev2gev) < 10 ) continue;
        if( fabs((*muon_itr)->eta()) > 2.4 ) continue;
        bool is_prompt = ((*muon_itr)->auxdata<unsigned int>("classifierParticleType") == 6);
        if(not is_prompt) continue;
        xAOD::TruthParticle* muon = new xAOD::TruthParticle();
        muon->makePrivateStore(**muon_itr);
        v_muon.push_back(muon);
    }
    sort(v_muon.begin(), v_muon.end(), byPt);

    // --------------------------------------
    //  Jets
    // --------------------------------------
    const xAOD::JetContainer* jets = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( jets, "AntiKt4TruthJets" ));

    std::vector<xAOD::Jet*> v_jet;

    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();
    for( ; jet_itr != jet_end; ++jet_itr) {
        if(( (*jet_itr)->pt() * mev2gev) < 20) continue;
        xAOD::Jet* jet = new xAOD::Jet();
        jet->makePrivateStore( **jet_itr );
        v_jet.push_back(jet);
    }

    // --------------------------------------
    //  Overlap
    // --------------------------------------

    // remove jets OR with electrons at dR < 0.2
    for(int i_jet = 0 ; i_jet < (int)v_jet.size() ; i_jet++) {
        for(int i_el = 0 ; i_el < (int)v_electron.size() ; i_el++) {
            if( v_jet.at(i_jet)->p4().DeltaR(v_electron.at(i_el)->p4()) < 0.2) {
                // remove jet element and move backward
                v_jet.erase( v_jet.begin() + i_jet );
                i_jet--;
                break;
            }
        } // i_el
    } // i_jet

    // remove electrons OR with jets at dR < 0.4
    for(int i_el = 0 ; i_el < (int)v_electron.size() ; i_el++) {
        for(int i_jet = 0 ; i_jet < (int)v_jet.size() ; i_jet++) {
            if( v_electron.at(i_el)->p4().DeltaR(v_jet.at(i_jet)->p4()) < 0.4 ) {
                // remove electron element and move backward
                v_electron.erase( v_electron.begin() + i_el );
                i_el--;
                break;
            }
        } // i_jet
    } // i_el

    // remove muons OR with jets at dR < 0.4
    for(int i_mu = 0 ; i_mu < (int)v_muon.size() ; i_mu++) {
        for(int i_jet = 0 ; i_jet < (int)v_jet.size() ; i_jet++) {
            if( v_muon.at(i_mu)->p4().DeltaR(v_jet.at(i_jet)->p4()) < 0.4) {
                // remove muon element and move backward
                v_muon.erase( v_muon.begin() + i_mu );
                i_mu--;
                break;
            }
        } // i_jet
    } // i_mu 

    // --------------------------------------
    //  Leptons
    // --------------------------------------
    std::vector<xAOD::TruthParticle*> v_leptons;

    for(unsigned int ie = 0; ie < v_electron.size(); ie++) {
        v_leptons.push_back(v_electron[ie]);
    }
    for(unsigned int im = 0; im < v_muon.size(); im++) {
        v_leptons.push_back(v_muon[im]);
    }
    sort(v_leptons.begin(), v_leptons.end(), byPt);

    // --------------------------------------
    //  MET
    // --------------------------------------
    const xAOD::MissingETContainer* met = 0;
    EL_RETURN_CHECK("execute()", event->retrieve( met, "MET_Truth" ));

    xAOD::MissingETContainer::const_iterator met_it = met->begin();

    JftMultiplicity_t truthJetFlavorMultiplicity;
    std::vector<xAOD::Jet*> v_bjet;

    for(const auto jet : v_jet) {
        int flavor = abs( jet->auxdata<int>("PartonTruthLabelID") );
        if(flavor == 5)
            v_bjet.push_back(jet);

        if(truthJetFlavorMultiplicity.find(flavor) != truthJetFlavorMultiplicity.end())
            truthJetFlavorMultiplicity[flavor] += 1;
        else
            truthJetFlavorMultiplicity[flavor] = 1;
    }

    ////////////////////////////////////////////////////
    // Fill Tree
    ////////////////////////////////////////////////////

    // Event vars
    m_br_runNumber      = runNumber;
    m_br_eventNumber    = eventNumber;
    m_br_eventWeight    = eventWeight;

    // Lepton kinematics
    for(const auto l : v_leptons) { m_br_lpt.push_back (l->pt() * mev2gev); }
    for(const auto l : v_leptons) { m_br_leta.push_back(l->eta()); }
    for(const auto l : v_leptons) { m_br_lphi.push_back(l->phi()); }


    // Fill the output tree
    m_outTree->Fill();
*/
    return EL::StatusCode::SUCCESS;
    
    
}

EL::StatusCode TruthSelector :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TruthSelector :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode TruthSelector :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}
