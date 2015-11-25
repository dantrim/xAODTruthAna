#include "xAODTruthAna/TruthSelector.h"

//std/stl
#include <iostream>
#include <cstdio> // printf
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
    cout << "TruthSelector::initializeOutputTree" << endl;
    m_outTreeFile = new TFile(m_outputFileName.c_str(), "recreate");
    m_outputTree = new TTree("superTruth", "superTruth");
    m_outputTree->SetAutoSave(10000000);
    m_outputTree->SetMaxTreeSize(3000000000u);
}
// ---------------------------------------------------------
void TruthSelector::initializeTreeBranches()
{
    cout << "TruthSelector::initializeTreeBranches" << endl;
    m_outputTree->Branch("runNumber",
                              &m_br_runNumber             );
    m_outputTree->Branch("eventNumber",
                              &m_br_eventNumber           );
    m_outputTree->Branch("eventweight",
                              &m_br_eventWeight           );
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
// ---------------------------------------------------------
Bool_t TruthSelector::Process(Long64_t entry)
{
    if(m_dbg > 10) cout << "TruthSelector::Process" << endl;
    static Long64_t chainEntry = -1;
    chainEntry++;

    m_event->getEntry(chainEntry);

    // ---------------------------------
    //  Event information
    // ---------------------------------
    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), m_event->retrieve(ei, "EventInfo") );

    unsigned long long runNumber = ei->runNumber();
    unsigned long long eventNumber = ei->eventNumber();
    double eventWeight = ei->mcEventWeight();

    if(m_dbg || chainEntry%5000==0){
        cout << "*************************************************************" << endl;
        cout << "**** Processing entry " << setw(6) << chainEntry
             << "  run " << setw(6) << runNumber
             << "  event " << setw(7) << eventNumber << " ****" << endl;
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


    m_br_runNumber = runNumber;
    m_br_eventNumber = eventNumber;
    m_br_eventWeight = eventWeight;
    m_outputTree->Fill();

    return kTRUE;
    
}
// ---------------------------------------------------------
void TruthSelector::Terminate()
{
    cout << "TruthSelector::Terminate" << endl;
    saveOutputTree();
}



        
