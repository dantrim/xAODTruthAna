#include "xAODTruthAna/TruthSelectorBase.h"

//std/stl
#include <iostream>
#include <cstdio> // printf
#include <sstream>
using namespace std;

//ROOT
#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TNamed.h"
#include "TString.h"
#include "TVector3.h"

//xAOD/EDM
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODEventInfo/EventInfo.h"

#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

//mystery ROOT macro
ClassImp(TruthSelectorBase)


//////////////////////////////////////////////////////////////////////////////
TruthSelectorBase::TruthSelectorBase() :
    n_evt_processed(0),
    n_evt_stored(0),
    m_dbg(0),
    m_validation(false),
    m_output_tree_file(NULL),
    m_output_tree(NULL),
    m_tree(NULL),
    m_event(new xAOD::TEvent(xAOD::TEvent::kClassAccess) ),
    m_input_samplename(""),
    m_output_filename("")
{
}
//////////////////////////////////////////////////////////////////////////////
TruthSelectorBase::~TruthSelectorBase()
{
    delete m_event;
}
//////////////////////////////////////////////////////////////////////////////
void TruthSelectorBase::Init(TTree *tree)
{
    string fn = "TruthSelectorBase::Init    ";
    cout << fn << "Input sample : " << m_input_samplename << endl;
    if(tree) {
        cout << fn << "Found tree with " << tree->GetEntries() << " entries" << endl;
    }
    else {
        cout << fn << "ERROR Input tree is null" << endl;
        exit(1);
    }
    m_tree = tree;
    return;
}
//////////////////////////////////////////////////////////////////////////////
Bool_t TruthSelectorBase::Notify()
{
    string fn = "TruthSelectorBase::Notify    ";
    if(!m_tree) {
        cout << fn << "ERROR Tree is null!" << endl;
        return kTRUE;
    }
    cout << fn << "Opening a new file" << endl;
    TFile* ifile = 0;
    TChain* chain = dynamic_cast<TChain*>(m_tree);
    if(chain) {
        ifile = chain->GetFile();
        cout << fn << "chain : " << ifile->GetName() << endl;
    }
    else {
        ifile = m_tree->GetCurrentFile();
        cout << fn << "not chain : " << ifile->GetName() << endl;
    }
    RETURN_CHECK(GetName(), m_event->readFrom(ifile, kFALSE) );
    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void TruthSelectorBase::SlaveBegin(TTree* /*tree*/)
{
    string fn = "TruthSelectorBase::SlaveBegin    ";
    if(m_dbg > 10) cout << fn << endl;
    m_timer.Start();
}
//////////////////////////////////////////////////////////////////////////////
void TruthSelectorBase::Begin(TTree* /*tree*/)
{
    string fn = "TruthSelectorBase::Begin    ";
    cout << fn << endl;
    return;
}
//////////////////////////////////////////////////////////////////////////////
Bool_t TruthSelectorBase::Process(Long64_t entry)
{
    string fn = "TruthSelectorBase::Process    ";
    if(m_dbg > 10) cout << fn << endl;
    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void TruthSelectorBase::Terminate()
{
    string fn = "TruthSelectorBase::Terminate    ";
    if(m_dbg > 10) cout << fn << endl;
    m_timer.Stop();
    cout << timer_summary() << endl;
}
//////////////////////////////////////////////////////////////////////////////
void TruthSelectorBase::save_output()
{
    string fn = "TruthSelectorBase::save_output";
    cout << fn << endl;

    m_output_tree_file = m_output_tree->GetCurrentFile();
    if(m_output_tree_file) {
        TDirectory* current_directory = gROOT->CurrentDirectory();
        m_output_tree_file->cd();
        TNamed n("inputSample", m_input_samplename.c_str());
        n.Write();
        current_directory->cd();
    }
    m_output_tree->Write(0, TObject::kOverwrite);
    cout << fn << "    Output tree saved to " << m_output_tree_file->GetName() << endl;
    m_output_tree_file->Close();
}
//////////////////////////////////////////////////////////////////////////////
std::string TruthSelectorBase::timer_summary()
{
    double real_time = m_timer.RealTime();
    double cpu_time = m_timer.CpuTime();
    int hours = int(real_time / 3600);
    real_time -= hours * 3600;
    int min = int(real_time / 60);
    real_time -= min * 60;
    int sec = int(real_time);

    float speed = n_evt_processed / m_timer.RealTime() / 1000;
    TString line1; line1.Form("Real %d:%02d:%02d, CPU %.3f", hours, min, sec, cpu_time);
    TString line2; line2.Form("[kHz]: %2.3f", speed);
    ostringstream oss;
    oss << "-------------------------------------------\n"
        << " Number of events processed : " << n_evt_processed << endl
        << " Number of events saved     : " << n_evt_stored << endl
        << "\t Analysis time      : " << line1 << endl
        << "\t Analysis speed     : " << line2 << endl
        << "-------------------------------------------" << endl;
    return oss.str();
}
