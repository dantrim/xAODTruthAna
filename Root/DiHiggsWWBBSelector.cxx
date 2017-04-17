#include "xAODTruthAna/DiHiggsWWBBSelector.h"


//std/stl
#include <iostream>
#include <cstdio> // printf
#include <sstream>
#include <vector>
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

// ana
#include "xAODTruthAna/xaod_utils.h"


//mystery ROOT macro
ClassImp(DiHiggsWWBBSelector)

std::string ana = "DiHiggsWWBBSelector";

struct IsPdg {
    int pdg;
    IsPdg(int p) : pdg(p) {}
    bool operator()(const xAOD::TruthParticle *l) { return (l && l->absPdgId()==pdg); }
};
IsPdg isElectron(11), isMuon(13), isTau(15);

bool is_ee(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
    { return isElectron(l0) and isElectron(l1); }

bool is_mm(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
    { return isMuon(l0) and isMuon(l1); }

bool is_em(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
    { return ((isElectron(l0) and isMuon(l1)) or (isMuon(l0) and isElectron(l1))); }

//////////////////////////////////////////////////////////////////////////////
DiHiggsWWBBSelector::DiHiggsWWBBSelector() :
    m_output_setup(false)
{
    cout << ana << endl;
}
//////////////////////////////////////////////////////////////////////////////
DiHiggsWWBBSelector::~DiHiggsWWBBSelector()
{
    cout << ana << "destructor"  << endl;
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsWWBBSelector::SlaveBegin(TTree* /*tree*/)
{
    string fn = ana + "::SlaveBegin    ";
    if(dbg() > 10) cout << fn << endl;
    timer()->Start();
}
//////////////////////////////////////////////////////////////////////////////
Bool_t DiHiggsWWBBSelector::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    n_evt_processed++;
    event()->getEntry(chain_entry);

    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo") );

    if(dbg() || chain_entry%1000==0) {
        cout << "*** Processing entry " << setw(6) << chain_entry
             << " run " << setw(6) << ei->runNumber()
             << " event " << setw(7) << ei->eventNumber() << " ***" << endl;
    }


    clear_variables();
    validation();


    output_tree()->Fill();



    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsWWBBSelector::Terminate()
{
    string fn = ana + "::Terminate    ";
    cout << fn << endl;
    timer()->Stop();
    cout << timer_summary() << endl;

    save_output();
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsWWBBSelector::initialize_output_tree(bool /*is_validation*/)
{
    string fn = ana + "::initialize_output_tree";
    if(dbg() > 10) cout << fn << endl; 
    //TFile* f = output_file();
    //TTree* t = output_tree();
    if(m_output_tree_file) delete m_output_tree_file;
    if(m_output_tree) delete m_output_tree;

    m_output_tree_file = new TFile(output_filename().c_str(), "recreate");
    m_output_tree = new TTree("superTruth", "superTruth");
    m_output_tree->SetAutoSave(10000000);
    m_output_tree->SetMaxTreeSize(3000000000u);

}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsWWBBSelector::initialize_tree_branches(bool is_validation)
{
    string fn = ana + "::initialize_tree_branches";
    if(!is_validation) {
        cout << fn << "    ERROR Tree branches undefined for non-validation!" << endl;
        exit(1);
    }
    else {
        TTree* tree = output_tree();
        if(!tree) { cout << fn << "    ERROR tree is null!" << endl; exit(1); }
        tree->Branch("isWWbb", &m_br_is_wwbb);
        tree->Branch("nLeptons", &m_br_n_leptons);
        tree->Branch("higgs_mass_w", &m_br_higgs_mass_w);
        tree->Branch("higgs_mass_b", &m_br_higgs_mass_b);

        tree->Branch("higgs_mass_lvlv", &m_br_higgs_mass_lvlv);
        tree->Branch("higgs_mass_bjets", &m_br_higgs_mass_bjets);
        tree->Branch("X_mass_lvlvb", &m_br_X_mass_lvlvb);

        tree->Branch("X_mass", &m_br_X_mass);
        tree->Branch("mW_0", &m_br_w0_mass);
        tree->Branch("mW_1", &m_br_w1_mass);

        tree->Branch("dileptonFlav", &m_br_dilepton_flavor);
        tree->Branch("l_pdgId", &m_br_l_pdgId);
        tree->Branch("dphi_ll", &m_br_dphi_ll);
        tree->Branch("pTll", &m_br_pTll);
        tree->Branch("l_pt", &m_br_l_pt);
        tree->Branch("l_eta", &m_br_l_eta);

        tree->Branch("dphi_bb", &m_br_dphi_bb);
        tree->Branch("b_pt", &m_br_b_pt);
        tree->Branch("b_eta", &m_br_b_eta);

        tree->Branch("dphi_ll_bb", &m_br_dphi_ll_bb);
    }

}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsWWBBSelector::clear_variables()
{
    m_br_is_wwbb = 0;
    m_br_n_leptons = -1;
    m_br_higgs_mass_w = -1;
    m_br_higgs_mass_b = -1;
    m_br_X_mass = -1;

    m_br_higgs_mass_lvlv = -1;
    m_br_higgs_mass_bjets = -1;
    m_br_X_mass_lvlvb = -1;

    m_br_w0_mass = -1;
    m_br_w1_mass = -1;

    m_br_dphi_ll = -5;
    m_br_pTll = -5;
    m_br_dilepton_flavor = -1;
    m_br_l_pdgId.clear();
    m_br_l_pt.clear();
    m_br_l_eta.clear();

    m_br_dphi_bb = -5;
    m_br_b_pt.clear();
    m_br_b_eta.clear();

    m_br_dphi_ll_bb = -5;
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsWWBBSelector::validation()
{
    if(!m_output_setup) {
        initialize_output_tree(true);
        initialize_tree_branches(true);
        m_output_setup = true;
    }
    

    //const xAOD::TruthVertexContainer* vertices = 0;
    auto vertices = Retrieve<xAOD::TruthVertexContainer>(event(), "TruthVertices");

    vector<const xAOD::TruthParticle*> higgses;
    vector<const xAOD::TruthParticle*> wbosons;
    vector<const xAOD::TruthParticle*> bquarks;

    bool is_ww_bb = false;
    

    // get the Higgs decay vertices
    int n_vert = 0;
    for(auto& vertex : *vertices) {
        if(vertex->nIncomingParticles() == 1 && vertex->nOutgoingParticles() == 2) { // h -> WW or h -> bb
            const xAOD::TruthParticle* mother = vertex->incomingParticle(0);
            if(mother->absPdgId() != 25) continue; // only higgs decays
            for(unsigned int i = 0; i < 2; i++) {
                const xAOD::TruthParticle *child = vertex->outgoingParticle(i);
                if(child->absPdgId()==5) bquarks.push_back(child);
                else if(child->isW()) wbosons.push_back(child);
            } // i
            higgses.push_back(mother);
        } // 1 -> 2
    } // vertex 

    if(higgses.size()==2 && wbosons.size()==2 && bquarks.size()==2) is_ww_bb = true;

    //cout << "validation    # higgs : " << higgses.size() << "  # W bosons : " << wbosons.size() << "  # b-quarks : " << bquarks.size() << "  WWbb ? " << (wbosons.size()==2 ? "YES" : "NO") << endl;

    if(is_ww_bb) {

        ///////////////////////////////////////////////////
        // get the masses of the higgses and X from each leg
        ///////////////////////////////////////////////////
        double hw_mass = (wbosons.at(0)->p4() + wbosons.at(1)->p4()).M() * mev2gev;
        double hb_mass = (bquarks.at(0)->p4() + bquarks.at(1)->p4()).M() * mev2gev;
        double x_mass = (wbosons.at(0)->p4() + wbosons.at(1)->p4()
                            + bquarks.at(0)->p4() + bquarks.at(1)->p4()).M() * mev2gev;

        ///////////////////////////////////////////////////
        // get W-bosons decay children
        ///////////////////////////////////////////////////
        vector<const xAOD::TruthParticle*> wleptons;
        vector<const xAOD::TruthParticle*> wneutrinos;
        const xAOD::TruthParticle *wleptons_arr[2] = { NULL };
        const xAOD::TruthParticle *wneutrinos_arr[2] = { NULL };
        for(unsigned int iw = 0; iw < 2; iw++) {
            bool self_decay = false;
            const xAOD::TruthParticle* wboson = wbosons.at(iw);
            do {
                for(unsigned int i = 0; i < wboson->nChildren(); i++) {
                    self_decay = false;
                    const xAOD::TruthParticle* child = wboson->child(i);
                    if(child->pdgId() == wboson->pdgId()) {
                        wboson = child;
                        self_decay = true;
                        break;
                    }
                    else if(child->isChLepton()) {
                        wleptons_arr[iw] = child;
                    }
                    else if(child->isNeutrino()) {
                        wneutrinos_arr[iw] = child;
                    }
                        
                } // i
            } while (self_decay);
        } // iw

        // only "light" leptons
        if(wleptons_arr[0] && (wleptons_arr[0]->absPdgId()!=15)) wleptons.push_back(wleptons_arr[0]);
        if(wleptons_arr[1] && (wleptons_arr[1]->absPdgId()!=15)) wleptons.push_back(wleptons_arr[1]);
        if(wneutrinos_arr[0]) wneutrinos.push_back(wneutrinos_arr[0]);
        if(wneutrinos_arr[1]) wneutrinos.push_back(wneutrinos_arr[1]);

        double w_mass_0 = -1;
        double w_mass_1 = -1;

        double h_mass_lvlv = -1;

        double dphi_ll = -5;

        double dphi_bb = bquarks.at(0)->p4().DeltaPhi(bquarks.at(1)->p4());

        double dphi_ll_bb = -5;

        int flavor = -1;
        if(wleptons.size()==2) {
            if(is_ee(wleptons.at(0), wleptons.at(1))) flavor = 0;
            else if(is_mm(wleptons.at(0), wleptons.at(1))) flavor = 1;
            else if(is_em(wleptons.at(0), wleptons.at(1))) flavor = 2;
            else { flavor = 3; } // tau stuff
            m_br_dilepton_flavor = flavor;
        } // == 2

        // dilepton
        if(wleptons.size()==2 && flavor>= 0 && flavor <3) { // only light leptons
            w_mass_0 = (wleptons.at(0)->p4() + wneutrinos.at(1)->p4()).M() * mev2gev;
            w_mass_1 = (wleptons.at(1)->p4() + wneutrinos.at(1)->p4()).M() * mev2gev;

            // higgs mass from leptons and neutrinos
            h_mass_lvlv = (wleptons.at(0)->p4() + wneutrinos.at(1)->p4()
                            + wleptons.at(1)->p4() + wneutrinos.at(1)->p4()).M() * mev2gev;

            ////////////////////////////////
            // dilepton
            ////////////////////////////////

            // angles between leptons
            dphi_ll = wleptons.at(0)->p4().DeltaPhi(wleptons.at(1)->p4());

            m_br_pTll = (wleptons.at(0)->p4() + wleptons.at(1)->p4()).Pt() * mev2gev;
            
            

            for(int i = 0; i < 2; i++) {
                m_br_l_pt.push_back(wleptons.at(i)->p4().Pt()*mev2gev);
                m_br_l_eta.push_back(wleptons.at(i)->p4().Eta());
                m_br_l_pdgId.push_back(wleptons.at(i)->absPdgId());
            }

            ///////////////////////////////
            // di-b system
            ///////////////////////////////

            // angles between dilepton and b-quark systems
            dphi_ll_bb = (wleptons.at(0)->p4() + wleptons.at(1)->p4()).DeltaPhi(
                            (bquarks.at(0)->p4() + bquarks.at(1)->p4()));

            for(int i = 0; i < 2; i++) {
                m_br_b_pt.push_back(bquarks.at(i)->p4().Pt()*mev2gev);
                m_br_b_eta.push_back(bquarks.at(i)->p4().Eta());
            }
        
        } // == 2 leptons

    


        ///////////////////////////////////////////////////
        // assign the variables to the branches
        ///////////////////////////////////////////////////
        m_br_is_wwbb = 1;
        m_br_n_leptons = wleptons.size();
        m_br_higgs_mass_w = hw_mass;
        m_br_higgs_mass_b = hb_mass;
        m_br_X_mass = x_mass;

        m_br_higgs_mass_lvlv = h_mass_lvlv;

        m_br_w0_mass = w_mass_0;
        m_br_w1_mass = w_mass_1;

        m_br_dphi_ll = dphi_ll;

        m_br_dphi_bb = dphi_bb;

        m_br_dphi_ll_bb = dphi_ll_bb;


    }

}
