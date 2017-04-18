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

typedef map<int, int> JetFlavMult; // first=flavor, second=mult

struct IsPdg {
    int pdg;
    IsPdg(int p) : pdg(p) {}
    bool operator()(const xAOD::TruthParticle *l) { return (l && l->absPdgId()==pdg); }
};
IsPdg isElectron(11), isMuon(13), isTau(15);

struct PtGreater {
    bool operator()(const xAOD::TruthParticle* a, const xAOD::TruthParticle* b)
        { return a->pt()>b->pt(); }
} byPt;

struct PtGreaterJet {
    bool operator()(const xAOD::Jet* a, const xAOD::Jet* b) { return a->pt()>b->pt(); }
} byPtJet;

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
        cout << "*** Processing entry " << setw(6) << chain_entry << " *** " << endl;
             //<< " run " << setw(6) << ei->runNumber()
             //<< " event " << setw(7) << ei->eventNumber() << " ***" << endl;
    }


    clear_variables();

    if(!m_output_setup) {
        initialize_output_tree(validation());
        initialize_tree_branches(validation());
        m_output_setup = true;
    }

    if(validation()) {
        validation_observables();
    }
    else {
        get_observables();
    }


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
        init_ana_tree();
    }
    else {
        init_validation_tree();
    }
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsWWBBSelector::init_ana_tree()
{
    string fn = ana + "::init_ana_tree";

    TTree* tree = output_tree();
    if(!tree) { cout << fn << "    ERROR tree is null!" << endl; exit(1); }

    tree->Branch("nLeptons", &m_br_n_leptons);
    tree->Branch("dileptonFlav", &m_br_dilepton_flavor);

    // leptons
    tree->Branch("l_pt", &m_br_l_pt);
    tree->Branch("l_eta", &m_br_l_eta);
    tree->Branch("l_pdgId", &m_br_l_pdgId);
    tree->Branch("l_q", &m_br_l_q);

    // dilepton
    tree->Branch("dphi_ll", &m_br_dphi_ll);
    tree->Branch("pTll", &m_br_pTll);
    tree->Branch("dRll", &m_br_dR_ll);
    tree->Branch("cosThetaB", &m_br_cosThetaB);
    tree->Branch("cosTheta1", &m_br_cosTheta1);
    tree->Branch("cosTheta2", &m_br_cosTheta2);
    tree->Branch("dphi_boost_ll", &m_br_dphi_boost_ll);

    // jets
    tree->Branch("n_jet", &m_br_jet_n);
    tree->Branch("n_ljet", &m_br_ljet_n);
    tree->Branch("n_bjet", &m_br_bjet_n);
    tree->Branch("n_cjet", &m_br_cjet_n);
    tree->Branch("n_hjet", &m_br_hjet_n);
    tree->Branch("n_gjet", &m_br_gjet_n);

    // jet kin
    tree->Branch("j_pt", &m_br_j_pt);
    tree->Branch("j_eta", &m_br_j_eta);

    // ljet kin
    tree->Branch("lj_pt", &m_br_lj_pt);
    tree->Branch("lj_eta", &m_br_lj_eta);

    // bjet kin
    tree->Branch("bj_pt", &m_br_bj_pt);
    tree->Branch("bj_eta", &m_br_bj_eta);
    tree->Branch("mass_bb", &m_br_mass_bb);

    // double b
    tree->Branch("dRbb", &m_br_dR_bb);
    tree->Branch("dphi_bb", &m_br_dphi_bb);

    // WWbb
    tree->Branch("dphi_ll_bb", &m_br_dphi_ll_bb);
    tree->Branch("dR_ll_bb", &m_br_dR_ll_bb);
    tree->Branch("dphi_WW_bb", &m_br_dphi_WW_bb);
    tree->Branch("mass_X", &m_br_mass_X);
    tree->Branch("mass_X_scaled", &m_br_mass_X_scaled);
    tree->Branch("CMboostMag", &m_br_CMboostMag);
    tree->Branch("dphi_boost_CM_ll_met", &m_br_dphi_boost_CM_ll_met);

    // met
    tree->Branch("met_et", &m_br_met_et);
    tree->Branch("met_phi", &m_br_met_phi);

    // lepton + met
    tree->Branch("dphi_met_ll", &m_br_dphi_met_ll);
    tree->Branch("mass_met_ll", &m_br_mass_met_ll);
    tree->Branch("met_pTll", &m_br_met_pTll);
    tree->Branch("dphi_boost_ll_met", &m_br_dphi_boost_ll_met);

    tree->Branch("HT2", &m_br_HT2);
    tree->Branch("HT2Ratio", &m_br_HT2Ratio);

    tree->Branch("HT2_noMET", &m_br_HT2_noMET);
    tree->Branch("HT2Ratio_noMET", &m_br_HT2Ratio_noMET);

}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsWWBBSelector::init_validation_tree()
{
    string fn = ana + "::init_validation_tree";

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
    tree->Branch("pTbb", &m_br_pTbb);
    tree->Branch("b_pt", &m_br_b_pt);
    tree->Branch("b_eta", &m_br_b_eta);
    
    tree->Branch("dphi_ll_bb", &m_br_dphi_ll_bb);
    tree->Branch("dR_bb", &m_br_dR_bb);

}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsWWBBSelector::clear_variables()
{

    // "validation"
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
    m_br_pTbb  = -1;
    m_br_b_pt.clear();
    m_br_b_eta.clear();

    m_br_dphi_ll_bb = -5;
    m_br_dR_bb = -10;

    // "non-validation"
    m_br_l_q.clear();
    m_br_jet_n = -1;
    m_br_ljet_n = -1;
    m_br_bjet_n = -1;
    m_br_cjet_n = -1;
    m_br_hjet_n = -1;
    m_br_gjet_n = -1;

    // jet kin
    m_br_j_pt.clear();
    m_br_j_eta.clear();

    // ljet kin
    m_br_lj_pt.clear();
    m_br_lj_eta.clear();

    // bjet kin
    m_br_bj_pt.clear();
    m_br_bj_eta.clear();
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsWWBBSelector::validation_observables()
{
    

    //const xAOD::TruthVertexContainer* vertices = 0;
    auto vertices = Retrieve<xAOD::TruthVertexContainer>(event(), "TruthVertices");

    vector<const xAOD::TruthParticle*> higgses;
    vector<const xAOD::TruthParticle*> wbosons;
    vector<const xAOD::TruthParticle*> bquarks;

    bool is_ww_bb = false;
    

    // get the Higgs decay vertices
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
        m_br_pTbb = (bquarks.at(0)->p4() + bquarks.at(1)->p4()).Pt() * mev2gev;

        double dphi_ll_bb = -5;

        const xAOD::JetContainer* jets = Retrieve<const xAOD::JetContainer>(event(), "AntiKt4TruthJets");
        std::vector<xAOD::Jet*> v_jets, v_bjets;

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
            // sort now that we don't need association
            sort(wleptons.begin(), wleptons.end(), byPt);

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

        
        } // == 2 leptons

        m_br_dR_bb = bquarks.at(0)->p4().DeltaR(bquarks.at(1)->p4());

        for(int i = 0; i < 2; i++) {
            m_br_b_pt.push_back(bquarks.at(i)->p4().Pt()*mev2gev);
            m_br_b_eta.push_back(bquarks.at(i)->p4().Eta());
        }

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

    // FILL THE OUTPUT TREE
    output_tree()->Fill();

}
//////////////////////////////////////////////////////////////////////////////
double DiHiggsWWBBSelector::cosThetaB(vector<xAOD::TruthParticle*> leptons)
{
    if(leptons.size()!=2) return -5;

    TLorentzVector lp, lm;
    for(int i = 0; i < 2; i++) {
        if(leptons.at(i)->charge() < 0) lm = leptons.at(i)->p4();
        else if(leptons.at(i)->charge() > 0) lp = leptons.at(i)->p4();
    }
    TLorentzVector ll = lp + lm;
    TVector3 boost = ll.BoostVector();
    lp.Boost(-boost);
    lm.Boost(-boost);

    double cos = tanh((lp.Eta()-lm.Eta())/2.);
    return cos;
}
//////////////////////////////////////////////////////////////////////////////
double DiHiggsWWBBSelector::cosTheta1(vector<xAOD::TruthParticle*> leptons)
{
    if(leptons.size()!=2) return -5;
    TLorentzVector lp, lm;
    for(int i = 0; i < 2; i++) {
        if(leptons.at(i)->charge() < 0) lm = leptons.at(i)->p4();
        else if(leptons.at(i)->charge() > 0) lp = leptons.at(i)->p4();
    }
    TLorentzVector ll = lp + lm;
    TVector3 boost = ll.BoostVector();
    lp.Boost(-boost);
    TVector3 lpboost = lp.BoostVector();

    double cos = boost.Unit().Dot(lpboost.Unit());
    return cos;

}
//////////////////////////////////////////////////////////////////////////////
double DiHiggsWWBBSelector::cosTheta2(vector<xAOD::TruthParticle*> leptons, TLorentzVector met)
{
    if(leptons.size()!=2) return -5;
    TLorentzVector lp, lm;
    for(int i = 0; i < 2; i++) {
        if(leptons.at(i)->charge() < 0) lm = leptons.at(i)->p4();
        else if(leptons.at(i)->charge() > 0) lp = leptons.at(i)->p4();
    }
    TLorentzVector ll = lp + lm;
    TVector3 boost = (ll + met).BoostVector();
    lp.Boost(-boost);
    lm.Boost(-boost);
    double cos = tanh((lp.Eta()-lm.Eta())/2.);
    return cos;
}
//////////////////////////////////////////////////////////////////////////////
double DiHiggsWWBBSelector::dphi_boost_ll(vector<xAOD::TruthParticle*> leptons)
{
    if(leptons.size()!=2) return -5;
    TLorentzVector lp, lm;
    for(int i = 0; i < 2; i++) {
        if(leptons.at(i)->charge() < 0) lm = leptons.at(i)->p4();
        else if(leptons.at(i)->charge() > 0) lp = leptons.at(i)->p4();
    }
    TLorentzVector ll = lp + lm;
    TVector3 boost = ll.BoostVector();
    lp.Boost(-boost);
    lm.Boost(-boost);
    double dphi = (lp + lm).Vect().DeltaPhi(boost);
    return dphi;

}
//////////////////////////////////////////////////////////////////////////////
double DiHiggsWWBBSelector::dphi_boost_ll_met(vector<xAOD::TruthParticle*> leptons, TLorentzVector met)
{
    if(leptons.size()!=2) return -5;
    TLorentzVector lp, lm;
    for(int i = 0; i < 2; i++) {
        if(leptons.at(i)->charge() < 0) lm = leptons.at(i)->p4();
        else if(leptons.at(i)->charge() > 0) lp = leptons.at(i)->p4();
    }
    TLorentzVector ll = lp + lm;
    TLorentzVector llmet = ll + met;
    TVector3 boost = llmet.BoostVector();
    lp.Boost(-boost);
    lm.Boost(-boost);
    double dphi = (lp + lm).Vect().DeltaPhi(boost);
    return dphi;
}
//////////////////////////////////////////////////////////////////////////////
double DiHiggsWWBBSelector::dphi_boost_CM_ll_met(vector<xAOD::TruthParticle*> leptons,
        TLorentzVector met, vector<xAOD::Jet*> bjets)
{
    TVector3 boost_to_cm = (leptons.at(0)->p4() + leptons.at(1)->p4() + met 
                + bjets.at(0)->p4() + bjets.at(1)->p4()).BoostVector();

    TLorentzVector l0 = leptons.at(0)->p4();
    TLorentzVector l1 = leptons.at(1)->p4();
    l0.Boost(-boost_to_cm);
    l1.Boost(-boost_to_cm);
    met.Boost(-boost_to_cm);

    TVector3 boost = (l0 + l1 + met).BoostVector();
    l0.Boost(-boost);
    l1.Boost(-boost);
    double dphi = (l0 + l1).Vect().DeltaPhi(boost);
    return dphi;

}
//////////////////////////////////////////////////////////////////////////////
Bool_t DiHiggsWWBBSelector::get_observables()
{
    string fn = ana + "::get_observables";

    ///////////////////////////////////////////////
    // collect objects
    ///////////////////////////////////////////////

    // jets
    const xAOD::JetContainer *jets = 0;
    RETURN_CHECK( GetName(), event()->retrieve( jets, "AntiKt4TruthJets" ));
    vector<xAOD::Jet*> v_jet;
    xAOD::JetContainer::const_iterator jet_itr = jets->begin();
    xAOD::JetContainer::const_iterator jet_end = jets->end();
    //for(const auto& j : *jets) {
    for( ; jet_itr != jet_end; ++jet_itr) {
        if( ((*jet_itr)->p4().Pt()*mev2gev < 20) ) continue;
        //if(j->p4().Pt()*mev2gev < 20) continue;
        xAOD::Jet* jet = new xAOD::Jet();
        jet->makePrivateStore(**jet_itr);
        v_jet.push_back(jet);
    }
    sort(v_jet.begin(), v_jet.end(), byPtJet);

    // electrons
    const xAOD::TruthParticleContainer* electrons = 0;
    RETURN_CHECK( GetName(), event()->retrieve( electrons, "TruthElectrons" ));
    vector<xAOD::TruthParticle*> v_electron;
    for(const auto& e : *electrons)  {
        if(e->pt()*mev2gev < 10) continue;
        if(std::abs(e->eta()) > 2.47) continue;
        bool is_prompt = (e->auxdata<unsigned int>( "classifierParticleType" ) == 2);
        if(!is_prompt) continue;
        xAOD::TruthParticle* electron = new xAOD::TruthParticle();
        electron->makePrivateStore(*e);
        v_electron.push_back(electron);
    }
    sort(v_electron.begin(), v_electron.end(), byPt);

    // muons
    const xAOD::TruthParticleContainer* muons = 0;
    RETURN_CHECK( GetName(), event()->retrieve( muons, "TruthMuons" ));
    vector<xAOD::TruthParticle*> v_muon;
    for(const auto& m : *muons) {
        if(m->pt()*mev2gev < 10) continue;
        if(std::abs(m->eta()) > 2.4) continue;
        bool is_prompt = (m->auxdata<unsigned int>( "classifierParticleType" ) == 6);
        if(!is_prompt) continue;
        xAOD::TruthParticle* muon = new xAOD::TruthParticle();
        muon->makePrivateStore(*m);
        v_muon.push_back(muon);
    }
    sort(v_muon.begin(), v_muon.end(), byPt);

    // MET
    const xAOD::MissingETContainer* met = 0;
    RETURN_CHECK( GetName(), event()->retrieve( met, "MET_Truth" ));

    ///////////////////////////////////////////////
    // overlap removal
    ///////////////////////////////////////////////

    float DR_J_E = 0.2; // remove j, keep e
    float DR_E_J = 0.4; // remove e, keep j
    float DR_M_J = 0.4; // remove m, keep j

    // jet removal
    for(int i_jet = 0; i_jet < (int)v_jet.size(); i_jet++) {
        for(int i_el = 0; i_el < (int)v_electron.size(); i_el++) {
            if(v_jet.at(i_jet)->p4().DeltaR(v_electron.at(i_el)->p4()) < DR_J_E) {
                if(dbg()>20) {
                    cout << fn << "    removing jet (" << v_jet.at(i_jet)->p4().Pt()*mev2gev
                               << ", " << v_jet.at(i_jet)->p4().Eta()
                               << ", " << v_jet.at(i_jet)->p4().Phi() << ") "
                               << "overlapping with electron ("
                               << v_electron.at(i_el)->p4().Pt()*mev2gev
                               << ", " << v_electron.at(i_el)->p4().Eta()
                               << ", " << v_electron.at(i_el)->p4().Phi() << ") : "
                               << v_jet.at(i_jet)->p4().DeltaR(v_electron.at(i_el)->p4()) << endl;
                }
                v_jet.erase(v_jet.begin()+i_jet);
                i_jet--;
                break;
            } // dR match
        } // i_el
    } // i_jet

    // electron removal
    for(int i_el = 0; i_el < (int)v_electron.size(); i_el++) {
        for(int i_jet = 0; i_jet < (int)v_jet.size(); i_jet++) {
            if(v_electron.at(i_el)->p4().DeltaR(v_jet.at(i_jet)->p4()) < DR_E_J) {
                if(dbg()>20) {
                    cout << fn << "    removing electron (" << v_electron.at(i_el)->p4().Pt()*mev2gev
                               << ", " << v_electron.at(i_el)->p4().Eta()
                               << ", " << v_electron.at(i_el)->p4().Phi() << ") "
                               << "overlapping with jet ("
                               << v_jet.at(i_jet)->p4().Pt()*mev2gev
                               << ", " << v_jet.at(i_jet)->p4().Eta()
                               << ", " << v_jet.at(i_jet)->p4().Phi() << ") : "
                               << v_jet.at(i_jet)->p4().DeltaR(v_electron.at(i_el)->p4()) << endl;
                }
                v_electron.erase(v_electron.begin() + i_el);
                i_el--;
                break;
            } // dR Match
        } // i_jet
    } // i_el

    // muon removal
    for(int i_mu = 0; i_mu < (int)v_muon.size(); i_mu++) {
        for(int i_jet = 0; i_jet < (int)v_jet.size(); i_jet++) {
            if(v_muon.at(i_mu)->p4().DeltaR(v_jet.at(i_jet)->p4()) < DR_M_J) {
                if(dbg()>20) {
                    cout << fn << "    removing muon (" << v_electron.at(i_mu)->p4().Pt()*mev2gev
                               << ", " << v_muon.at(i_mu)->p4().Eta()
                               << ", " << v_muon.at(i_mu)->p4().Phi() << ") "
                               << "overlapping with jet ("
                               << v_jet.at(i_jet)->p4().Pt()*mev2gev
                               << ", " << v_jet.at(i_jet)->p4().Eta()
                               << ", " << v_jet.at(i_jet)->p4().Phi() << ") : "
                               << v_jet.at(i_jet)->p4().DeltaR(v_muon.at(i_mu)->p4()) << endl;
                }
                v_muon.erase(v_muon.begin()+i_mu);
                i_mu--;
                break;
            } // dR match
        } // i_jet
    } // i_mu

    ///////////////////////////////////////////////
    // group leptons together
    ///////////////////////////////////////////////
    vector<xAOD::TruthParticle*> v_lepton;
    for(int i = 0; i < (int)v_electron.size(); i++)
        v_lepton.push_back(v_electron.at(i));
    for(int i = 0; i < (int)v_muon.size(); i++)
        v_lepton.push_back(v_muon.at(i));
    sort(v_lepton.begin(), v_lepton.end(), byPt);

    if( !(v_lepton.size()==2) ) return kTRUE; // require dilepton events

    ///////////////////////////////////////////////
    // jet flavors
    ///////////////////////////////////////////////
    JetFlavMult jet_flav_mult;
    vector<xAOD::Jet*> v_bjet, v_cjet, v_ljet, v_hjet, v_gjet;
    for(const auto j : v_jet) {
        int flavor = std::abs(j->auxdata<int>("PartonTruthLabelID"));
        if(flavor < 4)
            v_ljet.push_back(j);
        if(flavor==4)
            v_cjet.push_back(j);
        if(flavor==5)
            v_bjet.push_back(j);
        if(flavor==4 || flavor==5)
            v_hjet.push_back(j);
        if(flavor==21)
            v_gjet.push_back(j);

        if(jet_flav_mult.find(flavor)!=jet_flav_mult.end())
            jet_flav_mult[flavor] += 1;
        else
            jet_flav_mult[flavor] = 1;
    } // j

    if(dbg() > 5) {
        for(auto flav : jet_flav_mult)
            cout << "jet flavor: " << flav.first << "  multiplicity: " << flav.second << endl;
    }

    ///////////////////////////////////////////////
    // lepton kin
    ///////////////////////////////////////////////
    m_br_n_leptons = (int)v_lepton.size();
    for(const auto l : v_lepton) {
        m_br_l_pt.push_back     (l->pt()*mev2gev);
        m_br_l_eta.push_back    (l->eta());
        m_br_l_pdgId.push_back  (l->pdgId());
        m_br_l_q.push_back      (l->charge());
    } // l

    ///////////////////////////////////////////////
    // dilepton vars
    ///////////////////////////////////////////////
    m_br_pTll = (v_lepton.at(0)->p4() + v_lepton.at(1)->p4()).Pt()*mev2gev;
    m_br_dphi_ll = v_lepton.at(0)->p4().DeltaPhi(v_lepton.at(1)->p4());
    m_br_dR_ll = v_lepton.at(0)->p4().DeltaR(v_lepton.at(1)->p4());

    // ~decay angle
    m_br_cosThetaB = cosThetaB(v_lepton); // cosThetaB from 3-body
    m_br_cosTheta1 = cosTheta1(v_lepton);
    m_br_dphi_boost_ll = dphi_boost_ll(v_lepton);
    
    ///////////////////////////////////////////////
    // basic jet vars
    ///////////////////////////////////////////////
    m_br_jet_n = (int)v_jet.size();
    m_br_ljet_n = (int)v_ljet.size();
    m_br_cjet_n = (int)v_cjet.size();
    m_br_bjet_n = (int)v_bjet.size();
    m_br_hjet_n = (int)v_hjet.size();
    m_br_gjet_n = (int)v_gjet.size();
    
    // jet kin
    for(const auto j : v_jet) {
        m_br_j_pt.push_back     (j->pt()*mev2gev);
        m_br_j_eta.push_back    (j->eta());
    } // j

    // light-jet kin
    for(const auto j : v_ljet) {
        m_br_lj_pt.push_back    (j->pt()*mev2gev);
        m_br_lj_eta.push_back   (j->eta());
    }

    // b-jet kin
    for(const auto j : v_bjet) {
        m_br_bj_pt.push_back    (j->pt()*mev2gev);
        m_br_bj_eta.push_back   (j->eta());
    }

    if(v_bjet.size()>=2) {
        m_br_pTbb = (v_bjet.at(0)->p4() + v_bjet.at(1)->p4()).Pt()*mev2gev;
        m_br_dphi_bb = v_bjet.at(0)->p4().DeltaPhi(v_bjet.at(1)->p4());
        m_br_dR_bb = v_bjet.at(0)->p4().DeltaR(v_bjet.at(1)->p4());

        m_br_dphi_ll_bb = (v_lepton.at(0)->p4() + v_lepton.at(1)->p4()).
                            DeltaPhi((v_bjet.at(0)->p4() + v_bjet.at(1)->p4()));
        m_br_dR_ll_bb = (v_lepton.at(0)->p4() + v_lepton.at(1)->p4()).
                            DeltaR((v_bjet.at(0)->p4() + v_bjet.at(1)->p4()));
        m_br_mass_bb = (v_bjet.at(0)->p4() + v_bjet.at(1)->p4()).M()*mev2gev;
    }

    ///////////////////////////////////////////////
    // met related
    ///////////////////////////////////////////////

    TLorentzVector met_tlv;
    met_tlv.SetPxPyPzE((*met)["NonInt"]->mpx(),
                        (*met)["NonInt"]->mpy(),
                        0.,
                        (*met)["NonInt"]->met());
    m_br_met_et = met_tlv.Pt() * mev2gev;
    m_br_met_phi = met_tlv.Phi();


    // lepton + met
    m_br_dphi_met_ll = met_tlv.DeltaPhi((v_lepton.at(0)->p4() + v_lepton.at(1)->p4()));
    m_br_cosTheta2 = cosTheta2(v_lepton, met_tlv);
    m_br_mass_met_ll = (met_tlv + v_lepton.at(0)->p4() + v_lepton.at(1)->p4()).M()*mev2gev;
    m_br_met_pTll = (met_tlv + v_lepton.at(0)->p4() + v_lepton.at(1)->p4()).Pt()*mev2gev;
    m_br_dphi_boost_ll_met = dphi_boost_ll_met(v_lepton, met_tlv);

    if(v_bjet.size()>=2) {
        // lepton + b
        m_br_dphi_met_bb = met_tlv.DeltaPhi((v_bjet.at(0)->p4() + v_bjet.at(1)->p4()));

        // (lepton + met) + b
        m_br_dphi_WW_bb = (met_tlv + v_lepton.at(0)->p4() + v_lepton.at(1)->p4()).
                            DeltaPhi((v_bjet.at(0)->p4() + v_bjet.at(1)->p4()));

        m_br_mass_X = (met_tlv + v_lepton.at(0)->p4() + v_lepton.at(1)->p4()
                            + v_bjet.at(0)->p4() + v_bjet.at(1)->p4()).M()*mev2gev;

        TLorentzVector bbsyst = (v_bjet.at(0)->p4() + v_bjet.at(1)->p4());
        TLorentzVector bb;
        double scaling = 125.09 / (bbsyst.M() * mev2gev);
        bb.SetPtEtaPhiE(bbsyst.Pt() * scaling, bbsyst.Eta(), bbsyst.Phi(), bbsyst.E() * scaling);
        m_br_mass_X_scaled = (met_tlv + v_lepton.at(0)->p4() + v_lepton.at(1)->p4() + bb).M() * mev2gev;

        m_br_CMboostMag = (met_tlv + v_lepton.at(0)->p4() + v_lepton.at(1)->p4()
                            + v_bjet.at(0)->p4() + v_bjet.at(1)->p4()).BoostVector().Mag();

        m_br_dphi_boost_CM_ll_met = dphi_boost_CM_ll_met(v_lepton, met_tlv, v_bjet);

        m_br_HT2 = (v_bjet.at(0)->p4() + v_bjet.at(1)->p4()).Vect().Perp() * mev2gev
                        + (v_lepton.at(0)->p4() + v_lepton.at(1)->p4() + met_tlv).Vect().Perp() * mev2gev;
        double denom = v_bjet.at(0)->p4().Vect().Perp() + v_bjet.at(1)->p4().Vect().Perp()
                + v_lepton.at(0)->p4().Vect().Perp() + v_lepton.at(1)->p4().Vect().Perp()
                + met_tlv.Vect().Perp();
        denom = denom * mev2gev;

        //double denom = (v_bjet.at(0)->p4() + v_bjet.at(1)->p4() + v_lepton.at(0)->p4()
        //            + v_lepton.at(1)->p4() + met_tlv).Vect().Perp() * mev2gev;
        m_br_HT2Ratio = m_br_HT2 / denom;

        m_br_HT2_noMET = (v_bjet.at(0)->p4() + v_bjet.at(1)->p4()).Vect().Perp() * mev2gev
                        + (v_lepton.at(0)->p4() + v_lepton.at(1)->p4()).Vect().Perp() * mev2gev;
        denom = v_bjet.at(0)->p4().Vect().Perp() + v_bjet.at(1)->p4().Vect().Perp()
                + v_lepton.at(0)->p4().Vect().Perp() + v_lepton.at(1)->p4().Vect().Perp();
        denom = denom * mev2gev;

        m_br_HT2Ratio_noMET = m_br_HT2_noMET / denom;

    }


    // FILL THE OUTPUT TREE
    output_tree()->Fill();
    return kTRUE;
}
