#include "xAODTruthAna/DiHiggsResonance.h"

//std/stl
#include <iostream>
#include <cstdio> // printf
#include <sstream>
#include <vector>
#include <math.h> // sqrt, cos
#include <fstream>
using namespace std;


//ROOT
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"

//xAOD/EDM
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthEventContainer.h"

//ANA
#include "xAODTruthAna/xaod_utils.h"
#include "xAODTruthAna/MT2_ROOT.h"

//ClassImp(DiHiggsResonance)
//

namespace res {

std::string ANA = "DiHiggsResonance";

struct is_pdg {
    int pdg;
    is_pdg(int p) : pdg(p) {}
    bool operator()(const xAOD::TruthParticle *l) { return (l && l->absPdgId()==pdg); }
};
is_pdg is_ele(11);
is_pdg is_muo(13);
is_pdg is_tau(15);

struct pt_greater {
    bool operator()(const xAOD::TruthParticle* a, const xAOD::TruthParticle* b)
            { return a->pt() > b->pt(); }
} by_pt;

struct pt_greaterJet {
    bool operator()(const xAOD::Jet* a, const xAOD::Jet* b) { return a->pt() > b->pt(); }
} by_ptJet;

bool is_ee(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
    { return is_ele(l0) && is_ele(l1); }
bool is_mm(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
    { return is_muo(l0) && is_muo(l1); }
bool is_em(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
    { return ((is_ele(l0) && is_muo(l1)) || (is_muo(l0) && is_ele(l1))); }

//////////////////////////////////////////////////////////////////////////////
DiHiggsResonance::DiHiggsResonance() :
    m_weight(1.0),
    m_outfile(nullptr),
    file_setup(false),
    total_sumw(0.0),
    m_xsec(0.0),
    m_sumw(0.0),
    m_lumi(100.)
{
    cout << ANA << endl;


}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsResonance::set_x_mass(int mass)
{
    stringstream hypo;
    if(mass>0) hypo << mass << " GeV resonance";
    else { hypo << "hh non-resonant"; }

    cout << ANA << "    signal hypothesis : " << hypo.str() << endl;

    m_x_mass = mass;
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsResonance::set_lumi(float lumi)
{
    cout << ANA << "    Setting luminosity assumption to : " << lumi << " fb-1" << endl;
    m_lumi = lumi;
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsResonance::SlaveBegin(TTree* /*tree*/)
{
    timer()->Start();

}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsResonance::Terminate()
{
    timer()->Stop();
    cout << timer_summary() << endl;


    m_outfile->cd();
    h_mt1_scaled->Write();
    h_mt1->Write();
    m_outfile->Write();
    m_outfile->Close();

    cout << "DiHiggsResonance::Terminate    Histos stored to file: " << m_outfile->GetName() << endl;

    //m_outfile->Write();
    //m_outfile->Close();

}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsResonance::initialize_sumw_map()
{
    sumw_map.clear();
    string filename = "dihiggs_sumw_list.txt";
    std::ifstream ifs(filename.c_str());
    int xmass;
    double sumw;
    while(ifs >> xmass >> sumw) {
        sumw_map[xmass] = sumw; 
    }
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsResonance::initialize_xsec_map()
{
    dsid_map.clear();
    string filename = "dihiggs_dsid_map.txt";
    std::ifstream dfs(filename.c_str());
    int d;
    int x;
    while(dfs >> d >> x) {
        dsid_map[d] = x;
    }


    xsec_map.clear();
    filename = "dihiggs_xsec.txt";
    std::ifstream ifs(filename.c_str());
    int dsid;
    double xsec;
    while(ifs >> dsid >> xsec) {
        xsec_map[dsid_map[dsid]] = xsec;
    }
    for(auto x : xsec_map) {
        cout << "X = " << x.first << "  xsec = " << x.second << endl;
    }
    
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsResonance::initialize_histos()
{
    float low = 0.0;
    float high = 1.6*x_mass();
    int n_bins = 40;

    stringstream name;
    stringstream title;

    name << "h_mt1_" << x_mass();
    title << "MT_1;MT_1 [GeV];Entries";
    h_mt1 = new TH1F(name.str().c_str(), title.str().c_str(), n_bins, low, high);

    name.str("");
    title.str("");

    name << "h_mt1_scaled_" << x_mass();
    title << "MT_1_scaled;MT_1_scaled [GeV];Entries";
    h_mt1_scaled = new TH1F(name.str().c_str(), title.str().c_str(), n_bins, low, high);


}
//////////////////////////////////////////////////////////////////////////////
Bool_t DiHiggsResonance::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    n_evt_processed++;
    event()->getEntry(chain_entry);

    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo") );

    string fn = "DiHiggsResonance::Process    ";
    if(!file_setup) {
        stringstream fname;
        fname << "dihiggs_resonance_";
        string x = "hh";
        if(x_mass() > 0) {
            fname << "X" << x_mass();
        }
        else {
            fname << x;
        }
        fname << ".root";
            
        m_outfile = new TFile(fname.str().c_str(), "RECREATE");
        m_outfile->cd();
        file_setup = true;

        initialize_histos();

        initialize_sumw_map();
        initialize_xsec_map();
    }

    if(dbg() || chain_entry%5000==0) {
        cout << fn << " **** Processing entry " << setw(6) << chain_entry << " *** " << endl;
    }

    // get the MC event weight
    m_weight = ei->mcEventWeight();
    total_sumw += m_weight;

    m_xsec = xsec_map[x_mass()];
    m_sumw = sumw_map[x_mass()];

    /////////////////////////////////////////////
    //  get physics objects
    /////////////////////////////////////////////

    // jets
    vector<xAOD::Jet*> jets;
    const xAOD::JetContainer* xjets = 0;
    RETURN_CHECK( GetName(), event()->retrieve( xjets, "AntiKt4TruthJets" ));
    xAOD::JetContainer::const_iterator jet_itr = xjets->begin();
    xAOD::JetContainer::const_iterator jet_end = xjets->end();
    for( ; jet_itr != jet_end; ++jet_itr) {
        if(!( (*jet_itr)->p4().Pt() > 20.) ) continue;
        xAOD::Jet* j = new xAOD::Jet();
        j->makePrivateStore(**jet_itr);
        jets.push_back(j);
    }
    sort(jets.begin(), jets.end(), by_ptJet);

    // electrons
    vector<xAOD::TruthParticle*> electrons;
    const xAOD::TruthParticleContainer* xelectrons = 0;
    RETURN_CHECK( GetName(), event()->retrieve( xelectrons, "TruthElectrons" ));
    for(const auto& e : *xelectrons) {
        if( !(e->pt() > 10) ) continue;
        if(std::abs(e->eta()) > 2.47) continue;
        bool is_prompt = (e->auxdata<unsigned int>( "classifierParticleType" ) == 2);
        if(!is_prompt) continue;
        xAOD::TruthParticle* ele = new xAOD::TruthParticle();
        ele->makePrivateStore(*e);
        electrons.push_back(ele);
    }
    sort(electrons.begin(), electrons.end(), by_pt);

    // muons
    vector<xAOD::TruthParticle*> muons;
    const xAOD::TruthParticleContainer* xmuons = 0;
    RETURN_CHECK( GetName(), event()->retrieve( xmuons, "TruthMuons" ));
    for(const auto& m : *xmuons) {
        if( !(m->pt() > 10.) ) continue;
        if(std::abs(m->eta()) > 2.4) continue;
        bool is_prompt = (m->auxdata<unsigned int>( "classifierParticleType" ) == 6);
        if(!is_prompt) continue;
        xAOD::TruthParticle* mu = new xAOD::TruthParticle();
        mu->makePrivateStore(*m);
        muons.push_back(mu);
    }
    sort(muons.begin(), muons.end(), by_pt);

    // MET
    const xAOD::MissingETContainer* met = 0;
    RETURN_CHECK( GetName(), event()->retrieve( met, "MET_Truth" ));

    ///////////////////////////////////////////////
    // Overlap Removal
    ///////////////////////////////////////////////

    // TODO implement OR procedure closer to what we do in reco level

    float DR_J_E = 0.2; // remove j, keep e
    float DR_E_J = 0.4; // remove e, keep j
    float DR_M_J = 0.4; // remove m, keep j

    // jet removal
    for(int ijet = 0; ijet < (int)jets.size(); ijet++) {
        for(int iel = 0; iel < (int)electrons.size(); iel++) {
            if(jets.at(ijet)->p4().DeltaR(electrons.at(iel)->p4()) < DR_J_E) {
                jets.erase(jets.begin()+ijet);
                ijet--;
                break;
            } // dR match
        } // iel
    } // ijet

    // electron removal
    for(int iel = 0; iel < (int)electrons.size(); iel++) {
        for(int ijet = 0; ijet < (int)jets.size(); ijet++) {
            if(electrons.at(iel)->p4().DeltaR(jets.at(ijet)->p4()) < DR_E_J) {
                electrons.erase(electrons.begin() + iel);
                iel--;
                break;
            } // dR match
        } // ijet
    } // iel

    // muon removal
    for(int imu = 0; imu < (int)muons.size(); imu++) {
        for(int ijet = 0; ijet < (int)jets.size(); ijet++) {
            if(muons.at(imu)->p4().DeltaR(jets.at(ijet)->p4()) < DR_M_J) {
                muons.erase(muons.begin() + imu);
                imu--;
                break;
            } // dR match
        } // ijet
    } // imu

    ////////////////////////////////////////////////
    // group leptons
    ////////////////////////////////////////////////
    vector<xAOD::TruthParticle*> leptons;
    for(auto & el : electrons) leptons.push_back(el);
    for(auto & mu : muons) leptons.push_back(mu);
    sort(leptons.begin(), leptons.end(), by_pt);


    ///////////////////////////////////////////////
    // flavor "tagging"
    ///////////////////////////////////////////////

    // TODO implement "efficiencies" to choose c-jets as b-jets, etc..

    vector<xAOD::Jet*> bjets;
    for(const auto j : jets) {
        int flavor = std::abs(j->auxdata<int>("PartonTruthLabelID"));
        if(flavor==5) { bjets.push_back(j); }
    }

    /////////////////////////////////////////////
    //  selection
    /////////////////////////////////////////////
    // require dilepton events
    if(leptons.size() != 2) return kTRUE;
    if(!(bjets.size()>=2)) return kTRUE;
    if(x_mass() > 0) {
        look_at_resonance(leptons, bjets, met);
    }

    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsResonance::look_at_resonance(vector<xAOD::TruthParticle*> leptons,
        vector<xAOD::Jet*> bjets, const xAOD::MissingETContainer* met)
{
    float mt1 = get_MT_1(leptons, bjets, met, false);
    float mt1_scaled = get_MT_1(leptons, bjets, met, true);

    h_mt1_scaled->Fill(mt1_scaled);
    h_mt1->Fill(mt1);


}
//////////////////////////////////////////////////////////////////////////////
float DiHiggsResonance::get_MT_1(vector<xAOD::TruthParticle*> leptons,
        vector<xAOD::Jet*> bjets, const xAOD::MissingETContainer* met, bool do_bjet_rescaling)
{
    TLorentzVector l0;
    TLorentzVector l1;
    TLorentzVector b0;
    TLorentzVector b1;
    TLorentzVector MET;

    l0.SetPtEtaPhiM(leptons.at(0)->p4().Pt(), leptons.at(0)->p4().Eta(), leptons.at(0)->p4().Phi(),
            leptons.at(0)->p4().M());
    l1.SetPtEtaPhiM(leptons.at(1)->p4().Pt(), leptons.at(1)->p4().Eta(), leptons.at(1)->p4().Phi(),
            leptons.at(1)->p4().M());
    b0.SetPtEtaPhiM(bjets.at(0)->p4().Pt(), bjets.at(0)->p4().Eta(), bjets.at(0)->p4().Phi(),
            bjets.at(0)->p4().M());
    b1.SetPtEtaPhiM(bjets.at(1)->p4().Pt(), bjets.at(1)->p4().Eta(), bjets.at(1)->p4().Phi(),
            bjets.at(1)->p4().M());
    MET.SetPxPyPzE((*met)["NonInt"]->mpx(), (*met)["NonInt"]->mpy(), 0., (*met)["NonInt"]->met());

    TLorentzVector bjet_system = b0 + b1;
    if(do_bjet_rescaling) {
        //cout << "PERFORMING BJET SYSTEM MASS RESCALING" << endl;
        double mbb = bjet_system.M() * mev2gev;
        double scaling = 125.09/mbb;
        bjet_system.SetPtEtaPhiE(bjet_system.Pt() * scaling, bjet_system.Eta(),
            bjet_system.Phi(), bjet_system.E() * scaling);
    }

    TLorentzVector vis = (l0 + l1 + bjet_system);
    double pt_vis = vis.Pt();
    double m_vis = vis.M();
    double et_vis = sqrt(pt_vis*pt_vis + m_vis*m_vis);

    double x = (et_vis + MET.Pt());
    double y = (vis + MET).Pt();
    return sqrt( (x*x) - (y*y) ) * mev2gev;
}
//////////////////////////////////////////////////////////////////////////////
DiLepType DiHiggsResonance::get_lepton_type(const vector<xAOD::TruthParticle*> leptons)
{
    int n_e = 0;
    int n_m = 0;
    for(const auto& l : leptons) {
        if(l->absPdgId()==11) {
            n_e++;
        }
        else if(l->absPdgId()==13) {
            n_m++;
        }
    }

    if(n_e==2) {
        return DiLepType::EE;
    }
    else if(n_m==2) {
        return DiLepType::MM;
    }
    else if(n_e==1 && n_m==1) {
        return DiLepType::EM;
    }
    else {
        return DiLepType::Invalid;
    }
}
//////////////////////////////////////////////////////////////////////////////
}
