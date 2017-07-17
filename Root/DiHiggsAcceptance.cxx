#include "xAODTruthAna/DiHiggsAcceptance.h"

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

//ClassImp(DiHiggsAcceptance)

std::string ANA = "DiHiggsAcceptance";

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
DiHiggsAcceptance::DiHiggsAcceptance() :
    m_x_mass(0),
    m_weight(1.0),
    n_events_non_dilepton(0),
    n_events_non_dilepton_less(0),
    n_events_non_dilepton_more(0),
    m_outfile(nullptr),
    file_setup(false),
    total_sumw(0.0),
    m_xsec(0.0),
    m_sumw(0.0)
{
    cout << ANA << endl;


    for(unsigned int i = 0; i < DiLepType::Invalid; i++) {
        total_counts[i] = 0;
        passed_counts[i] = 0;
        total_w[i] = 0.0;
        passed_w[i] = 0.0;

        string histtype = "";
        if(i==DiLepType::EE) histtype = "EE";
        else if(i==DiLepType::MM) histtype = "MM";
        else if(i==DiLepType::EM) histtype = "EM";
        else if(i==DiLepType::ALL) histtype = "ALL";
        stringstream hname;
        stringstream htitle;
        hname << "h_total_w_" << histtype;
        htitle << "Total W (" << histtype << ");W;Entries";
        total_wh[i] = new TH1D(hname.str().c_str(), htitle.str().c_str(), 100, -0.4, 0.4);
        hname.str("");
        htitle.str("");
        hname << "h_passed_w_" << histtype;
        htitle << "Passed W (" << histtype << ");W;Entries";
        passed_wh[i] = new TH1D(hname.str().c_str(), htitle.str().c_str(), 100,  -0.4, 0.4);
    }
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsAcceptance::set_x_mass(int mass)
{
    stringstream hypo;
    if(mass>0) hypo << mass << " GeV resonance";
    else { hypo << "hh non-resonant"; }

    cout << ANA << "    signal hypothesis : " << hypo.str() << endl;

    m_x_mass = mass;
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsAcceptance::SlaveBegin(TTree* /*tree*/)
{
    timer()->Start();

}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsAcceptance::print_counts()
{
    stringstream sx;
    sx << "==================== COUNTS SUMMARY ======================\n";
    sx << "COUNTS " << x_mass() << "\n";
    sx << "COUNTS - - - - - - - - - - - - - - - - -\n";
    sx << "COUNTS EE TOTAL READ IN = " << total_w[DiLepType::EE] << "  [" << total_counts[DiLepType::EE] << "]\n";
    sx << "COUNTS EE TOTAL PASSED  = " << passed_w[DiLepType::EE] << "  [" << passed_counts[DiLepType::EE] << "]\n";
    if(total_counts[DiLepType::EE]>0)
    sx << "COUNTS EE ACCEPTANCE    = " << passed_w[DiLepType::EE]/total_w[DiLepType::EE] << "  [" << (float)passed_counts[DiLepType::EE]/(float)total_counts[DiLepType::EE] << "]\n";    
    else {
    sx << "COUNTS EE ACCEPTANCE    = -1\n"; 
    }
    sx << "COUNTS - - - - - - - - - - - - - - - - -\n";
    sx << "COUNTS MM TOTAL READ IN = " << total_w[DiLepType::MM] << "  [" << total_counts[DiLepType::MM] << "]\n";
    sx << "COUNTS MM TOTAL PASSED  = " << passed_w[DiLepType::MM] << "  [" << passed_counts[DiLepType::MM] << "]\n";
    if(total_counts[DiLepType::MM]>0)
    sx << "COUNTS MM ACCEPTANCE    = " << passed_w[DiLepType::MM] / total_w[DiLepType::MM] << "  [" << (float)passed_counts[DiLepType::MM]/(float)total_counts[DiLepType::MM] << "]\n";
    else
    sx << "COUNTS MM ACCEPTANCE    = -1\n";
    sx << "COUNTS - - - - - - - - - - - - - - - - -\n";
    sx << "COUNTS EM TOTAL READ IN = " << total_w[DiLepType::EM] << "  [" << total_counts[DiLepType::EM] << "]\n";
    sx << "COUNTS EM TOTAL PASSED  = " << passed_w[DiLepType::EM] << "  [" << passed_counts[DiLepType::EM] << "]\n";
    if(total_counts[DiLepType::EM]>0)
    sx << "COUNTS EM ACCEPTANCE    = " << passed_w[DiLepType::EM] / total_w[DiLepType::EM] << "  [" << (float)passed_counts[DiLepType::EM]/(float)total_counts[DiLepType::EM] << "]\n";
    else
    sx << "COUNTS EM ACCEPTANCe    = -1\n";
    sx << "COUNTS - - - - - - - - - - - - - - - - -\n";
    sx << "COUNTS ALL TOTAL READ IN= " << total_w[DiLepType::ALL] << "  [" << total_counts[DiLepType::ALL] << "]\n";
    sx << "COUNTS ALL PASSED       = " << passed_w[DiLepType::ALL] << "  [" << passed_counts[DiLepType::ALL] << "]\n";
    if(total_counts[DiLepType::ALL]>0)
    sx << "COUNTS ALL ACCEPTANCE   = " << passed_w[DiLepType::ALL] / total_w[DiLepType::ALL] << "  [" << (float)passed_counts[DiLepType::ALL]/(float)total_counts[DiLepType::ALL] << "]\n";
    else
    sx << "COUNTS ALL ACCEPTANCE   = -1\n";
    sx << "=============================================================\n";

    cout << sx.str() << endl;

}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsAcceptance::Terminate()
{
    print_counts();
    cout << "DiHiggsAcceptance::Terminate    There were " << n_events_non_dilepton << " events that were not dileptonic (< = " << n_events_non_dilepton_less << ", > = " << n_events_non_dilepton_more << ")" << endl;
    timer()->Stop();
    cout << timer_summary() << endl;
    for(int i = 0; i < DiLepType::Invalid; i++) {
        total_wh[i]->Write();
        passed_wh[i]->Write();
    }
    m_outfile->Write();
    m_outfile->Close();

    cout << "DiHiggsAcceptance::Terminate    TOTAL SUMW " << x_mass() << " " << total_sumw << endl;
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsAcceptance::initialize_sumw_map()
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
void DiHiggsAcceptance::initialize_xsec_map()
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
Bool_t DiHiggsAcceptance::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    n_evt_processed++;
    event()->getEntry(chain_entry);

    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo") );

    string fn = "DiHiggsAcceptance::Process    ";
    if(!file_setup) {
        stringstream fname;
        fname << "weights_";
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

    // incremenet all counts
    //total_counts[DiLepType::ALL]++;
    //total_w[DiLepType::ALL] += w();
    //total_wh[DiLepType::ALL]->Fill(w());

    // we should fill the "ALL" and denominator of each dilepton channel with the same
    // since we are calculating the acceptance of selecting (on top of the kinematic
    // selections) the given channel
    for(int i = 0; i < DiLepType::Invalid; i++) {
        total_counts[i]++;
        total_w[i] += w();
        total_wh[i]->Fill(w());
    }

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

    // require dilepton events
    if( !(leptons.size()==2) ) {
        n_events_non_dilepton++;
        if(leptons.size()>2) { n_events_non_dilepton_more++; }
        else { n_events_non_dilepton_less++; }
        return kTRUE; }

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
    if(x_mass() > 0) {
        resonance_acceptance(leptons, bjets, met);
    }
    else {
        non_resonant_acceptance();
    }


    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsAcceptance::resonance_acceptance(vector<xAOD::TruthParticle*> leptons,
            vector<xAOD::Jet*> bjets, const xAOD::MissingETContainer* met)
{
    string fn = "DiHiggsAcceptance::resonance_acceptance    ";
    float window_lower = 0.0;
    float window_upper = 0.0;
    x_mass_window_selection(window_lower, window_upper, 0.9, 1.08);


    DiLepType type = get_lepton_type(leptons);
    if(type==DiLepType::Invalid) {
        cout << fn << "Invalid dilepton type encountered, exiting" << endl;
        exit(1);
    }
    //total_counts[type]++;
    //total_w[type] += w();
    //total_wh[type]->Fill(w());

    // bjet multiplicity
    size_t n_bjets = bjets.size();
    if(!(n_bjets>=2)) return;

    // mll
    float mll = (leptons.at(0)->p4() + leptons.at(1)->p4()).M() * mev2gev;
    if(!(mll>20.)) return;

    // lepton pT
    float pt0 = leptons.at(0)->p4().Pt()*mev2gev;
    float pt1 = leptons.at(1)->p4().Pt()*mev2gev;
    if(!(pt0>25. && pt1>20.)) return;

    // SF Z-veto
    if(type==DiLepType::EE || type==DiLepType::MM) {
        if(!(fabs(mll-91.2)>10.)) return;
    }

    // dR ll
    float dRll = leptons.at(0)->p4().DeltaR(leptons.at(1)->p4());
    if(!(dRll<0.9)) return;
    //if(!(dRll<0.65)) return;

    // mt2_llbb
    float mt2_llbb = get_mt2_llbb(leptons, bjets, met);
    if(! (mt2_llbb > 90 && mt2_llbb < 140) ) return;

    // mbb
    float mbb = (bjets.at(0)->p4() + bjets.at(1)->p4()).M() * mev2gev;
    if(! (mbb > 100 && mbb < 140) ) return;

    // HT2Ratio
    float ht2ratio = get_ht2_ratio(leptons, bjets, met);
    if(! (ht2ratio>0.8) ) return;
    //if(! (ht2ratio>0.9) ) return;

    // MT_1
    float MT_1 = get_MT_1(leptons, bjets, met);
    if(! (MT_1 > window_lower && MT_1 < window_upper) ) return;

    passed_counts[DiLepType::ALL]++;
    passed_counts[type]++;
    passed_w[DiLepType::ALL] += w();
    passed_w[type] += w();
    passed_wh[DiLepType::ALL]->Fill(w());
    passed_wh[type]->Fill(w());

}
//////////////////////////////////////////////////////////////////////////////
float DiHiggsAcceptance::get_mt2_llbb(std::vector<xAOD::TruthParticle*> leptons,
        std::vector<xAOD::Jet*> bjets, const xAOD::MissingETContainer* met)
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

    TLorentzVector v0 = (l0 + l1);
    TLorentzVector v1 = (b0 + b1);

    ComputeMT2 calc = ComputeMT2(v0, v1, MET);
    return calc.Compute() * mev2gev;
}
//////////////////////////////////////////////////////////////////////////////
float DiHiggsAcceptance::get_ht2_ratio(vector<xAOD::TruthParticle*> leptons,
        vector<xAOD::Jet*> bjets, const xAOD::MissingETContainer* met)
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

    double numerator = ( (b0 + b1).Pt() + (l0 + l1 + MET).Pt() );
    double denominator = ( b0.Pt() + b1.Pt() + l0.Pt() + l1.Pt() + MET.Pt() );

    return (numerator / denominator);
}
//////////////////////////////////////////////////////////////////////////////
float DiHiggsAcceptance::get_MT_1(vector<xAOD::TruthParticle*> leptons,
        vector<xAOD::Jet*> bjets, const xAOD::MissingETContainer* met)
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

    TLorentzVector vis = (l0 + l1 + b0 + b1);
    double pt_vis = vis.Pt();
    double m_vis = vis.M();
    double et_vis = sqrt(pt_vis*pt_vis + m_vis*m_vis);

    double x = (et_vis + MET.Pt());
    double y = (vis + MET).Pt();
    return sqrt( (x*x) - (y*y) ) * mev2gev;
}
//////////////////////////////////////////////////////////////////////////////
void DiHiggsAcceptance::x_mass_window_selection(float& lower, float& upper, float dfactor, float ufactor)
{
    string fn = "DiHiggsAcceptance::x_mass_window_selection    ";
    lower = dfactor * x_mass();
    upper = ufactor * x_mass();
    //cout << fn << "Setting X mass window selection to (" << lower << "," << upper << ")" << endl; 
}
//////////////////////////////////////////////////////////////////////////////
DiLepType DiHiggsAcceptance::get_lepton_type(const vector<xAOD::TruthParticle*> leptons)
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
void DiHiggsAcceptance::non_resonant_acceptance()
{
    string fn = "DiHiggsAcceptance::non_resonant_acceptance    ";
    cout << fn << "Not yet implemented, exiting" << endl; 
    exit(1);
}
