#include "xAODTruthAna/WWbbTruth.h"

//std/stl
#include <iostream>
#include <cstdio> // printf
#include <sstream>
#include <math.h> // sqrt, cos
#include <fstream>
using namespace std;

//ROOT
#include "TH1F.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"

//xAOD/EDM
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthEventContainer.h"

//ANA
#include "xAODTruthAna/xaod_utils.h"
#include "xAODTruthAna/MT2_ROOT.h"

std::string ANANAME = "WWbbTruth";
const float GEV = 1e-3;

is_pdg IsElectron(11);
is_pdg IsMuon(13);
is_pdg IsTau(15);

pt_greater ByPt;
pt_greaterJet ByPtJet;

bool ee(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
{
    return IsElectron(l0) && IsElectron(l1);
}
bool mm(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
{
    return IsMuon(l0) && IsMuon(l1);
}
bool em(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
{
    return IsElectron(l0) && IsMuon(l1);
}
bool me(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
{
    return IsMuon(l0) && IsElectron(l1);
}
bool sf(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
{
    return ((ee(l0,l1)) || (mm(l0,l1)));
}
bool df(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
{
    return ((em(l0,l1)) || (me(l0,l1)));
}




//////////////////////////////////////////////////////////////////////////////
WWbbTruth::WWbbTruth() :
    m_outfilename(""),
    m_outfile_setup(false),
    m_do_sumw_calc(false),
    m_total_sumw(0.0),
    m_dsid(0),
    m_weight(0.0),
    m_eventweight(0.0),
    m_sumw(0.0),
    m_xsec(0.0),
    m_lumi(1.),
    m_rfile(0),
    m_tree(0)
{
    cout << "WWbbTruth()" << endl;
}
//////////////////////////////////////////////////////////////////////////////
void WWbbTruth::SlaveBegin(TTree* /*tree*/)
{
    timer()->Start();
}
//////////////////////////////////////////////////////////////////////////////
void WWbbTruth::initialize_xsec_map()
{
    xsec_map.clear();
    string filename = "xsec_wwbb.txt";
    std::ifstream dfs(filename.c_str());
    int d;
    double x;
    while(dfs >> d >> x) {
        xsec_map[d] = x;
    }

    cout << "------------------------------------------" << endl;
    cout << " XSEC LOADED " << endl;
    for(auto x : xsec_map) {
        cout << " > dsid " << x.first << "  xsec [pb] = " << x.second << endl;
    }
    cout << "------------------------------------------" << endl;
}
//////////////////////////////////////////////////////////////////////////////
void WWbbTruth::initialize_sumw_map()
{
    sumw_map.clear();

    string filename = "sumw_wwbb.txt";
    std::ifstream dfs(filename.c_str());
    if(!dfs.good()) {
        cout << "ERROR sumw map file could not be found!" << endl;
        exit(1);
    }
    int d;
    double s;
    while(dfs >> d >> s) {
        sumw_map[d] = s;
    }
    cout << "------------------------------------------" << endl;
    cout << " SUMW LOADED " << endl;
    for(auto x : sumw_map) {
        cout << " > dsid " << x.first << "  sumw = " << x.second << endl;
    }
    cout << "------------------------------------------" << endl;

    
}
//////////////////////////////////////////////////////////////////////////////
void WWbbTruth::setup_output_tree()
{
    stringstream ofn;
    ofn << "wwbb_truth_" << m_dsid; // << ".root";
    if(m_suffix != "") ofn << "_" << m_suffix;
    ofn << ".root";
    m_outfilename = ofn.str(); 
    m_rfile = new TFile(ofn.str().c_str(), "RECREATE");

    m_tree = new TTree("truth", "truth");


    // BRANCHES
    m_tree->Branch("dsid", &m_dsid);
    m_tree->Branch("w", &m_weight);
    m_tree->Branch("eventweight", &m_eventweight);
    m_tree->Branch("sumw", &m_sumw);

    // leptons
    m_tree->Branch("l_flav", &m_lepton_flavor);
    m_tree->Branch("l0_pt", &m_l0_pt);
    m_tree->Branch("l1_pt", &m_l1_pt);
    m_tree->Branch("l0_eta", &m_l0_eta);
    m_tree->Branch("l1_eta", &m_l1_eta);

    // jets
    m_tree->Branch("j0_pt", &m_j0_pt);
    m_tree->Branch("j1_pt", &m_j1_pt);
    m_tree->Branch("j2_pt", &m_j2_pt);
    m_tree->Branch("sj0_pt", &m_sj0_pt);
    m_tree->Branch("sj1_pt", &m_sj1_pt);
    m_tree->Branch("sj2_pt", &m_sj2_pt);
    m_tree->Branch("bj0_pt", &m_bj0_pt);
    m_tree->Branch("bj1_pt", &m_bj1_pt);
    m_tree->Branch("bj2_pt", &m_bj2_pt);
    m_tree->Branch("j0_eta", &m_j0_eta);
    m_tree->Branch("j1_eta", &m_j1_eta);
    m_tree->Branch("j2_eta", &m_j2_eta);
    m_tree->Branch("sj0_eta", &m_sj0_eta);
    m_tree->Branch("sj1_eta", &m_sj1_eta);
    m_tree->Branch("sj2_eta", &m_sj2_eta);
    m_tree->Branch("bj0_eta", &m_bj0_eta);
    m_tree->Branch("bj1_eta", &m_bj1_eta);
    m_tree->Branch("bj2_eta", &m_bj2_eta);

    m_tree->Branch("n_jets", &m_njets);
    m_tree->Branch("n_sjets", &m_nsjets);
    m_tree->Branch("n_bjets", &m_nbjets);

    // dilepton
    m_tree->Branch("mll", &m_mll);
    m_tree->Branch("ptll", &m_ptll);
    m_tree->Branch("dRll", &m_drll);
    m_tree->Branch("dphi_ll", &m_dphi_ll);

    // met
    m_tree->Branch("metphi", &m_metphi);
    m_tree->Branch("met", &m_met);

    // leptons + met
    m_tree->Branch("dr_llmet", &m_dr_llmet);
    m_tree->Branch("dphi_llmet", &m_dphi_llmet);
    m_tree->Branch("metptll", &m_metptll);

    // di-bjet system
    m_tree->Branch("mbb", &m_mbb);
    m_tree->Branch("dr_bb", &m_dr_bb);
    m_tree->Branch("dphi_bb", &m_dphi_bb);
    m_tree->Branch("ptbb", &m_ptbb);

    // lepton + di-bjet system
    m_tree->Branch("dr_llbb", &m_dr_llbb);
    m_tree->Branch("dphi_llbb", &m_dphi_llbb);
    m_tree->Branch("dphi_l0b0", &m_dphi_l0b0);

    // di-bjet system + met
    m_tree->Branch("dr_bbmet", &m_dr_bbmet);
    m_tree->Branch("dphi_bbmet", &m_dphi_bbmet);
    m_tree->Branch("pt_bbmet", &m_pt_bbmet);

    // total system
    m_tree->Branch("dphi_llmet_bb", &m_dphi_llmet_bb);
    m_tree->Branch("dr_llmet_bb", &m_dr_llmet_bb);
    m_tree->Branch("ht2", &m_ht2);
    m_tree->Branch("sumpt", &m_sumpt);
    m_tree->Branch("ht2ratio", &m_ht2ratio);
    m_tree->Branch("mt2_llbb", &m_mt2_llbb);
    m_tree->Branch("mt2_bb", &m_mt2_bb);

}
//////////////////////////////////////////////////////////////////////////////
Bool_t WWbbTruth::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    event()->getEntry(chain_entry);

    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo") );

    if(!m_outfile_setup && !m_do_sumw_calc) {
        cout << "setting up output ROOT file" << endl;

        initialize_sumw_map();
        initialize_xsec_map();

        m_xsec = xsec_map[m_dsid];
        m_sumw = sumw_map[m_dsid];

        setup_output_tree();

        m_outfile_setup = true;
    }

    string fn = "WWbbTruth::Process    ";
    if(dbg() || chain_entry%5000==0) {
        cout << fn << " **** Processing entry " << setw(6) << chain_entry << " **** " << endl;
    }

    if(m_do_sumw_calc) {
        m_total_sumw += ei->mcEventWeight();
    }
    else {
        m_weight = ei->mcEventWeight();
        m_eventweight = w();
        process_event();
    }

    return true;
}
//////////////////////////////////////////////////////////////////////////////
bool WWbbTruth::process_event()
{
    ///////////////////////////////////////////////////
    // get physics objects
    ///////////////////////////////////////////////////

    // jets
    vector<xAOD::Jet*> jets;
    const xAOD::JetContainer* xjets = 0;
    RETURN_CHECK( GetName(), event()->retrieve( xjets, "AntiKt4TruthJets" ));
    xAOD::JetContainer::const_iterator jet_itr = xjets->begin();
    xAOD::JetContainer::const_iterator jet_end = xjets->end();
    for( ; jet_itr != jet_end; ++jet_itr) {
        if(!( (*jet_itr)->p4().Pt() > 20. * GEV) ) continue;
        xAOD::Jet* j = new xAOD::Jet();
        j->makePrivateStore(**jet_itr);
        jets.push_back(j);
    } // jet_itr
    sort(jets.begin(), jets.end(), ByPtJet);

    // electrons
    vector<xAOD::TruthParticle*> electrons;
    const xAOD::TruthParticleContainer* xelectrons = 0;
    RETURN_CHECK( GetName(), event()->retrieve( xelectrons, "TruthElectrons" ));
    for(const auto& e : *xelectrons) {
        if( !((e->pt() * GEV) > 10.) ) continue;
        if(std::abs(e->eta()) > 2.47) continue;
        bool is_prompt = (e->auxdata<unsigned int>( "classifierParticleType" ) == 2);
        if(!is_prompt) continue;
        xAOD::TruthParticle* ele = new xAOD::TruthParticle();
        ele->makePrivateStore(*e);
        electrons.push_back(ele);
    }
    sort(electrons.begin(), electrons.end(), ByPt);

    // muons
    vector<xAOD::TruthParticle*> muons;
    const xAOD::TruthParticleContainer* xmuons = 0;
    RETURN_CHECK( GetName(), event()->retrieve( xmuons, "TruthMuons" ));
    for(const auto& m : *xmuons) {
        if( !( (m->pt() * GEV) > 10.) ) continue;
        if(std::abs(m->eta()) > 2.4) continue;
        bool is_prompt = (m->auxdata<unsigned int>("classifierParticleType" ) == 6);
        if(!is_prompt) continue;
        xAOD::TruthParticle* mu = new xAOD::TruthParticle();
        mu->makePrivateStore(*m);
        muons.push_back(mu);
    }
    sort(muons.begin(), muons.end(), ByPt);

    // MET
    const xAOD::MissingETContainer* met = 0;
    RETURN_CHECK( GetName(), event()->retrieve(met, "MET_Truth" ));


    //////////////////////////////////////////////////////////////////
    // Overlap Removal
    //////////////////////////////////////////////////////////////////

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
            } // overlap
        } // iel
    } // ijet

    // electron removal
    for(int iel = 0; iel < (int)electrons.size(); iel++) {
        for(int ijet = 0; ijet < (int)jets.size(); ijet++) {
            if(electrons.at(iel)->p4().DeltaR(jets.at(ijet)->p4()) < DR_E_J) {
                electrons.erase(electrons.begin() + iel);
                iel--;
                break;
            } // overlap
        } // ijet
    } // iel

    // muon removal
    for(int imu = 0; imu < (int)muons.size(); imu++) {
        for(int ijet = 0; ijet < (int)jets.size(); ijet++) {
            if(muons.at(imu)->p4().DeltaR(jets.at(ijet)->p4()) < DR_M_J) {
                muons.erase(muons.begin() + imu);
                imu--;
                break;
            } // overlap
        } // ijet
    } // imu

    /////////////////////////////////////////////////////
    // group leptons
    /////////////////////////////////////////////////////
    vector<xAOD::TruthParticle*> leptons;
    for(auto & el : electrons) leptons.push_back(el);
    for(auto & mu : muons) leptons.push_back(mu);
    sort(leptons.begin(), leptons.end(), ByPt);

    // enforce dilepton events
    if( !(leptons.size()==2) ) {
        return true;
    }

    /////////////////////////////////////////////////////
    // jet labelling
    /////////////////////////////////////////////////////
    vector<xAOD::Jet*> bjets;
    vector<xAOD::Jet*> sjets;
    for(const auto j : jets) {
        int flavor = std::abs(j->auxdata<int>("PartonTruthLabelID"));
        if(flavor==5) { bjets.push_back(j);}
        else { sjets.push_back(j); }
    }
    sort(bjets.begin(), bjets.end(), ByPtJet);
    sort(sjets.begin(), sjets.end(), ByPtJet);

    fill_tree(leptons, jets, sjets, bjets, met);


    return true;
}
//////////////////////////////////////////////////////////////////////////////
void WWbbTruth::fill_tree(vector<xAOD::TruthParticle*> leptons,
        vector<xAOD::Jet*> jets, vector<xAOD::Jet*> sjets, vector<xAOD::Jet*> bjets,
        const xAOD::MissingETContainer* met)
{
    // lepton 4vec
    m_l0_pt = leptons.at(0)->pt() * GEV;
    m_l1_pt = leptons.at(1)->pt() * GEV;
    m_l0_eta = leptons.at(0)->eta();
    m_l1_eta = leptons.at(1)->eta();

    // jet 4vec
    for(size_t ij = 0; ij < jets.size(); ij++) {
        if(ij>=3) break;
        if(ij==0) {
            m_j0_pt = jets.at(0)->pt() * GEV;
            m_j0_eta = jets.at(0)->eta();
        }
        if(ij==1) {
            m_j1_pt = jets.at(1)->pt() * GEV;
            m_j1_eta = jets.at(1)->eta();
        }
        if(ij==2) {
            m_j2_pt = jets.at(2)->pt() * GEV;
            m_j2_eta = jets.at(2)->eta();
        }
    }
    // sjet 4vec
    for(size_t ij = 0; ij < sjets.size(); ij++) {
        if(ij>=3) break;
        if(ij==0) {
            m_sj0_pt = sjets.at(0)->pt() * GEV;
            m_sj0_eta = sjets.at(0)->eta();
        }
        if(ij==1) {
            m_sj1_pt = sjets.at(1)->pt() * GEV;
            m_sj1_eta = sjets.at(1)->eta();
        }
        if(ij==2) {
            m_sj2_pt = sjets.at(2)->pt() * GEV;
            m_sj2_eta = sjets.at(2)->eta();
        }
    }
    // bjet 4vec
    for(size_t ij = 0; ij < bjets.size(); ij++) {
        if(ij>=3) break;
        if(ij==0) {
            m_bj0_pt = bjets.at(0)->pt() * GEV;
            m_bj0_eta = bjets.at(0)->eta();
        }
        if(ij==1) {
            m_bj1_pt = bjets.at(1)->pt() * GEV;
            m_bj1_eta = bjets.at(1)->eta();
        }
        if(ij==2) {
            m_bj2_pt = bjets.at(2)->pt() * GEV;
            m_bj2_eta = bjets.at(2)->eta();
        }
    }

    m_njets = jets.size();
    m_nsjets = sjets.size();
    m_nbjets = bjets.size();

    // dilepton
    m_mll = (leptons.at(0)->p4() + leptons.at(1)->p4()).M() * GEV;
    m_ptll = (leptons.at(0)->p4() + leptons.at(1)->p4()).Pt() * GEV;
    m_drll = (leptons.at(0)->p4().DeltaR(leptons.at(1)->p4()));
    m_dphi_ll = (leptons.at(0)->p4().DeltaPhi(leptons.at(1)->p4()));

    // met
    TLorentzVector metTLV;
    metTLV.SetPxPyPzE((*met)["NonInt"]->mpx(), (*met)["NonInt"]->mpy(), 0., (*met)["NonInt"]->met());
    m_metphi = metTLV.Phi();
    m_met = metTLV.Pt() * GEV;

    // leptons + met
    m_dr_llmet = (leptons.at(0)->p4() + leptons.at(1)->p4()).DeltaR(metTLV);
    m_dphi_llmet = (leptons.at(0)->p4() + leptons.at(1)->p4()).DeltaPhi(metTLV);
    m_metptll = (leptons.at(0)->p4() + leptons.at(1)->p4() + metTLV).Pt() * GEV;

    // di-bjet system
    if(bjets.size()>=2) {
        m_mbb = (bjets.at(0)->p4() + bjets.at(1)->p4()).M() * GEV;
        m_dr_bb = (bjets.at(0)->p4().DeltaR(bjets.at(1)->p4()));
        m_dphi_bb = (bjets.at(0)->p4().DeltaPhi(bjets.at(1)->p4()));
        m_ptbb = (bjets.at(0)->p4() + bjets.at(1)->p4()).Pt() * GEV;

        m_dr_llbb = (leptons.at(0)->p4() + leptons.at(1)->p4()).DeltaR( bjets.at(0)->p4() + bjets.at(1)->p4() );
        m_dphi_llbb = (leptons.at(0)->p4() + leptons.at(1)->p4()).DeltaPhi( bjets.at(0)->p4() + bjets.at(1)->p4() );

        m_dr_bbmet = (bjets.at(0)->p4() + bjets.at(1)->p4()).DeltaR(metTLV);
        m_dphi_bbmet = (bjets.at(0)->p4() + bjets.at(1)->p4()).DeltaPhi(metTLV);
        m_pt_bbmet = (bjets.at(0)->p4() + bjets.at(1)->p4() + metTLV).Pt() * GEV;
    }

    if(bjets.size()>0) {
        m_dphi_l0b0 = leptons.at(0)->p4().DeltaPhi(bjets.at(0)->p4());
    }

    // total system
    if(bjets.size()>=2) {
        m_dphi_llmet_bb = (leptons.at(0)->p4() + leptons.at(1)->p4() + metTLV).DeltaPhi(bjets.at(0)->p4() + bjets.at(1)->p4());
        m_dr_llmet_bb = (leptons.at(0)->p4() + leptons.at(1)->p4() + metTLV).DeltaR(bjets.at(0)->p4() + bjets.at(1)->p4());


        m_ht2 = m_ptbb + m_metptll;

        m_sumpt = (metTLV.Pt() + leptons.at(0)->p4().Pt() + leptons.at(1)->p4().Pt()
                    + bjets.at(0)->p4().Pt() + bjets.at(1)->p4().Pt()) * GEV;
        m_ht2ratio = (m_ht2 / m_sumpt);

        // mt2_llbb
        TLorentzVector v0 = (leptons.at(0)->p4() + leptons.at(1)->p4());
        TLorentzVector v1 = (bjets.at(0)->p4() + bjets.at(1)->p4());
        ComputeMT2 calc_mt2llbb = ComputeMT2(v0, v1, metTLV);
        m_mt2_llbb = calc_mt2llbb.Compute() * GEV;

        // mt2_bb
        ComputeMT2 calc_mt2bb = ComputeMT2(bjets.at(0)->p4(), bjets.at(1)->p4(), metTLV);
        m_mt2_bb = calc_mt2bb.Compute() * GEV;
    }


    ////////////////////////////////////////////////////////////////
    m_tree->Fill();
}
//////////////////////////////////////////////////////////////////////////////
void WWbbTruth::Terminate()
{
    cout << "WWbbTruth::Terminate" << endl;

    if(m_do_sumw_calc) {
        cout << "WWbbTruth::Terminate    TOTAL SUMW = " << m_total_sumw << endl;
    }
    else {
        m_rfile->cd();
        m_tree->Write();
        m_rfile->Write();

    }
}


