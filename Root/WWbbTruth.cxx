#include "xAODTruthAna/WWbbTruth.h"

//std/stl
#include <iostream>
#include <cstdio> // printf
#include <sstream>
#include <math.h> // sqrt, cos, tanh
#include <fstream>
using namespace std;

//ROOT
#include "TH1F.h"
#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"

//xAOD/EDM
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
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
    m_use_bjet_eff(false),
    m_is_hh_signal(false),
    m_rfile(0),
    m_tree(0),
    n_tree_fills(0)
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

    // weight tree
    h_cutflow_weighted = new TH1D("CutflowWeighted", "CutflowWeighted", 1000, 0., 1000.);
    h_cutflow_weighted->Sumw2();
    h_cutflow_weighted->SetDirectory(m_rfile);

    m_tree = new TTree("truth", "truth");

    // BRANCHES
    m_tree->Branch("dsid", &m_dsid);
    m_tree->Branch("w", &m_weight);
    m_tree->Branch("eventweight", &m_eventweight);
    m_tree->Branch("sumw", &m_sumw);
    m_tree->Branch("use_bjet_eff", &m_use_bjet_eff);
    m_tree->Branch("hh_signal", &m_is_hh_signal);
    m_tree->Branch("mcEventWeights", &m_mcEventWeights);

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
    m_tree->Branch("met_sumet", &m_met_sumet);
    m_tree->Branch("metphi_nonint", &m_metphi_nonint);
    m_tree->Branch("met_nonint", &m_met_nonint);

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
    m_tree->Branch("dr_l0b0", &m_dr_l0b0);
    m_tree->Branch("dphi_l0b1", &m_dphi_l0b1);
    m_tree->Branch("dr_l0b1", &m_dr_l0b1);

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
    m_tree->Branch("mt1", &m_mt1);
    m_tree->Branch("mt1_scaled", &m_mt1_scaled);

    // hh signal
    m_tree->Branch("mHH", &m_hh);
    m_tree->Branch("hhPt", &m_hh_pt);
    m_tree->Branch("dphi_hh", &m_dphi_hh);
    m_tree->Branch("drhh", &m_drhh);
    m_tree->Branch("h0_pt", &m_h0_pt);
    m_tree->Branch("h1_pt", &m_h1_pt);
    m_tree->Branch("h0_eta", &m_h0_eta);
    m_tree->Branch("h1_eta", &m_h1_eta);


    // THREE-BODY VARIABLES
    m_tree->Branch("nJets", &m_3b_njets);
    m_tree->Branch("nSJets", &m_3b_nsjets);
    m_tree->Branch("nBJets", &m_3b_nbjets);
    m_tree->Branch("nLeptons", &m_3b_nleptons);
    m_tree->Branch("nMuons", &m_3b_nmuons);
    m_tree->Branch("nElectrons", &m_3b_nelectrons);
    m_tree->Branch("MDR", &m_3b_mdr);
    m_tree->Branch("DPB_vSS", &m_3b_dpb);
    m_tree->Branch("cosThetaB", &m_3b_cosThetaB);
    m_tree->Branch("RPT", &m_3b_rpt);
    m_tree->Branch("gamInvRp1", &m_3b_gamInvRp1);
    m_tree->Branch("l_pt", &m_3b_lepPt);
    m_tree->Branch("l_q", &m_3b_lepQ);

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
        vector<float> all_weights = ei->mcEventWeights();
        for(unsigned int iw = 0; iw < all_weights.size(); ++iw) h_cutflow_weighted->Fill(iw, all_weights.at(iw));
        all_weights.clear();
        m_mcEventWeights = ei->mcEventWeights();
        m_eventweight = w();
        process_event();
    }

    //delete met;
    //store()->clear();

    return true;
}
//////////////////////////////////////////////////////////////////////////////
bool WWbbTruth::process_event()
{
    ///////////////////////////////////////////////////
    // get physics objects
    ///////////////////////////////////////////////////

    // jets
//    vector<xAOD::Jet*> jets;
    vector<xAOD::Jet> jets;
    const xAOD::JetContainer* xjets = 0;
    //RETURN_CHECK( GetName(), event()->retrieve( xjets, "AntiKt4TruthDressedWZJets" ));
    //RETURN_CHECK( GetName(), event()->retrieve( xjets, "ConeTruthLabelID" ));
    RETURN_CHECK( GetName(), event()->retrieve( xjets, "AntiKt4TruthJets" ));
    xAOD::JetContainer::const_iterator jet_itr = xjets->begin();
    xAOD::JetContainer::const_iterator jet_end = xjets->end();
    //for( ; jet_itr != jet_end; ++jet_itr) {
    for(const auto& j : *xjets) {
        if(!(j->p4().Pt() > 20. * GEV)) continue;
        //if(!( (*jet_itr)->p4().Pt() > 20. * GEV) ) continue;
        //xAOD::Jet* j = new xAOD::Jet();
        //j->makePrivateStore(**jet_itr);
        //jets.push_back(*j);
        //delete j;
        //jets.push_back(**jet_itr);
        jets.push_back(*j);
    } // jet_itr
    sort(jets.begin(), jets.end(), ByPtJet);

    // electrons
    //vector<xAOD::TruthParticle*> electrons;
    electrons.clear();
    const xAOD::TruthParticleContainer* xelectrons = 0;
    RETURN_CHECK( GetName(), event()->retrieve( xelectrons, "TruthElectrons" ));
    for(const auto& e : *xelectrons) {
        if( !((e->pt() * GEV) > 10.) ) continue;
        if(std::abs(e->eta()) > 2.47) continue;
        bool is_prompt = (e->auxdata<unsigned int>( "classifierParticleType" ) == 2);
        if(!is_prompt) continue;
        //xAOD::TruthParticle* ele = new xAOD::TruthParticle();
        //ele->makePrivateStore(*e);
        //electrons.push_back(ele);
        electrons.push_back(e);
    }
    sort(electrons.begin(), electrons.end(), ByPt);

    // muons
    //vector<xAOD::TruthParticle*> muons;
    muons.clear();
    const xAOD::TruthParticleContainer* xmuons = 0;
    RETURN_CHECK( GetName(), event()->retrieve( xmuons, "TruthMuons" ));
    for(const auto& m : *xmuons) {
        if( !( (m->pt() * GEV) > 10.) ) continue;
        if(std::abs(m->eta()) > 2.4) continue;
        bool is_prompt = (m->auxdata<unsigned int>("classifierParticleType" ) == 6);
        if(!is_prompt) continue;
        //xAOD::TruthParticle* mu = new xAOD::TruthParticle();
        //mu->makePrivateStore(*m);
        //muons.push_back(mu);
        muons.push_back(m);
    }
    sort(muons.begin(), muons.end(), ByPt);

    // MET
    const xAOD::MissingETContainer* met = 0;
    RETURN_CHECK( GetName(), event()->retrieve(met, "MET_Truth" ));

    // MET NonInt
    //const xAOD::MissingETContainer* met_nonint = 0;
    //RETURN_CHECK( GetName(), event()->retrieve(met_nonint, "MET_Truth_NonInt"));


    //////////////////////////////////////////////////////////////////
    // Overlap Removal
    //////////////////////////////////////////////////////////////////

    float DR_J_E = 0.2; // remove j, keep e
    float DR_E_J = 0.4; // remove e, keep j
    float DR_M_J = 0.4; // remove m, keep j

    // jet removal
    for(int ijet = 0; ijet < (int)jets.size(); ijet++) {
        for(int iel = 0; iel < (int)electrons.size(); iel++) {
            if(jets.at(ijet).p4().DeltaR(electrons.at(iel)->p4()) < DR_J_E) {
            //if(jets.at(ijet)->p4().DeltaR(electrons.at(iel)->p4()) < DR_J_E) {
                jets.erase(jets.begin()+ijet);
                ijet--;
                break;
            } // overlap
        } // iel
    } // ijet

    // electron removal
    for(int iel = 0; iel < (int)electrons.size(); iel++) {
        for(int ijet = 0; ijet < (int)jets.size(); ijet++) {
            if(electrons.at(iel)->p4().DeltaR(jets.at(ijet).p4()) < DR_E_J) {
            //if(electrons.at(iel)->p4().DeltaR(jets.at(ijet)->p4()) < DR_E_J) {
                electrons.erase(electrons.begin() + iel);
                iel--;
                break;
            } // overlap
        } // ijet
    } // iel

    // muon removal
    for(int imu = 0; imu < (int)muons.size(); imu++) {
        for(int ijet = 0; ijet < (int)jets.size(); ijet++) {
            if(muons.at(imu)->p4().DeltaR(jets.at(ijet).p4()) < DR_M_J) {
            //if(muons.at(imu)->p4().DeltaR(jets.at(ijet)->p4()) < DR_M_J) {
                muons.erase(muons.begin() + imu);
                imu--;
                break;
            } // overlap
        } // ijet
    } // imu

    /////////////////////////////////////////////////////
    // group leptons
    /////////////////////////////////////////////////////
    vector<const xAOD::TruthParticle*> leptons;
    for(auto & el : electrons) leptons.push_back(el);
    for(auto & mu : muons) leptons.push_back(mu);
    sort(leptons.begin(), leptons.end(), ByPt);

    // enforce dilepton events
    if( !(leptons.size()>=1) ) {
        return true;
    }

    /////////////////////////////////////////////////////
    // jet labelling
    /////////////////////////////////////////////////////
    //vector<xAOD::Jet*> bjets;
    //vector<xAOD::Jet*> sjets;
    vector<xAOD::Jet> bjets;
    vector<xAOD::Jet> sjets;
    for(const auto j : jets) {
        //int flavor = j.auxdata<int>("GhostBHadronsFinalCount");

        //int flavor = std::abs(j.auxdata<int>("ConeTruthLabelID"));
        int flavor = std::abs(j.auxdata<int>("PartonTruthLabelID"));
       // if(flavor>0) { bjets.push_back(j);}
       // else { sjets.push_back(j); }
        if(flavor==5) { bjets.push_back(j);}
        else { sjets.push_back(j); }
    }
    sort(bjets.begin(), bjets.end(), ByPtJet);
    sort(sjets.begin(), sjets.end(), ByPtJet);

    if(!(bjets.size()>=2)) {
        return true;
    }

    fill_tree(leptons, jets, sjets, bjets, met);
    //delete met;

    return true;
}
//////////////////////////////////////////////////////////////////////////////
//vector<double> WWbbTruth::super_razor(vector<xAOD::TruthParticle*> leptons,
//        vector<xAOD::Jet*> jets, const xAOD::MissingETContainer* met)
//{
//    vector<double> out; // [MDR, RPT, DPB, gamInvRp1]
//    double shatr;
//    TVector3 dummyvec;
//    double dpb;
//    double dummy;
//    double gamma;
//    double rpt;
//    double mdr;
//
//    vector<TLorentzVector> tlv_leptons;
//    for(auto & l : leptons) {
//        TLorentzVector tlv;
//        tlv.SetPtEtaPhiE(l->p4().Pt() * 1e-3, l->p4().Eta(), l->p4().Phi(), l->p4().E() * 1e-3); 
//        tlv_leptons.push_back(tlv);
//    }
//    TLorentzVector metTLV;
//    metTLV.SetPxPyPzE((*met)["NonInt"]->mpx() * 1e-3, (*met)["NonInt"]->mpy() * 1e-3, 0., (*met)["NonInt"]->met() * 1e-3);
//    m_metphi = metTLV.Phi();
//    m_met = metTLV.Pt() * GEV;
//    superRazor(tlv_leptons, metTLV, dummyvec, dummyvec, dummyvec, dummyvec,
//            shatr, dpb, dummy, gamma, dummy, mdr, rpt); 
//
//    out.push_back(mdr);
//    out.push_back(rpt);
//    out.push_back(dpb);
//    out.push_back(gamma);
//
//
//    return out;
//
//}
////////////////////////////////////////////////////////////////////////////////
//void WWbbTruth::superRazor(vector<TLorentzVector> leptons, TLorentzVector met, TVector3& vBETA_z, TVector3& pT_CM,
//                TVector3& vBETA_T_CMtoR, TVector3& vBETA_R,
//                double& SHATR, double& dphi_LL_vBETA_T, double& dphi_L1_L2,
//                double& gamma_R, double& dphi_vBETA_R_vBETA_T,
//                double& MDELTAR, double& rpt)
//{
//    //
//    // Code written by Christopher Rogan <crogan@cern.ch>, 04-23-13
//    // Details given in paper (http://arxiv.org/abs/1310.4827) written by
//    // Matthew R. Buckley, Joseph D. Lykken, Christopher Rogan, Maria Spiropulu
//    //
//    if (leptons.size() < 2) return;
//
//    // necessary variables
//    TLorentzVector l0 = leptons.at(0);
//    TLorentzVector l1 = leptons.at(1);
//
//    //
//    // Lab frame
//    //
//    //Longitudinal boost
//    vBETA_z = (1. / (l0.E() + l1.E()))*(l0 + l1).Vect();
//    vBETA_z.SetX(0.0);
//    vBETA_z.SetY(0.0);
//
//    l0.Boost(-vBETA_z);
//    l1.Boost(-vBETA_z);
//
//    //pT of CM frame
//    pT_CM = (l0 + l1).Vect() + met.Vect();
//    pT_CM.SetZ(0.0);
//
//    TLorentzVector ll = l0 + l1;
//    //invariant mass of the total event
//    SHATR = sqrt(2.*(ll.E()*ll.E() - ll.Vect().Dot(pT_CM)
//        + ll.E()*sqrt(ll.E()*ll.E() + pT_CM.Mag2() - 2.*ll.Vect().Dot(pT_CM))));
//
//    //rpt = pT_CM.Mag() / (pT_CM.Mag() + SHATR / 4.);
//    float jt = (leptons.at(0) + leptons.at(1) + met).Pt();
//    rpt = jt / ( jt + SHATR / 4. );
//
//    vBETA_T_CMtoR = (1. / sqrt(pT_CM.Mag2() + SHATR*SHATR))*pT_CM;
//
//    l0.Boost(-vBETA_T_CMtoR);
//    l1.Boost(-vBETA_T_CMtoR);
//    ll.Boost(-vBETA_T_CMtoR);
//
//    //
//    //R-frame
//    //
//    dphi_LL_vBETA_T = fabs((ll.Vect()).DeltaPhi(vBETA_T_CMtoR));
//
//    dphi_L1_L2 = fabs(l0.Vect().DeltaPhi(l1.Vect()));
//
//    vBETA_R = (1. / (l0.E() + l1.E()))*(l0.Vect() - l1.Vect());
//
//    gamma_R = 1. / sqrt(1. - vBETA_R.Mag2());
//
//    gamma_R = sqrt( 2 * (l0.Vect().Mag() * l1.Vect().Mag() + l0.Vect().Dot(l1.Vect())) ) / (l0.Vect().Mag() + l1.Vect().Mag());
//
//    dphi_vBETA_R_vBETA_T = fabs(vBETA_R.DeltaPhi(vBETA_T_CMtoR));
//
//    l0.Boost(-vBETA_R);
//    l1.Boost(vBETA_R);
//
//    //
//    //R+1 frame
//    //
//    MDELTAR = 2.*l0.E();
//    
//
//    return;
//}
//////////////////////////////////////////////////////////////////////////////
void WWbbTruth::fill_tree(vector<const xAOD::TruthParticle*> leptons,
        vector<xAOD::Jet> jets, vector<xAOD::Jet> sjets, vector<xAOD::Jet> bjets,
        const xAOD::MissingETContainer* met)
{
    // lepton 4vec
    m_l0_pt = leptons.at(0)->pt() * GEV;
    m_l0_eta = leptons.at(0)->eta();
    if(leptons.size()>1) {
        m_l1_pt = leptons.at(1)->pt() * GEV;
        m_l1_eta = leptons.at(1)->eta();
    }

    // jet 4vec
    for(size_t ij = 0; ij < jets.size(); ij++) {
        if(ij>=3) break;
        if(ij==0) {
            m_j0_pt = jets.at(0).pt() * GEV;
            m_j0_eta = jets.at(0).eta();
        }
        if(ij==1) {
            m_j1_pt = jets.at(1).pt() * GEV;
            m_j1_eta = jets.at(1).eta();
        }
        if(ij==2) {
            m_j2_pt = jets.at(2).pt() * GEV;
            m_j2_eta = jets.at(2).eta();
        }
    }
    // sjet 4vec
    for(size_t ij = 0; ij < sjets.size(); ij++) {
        if(ij>=3) break;
        if(ij==0) {
            m_sj0_pt = sjets.at(0).pt() * GEV;
            m_sj0_eta = sjets.at(0).eta();
        }
        if(ij==1) {
            m_sj1_pt = sjets.at(1).pt() * GEV;
            m_sj1_eta = sjets.at(1).eta();
        }
        if(ij==2) {
            m_sj2_pt = sjets.at(2).pt() * GEV;
            m_sj2_eta = sjets.at(2).eta();
        }
    }
    // bjet 4vec
    for(size_t ij = 0; ij < bjets.size(); ij++) {
        if(ij>=3) break;
        if(ij==0) {
            m_bj0_pt = bjets.at(0).pt() * GEV;
            m_bj0_eta = bjets.at(0).eta();
        }
        if(ij==1) {
            m_bj1_pt = bjets.at(1).pt() * GEV;
            m_bj1_eta = bjets.at(1).eta();
        }
        if(ij==2) {
            m_bj2_pt = bjets.at(2).pt() * GEV;
            m_bj2_eta = bjets.at(2).eta();
        }
    }

    // met
    TLorentzVector metTLV;
    metTLV.SetPxPyPzE((*met)["NonInt"]->mpx(), (*met)["NonInt"]->mpy(), 0., (*met)["NonInt"]->met());
    m_metphi = metTLV.Phi();
    m_met = metTLV.Pt() * GEV;
    m_met_sumet = (*met)["NonInt"]->sumet() * GEV;

    m_njets = jets.size();
    m_nsjets = sjets.size();
    m_nbjets = bjets.size();

    m_3b_nleptons = leptons.size();
    if(bjets.size()>=2) {
        m_mbb = (bjets.at(0).p4() + bjets.at(1).p4()).M() * GEV;
        m_dr_bb = (bjets.at(0).p4().DeltaR(bjets.at(1).p4()));
        m_dphi_bb = (bjets.at(0).p4().DeltaPhi(bjets.at(1).p4()));
        m_ptbb = (bjets.at(0).p4() + bjets.at(1).p4()).Pt() * GEV;

        m_dr_bbmet = (bjets.at(0).p4() + bjets.at(1).p4()).DeltaR(metTLV);
        m_dphi_bbmet = (bjets.at(0).p4() + bjets.at(1).p4()).DeltaPhi(metTLV);
        m_pt_bbmet = (bjets.at(0).p4() + bjets.at(1).p4() + metTLV).Pt() * GEV;

        // mt2_bb
        ComputeMT2 calc_mt2bb = ComputeMT2(bjets.at(0).p4(), bjets.at(1).p4(), metTLV);
        m_mt2_bb = calc_mt2bb.Compute() * GEV;
    }
    if(bjets.size()>0 && leptons.size()>0) {
        m_dphi_l0b0 = leptons.at(0)->p4().DeltaPhi(bjets.at(0).p4());
        m_dr_l0b0 = leptons.at(0)->p4().DeltaR(bjets.at(0).p4());
        if(bjets.size()>1) {
            m_dphi_l0b1 = leptons.at(0)->p4().DeltaPhi(bjets.at(1).p4());
            m_dr_l0b1 = leptons.at(0)->p4().DeltaR(bjets.at(1).p4());
        }
    }

    m_3b_njets = jets.size();
    m_3b_nsjets = sjets.size();
    m_3b_nbjets = bjets.size();

    if(!(leptons.size()>=2)) return;

    // dilepton
    m_mll = (leptons.at(0)->p4() + leptons.at(1)->p4()).M() * GEV;
    m_ptll = (leptons.at(0)->p4() + leptons.at(1)->p4()).Pt() * GEV;
    m_drll = (leptons.at(0)->p4().DeltaR(leptons.at(1)->p4()));
    m_dphi_ll = (leptons.at(0)->p4().DeltaPhi(leptons.at(1)->p4()));


    //TLorentzVector metTLV2;
    //metTLV2.SetPxPyPzE((*met)["Int"]->met()
    //m_met = (*met)["NonInt"]->met() * 1e-3;
    //m_met_nonint = (*met)["Int"]->met() * 1e-3;

    // leptons + met
    m_dr_llmet = (leptons.at(0)->p4() + leptons.at(1)->p4()).DeltaR(metTLV);
    m_dphi_llmet = (leptons.at(0)->p4() + leptons.at(1)->p4()).DeltaPhi(metTLV);
    m_metptll = (leptons.at(0)->p4() + leptons.at(1)->p4() + metTLV).Pt() * GEV;

    // di-bjet system
    if(bjets.size()>=2) {

        m_dr_llbb = (leptons.at(0)->p4() + leptons.at(1)->p4()).DeltaR( bjets.at(0).p4() + bjets.at(1).p4() );
        m_dphi_llbb = (leptons.at(0)->p4() + leptons.at(1)->p4()).DeltaPhi( bjets.at(0).p4() + bjets.at(1).p4() );

    }


    // total system
    if(bjets.size()>=2) {
        m_dphi_llmet_bb = (leptons.at(0)->p4() + leptons.at(1)->p4() + metTLV).DeltaPhi(bjets.at(0).p4() + bjets.at(1).p4());
        m_dr_llmet_bb = (leptons.at(0)->p4() + leptons.at(1)->p4() + metTLV).DeltaR(bjets.at(0).p4() + bjets.at(1).p4());

        m_ht2 = m_ptbb + m_metptll;

        m_sumpt = (metTLV.Pt() + leptons.at(0)->p4().Pt() + leptons.at(1)->p4().Pt()
                    + bjets.at(0).p4().Pt() + bjets.at(1).p4().Pt()) * GEV;
        m_ht2ratio = (m_ht2 / m_sumpt);

        // mt2_llbb
        TLorentzVector v0 = (leptons.at(0)->p4() + leptons.at(1)->p4());
        TLorentzVector v1 = (bjets.at(0).p4() + bjets.at(1).p4());
        ComputeMT2 calc_mt2llbb = ComputeMT2(v0, v1, metTLV);
        m_mt2_llbb = calc_mt2llbb.Compute() * GEV;


        // MT_1
        TLorentzVector bjet_system = (bjets.at(0).p4() + bjets.at(1).p4());
        float mbb_for_scaling = bjet_system.M();
        float scaling = 125.09 / mbb_for_scaling;

        TLorentzVector vis_system = (leptons.at(0)->p4() + leptons.at(1)->p4() + bjet_system); 
        float pt_vis = vis_system.Pt();
        float m_vis = vis_system.M();
        float et_vis = sqrt(pt_vis*pt_vis + m_vis*m_vis);
        float x = (et_vis + metTLV.Pt());
        float y = (vis_system + metTLV).Pt();
        m_mt1 = sqrt( (x*x) - (y*y) ) * GEV;

        // MT_1_scaled
        bjet_system.SetPtEtaPhiE(bjet_system.Pt() * scaling, bjet_system.Eta(),
                bjet_system.Phi(), bjet_system.E() * scaling);
        vis_system = (leptons.at(0)->p4() + leptons.at(1)->p4() + bjet_system);
        pt_vis = vis_system.Pt();
        m_vis = vis_system.M();
        et_vis = sqrt(pt_vis*pt_vis + m_vis*m_vis);
        x = (et_vis + metTLV.Pt());
        y = (vis_system + metTLV).Pt();
        m_mt1_scaled = sqrt( (x*x) - (y*y) ) * GEV;

    }

    #warning not yet considering HH signal variables
    m_hh = 0.;
    m_hh_pt = 0.;
    m_dphi_hh = -5.0;
    m_drhh = -1.0;
    m_h0_pt = 0.;
    m_h1_pt = 0.;
    m_h0_eta = 0.;
    m_h1_eta = 0.;

    // THREE-BODY STUFF
    int n_mu = 0;
    int n_el = 0;
    if(ee(leptons.at(0), leptons.at(1))) {
        n_mu = 0;
        n_el = 2;
    }
    else if(mm(leptons.at(0), leptons.at(1))) {
        n_mu = 2;
        n_el = 0;
    }
    else if(df(leptons.at(0), leptons.at(1))) {
        n_mu = 1;
        n_el = 1;
    }

    m_3b_nmuons = n_mu;
    m_3b_nelectrons = n_el;
    //vector<double> srazor = super_razor(leptons, bjets, met);
    //m_3b_mdr = srazor.at(0);
    //m_3b_dpb = srazor.at(2);
    //m_3b_rpt = srazor.at(1);
    //m_3b_gamInvRp1 = srazor.at(3);

    // cosThetaB
    TLorentzVector lp, lm; 
    for(auto & l : leptons) {
        if(l->charge() < 0) lm = l->p4(); 
        else if(l->charge() > 0) lp = l->p4();
    }
    TLorentzVector ll = lp + lm;
    TVector3 boost = ll.BoostVector();
    lp.Boost(-boost);
    lm.Boost(-boost);
    m_3b_cosThetaB = tanh((lp.Eta() - lm.Eta())/2.);

    //m_3b_lepPt.push_back(leptons.at(0)->pt() * 1e-3);
    //m_3b_lepPt.push_back(leptons.at(1)->pt() * 1e-3);

    //m_3b_lepQ.push_back(leptons.at(0)->charge());
    //m_3b_lepQ.push_back(leptons.at(1)->charge());
    


    ////////////////////////////////////////////////////////////////
    n_tree_fills++;
    m_tree->Fill();
}
//////////////////////////////////////////////////////////////////////////////
void WWbbTruth::Terminate()
{
    cout << "WWbbTruth::Terminate" << endl;

    if(m_do_sumw_calc) {
        cout << "WWbbTruth::Terminate    TOTAL SUMW = " << std::setprecision(6) << m_total_sumw << endl;
    }
    else {
        m_rfile->cd();
        m_tree->Write();
        m_rfile->Write();
    }
    cout << "WWbbTruth::Terminate    n_tree_fills = " << n_tree_fills << endl;
}


