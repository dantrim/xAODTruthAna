#include "xAODTruthAna/HHFilterVal.h"

//std/stl
#include <iostream>
#include <cstdio>
#include <sstream>
#include <math.h>
#include <fstream>
#include <map>
using namespace std;

//ROOT
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
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
#include "xAODTruthAna/helpers.h"

//std::string ANANAME = "HHFilterVal";

////////////////////////////////////////////////////////////////////////////
HHFilterVal::HHFilterVal() :
    m_outfilename(""),
    m_suffix(""),
    m_outfile_setup(false),
    m_do_sumw_calc(false),
    m_lumi(1.),
    m_dsid(0),
    n_bad_higgs_decays(0),
    n_bb_higgs_decays(0),
    n_ww_higgs_decays(0),
    m_rfile(nullptr)
{
    cout << "HHFilterVal::HHFilterVal" << endl;
    for(int i = 0; i < DileptonFlavor::Invalid; i++) n_dilepton_flavor[i] = 0;
}
////////////////////////////////////////////////////////////////////////////
void HHFilterVal::SlaveBegin(TTree* /*tree*/)
{
    timer()->Start();
}
////////////////////////////////////////////////////////////////////////////
void HHFilterVal::setup_output_tree()
{
    stringstream ofn;
    ofn << "hh_validation_" << m_dsid;
    if(m_suffix != "") ofn << "_" << m_suffix;
    ofn << ".root";
    m_outfilename = ofn.str(); 
    m_rfile = new TFile(ofn.str().c_str(), "RECREATE");
    if(m_rfile->IsZombie()) {
        cout << "HHFilterVal::setup_output_tree   ERROR Could not create output file, it is zombie!" << endl;
        exit(1);
    }
    m_rfile->cd();

    h_n_higgses = new TH1F("h_n_higgses", "Number of Truth Higgses in the Event;Number of Higgses; Events", 5, 0, 5);
    h2_n_b_by_status = new TH2F("h2_n_b_by_status", "Number of b quarks by Truth Status; Number of Higgses; Truth Status", 8, 0, 8, 50, 0, 50);
    h2_n_w_by_status = new TH2F("h2_n_w_by_status", "Number of W-bosons by Truth Status; Number of W's; Truth Status", 5, 0, 5, 50, 0, 50);

    h_noh_h0_pt = new TH1F("h_noh_h0_pt", "Leading Higgs pT in Events with no h #rightarrow WW; p_{T} [GeV]; Events", 100, 0, -1);
    h_noh_h1_pt = new TH1F("h_noh_h1_pt", "Sub-leading Higgs pT in Events with no h #rightarrow WW; p_{T} [GeV]; Events", 100, 0, -1);
    h_noh_h0_eta = new TH1F("h_noh_h0_eta", "Leading Higgs #eta in Events with no h #rightarrow WW; #eta;Events", 100, -6, 6);
    h_noh_h1_eta = new TH1F("h_noh_h1_eta", "Sub-leading Higgs #eta in Events with no h #rightarrow WW; #eta;Events", 100, -6, 6);
    h_noh_h0_m = new TH1F("h_noh_h0_m", "Leading Higgs mass in Events with no h #rightarrow WW; m_{h0} [GeV]; Events", 100, 0, -1);
    h_noh_h1_m = new TH1F("h_noh_h1_m", "Sub-leading Higgs mass in Events with no h #rightarrow WW; m_{h1} [GeV]; Events", 100, 0, -1);

    // histograms
    h_n_h = new TH1F("h_n_h", "Number of Higgses Per Event;Number;Events",4,0,4);
    h_mh0 = new TH1F("h_mh0", "Leading higgs Invariant Mass; m_{h0} [GeV]; a.u.", 100, 0, -1);
    h_mh1 = new TH1F("h_mh1", "Sub-leading higgs Invariant Mass; m_{h1} [GeV]; a.u.", 100, 0, -1);
    h_mhh = new TH1F("h_mhh", "Dihiggs Invariant Mass;m_{hh} [GeV];a.u.", 100, 0, -1);
    h_h0_pt = new TH1F("h_h0_pt", "Leading higgs p_{T};p_{T}^{h0} [GeV];a.u.", 100, 0, -1);
    h_h1_pt = new TH1F("h_h1_pt", "Sub-leading higgs p_{T};p_{T}^{h1} [GeV];a.u.", 100, 0, -1);
    h_hh_pt = new TH1F("h_hh_pt", "Dihiggs p_{T}; p_{T}^{hh} [GeV];a.u.", 100, 0, -1);

    h_n_ww = new TH1F("h_n_ww", "Event has h #rightarrow WW; Has h #rightarrow WW?;Events",2,0,2);
    h_w0_pt = new TH1F("h_w0_pt", "Leading W p_{T};p_{T}^{W0} [GeV];a.u.", 100, 0, -1);
    h_w1_pt = new TH1F("h_w1_pt", "Sub-leading W p_{T};p_{T}^{W1} [GeV];a.u.",100,0,-1);
    h_mw0 = new TH1F("h_mw0", "Leading W Invariant Mass;m_{W0} [GeV];a.u.", 100, 0, -1);
    h_mw1 = new TH1F("h_mw1", "Sub-leading W Invariant Mass; m_{W1} [GeV];a.u.", 100, 0, -1);
    h_mww = new TH1F("h_mww", "WW Invariant Mass;m_{WW} [GeV];a.u.", 100, 0, -1);
    h_mww_from_child = new TH1F("h_mww_from_child", "WW Invariant Mass (from l+#nu); m_{WW} [GeV]; a.u.", 100, 0, -1);

    h_n_bb = new TH1F("h_n_bb", "Event has h #rightarrow bb; Has h #rightarrow bb?;Events",2,0,2);
    h_b0_pt = new TH1F("h_b0_pt", "Leading b p_{T};p_{T}^{b0} [GeV];a.u.", 100, 0, -1);
    h_b1_pt = new TH1F("h_b1_pt", "Sub-leading b p_{T}; p_{T}^{b1} [GeV];a.u.", 100, 0, -1);
    h_mbb = new TH1F("h_mbb", "bb Invariant Mass;m_{bb} [GeV];a.u.", 100, 0, -1);

    h_lepton_flavor = new TH1F("h_lepton_flavor", "Dilepton Flavor of h #rightarrow WW;;a.u.", 9,0,9);
    vector<string> labels = { "EE", "MM", "TT", "EM", "ME", "TE", "ET", "TM", "MT", "SF", "DF" };
    for(int i = 0; i < DileptonFlavor::SF; i++) {
        h_lepton_flavor->GetXaxis()->SetBinLabel(i+1, labels[i].c_str()); 
    }
    h_lepton_flavor_sfdf = new TH1F("h_lepton_flavor_sfdf", "Dilepton Flavor of h #rightarrow WW;;a.u.", 2, 0, 2);
    for(int i = DileptonFlavor::SF; i < DileptonFlavor::Invalid; i++)  {
        h_lepton_flavor_sfdf->GetXaxis()->SetBinLabel( (i-DileptonFlavor::SF)+1, labels[i].c_str());
    }

    h_n_leptons = new TH1F("h_n_leptons", "Number of Charged Leptons; Number;Events",8,0,8);
    h_l0_pt = new TH1F("h_l0_pt", "Lead Lepton p_{T};p_{T}^{l0} [GeV];a.u.", 40, 0, 200);
    h_l1_pt = new TH1F("h_l1_pt", "Sub-lead Lepton p_{T};p_{T}^{l1} [GeV]; a.u.", 20, 0, 100); 
    h_l0_eta = new TH1F("h_l0_eta", "Lead Lepton #eta;#eta;a.u.", 60, -3,3);
    h_l1_eta = new TH1F("h_l1_eta", "Sub-lead Lepton #eta;#eta;a.u.", 60, -3,3);

    h_l0_pt_from_w = new TH1F("h_l0_pt_from_w", "Lead leptton (from W decay); p_{T}^{l0} [GeV]; a.u.", 40, 0, 200);
    h_l1_pt_from_w = new TH1F("h_l1_pt_from_w", "Sub-lead lepton (from W decay); p_{T}^{l1} [GeV]; a.u.", 20, 0, 100);

}
////////////////////////////////////////////////////////////////////////////
Bool_t HHFilterVal::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    event()->getEntry(chain_entry);

    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo"));

    if(!m_outfile_setup) {
        setup_output_tree();
        m_outfile_setup = true;
    }


    // Truth Vertices
    const xAOD::TruthVertexContainer* truthVertices = 0;
    RETURN_CHECK(GetName(), event()->retrieve(truthVertices, "TruthVertices") );

    const xAOD::TruthParticle *higgs_bosons[2] = { nullptr };

    vector<HiggsCandidate> higgses;
    HiggsCandidate higgs_ww("WW");
    HiggsCandidate higgs_bb("BB");

    vector<int> children;
    int n_higgses = 0;
    //cout << "--------------------------------" << endl;
    for(const auto & vertex : *truthVertices) {
        if(vertex->nIncomingParticles() == 1 && vertex->nOutgoingParticles() >= 2) {
            const xAOD::TruthParticle* mother = vertex->incomingParticle(0);
            if(mother && mother->absPdgId()==25) {
                n_higgses++;
                if(mother->nChildren()>=2) {

                    int n_b_children = 0;
                    int n_w_children = 0;

                    vector<int> child_idx;

                    const xAOD::TruthParticle* c0 = nullptr;
                    const xAOD::TruthParticle* c1 = nullptr;

                    int n_children = mother->nChildren();
                    for(int i = 0; i < n_children; i++) {
                        if(mother->child(i)->absPdgId()==24) {
                            n_w_children++;
                            child_idx.push_back(i);
                        }
                        else if(mother->child(i)->absPdgId()==5) {
                            n_b_children++;
                            child_idx.push_back(i);
                        }
                    }

                    bool is_bb = (n_b_children==2 && n_w_children==0);
                    bool is_ww = (n_b_children==0 && n_w_children==2);

                    //const xAOD::TruthParticle* c0 = mother->child(0);
                    //const xAOD::TruthParticle* c1 = mother->child(1);
                    //children.push_back(c0->pdgId());
                    //children.push_back(c1->pdgId());
                    //cout << "mother higgs status = " << mother->status() << "  [child0: pdgId=" << c0->pdgId() << ", status = " << c0->status() << "] [child1: pdgId=" << c1->pdgId() << ", status = " << c1->status() << "]" << endl;
                    //if(c0->absPdgId()==5 && c1->absPdgId()==5) {
                    if(is_bb) {
                        if(higgs_bb.filled) {
                            cout << "HHFilterVal::Process    WARNING We found another Hbb candidate in the event!" << endl;
                        }

                        c0 = mother->child(child_idx.at(0));
                        c1 = mother->child(child_idx.at(1));

                        higgs_bb.filled = true;
                        higgs_bb.parent = mother;
                        higgs_bb.c0 = c0;
                        higgs_bb.c1 = c1;

                        //const xAOD::TruthParticle* bq0 = c0;
                        //bool selfDecay = false;
                        ////do {
                        ////    for(unsigned int i = 0; i < bq0->nChildren(); i++) {
                        ////        selfDecay = false;
                        ////        const xAOD::TruthParticle* child = bq0->child(i);
                        ////        if(child->pdgId()==bq0->pdgId()) {
                        ////            bq0 = child;
                        ////            selfDecay = true;
                        ////            break;
                        ////        }
                        ////    } // i
                        ////} while(selfDecay);
                        //higgs_bb.c0 = bq0;

                        //const xAOD::TruthParticle* bq1 = c1;
                        //selfDecay = false;
                        ////do {
                        ////    for(unsigned int i = 0; i < bq1->nChildren(); i++) {
                        ////        selfDecay = false;
                        ////        const xAOD::TruthParticle* child = bq1->child(i);
                        ////        if(child->pdgId()==bq1->pdgId()) {
                        ////            bq1 = child;
                        ////            selfDecay = true;
                        ////            break;
                        ////        }
                        ////    } // i
                        ////} while(selfDecay);
                        //higgs_bb.c1 = bq1;

                        n_bb_higgs_decays++;

                    }
                    //else if(c0->absPdgId()==24 && c1->absPdgId()==24) {
                    else if(is_ww) {
                        if(higgs_ww.filled) {
                            cout << "HHFilterVal::Process    WARNING We found another HWW candidate in the event!" << endl;
                        }

                        c0 = mother->child(child_idx.at(0));
                        c1 = mother->child(child_idx.at(1));

                        higgs_ww.filled = true;
                        higgs_ww.parent = mother;
                        higgs_ww.c0 = c0;
                        higgs_ww.c1 = c1;

                       // const xAOD::TruthParticle* w0 = c0;
                       // bool selfDecay = false;
                       // //do {
                       // //    for(unsigned int i = 0; i < w0->nChildren(); i++) {
                       // //        selfDecay = false;
                       // //        const xAOD::TruthParticle* child = w0->child(i);
                       // //        if(child->pdgId()==w0->pdgId()) {
                       // //            w0 = child;
                       // //            selfDecay = true;
                       // //            break;
                       // //        }
                       // //    } // i
                       // //} while(selfDecay);
                       // higgs_ww.c0 = w0;

                       // const xAOD::TruthParticle* w1 = c1;
                       // selfDecay = false;
                       // //do {
                       // //    for(unsigned int i = 0; i < w1->nChildren(); i++) {
                       // //        selfDecay = false;
                       // //        const xAOD::TruthParticle* child = w1->child(i);
                       // //        if(child->pdgId()==w1->pdgId()) {
                       // //            w1 = child;
                       // //            selfDecay = true;
                       // //            break;
                       // //        }
                       // //    } // i
                       // //} while(selfDecay);
                       // higgs_ww.c1 = w1;

                        n_ww_higgs_decays++;
                    }
                    else {
                        cout << "HHFilterVal::Process    WARNING We found a non WW or bb Higgs decay: children are [" << c0->pdgId() << ", " << c1->pdgId() << "]" << endl;
                        n_bad_higgs_decays++;
                    }
                }
                else {
                    cout << "HHFilterVal::Process   WARNING Found a Higgs not with 2 children (n chidlren = " << mother->nChildren() << ")" << endl;

                }
            }
        }
    } // vertex

    h_n_h->Fill(n_higgses);
    if(higgs_ww.filled) {
        h_n_ww->Fill(1.0);
    }
    else h_n_ww->Fill(0.0);

    if(higgs_bb.filled) {
        h_n_bb->Fill(1.0);
    }
    else h_n_bb->Fill(0.0);

    //if(!(higgs_ww.filled &&  higgs_bb.filled)) {
    //    cout << "WW filed? " << higgs_ww.filled << "  BB filled? " << higgs_bb.filled << ", n_higgs = " << n_higgses << endl;
    //    cout << " --> all children: ";
    //    for(auto ch : children) cout << ch << " ";
    //    cout << endl;
    //}

    higgses.clear();
    higgses.push_back(higgs_bb);
    higgses.push_back(higgs_ww);

    WCandidate w0("W0");
    WCandidate w1("W1");
    if(higgs_ww.filled) {
        if(!get_w_decay(higgs_ww.c0, w0)) { cout << "WHOOPS W0" << endl; return true; }
        if(!get_w_decay(higgs_ww.c1, w1)) { cout << "WHOOPS W1" << endl; return true; }
    }
    //else {
    //    cout << "HHFilterVal::Process    Did not fill WW" << endl;
    //}

    vector<WCandidate> wbosons;
    wbosons.push_back(w0);
    wbosons.push_back(w1);

    DileptonFlavor dilepton_flavor = DileptonFlavor::Invalid;

    float w = ei->mcEventWeight();

    int n_leptonic_w = 0;
    if(higgs_ww.filled) {
        for(auto wb : wbosons) {
            if(wb.filled) {
            //if(is_leptonic_w(wb)) {
                n_leptonic_w++;
            }
            else {
                cout << "HHFilterVal::Process    Found non-leptonically decaying W!" << endl;
            }
        }
        if(n_leptonic_w == 2) {
            const xAOD::TruthParticle* lep0 = wbosons.at(0).lepton;
            const xAOD::TruthParticle* lep1 = wbosons.at(1).lepton;

            if(lep0->pt() > lep1->pt()) {
                if(!lep0->isTau())
                h_l0_pt_from_w->Fill(lep0->pt() * helpers::MEV2GEV, w);
                if(!lep1->isTau())
                h_l1_pt_from_w->Fill(lep1->pt() * helpers::MEV2GEV, w);
            }
            else if(lep0->pt() < lep1->pt()) {
                if(!lep0->isTau())
                h_l0_pt_from_w->Fill(lep1->pt() * helpers::MEV2GEV, w);
                if(!lep1->isTau())
                h_l1_pt_from_w->Fill(lep0->pt() * helpers::MEV2GEV, w);
            }

            dilepton_flavor = get_dilepton_flavor_from_leptons(lep0, lep1);

            h_lepton_flavor->Fill( (int)dilepton_flavor); //, w);

            n_dilepton_flavor[dilepton_flavor]++;
            if(dilepton_flavor==DileptonFlavor::EE || dilepton_flavor==DileptonFlavor::MM || dilepton_flavor==DileptonFlavor::TT) {
                n_dilepton_flavor[DileptonFlavor::SF]++;
                h_lepton_flavor_sfdf->Fill(0.0); //, w);
            }
            else if(dilepton_flavor==DileptonFlavor::EM || dilepton_flavor==DileptonFlavor::ME || dilepton_flavor==DileptonFlavor::ET || dilepton_flavor==DileptonFlavor::MT || dilepton_flavor==DileptonFlavor::TE || dilepton_flavor==DileptonFlavor::TM) {
                n_dilepton_flavor[DileptonFlavor::DF]++;
                h_lepton_flavor_sfdf->Fill(1.0); //, w);
            }
            
        }
        //cout << "HHFilterVall::Process    n_leptonic_w = " << n_leptonic_w << endl;
    }

    if(higgs_ww.filled && higgs_bb.filled) {
        fill_histograms(higgs_ww, higgs_bb, wbosons);
    }

    const xAOD::TruthParticleContainer* truthParticles = 0;
    RETURN_CHECK(GetName(), event()->retrieve(truthParticles, "TruthParticles") );

    int n_ele = 0;
    int n_muo = 0;
    int n_higgses_tpr = 0;
    int n_bs = 0;
    stringstream status_stream;
    vector<TLorentzVector> leptons;
    vector<TLorentzVector> initial_higgses;
    for(const auto & tpr : *truthParticles) {
        //if(tpr->status() != 1) continue;
        if(tpr->absPdgId()==25) {
            //n_higgses_tpr++;
            status_stream << " " << tpr->status() << "  (n parents = " << tpr->nParents() << ", n children = " << tpr->nChildren();// << ") ";
            if(tpr->nChildren()==3) {
                status_stream << "  [children: ";
                for(int i = 0; i < 3; i++) status_stream << tpr->child(i)->pdgId() << " ";
                status_stream <<"]";
            }
            status_stream << ") ";
            if(tpr->nParents()==2) { initial_higgses.push_back( tpr->p4() ); n_higgses_tpr++; }
            //cout << "FOUND HIGGS" << endl;
        }
        //if(tpr->nParents()>0 && tpr->parent()->absPdgId()==25) {
        //    cout << "FOUND HIGGS AS PARENT!" << endl;
        //}
        if(!tpr->isChLepton()) continue;
        if(tpr->isElectron()) n_ele++;
        else if(tpr->isMuon()) n_muo++;
        //if(tpr->status()==1) // && (tpr->isElectron() || tpr->isMuon()))
        if(tpr->isElectron() || tpr->isMuon())
        leptons.push_back(tpr->p4());
    }
    if(n_higgses !=2) 
    cout << " FOUND " << n_higgses_tpr << " status: " << status_stream.str() << endl;
    helpers::pt_greaterTLV sortTLV;
    std::sort(initial_higgses.begin(), initial_higgses.end(), sortTLV);
    h_n_higgses->Fill(n_higgses_tpr);

    if(n_higgses != 2) {
        if(initial_higgses.size()>0) {
            h_noh_h0_pt->Fill(initial_higgses.at(0).Pt() * helpers::MEV2GEV);
            h_noh_h0_eta->Fill(initial_higgses.at(0).Eta());
            h_noh_h0_m->Fill(initial_higgses.at(0).M() * helpers::MEV2GEV);
        }
        if(initial_higgses.size()>1) {
            h_noh_h1_pt->Fill(initial_higgses.at(1).Pt() * helpers::MEV2GEV);
            h_noh_h1_eta->Fill(initial_higgses.at(1).Eta());
            h_noh_h1_m->Fill(initial_higgses.at(1).M() * helpers::MEV2GEV);
        }

    }

    leptons.clear();
    
    const xAOD::TruthParticleContainer* x_elec = 0;
    RETURN_CHECK(GetName(), event()->retrieve(x_elec, "TruthElectrons") ); 
    const xAOD::TruthParticleContainer* x_muon = 0;
    RETURN_CHECK(GetName(), event()->retrieve(x_muon, "TruthMuons") );

    for(const auto & lep : *x_elec) {
        if(lep->status()==1) leptons.push_back(lep->p4());
    }
    for(const auto & lep : *x_muon) {
        if(lep->status()==1) leptons.push_back(lep->p4());
    }

//    h_n_leptons->Fill( (n_ele+n_muo) );
    h_n_leptons->Fill( leptons.size() );
    std::sort(leptons.begin(), leptons.end(), sortTLV);
    if(leptons.size()>0) {
        h_l0_pt->Fill(leptons.at(0).Pt() * helpers::MEV2GEV);//, w);
        h_l0_eta->Fill(leptons.at(0).Eta());//, w);
    }
    if(leptons.size()>1) {
        h_l1_pt->Fill(leptons.at(1).Pt() * helpers::MEV2GEV);//, w);
        h_l1_eta->Fill(leptons.at(1).Eta());//, w);
    }
    

    


    return true;
}
////////////////////////////////////////////////////////////////////////////
bool HHFilterVal::get_w_decay(const xAOD::TruthParticle* parent, HHFilterVal::WCandidate& wboson)
{

    if(parent==nullptr) {
        cout << "HHFilterVal::get_w_decay    ERROR Input parent W is null" << endl;
        return false;
    }
    if(!parent->isW()) {
        cout << "HHFilterVal::get_w_decay    ERROR Input parent W is not a W!" << endl;
        return false;
    }

    const xAOD::TruthParticle* wint = parent;
    const xAOD::TruthParticle* neutrino = nullptr;
    const xAOD::TruthParticle* lepton = nullptr;

    bool selfDecay = false;
    do  {
        for(unsigned int i = 0; i < wint->nChildren(); i++) {
            selfDecay = false;
            const xAOD::TruthParticle* child = wint->child(i);
            if(child->pdgId()==wint->pdgId()) {
                wint=child;
                selfDecay = true;
                break;
            }
            else if(child->isChLepton()) {
                lepton = child;
            }
            else if(child->isNeutrino()) {
                neutrino = child;
            }
        } // i
    } while (selfDecay);

    bool is_bad = false;
    if(neutrino==nullptr) {
        cout << "HHFilterVal::get_w_decay    WARNING No neutrino found (parent pt,eta = " << parent->pt() <<", " << parent->eta() << ")" << endl;
        is_bad = true;
    }
    if(lepton==nullptr) {
        cout << "HHFilterVal::get_w_decay    WARNING No lepton found (parent pt,eta= " << parent->pt() << ", " << parent->eta() << ")" << endl;
        is_bad = true;
    }
    if(is_bad) { wboson.filled = false; return false; }
    wboson.w = parent;
    wboson.lepton = lepton;
    wboson.neutrino = neutrino;
    wboson.filled = true;
    
    

//    cout << "HHFilterVal::get_w_decay   INPUT W PT="<<parent->pt() <<", WINT PT=" << wint->pt() << ", lepton PT=" << lepton->pt() << ", neutrino PT=" <<  neutrino->pt() << endl; //"  (lep+v) = " << (*(lepton->p4() + neutrino->p4())).Pt() << endl;


//    HiggsCandidate out("");
//    if(!parent->isW()) { out.filled = false; return out; }
//    //int n_children = parent->nChildren();
//    //if(n_children!=2) {
//    //    cout << "HHFilterVal::get_w_decay    W boson decay does not have 2 children! (has=" << n_children << ", ";
//    //    for(int i = 0; i < (int)n_children; i++) cout << " " << parent->child(i)->pdgId();
//    //    cout << ")" << endl;
//    //    out.filled = false;
//    //    return out;
//    //}
//
//    const xAOD::TruthParticle* w = parent;
//    const xAOD::TruthParticle* lepton = nullptr;
//    const xAOD::TruthParticle* neutrino = nullptr;
//    bool selfDecay = false;
//    do {
//        for(unsigned int i = 0; i < w->nChildren(); i++) {
//            selfDecay = false;
//            const xAOD::TruthParticle* child = w->child(i);
//            if(child->pdgId()==w->pdgId()) {
//                w = child;
//                selfDecay = true;
//                break;
//            }
//        } // i
//    } while(selfDecay);
//
//    int n_children = parent->nChildren();
//    if(n_children!=2) {
//        cout << "HHFilterVal::get_w_decay   W boson decay does not have 2 children! (has=" << n_children << ", ";
//        for(int i = 0; i < (int)n_children; i++) cout << " " << w->child(i)->pdgId();
//        cout << ")" << endl;
//        out.filled = false;
//        return out;
//    }
//    for(int i = 0; i < 2; i++) {
//        if(w->child(i)->isChLepton()) {
//            if(lepton) {
//                cout << "HHFilterVal::get_w_decay    W boson has multiple leptons!" << endl;
//            }
//            lepton = w->child(i);
//        }
//        else if(w->child(i)->isNeutrino()) {
//            if(neutrino) {
//                cout << "HHFilterVal::get_w_decay    W boson has multiple neutrinos!" << endl;
//            }
//            neutrino = w->child(i);
//        }
//    }
//    bool is_bad = false;
//    if(lepton == nullptr) {
//        cout << "HHFilterVal::get_w_decay    W boson does not have lepton decay!" << endl;
//        is_bad = true;
//    }
//    if(neutrino == nullptr) {
//        cout << "HHFilterVal::get_w_decay    W boson does not have neutrino decay!" << endl;
//        is_bad = true;
//    }
//    if(is_bad) {
//        out.c0 = nullptr;
//        out.c1 = nullptr;
//        out.filled = false;
//        return out;
//    }
//    out.c0 = lepton;
//    out.c1 = neutrino;
//    out.filled = true;
    return true;
}
////////////////////////////////////////////////////////////////////////////
bool HHFilterVal::is_leptonic_w(HiggsCandidate& wboson)
{
    if(wboson.c0->isChLepton() && wboson.c1->isNeutrino())
        return true;
    else if(wboson.c0->isNeutrino() && wboson.c1->isChLepton())
        return true;
    return false;
}
////////////////////////////////////////////////////////////////////////////
const xAOD::TruthParticle* HHFilterVal::lepton_from_w(HiggsCandidate w)
{
    const xAOD::TruthParticle* lepton = nullptr;
    if(w.c0->isChLepton()) {
        lepton = w.c0;
    }
    else if(w.c1->isChLepton()) {
        lepton = w.c1;
    }
    return lepton;
}
////////////////////////////////////////////////////////////////////////////
const xAOD::TruthParticle* HHFilterVal::neutrino_from_w(HiggsCandidate w)
{
    const xAOD::TruthParticle* neutrino = nullptr;
    if(w.c0->isNeutrino()) {
        neutrino = w.c0;
    }
    else if(w.c1->isNeutrino()) {
        neutrino = w.c1;
    }
    return neutrino;
}
////////////////////////////////////////////////////////////////////////////
HHFilterVal::DileptonFlavor HHFilterVal::get_dilepton_flavor_from_leptons(
        const xAOD::TruthParticle* part0, const xAOD::TruthParticle* part1)
{
    bool e0 = part0->isElectron();
    bool m0 = part0->isMuon();
    bool t0 = part0->isTau();
    bool e1 = part1->isElectron();
    bool m1 = part1->isMuon();
    bool t1 = part1->isTau();

    if(e0 && e1) return DileptonFlavor::EE;
    else if(e0 && m1) return DileptonFlavor::EM;
    else if(e0 && t1) return DileptonFlavor::ET;
    else if(m0 && m1) return DileptonFlavor::MM;
    else if(m0 && e1) return DileptonFlavor::ME;
    else if(m0 && t1) return DileptonFlavor::MT;
    else if(t0 && t1) return DileptonFlavor::TT;
    else if(t0 && e1) return DileptonFlavor::TE;
    else if(t0 && m1) return DileptonFlavor::TM;
    else
        return DileptonFlavor::Invalid;
}
////////////////////////////////////////////////////////////////////////////
HHFilterVal::DileptonFlavor HHFilterVal::get_dilepton_flavor(vector<HiggsCandidate> wbosons)
{
    if(wbosons.size()!=2) return DileptonFlavor::Invalid;

    const xAOD::TruthParticle* lepton0 = nullptr;
    const xAOD::TruthParticle* lepton1 = nullptr;

    if(wbosons.at(0).c0->isChLepton()) lepton0 = wbosons.at(0).c0;
    else if(wbosons.at(0).c1->isChLepton()) lepton0 = wbosons.at(0).c1;

    if(wbosons.at(1).c0->isChLepton()) lepton0 = wbosons.at(0).c0;
    else if(wbosons.at(1).c1->isChLepton()) lepton0 = wbosons.at(0).c1;

    return get_dilepton_flavor_from_leptons(lepton0, lepton1);
}
////////////////////////////////////////////////////////////////////////////
bool HHFilterVal::fill_histograms(const HiggsCandidate& hww, const HiggsCandidate& hbb,
        const vector<WCandidate>& wbosons)
{
    helpers::pt_greaterTLV sortTLV;

    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo"));
    double w = ei->mcEventWeight();

    vector<TLorentzVector> vec_hh_tlv = { hww.parent->p4(), hbb.parent->p4() };
    std::sort(vec_hh_tlv.begin(), vec_hh_tlv.end(), sortTLV);

    TLorentzVector hh_system;
    hh_system += hww.parent->p4();
    hh_system += hbb.parent->p4();

    w = 1;

    h_mhh->Fill(hh_system.M() * helpers::MEV2GEV, w);
    h_h0_pt->Fill(vec_hh_tlv.at(0).Pt() * helpers::MEV2GEV, w);
    h_h1_pt->Fill(vec_hh_tlv.at(1).Pt() * helpers::MEV2GEV, w);
    h_hh_pt->Fill(hh_system.Pt() * helpers::MEV2GEV, w);
    h_mh0->Fill(vec_hh_tlv.at(0).M() * helpers::MEV2GEV, w);
    h_mh1->Fill(vec_hh_tlv.at(1).M() * helpers::MEV2GEV, w);

    TLorentzVector ww_system;
    ww_system += hww.c0->p4();
    ww_system += hww.c1->p4();
    vector<TLorentzVector> vec_ww_tlv = { hww.c0->p4(), hww.c1->p4() } ;
    std::sort(vec_ww_tlv.begin(), vec_ww_tlv.end(), sortTLV);

    h_w0_pt->Fill(vec_ww_tlv.at(0).Pt() * helpers::MEV2GEV, w);
    h_w1_pt->Fill(vec_ww_tlv.at(1).Pt() * helpers::MEV2GEV, w);
    h_mw0->Fill(vec_ww_tlv.at(0).M() * helpers::MEV2GEV, w);
    h_mw1->Fill(vec_ww_tlv.at(1).M() * helpers::MEV2GEV, w);
    h_mww->Fill(ww_system.M() * helpers::MEV2GEV, w);

    TLorentzVector lvlv_system;
    lvlv_system += wbosons.at(0).lepton->p4();
    lvlv_system += wbosons.at(0).neutrino->p4();
    lvlv_system += wbosons.at(1).lepton->p4();
    lvlv_system += wbosons.at(1).neutrino->p4();
    h_mww_from_child->Fill(lvlv_system.M() * helpers::MEV2GEV, w);

    TLorentzVector bb_system;
    bb_system += hbb.c0->p4();
    bb_system += hbb.c1->p4();
    vector<TLorentzVector> vec_bb_tlv = { hbb.c0->p4(), hbb.c1->p4() };
    std::sort(vec_bb_tlv.begin(), vec_bb_tlv.end(), sortTLV);

    h_b0_pt->Fill(vec_bb_tlv.at(0).Pt() * helpers::MEV2GEV, w);
    h_b1_pt->Fill(vec_bb_tlv.at(1).Pt() * helpers::MEV2GEV, w);
    h_mbb->Fill(bb_system.M() * helpers::MEV2GEV, w);

    

    return true; 
}
////////////////////////////////////////////////////////////////////////////
void HHFilterVal::Terminate()
{
    cout << "HHFtilerVal::Terminate" << endl;

    cout << "HHFilterVal::Terminate    Total counts [hWW, hbb] = [" <<n_ww_higgs_decays <<", " << n_bb_higgs_decays << "], [hbb/hWW] = ["<<(1.0*n_bb_higgs_decays)/(1.0*n_ww_higgs_decays) << "]" << endl;
    cout << "HHFilterVal::Terminate    Total bad HH decays: " << n_bad_higgs_decays << endl;
    cout << "Dilepton flavor counters:" << endl;
    cout << "  EE: " << n_dilepton_flavor[DileptonFlavor::EE] << endl;
    cout << "  MM: " << n_dilepton_flavor[DileptonFlavor::MM] << endl;
    cout << "  TT: " << n_dilepton_flavor[DileptonFlavor::TT] << endl;
    cout << "  EM: " << n_dilepton_flavor[DileptonFlavor::EM] << endl;
    cout << "  ME: " << n_dilepton_flavor[DileptonFlavor::ME] << endl;
    cout << "  ET: " << n_dilepton_flavor[DileptonFlavor::ET] << endl;
    cout << "  TE: " << n_dilepton_flavor[DileptonFlavor::TE] << endl;
    cout << "  MT: " << n_dilepton_flavor[DileptonFlavor::MT] << endl;
    cout << "  TM: " << n_dilepton_flavor[DileptonFlavor::TM] << endl;
    cout << "  SF: " << n_dilepton_flavor[DileptonFlavor::SF] << endl;
    cout << "  DF: " << n_dilepton_flavor[DileptonFlavor::DF] << endl;

    save_histograms();
}
////////////////////////////////////////////////////////////////////////////
bool HHFilterVal::save_histograms()
{
    if(!m_rfile) {
        cout << "HHFilterVal::save_histograms    Output file is null, will not save histos" << endl;
        return false;
    }
//    m_rfile->cd();
    cout << "this is the fucking file: " << m_rfile << endl;

    h_n_higgses->Write();
    h2_n_b_by_status->Write();
    h2_n_w_by_status->Write();

    h_noh_h0_pt->Write();
    h_noh_h1_pt->Write();
    h_noh_h0_eta->Write();
    h_noh_h1_eta->Write();
    h_noh_h0_m->Write();
    h_noh_h1_m->Write();

    h_n_h->Write();
    h_mh0->Write();
    h_mh1->Write();
    h_mhh->Write(); 
    h_h0_pt->Write();
    h_h1_pt->Write();
    h_hh_pt->Write();

    h_n_ww->Write();
    h_w0_pt->Write();
    h_w1_pt->Write();
    h_mw0->Write();
    h_mw1->Write();
    h_mww->Write();
    h_mww_from_child->Write();

    h_n_bb->Write();
    h_b0_pt->Write();
    h_b1_pt->Write();
    h_mbb->Write();

    h_n_leptons->Write();
    h_lepton_flavor->Write();
    h_lepton_flavor_sfdf->Write();

    h_l0_pt->Write();
    h_l1_pt->Write();
    h_l0_eta->Write();
    h_l1_eta->Write();

    h_l0_pt_from_w->Write();
    h_l1_pt_from_w->Write();

    //m_rfile->Write();
    //m_rfile->Close();

    return true;
}
