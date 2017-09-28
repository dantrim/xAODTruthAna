#ifndef WWBBTRUTH_H
#define WWBBTRUTH_H

//TruthAnalysis
#include "TruthSelectorBase.h"

//xAOD
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

//std/stl
#include <vector>
#include <string>

//ROOT
class TH1F;
class TFile;
class TChain;
class TTree;
class TBranch;

//namespace wwbb {

    struct is_pdg {
        int pdg;
        is_pdg(int p) : pdg(p) {}
        bool operator()(const xAOD::TruthParticle *l) { return (l && l->absPdgId()==pdg); }
    };

    struct pt_greater {
        bool operator()(const xAOD::TruthParticle* a, const xAOD::TruthParticle* b)
                { return a->pt() > b->pt(); }
    };

    struct pt_greaterJet {
        bool operator()(const xAOD::Jet* a, const xAOD::Jet* b) { return a->pt() > b->pt(); }
    };

    bool ee(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1);
    bool mm(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1);
    bool em(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1);
    bool me(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1);
    bool sf(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1);
    bool df(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1);

//}

class WWbbTruth : public TruthSelectorBase
{

    public :
        WWbbTruth();
        virtual ~WWbbTruth(){};

        void do_sumw_calc(bool doit) { m_do_sumw_calc = doit; }
        void set_lumi(double l) { m_lumi = l; }
        void set_dsid(int d) { m_dsid = d; }
        void set_suffix(std::string suf) { m_suffix = suf; }

        // weight to luminosity (m_lumi is numerical value in fb-1, convert to pb-1)
        double w() { return ((m_weight * m_xsec * (m_lumi * 1000.)) / m_sumw); }

        void initialize_sumw_map();
        void initialize_xsec_map();

        void setup_output_tree();
        bool process_event();
        void fill_tree(std::vector<xAOD::TruthParticle*> leptons,
                std::vector<xAOD::Jet*> jets, std::vector<xAOD::Jet*> sjets, std::vector<xAOD::Jet*> bjets,
                const xAOD::MissingETContainer* met);

        // TSelector overrides
        virtual void SlaveBegin(TTree* tree);
        virtual void Terminate();
        virtual Bool_t Process(Long64_t entry);


        //enum DileptonType {
        //    EE=0,
        //    MM,
        //    EM,
        //    ME,
        //    ALL,
        //    Invalid
        //}

    private :

        std::string m_outfilename;
        std::string m_suffix;
        bool m_outfile_setup;
        bool m_do_sumw_calc;
        double m_total_sumw;
        int m_dsid;
        double m_weight;
        double m_eventweight;
        double m_sumw;
        double m_xsec;
        double m_lumi; // in fb-1

        std::map<int, double> sumw_map;
        std::map<int, double> xsec_map;

        // ROOT STUFF
        TFile* m_rfile;
        TTree* m_tree;

        ///////////////////////////////////////////
        // VARIABLES FOR OUTPUT TREE
        ///////////////////////////////////////////

        // leptons
        int m_lepton_flavor;
        float m_l0_pt;
        float m_l1_pt;
        float m_l0_eta;
        float m_l1_eta;

        // jets
        float m_j0_pt;
        float m_j1_pt;
        float m_j2_pt;
        float m_sj0_pt;
        float m_sj1_pt;
        float m_sj2_pt;
        float m_bj0_pt;
        float m_bj1_pt;
        float m_bj2_pt;

        float m_j0_eta;
        float m_j1_eta;
        float m_j2_eta;
        float m_sj0_eta;
        float m_sj1_eta;
        float m_sj2_eta;
        float m_bj0_eta;
        float m_bj1_eta;
        float m_bj2_eta;

        int m_njets;
        int m_nsjets;
        int m_nbjets;

        // dilepton
        float m_mll;
        float m_ptll;
        float m_drll;
        float m_dphi_ll;

        // met
        float m_metphi;
        float m_met;

        // leptons + met
        float m_dr_llmet;
        float m_dphi_llmet;
        float m_metptll;

        // di-bjet system
        float m_mbb;
        float m_dr_bb;
        float m_dphi_bb;
        float m_ptbb; 

        // lepton + di-bjet system
        float m_dr_llbb;
        float m_dphi_llbb;
        float m_dphi_l0b0;

        // di-bjet system + met
        float m_dr_bbmet;
        float m_dphi_bbmet;
        float m_pt_bbmet;

        // total system
        float m_dphi_llmet_bb;
        float m_dr_llmet_bb;

        float m_ht2;
        float m_sumpt;
        float m_ht2ratio;

        float m_mt2_llbb;
        float m_mt2_bb;

        // THREE-BODY VARIABLES
        int m_3b_njets;
        int m_3b_nsjets;
        int m_3b_nbjets;
        int m_3b_nleptons;
        float m_3b_mdr;
        float m_3b_dpb;
        float m_3b_cosThetaB;
        float m_3b_rpt;
        float m_3b_gamInvRp1;
        std::vector<float> m_3b_lepPt;
        



}; // class

#endif
