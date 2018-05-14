#ifndef WWBBTRUTH_H
#define WWBBTRUTH_H

//TruthAnalysis
#include "TruthSelectorBase.h"

//xAOD
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODTruth/TruthVertexContainer.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

//std/stl
#include <vector>
#include <string>

//ROOT
class TH1F;
class TH1D;
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
        bool operator()(const xAOD::Jet a, const xAOD::Jet b) { return a.pt() > b.pt(); }
    };

    struct pt_greaterTLV {
        bool operator()(const TLorentzVector& a, const TLorentzVector& b) { return a.Pt() > b.Pt(); }
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

        void set_use_bjet_eff(bool doit) { m_use_bjet_eff = doit; }
        void set_is_hh_signal(bool doit) { m_is_hh_signal = doit; }
        void set_skip_maps(bool doit) { m_skip_maps = doit; }
        void set_top_sample(bool doit) { m_top_sample = doit; }
        bool use_bjet_eff() { return m_use_bjet_eff; }
        bool is_hh_signal() { return m_is_hh_signal; }
        bool skip_maps() { return m_skip_maps; }
        bool top_sample() { return m_top_sample; }

        void initialize_sumw_map();
        void initialize_xsec_map();

        void setup_output_tree();
        bool process_event();
        //std::vector<double> super_razor(std::vector<const xAOD::TruthParticle*> leptons,
        //        std::vector<xAOD::Jet*> jets, const xAOD::MissingETContainer* met);
        //void superRazor(std::vector<TLorentzVector> leptons, TLorentzVector met,
        //        TVector3& vBETA_z, TVector3& pT_CM, TVector3& vBETA_T_CMtoR, TVector3& vBETA_R,
        //        double& shatr, double& DPB, double& dphi_l1_l2, double& gamma, double& dphi_blah,
        //        double& mdr, double& rpt);

        bool calculate_top_mass_variables(std::vector<const xAOD::TruthParticle*> leptons,
                std::vector<xAOD::Jet> sjets, std::vector<xAOD::Jet> bjets,
                const xAOD::MissingETContainer* met);
        bool calculate_top_variables(std::vector<const xAOD::TruthParticle*> leptons,
                std::vector<xAOD::Jet> sjets, std::vector<xAOD::Jet> bjets,
                const xAOD::MissingETContainer* met);
        void fill_tree(std::vector<const xAOD::TruthParticle*> leptons,
                std::vector<xAOD::Jet> jets, std::vector<xAOD::Jet> sjets, std::vector<xAOD::Jet> bjets,
                const xAOD::MissingETContainer* met);

        //double fracParticleInJet(const xAOD::TruthParticle* thePart, const xAOD::Jet* jet,
        //            bool DR, bool nparts);
        //void findJetConstituents(const xAOD::Jet* jet, std::set<const xAOD::TruthParticle*>& constituents, bool dR) const;
        

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
        bool m_use_bjet_eff;
        bool m_is_hh_signal;
        bool m_skip_maps;
        bool m_top_sample;

        double m_total_read_w;

        std::map<int, double> sumw_map;
        std::map<int, double> xsec_map;

        std::vector<const xAOD::TruthParticle*> electrons;
        std::vector<const xAOD::TruthParticle*> muons;
        xAOD::TruthParticleContainer* tp_electrons;

        // ROOT STUFF
        TFile* m_rfile;
        TTree* m_tree;

        ///////////////////////////////////////////
        // VARIABLES FOR OUTPUT TREE
        ///////////////////////////////////////////

        // weights
        std::vector<float> m_mcEventWeights;
        TH1D* h_cutflow_weighted;

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
        float m_met_sumet;
        float m_metphi_nonint;
        float m_met_nonint;

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
        float m_dphi_l0b1;
        float m_dr_l0b0;
        float m_dr_l0b1;

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

        float m_mt1;
        float m_mt1_scaled;

        float m_mt2_llbb;
        float m_mt2_bb;
        float m_mt2;

        // hh signal
        float m_hh;
        float m_hh_pt;
        float m_dphi_hh;
        float m_drhh;
        float m_h0_pt;
        float m_h1_pt;
        float m_h0_eta;
        float m_h1_eta;

        // THREE-BODY VARIABLES
        int m_3b_njets;
        int m_3b_nsjets;
        int m_3b_nbjets;
        int m_3b_nleptons;
        int m_3b_nmuons;
        int m_3b_nelectrons;
        float m_3b_mdr;
        float m_3b_dpb;
        float m_3b_cosThetaB;
        float m_3b_rpt;
        float m_3b_gamInvRp1;
        std::vector<float> m_3b_lepPt;
        std::vector<int> m_3b_lepQ;

        // top samples
        std::vector<float> m_truth_wpt; // pt of the l+v system(s)
        std::vector<float> m_truth_wpt_lv;
        std::vector<float> m_truth_wmass; // mass of the l+v system(s)
        std::vector<float> m_truth_wmass_lv; // mass of the l+v system(s)
        std::vector<int> m_truth_wq; // W charge
        std::vector<float> m_truth_topmass;
        std::vector<int> m_truth_topq;
        int m_n_wms;
        int m_n_wps;
        int m_n_ws;
        int m_n_tops;

        float m_ttbar_pt;
        float m_top_pt0;
        float m_top_pt1;
        float m_top_mass0;
        float m_top_mass1;
        float m_ttbar_mass;
        float m_ttbar_eta;
        float m_ttbar_dphi;
        float m_ttbar_dr;
        float m_ttbar_pt_out_01;
        float m_ttbar_pt_out_10;

        float m_dphi_ttbar_met;
        float m_dphi_ttbar_sj0;
        float m_dphi_ttbar_sj1;
        float m_pt_ratio_sj0_ttbar;
        float m_pt_ratio_sj1_ttbar;
        float m_pt_ratio_bj0_ttbar;
        float m_pt_ratio_bj1_ttbar;
        float m_pt_ratio_bb_ttbar;
        float m_pt_ratio_NonTTHT_ttbar_scalar;
        float m_pt_ratio_NonTTHT_ttbar;
        float m_dphi_NonTTHT_ttbar;
        

        // add'l ttbar jet activity
        float m_dphi_sj0_l0;
        float m_dphi_sj1_l0;
        float m_dphi_sj0_l1;
        float m_dphi_sj1_l1;
        float m_dphi_sj0_ll;
        float m_dphi_sj1_ll;


        int n_tree_fills;


}; // class

#endif
