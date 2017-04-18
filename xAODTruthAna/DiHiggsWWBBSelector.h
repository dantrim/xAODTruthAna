#ifndef DIHIGGSWWBBSELECTOR_H
#define DIHIGGSWWBBSELECTOR_H

//TruthAnalysis
#include "TruthSelectorBase.h"
#include "xAODJet/Jet.h"

class DiHiggsWWBBSelector : public TruthSelectorBase
{
    public :
        DiHiggsWWBBSelector();
        virtual ~DiHiggsWWBBSelector();

        // TSelector
        virtual void SlaveBegin(TTree *tree);
        virtual void Terminate();
        virtual Bool_t Process(Long64_t entry);

        void initialize_output_tree(bool is_validation = false);
        void initialize_tree_branches(bool is_validation = false);
        void init_ana_tree();
        void init_validation_tree();
        void clear_variables();

        // variable calculation
        double cosThetaB(std::vector<xAOD::TruthParticle*> leptons);
        double cosTheta1(std::vector<xAOD::TruthParticle*> leptons);
        double cosTheta2(std::vector<xAOD::TruthParticle*> leptons, TLorentzVector met);
        double dphi_boost_ll(std::vector<xAOD::TruthParticle*> leptons);
        double dphi_boost_ll_met(std::vector<xAOD::TruthParticle*> leptons, TLorentzVector met);
        double dphi_boost_CM_ll_met(std::vector<xAOD::TruthParticle*> leptons, TLorentzVector met,
                                        std::vector<xAOD::Jet*> bjets);
    

        // yep
        void validation_observables();
        Bool_t get_observables();


    private :
        bool m_output_setup;


        //////////////////////////////////////////////
        // variables for output ntuple
        //////////////////////////////////////////////

        // validation
        int m_br_is_wwbb;
        int m_br_n_leptons;
        double m_br_higgs_mass_w;
        double m_br_higgs_mass_b;
        std::vector<double> m_br_higgs_pt;
        double m_br_X_mass;
        double m_br_higgs_mass_lvlv;
        double m_br_higgs_mass_bjets;
        double m_br_X_mass_lvlvb;

        double m_br_w0_mass;
        double m_br_w1_mass;

        double m_br_dphi_ll;
        double m_br_pTll;
        int m_br_dilepton_flavor;
        std::vector<int> m_br_l_pdgId;
        std::vector<double> m_br_l_pt;
        std::vector<double> m_br_l_eta;

        double m_br_dphi_bb;
        double m_br_pTbb;
        std::vector<double> m_br_b_pt;
        std::vector<double> m_br_b_eta;

        double m_br_dphi_ll_bb;
        double m_br_dR_bb;


        //-----------------------------------//
        // non-validation
        //-----------------------------------//
        std::vector<int> m_br_l_q;

        // dilepton
        double m_br_dR_ll;
        double m_br_cosThetaB;
        double m_br_cosTheta1;
        double m_br_cosTheta2;
        double m_br_dphi_boost_ll;

        // jets
        int m_br_jet_n;
        int m_br_ljet_n;
        int m_br_bjet_n;
        int m_br_cjet_n;
        int m_br_hjet_n;
        int m_br_gjet_n;

        // jet kin
        std::vector<double> m_br_j_pt;
        std::vector<double> m_br_j_eta;

        // light-jets kin
        std::vector<double> m_br_lj_pt;
        std::vector<double> m_br_lj_eta;

        // b-jets kin
        std::vector<double> m_br_bj_pt;
        std::vector<double> m_br_bj_eta;
        double m_br_mass_bb;

        // met
        double m_br_met_et;
        double m_br_met_phi;

        // lepton + met
        double m_br_dphi_met_ll;
        double m_br_mass_met_ll;
        double m_br_met_pTll;
        double m_br_dphi_boost_ll_met;

        // lepton + b
        double m_br_dphi_met_bb;

        // WWbb
        double m_br_dphi_WW_bb;
        double m_br_dR_ll_bb;
        double m_br_mass_X;
        double m_br_mass_X_scaled;
        double m_br_CMboostMag;
        double m_br_dphi_boost_CM_ll_met;

        double m_br_HT2;
        double m_br_HT2Ratio;

        double m_br_HT2_noMET;
        double m_br_HT2Ratio_noMET;

    ClassDef(DiHiggsWWBBSelector, 0);

}; // DiHiggsWWBBSelector


#endif
