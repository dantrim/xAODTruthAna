#ifndef DIHIGGSACCEPTANCE_H
#define DIHIGGSACCEPTANCE_H


//TruthAnalysis
#include "TruthSelectorBase.h"

//xAOD
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"

//std/stl
#include <map>
#include <vector>

//ROOT
#include "TH1D.h"
#include "TFile.h"
#include "TRandom3.h"

typedef std::map<int, int> Counts;
typedef std::map<int, float> CountsW;
typedef std::map<int, TH1D*> CountsH;

enum DiLepType {
    EE=0,
    MM,
    EM,
    ALL,
    Invalid
};

class DiHiggsAcceptance : public TruthSelectorBase
{
    public :
        DiHiggsAcceptance();
        virtual ~DiHiggsAcceptance(){};

        void set_x_mass(int mass);
        int x_mass() { return m_x_mass; }

        void set_lumi(float lumi);
        float lumi() { return m_lumi; }

        void initialize_sumw_map();
        void initialize_xsec_map();
        void load_mt1_windows();

        // weight to luminosity (m_lumi is numerical value in fb-1, convert to pb-1) 
        double w() { return ((m_weight * m_xsec * (m_lumi * 1000.)) / m_sumw); }

        void resonance_acceptance(std::vector<xAOD::TruthParticle*> leptons, std::vector<xAOD::Jet*> bjets,
                const xAOD::MissingETContainer* met, bool is_bjet_eff = false);
        void non_resonant_acceptance(std::vector<xAOD::TruthParticle*> leptons, std::vector<xAOD::Jet*> bjets,
                const xAOD::MissingETContainer* met, bool is_bjet_eff = false);

        float get_mt2_llbb(std::vector<xAOD::TruthParticle*> leptons, std::vector<xAOD::Jet*> bjets,
                const xAOD::MissingETContainer* met);
        float get_ht2_ratio(std::vector<xAOD::TruthParticle*> leptons, std::vector<xAOD::Jet*> bjets,
                const xAOD::MissingETContainer* met);
        float get_MT_1(std::vector<xAOD::TruthParticle*> leptons, std::vector<xAOD::Jet*> bjets,
                const xAOD::MissingETContainer* met, bool do_bjet_rescaling);
        float get_mt2_bb(std::vector<xAOD::Jet*> bjets, const xAOD::MissingETContainer* met);

        float get_dr_ll_bb(std::vector<xAOD::TruthParticle*> leptons, std::vector<xAOD::Jet*> bjets);

        void x_mass_window_selection(float& lower, float& upper, float dfactor = -1, float ufactor = -1);
        DiLepType get_lepton_type(const std::vector<xAOD::TruthParticle*> leptons);

        bool isBJet(const double eta, const int label);


        // TSelector overrides
        virtual void SlaveBegin(TTree* tree);
        virtual void Terminate();
        virtual Bool_t Process(Long64_t entry);

        void print_counts();

    private :
        int m_x_mass;
        float m_weight;
        float m_xsec;
        double m_sumw;
        double m_lumi; // luminosity in fb-1
        int n_events_non_dilepton;
        int n_events_non_dilepton_less;
        int n_events_non_dilepton_more;
        Counts total_counts;
        Counts passed_counts;
        CountsW total_w;
        CountsW passed_w;
        CountsH total_wh;
        CountsH passed_wh;

        Counts total_counts_beff;
        Counts passed_counts_beff;
        CountsW total_w_beff;
        CountsW passed_w_beff;
        CountsH total_wh_beff;
        CountsH passed_wh_beff;

        TFile* m_outfile;
        TRandom3* m_random;
        bool file_setup;
        double total_sumw;
        std::map<int, double> sumw_map;
        std::map<int, double> xsec_map;
        std::map<int, int> dsid_map;
        std::map<int, float> window_map_left;
        std::map<int, float> window_map_right;


    //ClassDef(DiHiggsAcceptance, 0);

};


#endif
