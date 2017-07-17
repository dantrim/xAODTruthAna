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

        void initialize_sumw_map();
        void initialize_xsec_map();

        // weight normalized to 100/fb (xsec is in pb)
        double w() { return ((m_weight * m_xsec * 100000.) / m_sumw); }

        void resonance_acceptance(std::vector<xAOD::TruthParticle*> leptons, std::vector<xAOD::Jet*> bjets,
                const xAOD::MissingETContainer* met);
        void non_resonant_acceptance();

        float get_mt2_llbb(std::vector<xAOD::TruthParticle*> leptons, std::vector<xAOD::Jet*> bjets,
                const xAOD::MissingETContainer* met);
        float get_ht2_ratio(std::vector<xAOD::TruthParticle*> leptons, std::vector<xAOD::Jet*> bjets,
                const xAOD::MissingETContainer* met);
        float get_MT_1(std::vector<xAOD::TruthParticle*> leptons, std::vector<xAOD::Jet*> bjets,
                const xAOD::MissingETContainer* met);

        void x_mass_window_selection(float& lower, float& upper, float dfactor = 0.9, float ufactor = 1.1);
        DiLepType get_lepton_type(const std::vector<xAOD::TruthParticle*> leptons);


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
        int n_events_non_dilepton;
        int n_events_non_dilepton_less;
        int n_events_non_dilepton_more;
        Counts total_counts;
        Counts passed_counts;
        CountsW total_w;
        CountsW passed_w;
        TFile* m_outfile;
        CountsH total_wh;
        CountsH passed_wh;
        bool file_setup;
        double total_sumw;
        std::map<int, double> sumw_map;
        std::map<int, double> xsec_map;
        std::map<int, int> dsid_map;


    //ClassDef(DiHiggsAcceptance, 0);

};


#endif
