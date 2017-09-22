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

namespace res {

enum DiLepType {
    EE=0,
    MM,
    EM,
    ALL,
    Invalid
};

class DiHiggsResonance : public TruthSelectorBase
{
    public :
        DiHiggsResonance();
        virtual ~DiHiggsResonance(){};

        void set_x_mass(int mass);
        int x_mass() { return m_x_mass; }

        void set_lumi(float lumi);
        float lumi() { return m_lumi; }

        void initialize_sumw_map();
        void initialize_xsec_map();
        void initialize_histos();

        void look_at_resonance(std::vector<xAOD::TruthParticle*> leptons, std::vector<xAOD::Jet*> bjets,
                const xAOD::MissingETContainer* met);

        // weight to luminosity (m_lumi is numerical value in fb-1, convert to pb-1) 
        double w() { return ((m_weight * m_xsec * (m_lumi * 1000.)) / m_sumw); }

        float get_MT_1(std::vector<xAOD::TruthParticle*> leptons, std::vector<xAOD::Jet*> bjets,
            const xAOD::MissingETContainer* met, bool do_bjet_rescaling);

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
        double m_lumi; // luminosity in fb-1
        TFile* m_outfile;
        bool file_setup;
        double total_sumw;
        std::map<int, double> sumw_map;
        std::map<int, double> xsec_map;
        std::map<int, int> dsid_map;

        TH1F* h_mt1_scaled;
        TH1F* h_mt1;


    //ClassDef(DiHiggsResonance, 0);

};
};

#endif
