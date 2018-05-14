#ifndef HHFilterVal_H
#define HHFilterVal_H

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

class TFile;
class TH1F;


class HHFilterVal : public TruthSelectorBase
{

    public :
        HHFilterVal();
        virtual ~HHFilterVal(){};

        void  do_sumw_calc(bool doit) { m_do_sumw_calc = doit; }
        void set_lumi(double l) { m_lumi = l; }
        void set_dsid(int d) { m_dsid = d; }
        void set_suffix(std::string suffix) { m_suffix = suffix; }

        void setup_output_tree();

        // TSelector Overrides
        virtual void SlaveBegin(TTree* tree);
        virtual void Terminate();
        virtual Bool_t Process(Long64_t entry);

        struct HiggsCandidate {
            HiggsCandidate(std::string in_name) : name(in_name), filled(false) {
                parent = nullptr;
                c0 = nullptr;
                c1 = nullptr;
            }
            HiggsCandidate() : name(""), filled(false) {
                parent = nullptr;
                c0 = nullptr;
                c1 = nullptr;
            }
            std::string name;
            bool filled;
            const xAOD::TruthParticle* parent;
            const xAOD::TruthParticle* c0;
            const xAOD::TruthParticle* c1;
        };

        struct WCandidate {
            WCandidate(std::string in_name) : name(in_name), filled(false) {
                w = nullptr;
                lepton = nullptr;
                neutrino = nullptr;
            }
            WCandidate() : name(""), filled(false) {
                w = nullptr;
                lepton = nullptr;
                neutrino = nullptr;
            }
            std::string name;
            bool filled;
            const xAOD::TruthParticle* w;
            const xAOD::TruthParticle* lepton;
            const xAOD::TruthParticle* neutrino;
        };

        enum DileptonFlavor {
            EE=0,
            MM,
            TT,
            EM,
            ME,
            TE,
            ET,
            TM,
            MT,
            SF,
            DF,
            Invalid
        };

        bool get_w_decay(const xAOD::TruthParticle* parent, WCandidate& wboson);
        bool is_leptonic_w(HiggsCandidate& hc);
        DileptonFlavor get_dilepton_flavor(std::vector<HiggsCandidate> wbosons);
        const xAOD::TruthParticle* lepton_from_w(HiggsCandidate w);
        const xAOD::TruthParticle* neutrino_from_w(HiggsCandidate w);
        DileptonFlavor get_dilepton_flavor_from_leptons(const xAOD::TruthParticle* part0,
                const xAOD::TruthParticle* part1);

        bool fill_histograms(const HiggsCandidate& hww, const HiggsCandidate& hbb,
            const std::vector<WCandidate>& wbosons);
        bool save_histograms();
    private :
        std::string m_outfilename;
        std::string m_suffix;
        bool m_outfile_setup;

        bool m_do_sumw_calc;
        double m_lumi; // luminosity if fb-1
        int m_dsid;
        int n_bad_higgs_decays; // unexpected higgs decays
        int n_bb_higgs_decays;
        int n_ww_higgs_decays; 

        int n_dilepton_flavor [DileptonFlavor::Invalid];

        // ROOT STUFF
        TFile* m_rfile;

        // histograms
        TH1F* h_n_h;
        TH1F* h_mh0;
        TH1F* h_mh1;
        TH1F* h_mhh;
        TH1F* h_h0_pt;
        TH1F* h_h1_pt;
        TH1F* h_hh_pt;

        TH1F* h_n_ww;
        TH1F* h_n_w;
        TH1F* h_w0_pt;
        TH1F* h_w1_pt;
        TH1F* h_mw0;
        TH1F* h_mw1;
        TH1F* h_mww;

        TH1F* h_n_bb;
        TH1F* h_b0_pt;
        TH1F* h_b1_pt;
        TH1F* h_mbb;

        TH1F* h_n_leptons;
        TH1F* h_lepton_flavor;
        TH1F* h_lepton_flavor_sfdf;
        TH1F* h_l0_pt;
        TH1F* h_l1_pt;
        TH1F* h_l0_eta;
        TH1F* h_l1_eta;




}; // class HHFilterVal

#endif
