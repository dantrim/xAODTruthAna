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

class WWbbTruth : public TruthSelectorBase
{

    public :
        WWbbTruth();
        virtual ~WWbbTruth(){};

        void do_sumw_calc(bool doit) { m_do_sumw_calc = doit; }
        void set_lumi(double l) { m_lumi = l; }
        void set_dsid(int d) { m_dsid = d; }

        // weight to luminosity (m_lumi is numerical value in fb-1, convert to pb-1)
        double w() { return ((m_weight * m_xsec * (m_lumi * 1000.)) / m_sumw); }

        void initialize_sumw_map();
        void initialize_xsec_map();


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
        bool m_outfile_setup;
        bool m_do_sumw_calc;
        double m_total_sumw;
        int m_dsid;
        double m_weight;
        double m_sumw;
        double m_xsec;
        double m_lumi; // in fb-1

        std::map<int, double> sumw_map;
        std::map<int, double> xsec_map;


}; // class

#endif
