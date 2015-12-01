#ifndef xAODTruthAna_TruthSelector_h
#define xAODTruthAna_TruthSelector_h

//ROOT
#include "TSelector.h"
//#include "TTree.h"
//#include "TChain.h"
//#include "TTreeFormula.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"

//std/stl
#include <string>
#include <vector>
#include "xAODTruth/TruthParticle.h"

class TTree;
namespace xAOD {
    class TEvent;
    //class TruthParticle;
}

class TruthSelector : public TSelector
{

    public :
        TruthSelector();
        virtual ~TruthSelector();

        // configure the run
        void initializeOutputTree();
        void initializeTreeBranches();
        void clearContainers();
        void saveOutputTree();
        void setDebug(int dbgLevel) { m_dbg = dbgLevel; }
        void setWZ2WW(bool doit) { m_wz2ww = doit; }

        // TSelector methods
        virtual Int_t Version() const { return 2; }
        virtual void Init(TTree *tree);             
        virtual Bool_t Notify();
        virtual void Begin(TTree* tree);                       
        virtual void SlaveBegin(TTree *tree);
        virtual void Terminate();                   
        virtual Bool_t Process(Long64_t entry); 

        xAOD::TruthParticle* getLeptonForMet(std::vector<xAOD::TruthParticle*>& v_mu,
                                                std::vector<xAOD::TruthParticle*>& v_el);


        std::string m_inputSampleName;
        std::string m_outputFileName;

        std::string timerSummary();

        
        //----------------------------------------
        // Variables for the output tree
        //----------------------------------------

        // event vars
        unsigned long long m_br_runNumber;
        unsigned long long m_br_eventNumber;
        double m_br_eventWeight;
        std::vector<double> m_br_eventWeights;

        // lepton vars
        int m_br_lepton_n;
        std::vector<double> m_br_lepton_pt;
        std::vector<double> m_br_lepton_eta;
        std::vector<double> m_br_lepton_phi;
        std::vector<int> m_br_lepton_pdgId;
        std::vector<int> m_br_lepton_charge;
        std::vector<int> m_br_lepton_type;
        std::vector<int> m_br_lepton_origin;
        int m_br_isSF;
        int m_br_isDF;
        int m_br_isEE;
        int m_br_isMM;
        int m_br_isEM;

        // jet vars
        int m_br_jet_n;
        std::vector<double> m_br_jet_pt;
        std::vector<double> m_br_jet_eta;
        std::vector<double> m_br_jet_phi;
        std::vector<double> m_br_jet_m;
        std::vector<int>    m_br_jet_flav;

        // light jet vars
        int m_br_ljet_n;
        std::vector<double> m_br_ljet_pt;
        std::vector<double> m_br_ljet_eta;
        std::vector<double> m_br_ljet_phi;
        std::vector<double> m_br_ljet_m;
        std::vector<int>    m_br_ljet_flav;

        // c jet vars
        int m_br_cjet_n;
        std::vector<double> m_br_cjet_pt;
        std::vector<double> m_br_cjet_eta;
        std::vector<double> m_br_cjet_phi;
        std::vector<double> m_br_cjet_m;
        std::vector<int>    m_br_cjet_flav;

        // b jet vars
        int m_br_bjet_n;
        std::vector<double> m_br_bjet_pt;
        std::vector<double> m_br_bjet_eta;
        std::vector<double> m_br_bjet_phi;
        std::vector<double> m_br_bjet_m;
        std::vector<int>    m_br_bjet_flav;

        // heavy-flavor jet vars
        int m_br_hjet_n;
        std::vector<double> m_br_hjet_pt;
        std::vector<double> m_br_hjet_eta;
        std::vector<double> m_br_hjet_phi;
        std::vector<double> m_br_hjet_m;
        std::vector<int>    m_br_hjet_flav;

        // gluon jet vars
        int m_br_gjet_n;
        std::vector<double> m_br_gjet_pt;
        std::vector<double> m_br_gjet_eta;
        std::vector<double> m_br_gjet_phi;
        std::vector<double> m_br_gjet_m;
        std::vector<int>    m_br_gjet_flav;

        // met vars
        double m_br_met_et;
        double m_br_met_phi;

        // mt2
        double m_br_mt2;

        // meff
        double m_br_meff;
        double m_br_meff_S2L;

        // dilepton
        double m_br_dphiLL;
        double m_br_ptLL;
        double m_br_drLL;
        double m_br_mll;

        // pbLL
        double m_br_pbLL;
        double m_br_dphi_met_pbLL;

        // r1
        double m_br_r1;
        double m_br_r1_S2L;

        // r2
        double m_br_r2;

        // jet lepton vars
        double m_br_dr_l0_j0;
        double m_br_dr_l1_j0;
        double m_br_dr_l0_bj0;
        double m_br_dr_l0_lj0;

        double m_br_dphi_l0_j0;
        double m_br_dphi_l1_j0;
        double m_br_dphi_l0_bj0;
        double m_br_dphi_l0_lj0;

        double m_br_dphi_ptLL_j0;
        double m_br_dphi_ptLL_bj0;
        double m_br_dphi_ptLL_lj0;

        // met lepton vars
        double m_br_dphi_l0_met;
        double m_br_dphi_l1_met;

        double m_br_dphi_ptLL_met;

        double m_br_Rsib;

        // met jet vars
        double m_br_dphi_met_j0;
        double m_br_dphi_met_bj0;
        double m_br_dphi_met_lj0;

        // WW-like vars
        double m_br_deltaX;
        double m_br_cosThetaB;

        // super-razor
        double m_br_shatr;
        double m_br_DPB;
        double m_br_dphi_l1_l2;
        double m_br_MDR;
        double m_br_cosThetaRp1;





    protected :
        TFile *m_outTreeFile;
        TTree *m_outputTree;
        TTree *m_tree;
        xAOD::TEvent* m_event;
        int m_dbg;

        bool m_wz2ww;

        int n_evtProcessed;
        int n_evtStored;

        TStopwatch m_timer;

    ClassDef(TruthSelector, 0);

}; // TruthSelector

#endif
