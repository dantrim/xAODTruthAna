#ifndef xAODTruthAna_TruthSelector_h
#define xAODTruthAna_TruthSelector_h

//infrastructure
//#include "xAODRootAccess/Init.h"
//#include "xAODRootAccess/TEvent.h"
//#include "xAODRootAccess/TStore.h"
//#include "xAODEventInfo/EventInfo.h"
//#include "xAODEventFormat/EventFormat.h"
//#include "xAODCore/ShallowCopy.h"

//ROOT
#include "TSelector.h"
//#include "TTree.h"
//#include "TChain.h"
//#include "TTreeFormula.h"
#include "TFile.h"

//std/stl
#include <string>
#include <vector>

class TTree;
namespace xAOD {
    class TEvent;
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

        // TSelector methods
        virtual Int_t Version() const { return 2; }
        virtual void Init(TTree *tree);             
        virtual Bool_t Notify();
        virtual void Begin(TTree* tree);                       
        virtual void SlaveBegin(TTree *tree);
        virtual void Terminate();                   
        virtual Bool_t Process(Long64_t entry);     
        

        std::string m_inputSampleName;
        std::string m_outputFileName;

        
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


    protected :
        TFile *m_outTreeFile;
        TTree *m_outputTree;
        TTree *m_tree;
        xAOD::TEvent* m_event;
        int m_dbg;

    ClassDef(TruthSelector, 0);

}; // TruthSelector

#endif
