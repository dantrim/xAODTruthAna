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

        // output tree
        unsigned long long m_br_runNumber;
        unsigned long long m_br_eventNumber;
        double m_br_eventWeight;


    protected :
        TFile *m_outTreeFile;
        TTree *m_outputTree;
        TTree *m_tree;
        xAOD::TEvent* m_event;
        int m_dbg;

    ClassDef(TruthSelector, 0);

}; // TruthSelector

#endif
