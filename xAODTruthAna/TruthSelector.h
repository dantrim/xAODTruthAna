#ifndef xAODTruthAna_TruthSelector_H
#define xAODTruthAna_TruthSelector_H

//EventLoop
#include <EventLoop/Algorithm.h>

//ROOT
#include "TTree.h"

//std/stl
#include <string>
#include <vector>

class TruthSelector : public EL::Algorithm
{
    public :
        // variables that don't get filled at submission time should
        // be protected from being sent from the submission node to
        // the worker node (by using //!)

        bool verbose; //!

        // standard constructor
        TruthSelector();

        // Output Tree
        TTree *m_outTree; //!
        std::string *m_outFileName;

        unsigned long long m_br_runNumber; //!
        unsigned long long m_br_eventNumber; //!
        double m_br_eventWeight; //!
        std::vector<float> m_br_lpt; //!
        std::vector<float> m_br_leta; //!
        std::vector<float> m_br_lphi; //!

        // Algorithm
        virtual EL::StatusCode setupJob(EL::Job& job);
        virtual EL::StatusCode fileExecute();
        virtual EL::StatusCode histInitialize();
        virtual EL::StatusCode initialize();
        virtual EL::StatusCode execute();
        virtual EL::StatusCode postExecute();
        virtual EL::StatusCode finalize();
        virtual EL::StatusCode histFinalize();
      
    // ClassDef needed to distribute the algorithm to the workers
    ClassDef(TruthSelector, 1); 

}; // class 


#endif
