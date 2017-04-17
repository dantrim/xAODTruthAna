#ifndef DIHIGGSWWBBSELECTOR_H
#define DIHIGGSWWBBSELECTOR_H

//TruthAnalysis
#include "TruthSelectorBase.h"

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
        void clear_variables();

        // yep
        void validation();

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
        double m_br_X_mass;
        double m_br_higgs_mass_lvlv;
        double m_br_higgs_mass_bjets;
        double m_br_X_mass_lvlvb;

        double m_br_w0_mass;
        double m_br_w1_mass;


    ClassDef(DiHiggsWWBBSelector, 0);

}; // DiHiggsWWBBSelector


#endif
