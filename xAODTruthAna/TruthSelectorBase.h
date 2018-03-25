#ifndef TRUTHSELECTOR_BASE_H
#define TRUTHSELECTOR_BASE_H

//ROOT
#include "TSelector.h"
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
    class TStore;
}

const double mev2gev = 1.0e-3;

#undef CHECK
#define CHECK( ARG ) \
    do { \
        const bool result = ARG; \
        if( ! result ) { \
            ::Error("TruthSelectorBase",  "Failed to execute: \"%s\"", \
                #ARG ); \
            exit(1); \
        } \
    } while (false)

//struct IsPdg {
//    int pdg;
//    IsPdg(int p) : pdg(p) {}
//    bool operator()(const xAOD::TruthParticle *l) { return (l && l->absPdgId()==pdg); }
//};
//IsPdg isElectron(11), isMuon(13), isTau(15);
//
//bool is_ee(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
//    { return isElectron(l0) and isElectron(l1); }
//
//bool is_mm(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
//    { return isMuon(l0) and isMuon(l1); }
//
//bool is_em(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
//    { return ((isElectron(l0) and isMuon(l1)) or (isMuon(l0) and isElectron(l1))); }


class TruthSelectorBase : public TSelector
{

    public :
        TruthSelectorBase();
        virtual ~TruthSelectorBase();

        virtual void set_validation(bool doit) { m_validation = doit; }
        bool validation() { return m_validation; }

        // TSelector
        virtual Int_t Version() const { return 2; }
        virtual void Init(TTree *tree);
        virtual Bool_t Notify();
        virtual void Begin(TTree *tree);
        virtual void SlaveBegin(TTree *tree);
        virtual void Terminate();
        virtual Bool_t Process(Long64_t entry);

        // misc
        void set_debug(int dbg_level) { m_dbg = dbg_level; }
        int dbg() { return m_dbg; }

        void set_input_samplename(std::string name) { m_input_samplename = name; }
        void set_output_filename(std::string name) { m_output_filename = name; }
        std::string input_samplename() { return m_input_samplename; }
        std::string output_filename() { return m_output_filename; }

        TStopwatch* timer() { return &m_timer; }
        std::string timer_summary();

        xAOD::TEvent* event() { return m_event; }
        xAOD::TStore* store() { return m_tstore; }

        TFile *m_output_tree_file;
        TTree *m_output_tree;
        TFile* output_file() { return m_output_tree_file; }
        TTree* output_tree() { return m_output_tree; }

        // counters
        int n_evt_processed;
        int n_evt_stored;

        void save_output();

    private :

        int m_dbg;
        bool m_validation;

        //TFile *m_output_tree_file;
        //TTree *m_output_tree;
        TTree *m_tree;
        xAOD::TEvent* m_event;
        xAOD::TStore* m_tstore;


        TStopwatch m_timer;


        std::string m_input_samplename;
        std::string m_output_filename;

    ClassDef(TruthSelectorBase, 0);

}; // TruthSelectorBase

#endif
