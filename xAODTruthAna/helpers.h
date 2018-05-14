#ifndef TRUTH_HELPERS_H
#define TRUTH_HELPERS_H

//xAOD
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODJet/Jet.h"

namespace helpers {

    const float MEV2GEV = 1e-3;

    struct is_pdg {
        int pdg;
        is_pdg(int p) : pdg(p) {}
        bool operator()(const xAOD::TruthParticle *l) { return (l && l->absPdgId() == pdg) ;}
    };

    struct pt_greater {
        bool operator()(const xAOD::TruthParticle* a, const xAOD::TruthParticle* b)
                { return a->pt() > b->pt(); }
    };

    struct pt_greaterJet {
        bool operator()(const xAOD::Jet* a, const xAOD::Jet* b) { return a->pt() > b->pt(); }
    };

    struct pt_greaterTLV {
        bool operator()(const TLorentzVector& a, const TLorentzVector& b) { return a.Pt() > b.Pt() ; }
    };

    bool ee(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1);
    bool mm(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1);
    bool em(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1);
    bool me(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1);
    bool sf(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1);
    bool df(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1);

};


#endif
