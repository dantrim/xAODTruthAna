#include "xAODTruthAna/helpers.h"

bool helpers::ee(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
{
    helpers::is_pdg is_ele(11);
    return ( is_ele(l0) && is_ele(l1) );
}

bool helpers::mm(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
{
    helpers::is_pdg is_muo(13);
    return ( is_muo(l0) && is_muo(l1) );
}

bool helpers::em(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
{
    helpers::is_pdg is_ele(11);
    helpers::is_pdg is_muo(13);
    return ( is_ele(l0) && is_muo(l1) );
}

bool helpers::me(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
{
    helpers::is_pdg is_ele(11);
    helpers::is_pdg is_muo(13); 
    return ( is_muo(l0) && is_ele(l1) );
}

bool helpers::sf(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
{
    helpers::is_pdg is_ele(11);
    helpers::is_pdg is_muo(13);
    return ( (is_ele(l0) && is_ele(l1)) || (is_muo(l0) && is_muo(l1)) );
}

bool helpers::df(const xAOD::TruthParticle* l0, const xAOD::TruthParticle* l1)
{
    helpers::is_pdg is_ele(11);
    helpers::is_pdg is_muo(13);
    return ( (is_ele(l0) && is_muo(l1)) || (is_muo(l0) && is_ele(l1)) );
}
