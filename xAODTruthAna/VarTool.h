#ifndef xAODTruthAna_VarTool_h
#define xAODTruthAna_VarTool_h

#include "TLorentzVector.h"

void SuperRazor(const TLorentzVector& l0, const TLorentzVector& l1, const TLorentzVector& met,
                double& shatr, double& DPB, double& dphi_l1_l2,
                double& MDR, double& cosThetaRp1);


#endif
