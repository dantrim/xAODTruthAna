// std
#include <iostream>
using namespace std;

#include "xAODTruthAna/VarTool.h"

// ROOT
#include "TVector3.h"

void SuperRazor(const TLorentzVector& l0_, const TLorentzVector& l1_, const TLorentzVector& met_,
                double& shatr, double& DPB, double& dphi_l1_l2,
                double& MDR, double& cosThetaRp1)
{
    TVector3 vBETA_z;
    TVector3 pT_CM;
    TVector3 vBETA_T_CMtoR;
    TVector3 vBETA_R;

    TLorentzVector l0 = l0_;
    TLorentzVector l1 = l1_;
    TLorentzVector met = met_;

    //
    // lab frame
    //
    // longitudinal boost
    vBETA_z = (1. / (l0.E() + l1.E())) * (l0 + l1).Vect();
    vBETA_z.SetX(0.0);
    vBETA_z.SetY(0.0);

    l0.Boost(-vBETA_z);
    l1.Boost(-vBETA_z);

    //pT of CM frame
    pT_CM = (l0 + l1).Vect() + met.Vect();
    pT_CM.SetZ(0.0);

    TLorentzVector ll = l0 + l1;

    shatr = sqrt(2.*(ll.E()*ll.E() - ll.Vect().Dot(pT_CM)
            + ll.E()*sqrt(ll.E()*ll.E() + pT_CM.Mag2() - 2.*ll.Vect().Dot(pT_CM))));

    vBETA_T_CMtoR = (1. / sqrt(pT_CM.Mag2() + shatr * shatr)) * pT_CM;

    l0.Boost(-vBETA_T_CMtoR);
    l1.Boost(-vBETA_T_CMtoR);
    ll.Boost(-vBETA_T_CMtoR);

    //
    //R-Frame
    //
    DPB = fabs((ll.Vect()).DeltaPhi(vBETA_T_CMtoR));
    dphi_l1_l2 = fabs(l0.Vect().DeltaPhi(l1.Vect()));
    vBETA_R = (1./(l0.E() + l1.E()))*(l0.Vect()-l1.Vect());

    l0.Boost(-vBETA_R);
    l1.Boost(vBETA_R);

    //
    //R+1 frame
    //
    MDR = 2.*l0.E();
    cosThetaRp1 = l0.Vect().Dot(vBETA_R) / (l0.Vect().Mag()*vBETA_R.Mag());

    return;
}
