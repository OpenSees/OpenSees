/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// Original implementation: José Abell (UANDES), Massimo Petracca (ASDEA)
//
// ASDPlasticMaterial
//
// Fully general templated material class for plasticity modeling

#include "DruckerPragerVonMisesLinearHardening.h"

//First constructor, creates a material at its "ground state" from its parameters.
DruckerPragerVonMisesLinearHardening::DruckerPragerVonMisesLinearHardening(int tag_in, double k0_in, double H_alpha, double H_k, double E, double nu, double rho_, double p0) :
    DPVMLHBase::ASDPlasticMaterial(tag_in, rho_, p0,
            DPLH_YFType(alpha, k),       // Point YF to internal variables
            LinearIsotropic3D_EL(E, nu), // Create Elasticity
            DPDLHLH_PFType(alpha, k),       // Point PF to the internal variables
            DPLHVarsType(alpha, k)),     // Declare the list of internal variables
    alpha(H_alpha),
    k(H_k, k0_in)
{

}

// Second constructor is not called by the user, instead it is called when creating a copy of the
// material. This must provide an initialization for the state variables and link the components
// to these variables appropriately.
DruckerPragerVonMisesLinearHardening::DruckerPragerVonMisesLinearHardening(int tag_in, double rho, double p0,
        DPLH_YFType &yf,
        LinearIsotropic3D_EL &el,
        DPDLHLH_PFType &pf,
        DPLHVarsType &vars) :
    DPVMLHBase::ASDPlasticMaterial(tag_in, this->getRho(),
            p0,     //Sets p0
            DPLH_YFType(alpha, k),    // Point YF to new internal variables
            LinearIsotropic3D_EL(el), // Create Elasticity -- use copy constructor here
            DPDLHLH_PFType(alpha, k),    // Point PF to the internal variables
            DPLHVarsType(alpha, k)),   // Declare the list of internal variables
    alpha(0),
    k(0, 0)
{

}

DruckerPragerVonMisesLinearHardening::DruckerPragerVonMisesLinearHardening() :
    DPVMLHBase::ASDPlasticMaterial(0, 0, 0,
            DPLH_YFType(alpha, k),       // Point YF to internal variables
            LinearIsotropic3D_EL(0, 0), // Create Elasticity
            DPDLHLH_PFType(alpha, k),       // Point PF to the internal variables
            DPLHVarsType(alpha, k)),     // Declare the list of internal variables
    alpha(0),
    k(0, 0)
{

}

//Checks whether predicted stress is less than zero, in which case sets stress to low confinement
//value and gives a reduced stiffness.
int DruckerPragerVonMisesLinearHardening::pre_integration_callback(const VoigtVector &depsilon,
        const VoigtVector &dsigma,
        const VoigtVector &TrialStress,
        const VoigtMatrix &Stiffness,
        double yf1,
        double yf2,
        bool & returns)
{
    using namespace ASDPlasticMaterialGlobals;
    static VoigtVector str(3, 3, 0);
    static VoigtMatrix stiff(3, 3, 3, 3, 0);
    double p = -(TrialStress(0, 0) + TrialStress(1, 1) + TrialStress(2, 2)) / 3;
    if (p < 0)
    {
        str *= 0;
        stiff *= 0;
        str(0, 0) = -0.1;
        str(1, 1) = -0.1;
        str(2, 2) = -0.1;
        stiff = Stiffness;// / 10000;
        stiff *= 1. / 10000.;
        this->setTrialStress(str);
        this->setStiffness(stiff);

        returns = true;
    }
    else
    {
        returns = false;
    }
    return 0;
}

void
DruckerPragerVonMisesLinearHardening::Print(ostream& s, int flag)
{
    s << "DruckerPragerVonMisesLinearHardening" << endln;
    s << "\tTag: " << this->getTag() << endln;
    s << " Please Implement Me !!! " << endl;
}