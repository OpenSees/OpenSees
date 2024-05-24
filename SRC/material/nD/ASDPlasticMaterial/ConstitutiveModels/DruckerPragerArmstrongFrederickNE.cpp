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

#include "DruckerPragerArmstrongFrederickNE.h"

//First constructor, creates a material at its "ground state" from its parameters.
DruckerPragerArmstrongFrederickNE::DruckerPragerArmstrongFrederickNE(int tag_in, double k0_in, double ha_alpha, double cr_alpha, double H_k, double K_in, double pa_in, double n_in, double sigma3_max_in, double nu_in, double rho_, double p0) :
    DPAFNEBase::ASDPlasticMaterial(tag_in, rho_, p0,
            DPAFNE_YFType(alpha, k),       // Point YF to internal variables
            DuncanChang_EL(K_in, pa_in, n_in, nu_in, sigma3_max_in), // Create Elasticity
            DPAFNE_PFType(alpha, k),       // Point PF to the internal variables
            DPAFNEVarsType(alpha, k)),     // Declare the list of internal variables
    alpha(ha_alpha, cr_alpha),
    k(H_k, k0_in)
{

}

// Second constructor is not called by the user, instead it is called when creating a copy of the
// material. This must provide an initialization for the state variables and link the components
// to these variables appropriately.
DruckerPragerArmstrongFrederickNE::DruckerPragerArmstrongFrederickNE(int tag_in, double rho, double p0, DPAFNE_YFType &yf,
        DuncanChang_EL &el,
        DPAFNE_PFType &pf,
        DPAFNEVarsType &vars) :
    DPAFNEBase::ASDPlasticMaterial(tag_in, this->getRho(),
            p0,     //Sets p0
            DPAFNE_YFType(alpha, k),    // Point YF to new internal variables
            DuncanChang_EL(el), // Create Elasticity -- use copy constructor here
            DPAFNE_PFType(alpha, k),    // Point PF to the internal variables
            DPAFNEVarsType(alpha, k)),   // Declare the list of internal variables
    alpha(0, 0.0),
    k(0, 0)
{

}

DruckerPragerArmstrongFrederickNE::DruckerPragerArmstrongFrederickNE() :
    DPAFNEBase::ASDPlasticMaterial(0, 0, 0,
            DPAFNE_YFType(alpha, k),       // Point YF to internal variables
            DuncanChang_EL(0, 0, 0, 0, 0), // Create Elasticity
            DPAFNE_PFType(alpha, k),       // Point PF to the internal variables
            DPAFNEVarsType(alpha, k)),     // Declare the list of internal variables
    alpha(0, 0),
    k(0, 0)
{

}

//Checks whether predicted stress is less than zero, in which case sets stress to low confinement
//value and gives a reduced stiffness.
int DruckerPragerArmstrongFrederickNE::pre_integration_callback(const VoigtVector &depsilon,
        const VoigtVector &dsigma,
        const VoigtVector &TrialStress,
        const VoigtMatrix &Stiffness,
        double yf1,
        double yf2,
        bool & returns)
{
    using namespace ASDPlasticMaterialGlobals;
    static VoigtVector stress(3, 3, 0);
    static VoigtVector plasticstrain(3, 3, 0);
    static VoigtMatrix stiff(3, 3, 3, 3, 0);
    double p = -(TrialStress(0, 0) + TrialStress(1, 1) + TrialStress(2, 2)) / 3;
    if (p < 0)
    {
        stress *= 0;
        stiff *= 0;
        stress(0, 0) = -0.1;
        stress(1, 1) = -0.1;
        stress(2, 2) = -0.1;
        stiff = Stiffness;// / 10000;
        double mu = stiff(0, 1, 0, 1);
        double mu_reduced = mu / 10000;

        //Reduce only the shear components of the stiffness tensor.
        stiff( 0, 0, 0, 0 ) += -2 * mu + 2 * mu_reduced; // lambda + 2 * mu;
        // stiff( 0, 0, 1, 1 ) += -mu + mu_reduced;// lambda;
        // stiff( 0, 0, 2, 2 ) += -mu + mu_reduced;// lambda;
        stiff( 0, 1, 0, 1 ) += -mu + mu_reduced;// mu;
        stiff( 0, 1, 1, 0 ) += -mu + mu_reduced;// mu;
        stiff( 0, 2, 0, 2 ) += -mu + mu_reduced;// mu;
        stiff( 0, 2, 2, 0 ) += -mu + mu_reduced;// mu;
        stiff( 1, 0, 0, 1 ) += -mu + mu_reduced;// mu;
        stiff( 1, 0, 1, 0 ) += -mu + mu_reduced;// mu;
        // stiff( 1, 1, 0, 0 ) += -mu + mu_reduced;// lambda;
        stiff( 1, 1, 1, 1 ) += -2 * mu + 2 * mu_reduced; // lambda + 2 * mu;
        // stiff( 1, 1, 2, 2 ) += -mu + mu_reduced;// lambda;
        stiff( 1, 2, 1, 2 ) += -mu + mu_reduced;// mu;
        stiff( 1, 2, 2, 1 ) += -mu + mu_reduced;// mu;
        stiff( 2, 0, 0, 2 ) += -mu + mu_reduced;// mu;
        stiff( 2, 0, 2, 0 ) += -mu + mu_reduced;// mu;
        stiff( 2, 1, 1, 2 ) += -mu + mu_reduced;// mu;
        stiff( 2, 1, 2, 1 ) += -mu + mu_reduced;// mu;
        // stiff( 2, 2, 0, 0 ) += -mu + mu_reduced;// lambda;
        // stiff( 2, 2, 1, 1 ) += -mu + mu_reduced;// lambda;
        stiff( 2, 2, 2, 2 ) += -2 * mu + 2 * mu_reduced; // lambda + 2 * mu;


        plasticstrain = this->getCommittedPlasticStrainTensor();
        plasticstrain(i, j) = plasticstrain(i, j) + depsilon(i, j);

        this->setTrialStress(stress);
        this->setStiffness(stiff);
        this->setTrialPlastic_Strain(depsilon);

        //Reuse the 'stress' variable now to reset the backstress to zero
        stress *= 0;
        alpha.setVar(stress);
        alpha.commit_tmp();

        // cout << "pre_integration_callback acted!\n";

        returns = true;
    }
    else
    {
        returns = false;
    }
    return 0;
}

void
DruckerPragerArmstrongFrederickNE::Print(ostream& s, int flag)
{
    s << "DruckerPragerArmstrongFrederickNE" << endln;
    s << "\tTag: " << this->getTag() << endln;
    s << " Please Implement Me !!! " << endl;
}