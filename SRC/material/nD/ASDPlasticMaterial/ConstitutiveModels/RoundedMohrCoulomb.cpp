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

#include "RoundedMohrCoulomb.h"

//First constructor, creates a material at its "ground state" from its parameters.
RoundedMohrCoulomb::RoundedMohrCoulomb(int tag_in, double rho_in, double E_in, double nu_in, double m_in, double qa_in, double pc_in, double e_in, double H_eta_in, double eta0_in, double p0) :
    RMCBase::ASDPlasticMaterial(tag_in, rho_in, p0,
                                          RMC_YFType(m_in, qa_in, pc_in, e_in, eta),    // Point YF to internal variables
                                          LinearIsotropic3D_EL(E_in, nu_in), // Create Elasticity
                                          RMC_PFType(alpha, eta),       // Point PF to the internal variables
                                          RMCVarsType(alpha, eta)),     // Declare the list of internal variables
    alpha(0.), // Set not to evolve. Dummy variable. But needed for the VM PF type which has two scalars
    eta(H_eta_in, eta0_in)
{

}

// Second constructor is not called by the user, instead it is called when creating a copy of the
// material. This must provide an initialization for the state variables and link the components
// to these variables appropriately.
RoundedMohrCoulomb::RoundedMohrCoulomb(int tag_in, double rho, double p0, RMC_YFType &yf,
                                       LinearIsotropic3D_EL &el,
                                       RMC_PFType &pf,
                                       RMCVarsType &vars) :
    RMCBase::ASDPlasticMaterial(tag_in, this->getRho(),
                                          p0,     //Sets p0
                                          RMC_YFType(0, 0, 0, 0, eta), // Point YF to new internal variables
                                          LinearIsotropic3D_EL(el), // Create Elasticity -- use copy constructor here
                                          RMC_PFType(alpha, eta),    // Point PF to the internal variables
                                          RMCVarsType(alpha, eta)),   // Declare the list of internal variables
    alpha(0.),
    eta(0.)
{

}

RoundedMohrCoulomb::RoundedMohrCoulomb() :
    RMCBase::ASDPlasticMaterial(0, 0, 0,
                                          RMC_YFType(0, 0, 0, 0, eta),   // Point YF to internal variables
                                          LinearIsotropic3D_EL(0, 0), // Create Elasticity
                                          RMC_PFType(alpha, eta),       // Point PF to the internal variables
                                          RMCVarsType(alpha, eta)),     // Declare the list of internal variables
    alpha(0),
    eta(0)
{

}


void
RoundedMohrCoulomb::Print(ostream& s, int flag)
{
    s << "RoundedMohrCoulomb" << endln;
    s << "\tTag: " << this->getTag() << endln;
    s << " eta =  " <<  eta << endl;
    RMCBase::Print(s);
}