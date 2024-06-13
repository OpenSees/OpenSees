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

#include "VonMisesLinearHardening.h"
#include "NDMaterial.h"
#include <iostream>
#include "../../../ltensor/LTensor.h"
#include "../ASDPlasticMaterialGlobals.h"
//First constructor, creates a material at its "ground state" from its parameters.
VonMisesLinearHardening::VonMisesLinearHardening(int tag_in, double k0_in, double H_alpha, double H_k, double E, double nu, double rho_) :
    VMLHBase::ASDPlasticMaterial(tag_in, rho_, 0.0, //Initial confinement can be 0 for this model
                                           VMLH_YFType(alpha, k),       // Point YF to internal variables
                                           LinearIsotropic3D_EL(E, nu), // Create Elasticity
                                           VMLH_PFType(alpha, k),       // Point PF to the internal variables
                                           VMLHVarsType(alpha, k)),     // Declare the list of internal variables
    alpha(H_alpha),
    k(H_k, k0_in)
{
    // cout << "First VMLH ctor" << endl;
    // cout << "tag_in  = " << tag_in << endl;
    // cout << "k0_in  = " << k0_in << endl;
    // cout << "H_alpha  = " << H_alpha << endl;
    // cout << "H_k  = " << H_k << endl;
    // cout << "E  = " << E << endl;
    // cout << "nu  = " << nu << endl;
    // cout << "rho_  = " << rho_ << endl;
    // cout << "yf  = " << yf(getStressTensor()) << endl;

}

// Second constructor is not called by the user, instead it is called when creating a copy of the
// material. This must provide an initialization for the state variables and link the components
// to these variables appropriately.
VonMisesLinearHardening::VonMisesLinearHardening(int tag_in, double p0, double rho, VMLH_YFType &yf,
        LinearIsotropic3D_EL &el,
        VMLH_PFType &pf,
        VMLHVarsType &vars) :
    VMLHBase::ASDPlasticMaterial(tag_in, this->getRho(), p0, // Initial confinement can be 0 for this model
                                           VMLH_YFType(alpha, k),     // Point YF to new internal variables
                                           LinearIsotropic3D_EL(el),  // Create Elasticity -- use copy constructor here
                                           VMLH_PFType(alpha, k),     // Point PF to the internal variables
                                           VMLHVarsType(alpha, k)),   // Declare the list of internal variables
    alpha(0),
    k(0, 0)
{
}

VonMisesLinearHardening::VonMisesLinearHardening() :
    VMLHBase::ASDPlasticMaterial(0, 0, 0.0, //Initial confinement can be 0 for this model
                                           VMLH_YFType(alpha, k),       // Point YF to internal variables
                                           LinearIsotropic3D_EL(0, 0),  // Create Elasticity
                                           VMLH_PFType(alpha, k),       // Point PF to the internal variables
                                           VMLHVarsType(alpha, k)),     // Declare the list of internal variables
    alpha(0),
    k(0, 0)
{}


void
VonMisesLinearHardening::Print(ostream& s, int flag)
{
    s << "VonMisesLinearHardening::" << endln;
    s << "\tTag: " << this->getTag() << endln;
    s << " Please Implement Me !!! " << endl;
    // s << "\tElastic_Modulus: " << E << endln;
    // s << "\tPoissons_Ratio: " << v << endln;
    // s << "\tDensity: " << rho << endln;
    // s << "\tVon_Mises_radius: "
    // s << "\tKinematic_hardening_rate: "
    // s << "\tIsotropic_hardening_rate: "
}

