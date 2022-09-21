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

#include "CamClayLT.h"
#include "NDMaterial.h"
#include <iostream>
#include "../../../ltensor/LTensor.h"
#include "../ASDPlasticMaterialGlobals.h"

//First constructor, creates a material at its "ground state" from its parameters.
CamClayLT::CamClayLT(int tag_in,
                     double M_,
                     double lambda_,
                     double kappa_,
                     double e0_,
                     double p0_,
                     double nu_,
                     double initial_confinement,
                     double rho_) :
    CCBase::ASDPlasticMaterial(tag_in, rho_, initial_confinement, //Initial confinement
                                         CC_YFType(M_, p0),       // Point YF to internal variables
                                         CamClay_EL(e0_, kappa_, nu_), // Create Elasticity
                                         CC_PFType(M_, p0),       // Point PF to the internal variables
                                         CCVarsType(p0)),     // Declare the list of internal variables
    p0( M_,  lambda_,  kappa_,  e0_,  p0_)
{


}

// Second constructor is not called by the user, instead it is called when creating a copy of the
// material. This must provide an initialization for the state variables and link the components
// to these variables appropriately.
CamClayLT::CamClayLT(int tag_in, double rho,
                     double pressure,
                     CC_YFType &yf,
                     CamClay_EL &el,
                     CC_PFType &pf,
                     CCVarsType &vars) :
    CCBase::ASDPlasticMaterial(tag_in, this->getRho(), pressure, // Initial confinement
                                         CC_YFType(0, p0),       // Point YF to internal variables
                                         CamClay_EL(0, 0, 0), // Create Elasticity
                                         CC_PFType(0, p0),       // Point PF to the internal variables
                                         CCVarsType(p0)),     // Declare the list of internal variables
    p0( 0,  0,  0,  0,  pressure)
{
}

CamClayLT::CamClayLT() :
    CCBase::ASDPlasticMaterial(0, 0, 0.0, //Initial confinement
                                         CC_YFType(0, p0),       // Point YF to internal variables
                                         CamClay_EL(0, 0, 0),  // Create Elasticity
                                         CC_PFType(0, p0),       // Point PF to the internal variables
                                         CCVarsType(p0)),     // Declare the list of internal variables
    p0(0, 0, 0, 0, 0)
{}


void
CamClayLT::Print(ostream& s, int flag)
{
    const double &p0_ = p0.getVariableConstReference();
    s << "CamClayLT:  p0 = " << p0_ << endl;
    // s << "CamClayLT::" << endln;
    // s << "\tTag: " << this->getTag() << endln;
    // s << " Please Implement Me !!! " << endl;
    // s << "\tElastic_Modulus: " << E << endln;
    // s << "\tPoissons_Ratio: " << v << endln;
    // s << "\tDensity: " << rho << endln;
    // s << "\tVon_Mises_radius: "
    // s << "\tKinematic_hardening_rate: "
    // s << "\tIsotropic_hardening_rate: "
}

