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


#include "../ASDPlasticMaterialGlobals.h"
#include "../ASDPlasticMaterial.h"
#include "../MaterialInternalVariables.h"

//Yield Functions
#include "../YieldFunctions/CamClay_YF.h"

//Plastic flow directions
#include "../PlasticFlowDirections/CamClay_PF.h"

//Elasticity Models
#include "../ElasticityModels/CamClay_EL.h"

//Evolving variables
#include "../EvolvingVariables/CamClay_p0_EV.h"



#include <classTags.h>
// New materials are created by subclassing instances of the ASDPlasticMaterial<.,.,.,.,>
// template class, with the appropriate components as template parameters.
// Heavy use of templating is made, therefore typedeffing is a friend in helping clear up the mess.

//Von Mises Model with linear hardening (CC)
class CamClayLT;  //This model we will define

//Typedefs for internal variables list, yield function, and plastic flow function
typedef MaterialInternalVariables < CamClay_p0_EV> CCVarsType;
typedef CamClay_YF < CamClay_p0_EV> CC_YFType;
typedef CamClay_PF < CamClay_p0_EV> CC_PFType;

//Create a helpful typedef for the base class from which we will inherit to create the new material.
typedef ASDPlasticMaterial <CamClay_EL,
        CC_YFType,
        CC_PFType,
        CCVarsType,
        ND_TAG_CEM_CamClayLT,
        CamClayLT > CCBase;

//Define the new class. We must provide two constructor and the evolving variables as data-members.
class CamClayLT : public CCBase
{
public:

    //First constructor, creates a material at its "ground state" from its parameters.
    CamClayLT(int tag_in,
              double M_,
              double lambda_,
              double kappa_,
              double e0_,
              double p0_,
              double nu_,
              double initial_confinement,
              double rho_);






    // Second constructor is not called by the user, instead it is called when creating a copy of the
    // material. This must provide an initialization for the state variables and link the components
    // to these variables appropriately.
    CamClayLT(int tag_in, double rho,
              double pressure,
              CC_YFType &yf,
              CamClay_EL &el,
              CC_PFType &pf,
              CCVarsType &vars);

    CamClayLT();
    void Print(ostream& s, int flag);
    //The state variables.

private:
    CamClay_p0_EV p0;
};

