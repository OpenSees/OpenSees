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
#include "../YieldFunctions/VonMises_YF.h"

//Plastic flow directions
#include "../PlasticFlowDirections/VonMises_PF.h"

//Elasticity Models
#include "../ElasticityModels/LinearIsotropic3D_EL.h"

//Evolving variables
#include "../EvolvingVariables/ArmstrongFrederickTensor_EV.h"
#include "../EvolvingVariables/LinearHardeningScalar_EV.h"



#include <classTags.h>
// New materials are created by subclassing instances of the ASDPlasticMaterial<.,.,.,.,>
// template class, with the appropriate components as template parameters.
// Heavy use of templating is made, therefore typedeffing is a friend in helping clear up the mess.

//Von Mises Model with linear hardening (VMAF)
class VonMisesArmstrongFrederick;  //This model we will define

//Typedefs for internal variables list, yield function, and plastic flow function
typedef MaterialInternalVariables < ArmstrongFrederickTensor_EV, LinearHardeningScalar_EV> VMAFVarsType;
typedef VonMises_YF < ArmstrongFrederickTensor_EV, LinearHardeningScalar_EV> VMAF_YFType;
typedef VonMises_PF < ArmstrongFrederickTensor_EV, LinearHardeningScalar_EV> VMAF_PFType;

//Create a helpful typedef for the base class from which we will inherit to create the new material.
typedef ASDPlasticMaterial <LinearIsotropic3D_EL,
        VMAF_YFType,
        VMAF_PFType,
        VMAFVarsType,
        ND_TAG_CEM_VonMisesArmstrongFrederick,
        VonMisesArmstrongFrederick > VMAFBase;

//Define the new class. We must provide two constructor and the evolving variables as data-members.
class VonMisesArmstrongFrederick : public VMAFBase
{
public:

    //First constructor, creates a material at its "ground state" from its parameters.
    VonMisesArmstrongFrederick(int tag_in, double k0_in, double ha_alpha, double cr_alpha, double H_k, double E, double nu, double rho_);

    // Second constructor is not called by the user, instead it is called when creating a copy of the
    // material. This must provide an initialization for the state variables and link the components
    // to these variables appropriately.
    VonMisesArmstrongFrederick(int tag_in, double rho, double p0, VMAF_YFType &yf,
                            LinearIsotropic3D_EL &el,
                            VMAF_PFType &pf,
                            VMAFVarsType &vars);

    VonMisesArmstrongFrederick();
    //The state variables.

private:
    ArmstrongFrederickTensor_EV alpha; // Backstress
    LinearHardeningScalar_EV k;     // Critical stress ratio (k = M under this formulation)

};

