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
#include "../YieldFunctions/DruckerPrager_YF.h"

//Plastic flow directions
#include "../PlasticFlowDirections/DruckerPrager_PF.h"

//Elasticity Models
#include "../ElasticityModels/NoTensionLinearIsotropic3D_EL.h"

//Evolving variables
#include "../EvolvingVariables/LinearHardeningTensor_EV.h"
#include "../EvolvingVariables/LinearHardeningScalar_EV.h"



#include <classTags.h>



// New materials are created by subclassing instances of the ASDPlasticMaterial<.,.,.,.,>
// template class, with the appropriate components as template parameters.
// Heavy use of templating is made, therefore typedeffing is a friend in helping clear up the mess.

//Drucker Prager Model with linear hardening (DPLH)
class DruckerPragerLinearHardening;  //This model we will define

//Activate pre_integration_callback to handle tension case
template< >
struct supports_pre_integration_callback<DruckerPragerLinearHardening>
{
    static const bool value = true;
};

//Typedefs for internal variables list, yield function, and plastic flow function
typedef MaterialInternalVariables < LinearHardeningTensor_EV, LinearHardeningScalar_EV> DPLHVarsType;
typedef DruckerPrager_YF < LinearHardeningTensor_EV, LinearHardeningScalar_EV> DPLH_YFType;
typedef DruckerPrager_PF < LinearHardeningTensor_EV, LinearHardeningScalar_EV> DPLH_PFType;

//Create a helpful typedef for the base class from which we will inherit to create the new material.
typedef ASDPlasticMaterial <NoTensionLinearIsotropic3D_EL,
        DPLH_YFType,
        DPLH_PFType,
        DPLHVarsType,
        ND_TAG_CEM_DruckerPragerLinearHardening,
        DruckerPragerLinearHardening > DPLHBase;

//Define the new class. We must provide two constructor and the evolving variables as data-members.
class DruckerPragerLinearHardening : public DPLHBase
{
public:

    //First constructor, creates a material at its "ground state" from its parameters.
    DruckerPragerLinearHardening(int tag_in, double k0_in, double H_alpha, double H_k, double E, double nu, double rho_, double p0) ;

    // Second constructor is not called by the user, instead it is called when creating a copy of the
    // material. This must provide an initialization for the state variables and link the components
    // to these variables appropriately.
    DruckerPragerLinearHardening(int tag_in, double rho, double p0, DPLH_YFType & yf,
                                 NoTensionLinearIsotropic3D_EL & el,
                                 DPLH_PFType & pf,
                                 DPLHVarsType & vars);

    // Empty constructor for parallel
    DruckerPragerLinearHardening() ;

    int pre_integration_callback(const VoigtVector&, const VoigtVector&, const VoigtVector&, const VoigtMatrix&, double, double, bool&);

    void Print(ostream& s, int flag);


    //The state variables.
private:
    LinearHardeningTensor_EV alpha; // Backstress
    LinearHardeningScalar_EV k;     // Critical stress ratio (k = M under this formulation)
    // int HARDENING_TYPE; //Determine the hardening type according to hardening rate. This is useful for consistent stiffness tensor. 
};

