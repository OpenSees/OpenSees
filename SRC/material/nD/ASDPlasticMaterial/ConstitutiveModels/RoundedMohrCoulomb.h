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
#include "../YieldFunctions/RoundedMohrCoulomb_YF.h"

//Plastic flow directions
#include "../PlasticFlowDirections/VonMises_PF.h"

//Elasticity
#include "../ElasticityModels/LinearIsotropic3D_EL.h"

//Evolving variables
#include "../EvolvingVariables/LinearHardeningScalar_EV.h"
#include "../EvolvingVariables/LinearHardeningTensor_EV.h"



#include <classTags.h>


//Rounded Mohr-Coulomb plasticity with non-associated (VM) flow and isotropic hardening
class RoundedMohrCoulomb;  //This model we will define

//Typedefs for internal variables list, yield function, and plastic flow function
typedef MaterialInternalVariables < LinearHardeningTensor_EV, LinearHardeningScalar_EV> RMCVarsType;
typedef RoundedMohrCoulomb_YF <  LinearHardeningScalar_EV> RMC_YFType;
typedef VonMises_PF < LinearHardeningTensor_EV, LinearHardeningScalar_EV> RMC_PFType;

//Create a helpful typedef for the base class from which we will inherit to create the n
typedef ASDPlasticMaterial <LinearIsotropic3D_EL,
        RMC_YFType,
        RMC_PFType,
        RMCVarsType,
        ND_TAG_CEM_RoundedMohrCoulomb,
        RoundedMohrCoulomb
        > RMCBase;

//Define the new class. We must provide two constructor and the evolving variables as data-members.
class RoundedMohrCoulomb : public RMCBase
{
public:

    //First constructor, creates a material at its "ground state" from its parameters.
    RoundedMohrCoulomb(int tag_in, double rho_in, double E_in, double nu_in, double m_in, double qa_in, double pc_in, double e_in, double H_k, double eta0_in, double p0) ;




    // Second constructor is not called by the user, instead it is called when creating a copy of the
    // material. This must provide an initialization for the state variables and link the components
    // to these variables appropriately.
    RoundedMohrCoulomb(int tag_in, double rho, double p0, RMC_YFType & yf,
                       LinearIsotropic3D_EL & el,
                       RMC_PFType & pf,
                       RMCVarsType & vars);

    // Empty constructor for parallel
    RoundedMohrCoulomb() ;

    void Print(ostream& s, int flag);


    //The state variables.
private:
    LinearHardeningTensor_EV alpha; // Backstress
    LinearHardeningScalar_EV eta;     // RMC Yf opening
};

