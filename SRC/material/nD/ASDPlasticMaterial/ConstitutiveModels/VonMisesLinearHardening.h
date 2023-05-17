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


#include "../ASDPlasticMaterial.h"



#include <classTags.h>
// New materials are created by subclassing instances of the ASDPlasticMaterial<.,.,.,.,>
// template class, with the appropriate components as template parameters.
// Heavy use of templating is made, therefore typedeffing is a friend in helping clear up the mess.



//Von Mises Model with linear hardening (VMLH)
class VonMisesLinearHardening;  //This model we will define

//Typedefs for internal variables list, yield function, and plastic flow function
// typedef MaterialInternalVariables < LinearHardeningTensor_EV, LinearHardeningScalar_EV> VMLHVarsType;
// typedef VonMises_YF < LinearHardeningTensor_EV, LinearHardeningScalar_EV> VMLH_YFType;
// typedef VonMises_PF < LinearHardeningTensor_EV, LinearHardeningScalar_EV> VMLH_PFType;


using VMSL = VonMisesRadiusIV<ScalarLinearHardeningFunction>;
using BSTL = BackStressIV<TensorLinearHardeningFunction>;
using VMLH_YFType = VonMises_YF<BSTL, VMSL>;
using VMLH_PFType = VonMises_PF<BSTL, VMSL>;

using VMLHBase = ASDPlasticMaterial <LinearIsotropic3D_EL,
        VMLH_YFType,
        VMLH_PFType,
        ND_TAG_ASDPlasticMaterial,
        VonMisesLinearHardening > ;




//Define the new class. We must provide two constructor and the evolving variables as data-members.
class VonMisesLinearHardening : public VMLHBase
{
public:
    VonMisesLinearHardening(int tag_in, double rho_) :
        VMLHBase::ASDPlasticMaterial(tag_in,rho_) 
        {
        	VMSL vmsl(1);
        	BSTL bstl(VoigtVector(0,0,0,0,0,0)); 

        	iv_storage.set(vmsl);
        	iv_storage.set(bstl);

		    YoungsModulus E(1);
		    PoissonsRatio nu(0.);
			ScalarLinearHardeningParameter HS(0.0);
			TensorLinearHardeningParameter HT(0.0);
			
		    parameters_storage.set(E);
		    parameters_storage.set(nu);
		    parameters_storage.set(HT);
		    parameters_storage.set(HS);

		    cout << "ASDPlasticMaterial" << endl;
	        cout << "  Yield Function          : " << yf.NAME << endl;
	        cout << "  Plastic flow direction  : " << pf.NAME << endl;
	        cout << "  Elasticity              : " << et.NAME << endl;
	        cout << "  # of Internal variables :" << iv_storage.size() <<  endl;
	        iv_storage.print_components();
	        cout << "  # of Parameters         :" << parameters_storage.size() <<  endl;
	        parameters_storage.print_components();
        }

};

