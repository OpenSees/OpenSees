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



// New materials are created by subclassing instances of the ASDPlasticMaterial<.,.,.,.,>
// template class, with the appropriate components as template parameters.
// Heavy use of templating is made, therefore typedeffing is a friend in helping clear up the mess.



//Von Mises Model with linear hardening (VMLH)
class VonMisesLinearHardening;  //This model we will define

//Model internal variables
using VMSL = VonMisesRadiusIV<ScalarLinearHardeningFunction>;
using BSTL = BackStressIV<TensorLinearHardeningFunction>;

//Select elasticity model
using EL = LinearIsotropic3D_EL;

//Select Yield-function model, using the internal variables
using YF = VonMises_YF<BSTL, VMSL>;

//Select Plastic-flow model, using the internal variables
using PF = VonMises_PF<BSTL, VMSL>;

using VMLHBase = ASDPlasticMaterial <EL,
        YF,
        PF,
        ND_TAG_ASDPlasticMaterial,
        VonMisesLinearHardening > ;


//Define the new class. We must provide two constructor and the evolving variables as data-members.
class VonMisesLinearHardening : public VMLHBase
{
public:
    VonMisesLinearHardening(int tag_in) :
        VMLHBase::ASDPlasticMaterial(tag_in) 
        {

        }
};

