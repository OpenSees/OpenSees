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
// template<typename Halpha, typename Hk>
// class VonMises;  //This model we will define

//Model internal variables
// template<typename Halpha>
// using alpha = BackStressIV<Halpha>;

// template<typename Hk>
// using k = VonMisesRadiusIV<Hk>;

// //Select elasticity model
// using EL = LinearIsotropic3D_EL;

// //Select Yield-function model, using the internal variables
// template<typename Halpha, typename Hk>
// using YF = VonMises_YF<alpha<Halpha>, k<Hk>>;

// //Select Plastic-flow model, using the internal variables
// template<typename Halpha, typename Hk>
// using PF = VonMises_PF<alpha<Halpha>, k<Hk>>;

// template<typename Halpha, typename Hk>
// using VMLHBase = ASDPlasticMaterial <EL,
//         YF<Halpha, Hk>,
//         PF<Halpha, Hk>,
//         ND_TAG_ASDPlasticMaterial
//         >;

// template<typename Halpha, typename Hk>
// using VonMises = VMLHBase<Halpha, Hk>;

// using VonMisesLinearHardening = VonMises<TensorLinearHardeningFunction, ScalarLinearHardeningFunction>;

// using VonMisesLinearHardening =
//     ASDPlasticMaterial <LinearIsotropic3D_EL,
//         VonMises_YF<
//             BackStressIV<TensorLinearHardeningFunction>, 
//             VonMisesRadiusIV<ScalarLinearHardeningFunction>
//             >,
//         VonMises_PF<
//             BackStressIV<TensorLinearHardeningFunction>, 
//             VonMisesRadiusIV<ScalarLinearHardeningFunction>
//             >,
//         ND_TAG_ASDPlasticMaterial
//         >;