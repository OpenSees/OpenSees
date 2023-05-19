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

// #include "ConstitutiveModels/VonMisesLinearHardening.h"
// #include "ConstitutiveModels/VonMisesFamily.h"
#include <classTags.h>
#include "ASDPlasticMaterial.h"



// using VonMisesLinearHardening =
//     ASDPlasticMaterial <LinearIsotropic3D_EL,
//         VonMises_YF<
//             BackStressIV<TensorLinearHardeningFunction>,
//             VonMisesRadiusIV<ScalarLinearHardeningFunction>
//         >,
//         VonMises_PF<
//             BackStressIV<TensorLinearHardeningFunction>,
//             VonMisesRadiusIV<ScalarLinearHardeningFunction>
//         >,
//     ND_TAG_ASDPlasticMaterial
//     >;

// using proto_map_key = std::tuple <
//                       const char * /* elasticity*/,
//                       const char * /* yield surface*/,
//                       const char * /* plastic potential*/,
//                       const char * /* internal-variables-and-hardening-laws*/
//                        >;
// using proto_map_value = NDMaterial*;
// using proto_map_t = std::map<proto_map_key, proto_map_value>;
// inline proto_map_t make_prototypes() {
//     proto_map_t m;

//     {    // a proto of VM + VM
//         auto iproto = new VonMisesLinearHardening(0);
//         std::string elname = iproto->getELName();
//         cout << "elname = " << elname.c_str() << endl;
//         std::string yfname = iproto->getYFName();
//         cout << "yfname = " << yfname.c_str() << endl;
//         std::string pfname = iproto->getPFName();
//         cout << "pfname = " << pfname.c_str() << endl;
//         std::string ivname = iproto->getIVName();
//         cout << "ivname = " << ivname.c_str() << endl;
//         m[std::make_tuple(elname.c_str(),yfname.c_str(),pfname.c_str(),ivname.c_str())] = iproto;
//     }

//     return m;
// }


// template<>
//     "VonMises_YF",
//     "VonMises_PF",
//     "LinearIsotropic3D",
//     "BackStress(LinearHardeningForTensor):VonMisesRadius(LinearHardeningForScalar):"

