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

#include "ConstitutiveModels/VonMisesLinearHardening.h"



// ASDPlasticMaterial * get_instace_of_ASDPlasticMaterial(const char * name)
// {
	
// }


    // using proto_map_key = std::tuple<
    //     const char* /* yield surface*/, 
    //     const char* /* plastic potential*/>;
    // using proto_map_value NDMaterial*;
    // using proto_map_t = std::map<proto_map_key, proto_map_value>;
    // inline make_prototypes() {
    //     proto_map_t m;
        
    //     { // a proto of VM + VM
    //         auto iproto = new ASDPlastic<VMYS, VMPP>();
    //         m[std::make_tuple(VMYS::NAME, VMPP::NAME)] = iprotp;
    //     }
        
    //     return m;
    // }
    // static proto_map_t TheAvailableProtos = make_prototypes();
    



// #include "ConstitutiveModels/DruckerPragerLinearHardening.h"
// #include "ConstitutiveModels/DruckerPragerVonMisesLinearHardening.h"
// #include "ConstitutiveModels/DruckerPragerArmstrongFrederickLE.h"
// #include "ConstitutiveModels/DruckerPragerArmstrongFrederickNE.h"
// #include "ConstitutiveModels/DruckerPragerNonAssociateLinearHardening.h"
// #include "ConstitutiveModels/DruckerPragerNonAssociateArmstrongFrederick.h"
// #include "ConstitutiveModels/VonMisesArmstrongFrederick.h"
// #include "ConstitutiveModels/CamClayLT.h"
// #include "ConstitutiveModels/RoundedMohrCoulomb.h"
// #include "ConstitutiveModels/sanisand2004.h"