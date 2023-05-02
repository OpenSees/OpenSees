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

#ifndef VonMises_PF_H
#define VonMises_PF_H

#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterialGlobals.h"

#include <cmath>
#include <typeinfo>


using namespace std;

template<class AlphaHardeningType, class KHardeningType>
class VonMises_PF : public PlasticFlowBase<VonMises_PF<AlphaHardeningType, KHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "VonMises_PF";

    VonMises_PF():
        PlasticFlowBase<VonMises_PF<AlphaHardeningType , KHardeningType >>::PlasticFlowBase()  // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                { }

    template <typename StorageType, typename ParameterStorageType>
    const VoigtVector& operator()(
    	const VoigtVector &depsilon, 
    	const VoigtVector& sigma, 
    	const StorageType& internal_variables_storage, 
    	const ParameterStorageType& parameters_storage)
    {
        //Identical to derivative of VonMises_YF wrt sigma (a.k.a nij)
        const AlphaHardeningType& AHT = 
        internal_variables_storage.template get<AlphaHardeningType> ();

        VoigtVector alpha = AHT.trial_value;

        result = sigma.deviator() - alpha;

        double den = sqrt(result.dot(result));
        if (abs(den) > 100*ASDPlasticMaterialGlobals::MACHINE_EPSILON)
            result = result / den;

        return result;
    }

    using internal_variables_t = std::tuple<AlphaHardeningType, KHardeningType>;
    using parameters_t = std::tuple<>;

private:

    static VoigtVector result; //For returning VoigtVectors
};

template<class AlphaHardeningType, class KHardeningType>
VoigtVector VonMises_PF<AlphaHardeningType , KHardeningType >::result;


#endif