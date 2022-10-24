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

#include "../EvolvingVariable.h"
#include <cmath>



template<class AlphaHardeningType, class KHardeningType>
class VonMises_PF : public PlasticFlowBase<VonMises_PF<AlphaHardeningType, KHardeningType>> // CRTP
{
public:

    typedef EvolvingVariable<VoigtVector, AlphaHardeningType> AlphaType;
    typedef EvolvingVariable<VoigtScalar, KHardeningType> KType;

    // PlasticFlowBase<VonMises_PF<HardeningType>>::PlasticFlowBase(), // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
    VonMises_PF( AlphaType &alpha_in, KType &k_in):
        PlasticFlowBase<VonMises_PF<AlphaHardeningType , KHardeningType >>::PlasticFlowBase(), // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                alpha_(alpha_in), k_(k_in) { }

    const VoigtVector& operator()(const VoigtVector &depsilon, const VoigtVector& sigma)
    {
        //Identical to derivative of VonMises_YF wrt sigma (a.k.a nij)
        const VoigtVector &alpha = alpha_.getVariableConstReference();

        result = sigma.deviator() - alpha;

        double den = sqrt(result.dot(result));
        if (abs(den) > 100*ASDPlasticMaterialGlobals::MACHINE_EPSILON)
            result = result / den;

        return result;
    }

private:
    AlphaType &alpha_;
    KType &k_;

    static VoigtVector result; //For returning VoigtVectors
};

template<class AlphaHardeningType, class KHardeningType>
VoigtVector VonMises_PF<AlphaHardeningType , KHardeningType >::result;


#endif