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
// ASDPlasticMaterial3D
//
// Fully general templated material class for plasticity modeling

#ifndef ConstantDilatancy_PF_H
#define ConstantDilatancy_PF_H

#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"
#include "../AllASDModelParameterTypes.h"

#include <cmath>
#include <typeinfo>


template<class AlphaHardeningType>
class ConstantDilatancy_PF : public PlasticFlowBase<ConstantDilatancy_PF<AlphaHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "ConstantDilatancy_PF";


    ConstantDilatancy_PF( ):
        PlasticFlowBase<ConstantDilatancy_PF<AlphaHardeningType >>::PlasticFlowBase()  // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                { }

    PLASTIC_FLOW_DIRECTION
    {
        double p = -sigma.meanStress();
        auto s = sigma.deviator();
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        double D = GET_PARAMETER_VALUE(Dilatancy);

        vv_out = s - alpha;
        double den = sqrt(vv_out.dot(vv_out));
        if (abs(den) > sqrt(s.dot(s))*ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
            vv_out = vv_out / den;

        vv_out -= D * kronecker_delta() / 3;


        return vv_out;
    }

    using internal_variables_t = std::tuple<AlphaHardeningType>;
    using parameters_t = std::tuple<Dilatancy>;

private:

    static VoigtVector vv_out; 
};


template<class AlphaHardeningType>
VoigtVector ConstantDilatancy_PF<AlphaHardeningType  >::vv_out;

#endif
