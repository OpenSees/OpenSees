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

#ifndef DruckerPrager_PF_H
#define DruckerPrager_PF_H

#include "../PlasticFlowBase.h"
#include "../ASDPlasticMaterial3DGlobals.h"
using namespace ASDPlasticMaterial3DGlobals;

#include <cmath>
#include <typeinfo>


template<class AlphaHardeningType, class EtaHardeningType>
class DruckerPrager_PF : public PlasticFlowBase<DruckerPrager_PF<AlphaHardeningType, EtaHardeningType>> // CRTP
{
public:

    static constexpr const char* NAME = "DruckerPrager_PF";


    DruckerPrager_PF( ):
        PlasticFlowBase<DruckerPrager_PF<AlphaHardeningType, EtaHardeningType >>::PlasticFlowBase()  // Note here that we need to fully-qualify the type of YieldFunctionBase, e.g. use scope resolution :: to tell compiler which instance of YieldFunctionBase will be used :/
                { }

    PLASTIC_FLOW_DIRECTION
    {
        auto s = sigma.deviator();
        auto alpha = GET_TRIAL_INTERNAL_VARIABLE(AlphaHardeningType);
        
        // Get dilation parameter (etabar) from parameters
        double etabar = GET_PARAMETER_VALUE(DP_etabar);  // dilation parameter (controls plastic volume change)
        
        // Compute deviatoric part
        VoigtVector dev_part = s - alpha;
        double den = sqrt(0.5*tensor_dot_stress_like(dev_part, dev_part));
        
        if (abs(den) > sqrt(0.5*tensor_dot_stress_like(s, s))*ASDPlasticMaterial3DGlobals::MACHINE_EPSILON)
            dev_part = dev_part / den;
        else
            dev_part *= 0.0;  // Zero out if denominator too small
            
        // Add pressure-dependent part: etabar * dp/dsigma = etabar/3 * I
        VoigtVector pressure_part;
        pressure_part *= 0.0;  // Initialize to zero
        pressure_part(0) = etabar / 3.0;  // sigma_xx component
        pressure_part(1) = etabar / 3.0;  // sigma_yy component  
        pressure_part(2) = etabar / 3.0;  // sigma_zz component
        // shear components remain zero
        
        vv_out = dev_part + pressure_part;
        
        return vv_out;
    }



    using internal_variables_t = std::tuple<AlphaHardeningType, EtaHardeningType>;
    using parameters_t = std::tuple<DP_etabar>;

private:

    static VoigtVector vv_out; 

};


template<class AlphaHardeningType, class EtaHardeningType>
VoigtVector DruckerPrager_PF<AlphaHardeningType, EtaHardeningType  >::vv_out;

#endif
