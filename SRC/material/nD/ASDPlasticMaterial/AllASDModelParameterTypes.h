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

#ifndef _AllASDModelParametersType
#define _AllASDModelParametersType
#include "ModelParameterType.h"

// ============================================================================
// Model Parameters associated with Elasticity
struct YoungsModulusName { static constexpr const char* name = "YoungsModulus";};
using YoungsModulus = ModelParameterType<double, YoungsModulusName>;
struct PoissonsRatioName { static constexpr const char* name = "PoissonsRatio";};
using PoissonsRatio = ModelParameterType<double, PoissonsRatioName>;

struct ReferenceYoungsModulusName { static constexpr const char* name = "ReferenceYoungsModulus";};  // Bulk modulus at reference pressure
using ReferenceYoungsModulus = ModelParameterType<double, YoungsModulusName>;
struct ReferencePressureName { static constexpr const char* name = "ReferencePressure";};       // The reference pressure
using ReferencePressure = ModelParameterType<double, YoungsModulusName>;
struct DuncanChang_MaxSigma3Name { static constexpr const char* name = "DuncanChang_MaxSigma3";};       // The reference pressure
using DuncanChang_MaxSigma3 = ModelParameterType<double, YoungsModulusName>;
struct DuncanChang_nName { static constexpr const char* name = "DuncanChang_n";};       // The reference pressure
using DuncanChang_n = ModelParameterType<double, YoungsModulusName>;


// ============================================================================
// Model Parameters associated with Yield Function
struct InternalFrictionAngleName { static constexpr const char* name = "InternalFrictionAngle";};
using InternalFrictionAngle = ModelParameterType<double, InternalFrictionAngleName>;

//RMC: Rounded-Mohr-Coulomb
struct RMC_m_Name { static constexpr const char* name = "RMC_m";};
using RMC_m = ModelParameterType<double, RMC_m_Name>;
struct RMC_qa_Name { static constexpr const char* name = "RMC_qa";};
using RMC_qa = ModelParameterType<double, RMC_qa_Name>;
struct RMC_pc_Name { static constexpr const char* name = "RMC_pc";};
using RMC_pc = ModelParameterType<double, RMC_pc_Name>;
struct RMC_e_Name { static constexpr const char* name = "RMC_e";};
using RMC_e = ModelParameterType<double, RMC_e_Name>;



// ============================================================================
// Model Parameters associated with Plastic flow

// ============================================================================
// Model Parameters associated with Evolving Variables
struct ScalarLinearHardeningParameter_Name { static constexpr const char* name = "ScalarLinearHardeningParameter";};
using ScalarLinearHardeningParameter = ModelParameterType<double, ScalarLinearHardeningParameter_Name>;
struct TensorLinearHardeningParameter_Name { static constexpr const char* name = "TensorLinearHardeningParameter";};
using TensorLinearHardeningParameter = ModelParameterType<double, TensorLinearHardeningParameter_Name>;


//For Armstrong-Frederick Hardening
struct AF_ha_Name { static constexpr const char* name = "AF_ha";};
using AF_ha = ModelParameterType<double, AF_ha_Name>;
struct AF_cr_Name { static constexpr const char* name = "AF_cr";};
using AF_cr = ModelParameterType<double, AF_cr_Name>;

//For Constant-Dilatancy
struct Dilatancy_Name { static constexpr const char* name = "Dilatancy";};
using Dilatancy = ModelParameterType<double, Dilatancy_Name>;

//For ScalarExponentialLinear hardening
struct ScalarExponentialLinear_Sigma0_Name { static constexpr const char* name = "ScalarExponentialLinear_Sigma0";};
using ScalarExponentialLinear_Sigma0 = ModelParameterType<double, ScalarExponentialLinear_Sigma0_Name>;
struct ScalarExponentialLinear_SigmaInf_Name { static constexpr const char* name = "ScalarExponentialLinear_SigmaInf";};
using ScalarExponentialLinear_SigmaInf = ModelParameterType<double, ScalarExponentialLinear_SigmaInf_Name>;
struct ScalarExponentialLinear_delta_Name { static constexpr const char* name = "ScalarExponentialLinear_delta";};
using ScalarExponentialLinear_delta = ModelParameterType<double, ScalarExponentialLinear_delta_Name>;




// ============================================================================
// Other Model Parameters 
struct MassDensity_Name { static constexpr const char* name = "MassDensity";};
using MassDensity = ModelParameterType<double, MassDensity_Name>;
struct InitialP0_Name { static constexpr const char* name = "InitialP0";};
using InitialP0 = ModelParameterType<double, InitialP0_Name>;


#endif //not defined _AllASDModelParametersType
