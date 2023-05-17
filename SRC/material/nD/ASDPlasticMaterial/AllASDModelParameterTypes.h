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

// ============================================================================
// Internal variables base template struct
template<typename T, class NAMER>
struct ModelParameterType {

    static constexpr const char* NAME = NAMER::name;
    T value;
    ModelParameterType() = default;
    ModelParameterType(T x) : value(x) {}
    inline const char* getName() const { return NAME; }
};

//Templated stream insertion
template<typename T, class NAMER>
std::ostream& operator<<(std::ostream& os, const ModelParameterType<T, NAMER>& param) {
    os << param.NAME << " = " << param.value;
    return os;
}

// ============================================================================
// Internal variables associated with Elasticity
struct YoungsModulusName { static constexpr const char* name = "YoungsModulus";};
using YoungsModulus = ModelParameterType<double, YoungsModulusName>;
struct PoissonsRatioName { static constexpr const char* name = "PoissonsRatio";};
using PoissonsRatio = ModelParameterType<double, PoissonsRatioName>;

// ============================================================================
// Internal variables associated with Yield Function
struct InternalFrictionAngleName { static constexpr const char* name = "InternalFrictionAngle";};
using InternalFrictionAngle = ModelParameterType<double, InternalFrictionAngleName>;

// ============================================================================
// Internal variables associated with Plastic flow

// ============================================================================
// Internal variables associated with Evolving Variables
struct ScalarLinearHardeningParameter_Name { static constexpr const char* name = "ScalarLinearHardeningParameter";};
using ScalarLinearHardeningParameter = ModelParameterType<double, ScalarLinearHardeningParameter_Name>;
struct TensorLinearHardeningParameter_Name { static constexpr const char* name = "TensorLinearHardeningParameter";};
using TensorLinearHardeningParameter = ModelParameterType<double, TensorLinearHardeningParameter_Name>;

#endif //not defined _AllASDModelParametersType