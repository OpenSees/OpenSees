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

#include "ASDPlasticMaterialGlobals.h"

//Base template struct for all internal variables
template <class EvolvingVariableType, class HardeningType, class NAMER>
struct InternalVariableType {
    static constexpr const char* NAME = NAMER::name;
    EvolvingVariableType trial_value;
    EvolvingVariableType committed_value;
    InternalVariableType() = default;
    InternalVariableType(EvolvingVariableType x) : trial_value(x), committed_value(x) {}
    void commit() {committed_value = trial_value;};
    void revert() {trial_value = committed_value;};
    inline const char* getName() const { return NAME; }

    template <class ParameterStorageType>
    EvolvingVariableType hardening_function(
                const VoigtVector &depsilon,
                const VoigtVector &m,
                const VoigtVector& sigma,
                const ParameterStorageType& parameters) const
    {
        return HardeningType::f(depsilon, m, sigma, parameters);
    }

    using parameters_t = typename HardeningType::parameters_t;

    // template <class ParameterStorageType>
    // using parameters_t = typename HardeningType::template parameters_t<ParameterStorageType>;

};

//Stream operator for internal variables
template <class EvolvingVariableType, class HardeningType, class NAMER>
std::ostream& operator<<(std::ostream& os, const InternalVariableType<EvolvingVariableType, HardeningType, NAMER>& param) {
    os << param.getName() << ": [Trial: " << param.trial_value.transpose() << ", Committed: " << param.committed_value.transpose() << "]";
    os << endl;
    os << "   HardeningType::parameters --> " << typeid(typename InternalVariableType<EvolvingVariableType, HardeningType, NAMER>::parameters_t).name() << endl;
    return os;
}

//Definitions of possible internal variables
struct BackStressName { static constexpr const char* name = "BackStress";};
template <class HardeningType>
using BackStressIV = InternalVariableType<VoigtVector, HardeningType, BackStressName>;

struct VonMisesRadiusName { static constexpr const char* name = "VonMisesRadius";};
template <class HardeningType>
using VonMisesRadiusIV = InternalVariableType<VoigtScalar, HardeningType, VonMisesRadiusName>;

// template<class T, class EVT> 
// T& operator << (T& stream, const BackStress<EVT>& x) {
//     return stream << x.value;
// }

//VonMises Cylinder radius type
// template <class EvolvingVariableType>
// struct VonMises_m {
//     static constexpr const char* NAME = "VonMises_m";
//     EvolvingVariableType value;
//     VonMises_m() = default;
//     VonMises_m(EvolvingVariableType x) : value(x) {}
//     inline const char* getName() const { return NAME; }
// };
// template<class T, class EVT> 
// T& operator << (T& stream, const VonMises_m<EVT>& x) {
//     return stream << x.value;
// }




// Base template struct
// template<typename Derived, typename EvolvingVariableType>
// struct EvolvingVariable {
//     static constexpr const char* NAME = Derived::NAME;
//     EvolvingVariableType value;
//     EvolvingVariable() = default;
//     EvolvingVariable(EvolvingVariableType x) : value(x) {}
//     inline const char* getName() const { return NAME; }
// };
// template<typename Derived, typename EvolvingVariableType>
// std::ostream& operator<<(std::ostream& os, const EvolvingVariable<Derived, EvolvingVariableType>& param) {
//     // os << param.getName() << ": " << param.value;
//     os << param.value;
//     return os;
// }

// // Backstress variable type
// struct BackStressTag { static constexpr const char* NAME = "BackStress"; };
// template<typename EvolvingVariableType>
// struct BackStress : public EvolvingVariable<BackStress<EvolvingVariableType>, EvolvingVariableType> {};

// // VonMises Cylinder radius type
// struct VonMises_mTag { static constexpr const char* NAME = "VonMises_m"; };
// template<typename EvolvingVariableType>
// struct VonMises_m : public EvolvingVariable<VonMises_m<EvolvingVariableType>, EvolvingVariableType> {};


// //Backstress variable type
// template <class EvolvingVariableType>
// struct BackStress {
//     static constexpr const char* NAME = "BackStress";
//     EvolvingVariableType value;
//     BackStress() = default;
//     BackStress(EvolvingVariableType x) : value(x) {}
//     inline const char* getName() const { return NAME; }
// };
// template<class T, class EVT> 
// T& operator << (T& stream, const BackStress<EVT>& x) {
//     return stream << x.value;
// }

// //VonMises Cylinder radius type
// template <class EvolvingVariableType>
// struct VonMises_m {
//     static constexpr const char* NAME = "VonMises_m";
//     EvolvingVariableType value;
//     VonMises_m() = default;
//     VonMises_m(EvolvingVariableType x) : value(x) {}
//     inline const char* getName() const { return NAME; }
// };
// template<class T, class EVT> 
// T& operator << (T& stream, const VonMises_m<EVT>& x) {
//     return stream << x.value;
// }

