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
    int size() const {return trial_value.size();}

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
