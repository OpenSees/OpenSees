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


#ifndef LinearHardeningTensor_EV_H
#define LinearHardeningTensor_EV_H

#include "../EvolvingVariable.h"
#include "../ASDPlasticMaterialGlobals.h" 
#include "../AllASDModelParameterTypes.h"



class LinearHardeningTensor_EV : public EvolvingVariable<VoigtVector, LinearHardeningTensor_EV> //CRTP on LinearHardeningTensor_EV
{
public:

    LinearHardeningTensor_EV() : EvolvingVariable(VoigtVector(0,0,0,0,0,0)) {};
    LinearHardeningTensor_EV(VoigtVector v0) : EvolvingVariable(v0) {};

    template<class ParameterStorageType>
    const VoigtVector& getDerivative(const VoigtVector &depsilon,
                                     const VoigtVector &m,
                                     const VoigtVector& stress,
        const ParameterStorageType& parameters_storage) const
    {
        double H = parameters_storage.template get<LinearHardeningForTensor> ().value;
        derivative = H * m.deviator();
        return derivative;
    }

    // int sendSelf(int commitTag, Channel &theChannel) {return 0;}
    // int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {return 0;}

    using parameters_t = std::tuple<LinearHardeningForTensor>;

private:
    static VoigtVector derivative;     // Needs to be static so multiple instances only do one malloc call and we can return a const-reference
};


VoigtVector LinearHardeningTensor_EV::derivative;


#endif