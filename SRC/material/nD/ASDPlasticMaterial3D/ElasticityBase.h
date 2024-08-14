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

#ifndef ElasticityBase_H
#define ElasticityBase_H

#include "EigenAPI.h"
#include <Channel.h>

// Helper template to check if a class has a parameters_t type alias
template <typename T, typename = void>
struct has_parameters_t : std::false_type {};

template <typename T>
struct has_parameters_t<T, typename std::enable_if<!std::is_same<typename T::parameters_t, void>::value>::type> : std::true_type {};

#define ELASTICITY_MATRIX template<class ParameterStorageType> \
    const VoigtMatrix& operator()(const VoigtVector& stress, \
        const ParameterStorageType& parameters_storage) const 

#define GET_PARAMETER_VALUE(type) parameters_storage.template get<type> ().value


template <class T>
class ElasticityBase
{
public:

    ElasticityBase() {
        static_assert(has_parameters_t<T>::value, "Derived class must have a 'parameters_t' type alias.");
    }
    
    ELASTICITY_MATRIX
    {
        return static_cast<T*>(this)->operator()(stress, parameters_storage);
    }

protected:

    static VoigtMatrix EE_MATRIX; 
};

template <class T>
VoigtMatrix ElasticityBase<T>::EE_MATRIX;


#endif
