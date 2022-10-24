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

#ifndef ElasticityBase_H
#define ElasticityBase_H

#include "EigenAPI.h"
#include "EvolvingVariable.h"
#include <Channel.h>


template <class T>
class ElasticityBase
{
public:

    ElasticityBase()
    {
        // Derived classes will have custom constructors.
    }

    // Note the use of the Curiously-recurring-template-pattern
    // Operator () retunrs a const-reference. Thereforce, implementation must provide
    // a reference to a persistent data member. Usually this one is declared static
    // so that storage and be reused across instances, and we avoid calls to malloc.
    const VoigtMatrix& operator()(const VoigtVector& stress) 
    {
        return static_cast<T*>(this)->operator()(stress);
    }

    int sendSelf(int commitTag, Channel &theChannel)
    {
        return static_cast<T*>(this)->sendSelf(commitTag, theChannel);
    }

    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
    {
        return static_cast<T*>(this)->recvSelf( commitTag, theChannel, theBroker);
    }

private:
    // Derived classes will store their parameters and such.
};




#endif