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

#ifndef MaterialInternalVariables_H
#define MaterialInternalVariables_H

#include <iostream>
#include "../../ltensor/LTensor.h"
#include <Channel.h>


template <class... Qs> struct MaterialInternalVariables
{
    void evolve(double dlambda,
                const DTensor2& depsilon,
                const DTensor2& m,
                const DTensor2& sigma)
    {
        // Last guy does nothing, as it holds nothing.
        // Probably gets optimized out. (Empty struct optimization)
        // Needed to stop template recursion.
    }

    void setVars(const MaterialInternalVariables<>& vars) {}

    int sendSelf(int commitTag, Channel &theChannel)
    {
        return 0;
    }

    int receiveSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
    {
        return 0;
    }

    void commit() { }
    void revert() { }
    void commit_tmp() { }
    void revert_tmp() { }
    void print(std::ostream &s) { }

};

template <class Q, class... Qs>
struct MaterialInternalVariables<Q, Qs...> : MaterialInternalVariables<Qs...>
{
    MaterialInternalVariables(Q &q, Qs&... qs) : MaterialInternalVariables<Qs...>(qs...), q_i(q) {}

    void setVars(const MaterialInternalVariables<Q, Qs...>& vars)
    {
        q_i = vars.q_i;
        const MaterialInternalVariables<Qs...>* morevars = static_cast<const MaterialInternalVariables<Qs...>*>(&vars);

        static_cast<MaterialInternalVariables<Qs...>*>(this)->setVars(*morevars);
    }

    //Recurse, calling evolve on each internal variable.
    void evolve(double dlambda,
                const DTensor2& depsilon,
                const DTensor2& m,
                const DTensor2& sigma)
    {
        q_i.evolve( dlambda,  depsilon,  m,  sigma);
        static_cast<MaterialInternalVariables<Qs...>*>(this)->evolve(dlambda,  depsilon,  m,  sigma);
    }

    int sendSelf(int commitTag, Channel &theChannel)
    {
        q_i.sendSelf(commitTag, theChannel);
        static_cast<MaterialInternalVariables<Qs...>*>(this)->sendSelf(commitTag, theChannel);
        return 0;
    }

    int receiveSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
    {
        q_i.receiveSelf(commitTag, theChannel, theBroker);
        static_cast<MaterialInternalVariables<Qs...>*>(this)->receiveSelf(commitTag, theChannel, theBroker);
        return 0;
    }

    void commit()
    {
        q_i.commit();
        static_cast<MaterialInternalVariables<Qs...>*>(this)->commit();
    }

    void revert()
    {
        q_i.revert();
        static_cast<MaterialInternalVariables<Qs...>*>(this)->revert();
    }

    void revert_to_start()
    {
        q_i.revert_to_start();
        static_cast<MaterialInternalVariables<Qs...>*>(this)->revert_to_start();
    }

    void commit_tmp()
    {
        q_i.commit_tmp();
        static_cast<MaterialInternalVariables<Qs...>*>(this)->commit_tmp();
    }

    void revert_tmp()
    {
        q_i.revert_tmp();
        static_cast<MaterialInternalVariables<Qs...>*>(this)->revert_tmp();
    }

    void print(std::ostream &s)
    {
        q_i.print( s);
        static_cast<MaterialInternalVariables<Qs...>*>(this)->print(s);
    }


    Q& q_i;
};


#endif
