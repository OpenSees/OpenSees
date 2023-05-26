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

#ifndef EvolvingVariable_H
#define EvolvingVariable_H

#include <EigenAPI.h>
#include <Channel.h>
#include <iostream>
#include "ASDPlasticMaterialTraits.h"
using namespace std;


//Forward declaration needed for declaring the operator overload << (friend)
template<class VarType, class T>
class EvolvingVariable;

//Forward declaration
template<class VarType, class T>
std::ostream& operator<<(std::ostream& os, const EvolvingVariable<VarType, T>& obj);



template<class VarType, class T>
class EvolvingVariable
{
public:

    EvolvingVariable(VarType a_): a(a_), a_committed(a_), a_tmp(a_) { }

    const VarType &getDerivative(const VoigtVector &depsilon,
                                 const VoigtVector &m,
                                 const VoigtVector& sigma) const
    {
        return static_cast<const T*>(this)->getDerivative(depsilon,  m,  sigma);
    };

    EvolvingVariable<VarType, T> & operator= ( const EvolvingVariable<VarType, T> & other)
    {
        //Check self-assignment
        if (&other == this)
        {
            return *this;
        }

        a = other.a;
        a_tmp = other.a_tmp;
        a_committed = other.a_committed;

        return *this;
    }

    template <typename U = T>
    typename std::enable_if < !evolving_variable_implements_custom_evolve_function<U>::value, void >::type
    evolve(double dlambda,
           const VoigtVector& depsilon,
           const VoigtVector& m,
           const VoigtVector& sigma)
    {
        const VarType& h = getDerivative(depsilon, m, sigma);
        static VarType aux(a);
        aux = h;
        aux *= dlambda;
        a +=  aux;
    }

    template <typename U = T>
    typename std::enable_if < evolving_variable_implements_custom_evolve_function<U>::value, void >::type
    evolve(double dlambda,
           const VoigtVector& depsilon,
           const VoigtVector& m,
           const VoigtVector& sigma)
    {
        static_cast<U*>(this)->evolve(dlambda, depsilon, m, sigma);
    }


    const VarType& getVariableConstReference() const
    {
        return a;
    }

    VarType getVariable() const
    {
        return a;
    }


    const VarType& getCommittedVariableConstReference() const
    {
        return a_committed;
    }

    VarType getCommittedVariable() const
    {
        return a_committed;
    }

    void setVar(const VarType& v)
    {
        a = v;
    }

    void setCommittedVar(const VarType& v)
    {
        a_committed = v;
    }

    void commit()
    {
        // cout << "Final commit state from " << a_committed << " to " << a << endl;
        a_committed = a;
    }

    void revert()
    {
        a = a_committed;
    }

    void commit_tmp()
    {
        // cout << "Commiting from " << a_tmp << " to " << a << endl;
        a_tmp = a;
    }

    void revert_tmp()
    {
        // cout << "Reverting from " << a << " to " << a_tmp << endl;
        a = a_tmp;
    }

    void revert_to_start()
    {
        // cout << "Reverting from " << a << " to " << a_tmp << endl;
        a = 0 * a_tmp;
    }

    void print(std::ostream &s)
    {
        s << "     > a           = " << a << endl;
        s << "     > a_tmp       = " << a << endl;
        s << "     > a_committed = " << a << endl;
    }

    //Overloaded operators.
    operator VarType ()
    {
        return a;    // Convert into variable type
    }
    friend std::ostream& operator<< <>(std::ostream& os, const EvolvingVariable<VarType, T>& obj);

    int sendSelf(int commitTag, Channel &theChannel)
    {
        return static_cast<T*>(this)->sendSelf(commitTag, theChannel);
    }

    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
    {
        return static_cast<T*>(this)->recvSelf(commitTag, theChannel, theBroker);
    }

    // ///////////////////////////////////////////////////////////////////////////////////
    // Set the limit for the hardening saturation.
    // With small strain increment, the hardening saturation is implemented smoothly.
    // With great strain increment, the overshooting creates problems. So set limitation.
    // //////////////////////////////////////////////////////////////////////////////////
    template <typename U = T>
    typename std::enable_if < !requires_hardening_saturation_limit_check<U>::requires >::type
    check_hardening_saturation_limit_(VarType& a, VoigtVector const& plasticFlow_m) {}

    template <typename U = T>
    typename std::enable_if<requires_hardening_saturation_limit_check<U>::requires>::type
    check_hardening_saturation_limit_(VarType& a, VoigtVector const& plasticFlow_m)
    {
        static_cast<U*>(this)->check_hardening_saturation_limit(a, plasticFlow_m);
    }

    // ///////////////////////////////////////////////////////////////////////////////////
private:
    VarType a;
    VarType a_committed;
    VarType a_tmp;
    static VarType da_1;
    static VarType da_2;
};





// Forward stream operators to underlying class
template<class VarType, class T>
std::ostream& operator<<(std::ostream& os, const EvolvingVariable<VarType, T>& obj)
{
    os << obj.a;
    return os;
}

// template<class VarType, class T>
// VarType EvolvingVariable<VarType, T>::a_tmp;

#endif //EvolvingVariable_H