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

//Some basic trait classes to aid with implementation of CEP
#ifndef ASDPlasticityTraits_H
#define ASDPlasticityTraits_H

template< typename T >
struct supports_pre_integration_callback
{
    static const bool value = false;
};

// struct defines_pre_integration_callback
// {};

// template<>
// struct supports_pre_integration_callback<defines_pre_integration_callback>
// {
//     static const bool value = true;
// };


template< typename T >
struct supports_post_integration_callback
{
    static const bool value = false;
};

template< typename T >
struct supports_implicit_integration
{
    static const bool value = false;
};

template< typename VarType >
struct requires_hardening_saturation_limit_check
{
    static const bool value = false;
};


template< typename VarType >
struct yield_function_has_apex
{
    static const bool value = false;
};



#endif