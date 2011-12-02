/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           General platform for elaso-plastic constitutive model       #
#                    implementation                                              #
# CLASS:             ConstitutiveDriver (the class containing all constitutive   #
#                                        drivers--integration algorithm)         #
#                                                                                #
# VERSION:                                                                       #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )    #
# TARGET OS:         DOS || UNIX || . . .                                        #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                 #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                                 #
#                                                                                #
#                                                                                #
# DATE:              08-03-2000                                                  #
# UPDATE HISTORY:                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
# SHORT EXPLANATION: Currently two explict algorithms, i.e. ForwardEulerEPState  #
#                    and FESubIncrementation, and two implicit algorithms, i.e.  #                                                                              #
#                    BackwardEulerEPState and BESubIncrementation are 		 #
#                    implemented.                                                 #
//================================================================================
*/

#ifndef ConDriver_H
#define ConDriver_H

#include <stresst.h>
#include <straint.h>
#include <BJtensor.h>
#include <Template3Dep.h>

class ConstitutiveDriver
{
  public:

     ConstitutiveDriver() {};

     EPState ForwardEulerEPState(straintensor &strain_increment, 
                                 Template3Dep &material_point);
     
     EPState SemiBackwardEulerEPState( straintensor &strain_increment, 
                                       Template3Dep &material_point);

     EPState FESubIncrementation( straintensor &strain_increment,
                                  Template3Dep &material_point,
                                  int number_of_subincrements);

     EPState BackwardEulerEPState( straintensor &strain_increment, 
                                   Template3Dep &material_point);

     EPState BESubIncrementation( straintensor & strain_increment,
                                  Template3Dep & material_point,
                                  int number_of_subincrements);                                                 
						                                                  
    //================================================================================
    // this one is intended to shell the previous three and to decide 
    // ( according to the data stored in Material_Model object ) 
    // which constitutive tensor to return ( forward ( non-constistent
    // or backward ( consistent ) or . . . 
    
    //virtual tensor ConstitutiveTensor(stresstensor   & final_stress, 
    //                                 stresstensor   & start_stress,
    //                                 straintensor   & strain_increment,
    //                                 Material_Model & Criterion,
    //                                 double           just_this_PP );
    

    //================================================================================
    // trying to find intersection point
    // according to M. Crisfield's book
    // "Non-linear Finite Element Analysis of Solids and Structures "
    // chapter 6.6.1 page 168.
    //================================================================================
  private:
    EPState PredictorEPState(straintensor & strain_increment,
                              Template3Dep & material_point);

    stresstensor yield_surface_cross(const stresstensor & start_stress,
                                     const stresstensor & end_stress,
                                     const Template3Dep & material_point);

    double zbrentstress(const stresstensor & start_stress,
                        const stresstensor & end_stress,
                        const Template3Dep & material_point,
                        double x1, double x2, double tol);

    double func( const stresstensor &start_stress,
                 const stresstensor &end_stress,
                 const Template3Dep &material_point,
                 double alfa );

    //================================================================================
    // Routines to generate elastic constitutive tensors
    //================================================================================
    tensor ElasticStiffnessTensor(double Ey, double nu) const;
    //tensor ElasticStiffnessTensor(const EPState &EPS) const;

    tensor ElasticComplianceTensor(const double Ey, double nu) const;
    //tensor ElasticComplianceTensor(const EPState &EPS) const;

};


#endif

