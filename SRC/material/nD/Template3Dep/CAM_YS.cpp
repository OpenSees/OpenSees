
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Cam Clay yield criterion                          	 #
//#                      (Ref. Wood p113)                    			 #
//# CLASS:             CAMYieldSurface(for monotonic loading)                    #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, ZHaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              Mar. 28, 2001                                             #
//# UPDATE HISTORY:                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//================================================================================


#ifndef CAM_YS_CPP
#define CAM_YS_CPP

#include "CAM_YS.h"
#include <basics.h>


//================================================================================
// Normal constructor
//================================================================================

CAMYieldSurface::CAMYieldSurface(double Mp )  
{
   M = Mp;
} 


//================================================================================
//create a colne of itself
//================================================================================

YieldSurface * CAMYieldSurface::newObj() {  

     YieldSurface  *new_YS = new CAMYieldSurface(M);
     return new_YS;

}

//================================================================================
// Copy constrstructor
//================================================================================
//
//CAMYieldSurface::CAMYieldSurface(CAMYieldSurface &CAMYS ) { }
//


//================================================================================
//  Yield criterion evaluation function F(EPState)
//================================================================================

double CAMYieldSurface::f(const EPState *EPS) const
{
  double p =EPS->getStress().p_hydrostatic();
  double q =EPS->getStress().q_deviatoric();
  //cout << "p " << p;
	   
  double po = EPS->getScalarVar(1);
  double temp3 = q*q - M*M*p*(po - p);

  //printf("\n========Inside f = %.4f \n ", temp3);

  return temp3;
}


//================================================================================
// tensor dF/dsigma_ij
//================================================================================

tensor CAMYieldSurface::dFods(const EPState *EPS) const {
  
  tensor dFoverds( 2, def_dim_2, 0.0);
  tensor I2("I", 2, def_dim_2);

  double p = EPS->getStress().p_hydrostatic();
  double q = EPS->getStress().q_deviatoric();
  double po = EPS->getScalarVar( 1 );
	 tensor DpoDs = EPS->getStress().dpoverds();
	 tensor DqoDs = EPS->getStress().dqoverds();

  double dFoverdp = -1.0*M*M*( po - 2.0*p );
  double dFoverdq = 2.0*q;

  dFoverds = DpoDs  * dFoverdp +
             DqoDs  * dFoverdq;

  return dFoverds;

}


//================================================================================
// double xi_s1 = dF/dpo = -M*M*p  Derivative in terms of first scalar var 
//================================================================================

double CAMYieldSurface::xi_s1( const EPState *EPS ) const {

    //double p = EPS->getStress().Iinvariant1()/3.0;
    double p = EPS->getStress().p_hydrostatic();
    return -1.0 * M * M * p;

}


// moved to stresstensor Boris Jeremic@ucdavis.edu 21Aug2001
// // I took these functions from mmodel.cpp, programmed by Dr. Jeremic
// //#############################################################################
// tensor CAMYieldSurface::dpoverds( ) const
//   {
//     tensor ret(2, def_dim_2, 0.0);
//     tensor I2("I", 2, def_dim_2);
//     ret = I2*(-1.0/3.0);
//     ret.null_indices();
//     return ret;
//   }
// 
// //#############################################################################
// tensor CAMYieldSurface::dqoverds(const EPState *EPS) const
//   {
//     
//     stresstensor stress = EPS->getStress();
//     
//     tensor ret(2, def_dim_2, 0.0);
//     double q = stress.q_deviatoric();
//     stresstensor s( 0.0);
// //...    double J2D = stress.Jinvariant2();
// //...    double temp1 = sqrt(J2D);
// //...    double temp2 = sqrt(3.0)/(2.0*temp1);
//     double temp2 = (3.0/2.0)*(1/q);
//     s = stress.deviator();
//     ret = s*temp2;
//     ret.null_indices();
//     return ret;
//   }
// 
// //#############################################################################
// tensor CAMYieldSurface::dthetaoverds(const EPState *EPS) const
//   {
//     stresstensor stress = EPS->getStress();
// 
//     tensor ret(2, def_dim_2, 0.0);
//     stresstensor s( 0.0);
//     stresstensor t( 0.0);
//     tensor I2("I", 2, def_dim_2);
// 
// //    double EPS = pow(d_macheps(),(1./3.));
//     double J2D = stress.Jinvariant2();
//     double q     = stress.q_deviatoric();
//     double theta = stress.theta();
// 
// //out    while ( theta >= 2.0*PI )
// //out      theta = theta - 2.0*PI; // if bigger than full cycle
// //out    while ( theta >= 4.0*PI/3.0 )
// //out      theta = theta - 4.0*PI/3.0; // if bigger than four thirds of half cycle
// //out    while ( theta >= 2.0*PI/3.0 )
// //out      theta = theta - 2.0*PI/3.0; // if bigger than two third of half cycle
// //out    while ( theta >= PI/3.0 )
// //out      theta = 2.0*PI/3.0 - theta; // if bigger than one third of half cycle
// //out
// //out    if ( theta < 0.0001 )
// //out      {
// //out        ::printf("theta = %.10e in CAMYieldSurface::dthetaoverds(stress)\n",
// //out                           theta);
// //out//>><<>><<>><<        return ret;
// //out      }
// 
//     double c3t = cos(3.0*theta);
//     double s3t = sin(3.0*theta);
// 
//     double tempS = (3.0/2.0)*(c3t/(q*q*s3t));
//     double tempT = (9.0/2.0)*(1.0/(q*q*q*s3t));
// 
//     s = stress.deviator();
//     t = s("qk")*s("kp") - I2*(J2D*(2.0/3.0));
// 
//     s.null_indices();
//     t.null_indices();
// 
//     ret = s*tempS - t*tempT;
//     ret.null_indices();
//     return ret;
// }

//================================================================================
double CAMYieldSurface::getM() const
{       
    return M;
}


OPS_Stream& operator<< (OPS_Stream& os, const CAMYieldSurface & YS)
{
   os << "Cam Clay Yield Surface Parameters: " << endln;
   os << "M = " << YS.M << endln;
   return os;
}

#endif

