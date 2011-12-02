//================================================================================
//# COPY LEFT and RIGHT:                                                         #
//# Commercial    use    of    this  program without express permission of the   #
//# University  of  California, is strictly encouraged. Copyright and Copyleft   #
//# are covered by the following clause:                                         #
//#                                                                              #
//# Woody's license:                                                             #
//# ``This    source    code is Copyrighted in U.S., by the The Regents of the   #
//# University  of  California,  for  an indefinite period, and anybody caught   #
//# using  it  without  our  permission,  will be mighty good friends of ourn,   #
//# cause  we  don't give a darn. Hack it. Compile it. Debug it. Run it. Yodel   #
//# it. Enjoy it. We wrote it, that's all we wanted to do.'' bj                  #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Rounded Mohr Coulomb Potential Surface                    #
//# CLASS:                                                                       #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++                                                       #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic jeremic@ucdavis.edu                         #
//#                    Zhao Cheng,                                               #
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic                                 #
//#                                                                              #
//#                                                                              #
//# DATE:              12 Feb. 2003                                              #
//# UPDATE HISTORY:    Feb 25 2003                                               #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/

#ifndef RMC01_PS_CPP
#define RMC01_PS_CPP

#include "RMC01_PS.h"

//================================================================================
//create a clone of itself
//================================================================================

PotentialSurface * RMC01PotentialSurface::newObj()
  {
    PotentialSurface  *new_YS = new RMC01PotentialSurface();
    return new_YS;
  }



//================================================================================
// tensor dQ/dsigma_ij
//================================================================================

tensor RMC01PotentialSurface::dQods(const EPState *EPS) const
{
  tensor dQoverds( 2, def_dim_2, 0.0);
//  double p = EPS->getStress().p_hydrostatic(); // p
  double q = EPS->getStress().q_deviatoric(); // q
  double theta = EPS->getStress().theta(); // theta
//  double temp_phi = EPS->getScalarVar(1)*3.14159265358979/180.0; // frictional angle
//  double temp_cohesive = EPS->getScalarVar(2); // cohesion
  tensor DpoDs = EPS->getStress().dpoverds(); // dp/ds
  tensor DqoDs = EPS->getStress().dqoverds(); // dq/ds
  tensor DthetaoDs = EPS->getStress().dthetaoverds(); // d(theta)/ds
  double alfa = EPS->getScalarVar(1); // Take alfa & k as internal variables
//  double k = EPS->getScalarVar(2);    // instead of phi & conhesive
  double a1 = (3.0*1.7320508076*alfa) / (2.0+1.7320508076*alfa);  
//  double a1 = -6*sin(temp_phi)/(3.0-sin(temp_phi));
//  double a2 = -6*temp_cohesive*cos(temp_phi)/(3.0-sin(temp_phi));
  double e = (3.0-a1)/(3.0+a1);
  double Frou = g_0(theta, e); // r(theta)
  double Frou_prime = g_prime(theta, e); // r'(theta)
//  double dQoverdp = a1; // dQ/dp
//  double dQoverdq = Frou; // dQ/dq
//  double dQoverdtheta = q*Frou_prime; // dQ/d(theta)
  double dQoverdp = alfa*(-3.0); // dQ/dp
  double dQoverdq = Frou/1.7320508076; // dQ/dq
  double dQoverdtheta = q*Frou_prime/1.7320508076; // dQ/d(theta)
  
  dQoverds = DpoDs  * dQoverdp +
             DqoDs  * dQoverdq +
             DthetaoDs * dQoverdtheta;

  return dQoverds;

}


//================================================================================
// tensor d2Q/dsigma_ij_2
//================================================================================

tensor RMC01PotentialSurface::d2Qods2( const EPState *EPS ) const
{
  tensor d2Qoverds2( 4, def_dim_4, 0.0); // d2Q/ds2(pqmn)

//  double p = EPS->getStress().p_hydrostatic(); // p
  double q = EPS->getStress().q_deviatoric(); // q
  double theta = EPS->getStress().theta(); // theta
//  double temp_phi = EPS->getScalarVar(1)*3.14159265358979/180.0; // frictional angle
//  double temp_cohesive = EPS->getScalarVar(2); // conhesion
  double alfa = EPS->getScalarVar(1);
//  double k = EPS->getScalarVar(2);
  double a1 = (3.0*1.7320508076*alfa) / (2.0+1.7320508076*alfa);
  tensor DpoDs = EPS->getStress().dpoverds(); // dp/ds
  tensor DqoDs = EPS->getStress().dqoverds(); // dq/ds
  tensor DthetaoDs = EPS->getStress().dthetaoverds(); // d(theta)/ds
  tensor D2poDs2 = EPS->getStress().d2poverds2(); // d2p/ds2
  tensor D2qoDs2 = EPS->getStress().d2qoverds2(); // d2q/ds2
  tensor D2thetaoDs2 = EPS->getStress().d2thetaoverds2(); // d2(theta)/ds2
//  double a1 = -6*sin(temp_phi)/(3.0-sin(temp_phi));
//  double a2 = -6*temp_cohesive*cos(temp_phi)/(3.0-sin(temp_phi));
  double e = (3.0-a1)/(3.0+a1);
  double Frou = g_0(theta, e); // r(theta)
  double Frou_prime = g_prime(theta, e); // r'(theta)
  double Frou_second = g_second(theta, e); // r"(theta)
//  double dQoverdp = a1; // dQ/dp
//  double dQoverdq = Frou; // dQ/dq
//  double dQoverdtheta = q*Frou_prime; // dQ/d(theta)
  double dQoverdp = alfa*(-3.0); // dQ/dp
  double dQoverdq = Frou/1.7320508076; // dQ/dq
  double dQoverdtheta = q*Frou_prime/1.7320508076; // dQ/d(theta)
//  double a23 = Frou_prime; // d2Q/dqd(theta)
//  double a32 = a23; // d2Q/d(theta)dq
//  double a33 = q * Frou_second; // d2Q/d(theta)2
  double a23 = Frou_prime/1.7320508076; // d2Q/dqd(theta)
  double a32 = a23/1.7320508076; // d2Q/d(theta)dq
  double a33 = q * Frou_second/1.7320508076; // d2Q/d(theta)2  
  
//  d2Qoverds2 =  DthetaoDs("mn") * DqoDs("pq")     *  a23 +
//                DqoDs("mn") * DthetaoDs("pq")     *  a32 +
//                DthetaoDs("mn") * DthetaoDs("pq") *  a33 +
//                D2poDs2("pqmn")  * dQoverdp +
//                D2qoDs2("pqmn")  * dQoverdq +
//                D2thetaoDs2("pqmn") * dQoverdtheta;
  d2Qoverds2 =  DqoDs("pq") * DthetaoDs("mn")      *  a23 +
                DthetaoDs("pq") * DqoDs("mn")      *  a32 +
                DthetaoDs("pq") * DthetaoDs("mn")  *  a33 +
                D2poDs2("pqmn")  * dQoverdp +
                D2qoDs2("pqmn")  * dQoverdq +
                D2thetaoDs2("pqmn") * dQoverdtheta;
  return d2Qoverds2;

}

// For Consistent Algorithm, Z Cheng, Jan 2004
tensor RMC01PotentialSurface::d2Qodsds1(const EPState *EPS) const 
{  
  tensor I("I", 2, def_dim_2);
  tensor d2Qoverdsds1 = I;
  return d2Qoverdsds1;
}

//================================================================================
// friend ostream functions for output
//================================================================================

OPS_Stream & operator<< (OPS_Stream & os, const RMC01PotentialSurface & PS)
{
    os << "ROUNDED MC Potential Surface Parameters: " << endln;
    return os;
}


#endif

