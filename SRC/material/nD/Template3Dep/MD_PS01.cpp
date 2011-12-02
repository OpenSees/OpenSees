
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Manzari-Dafalias  potential surface 01(with Pc)           #
//#                    (Ref. Geotechnique v.47 No.2 255-272, 1997)               #
//# CLASS:             MDPotentialSurface01                                      #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, ZHaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              August 03 '93                                             #
//# UPDATE HISTORY:    August 08 '00                                             #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//================================================================================


#ifndef MD_PS01_CPP
#define MD_PS01_CPP

#include "MD_PS01.h"
#include <basics.h>
#include <math.h>

//================================================================================
// Normal constructor
//================================================================================

MDPotentialSurface01::MDPotentialSurface01(double pc ) 
{
  Pc = pc;
} 


//================================================================================
//create a colne of itself
//================================================================================

PotentialSurface * MDPotentialSurface01::newObj() {  

     PotentialSurface  *new_PS = new MDPotentialSurface01(Pc);
     return new_PS;

}

//================================================================================
//  tensor dQ/dsigma_ij: the normal to the potential surface
//================================================================================

tensor MDPotentialSurface01::dQods(const EPState *EPS) const { 

  tensor dQoverds( 2, def_dim_2, 0.0);
  tensor I2("I", 2, def_dim_2);

  tensor S = EPS->getStress().deviator();
  double p = EPS->getStress().p_hydrostatic();
  p = p -  Pc;

  tensor alpha = EPS->getTensorVar( 1 ); //the first tensor hardening var is alpha_ij

  //double m = EPS->getScalarVar( 1 );	 //the first scalar hardening var is m
  double D = EPS->getScalarVar( 2 );	 // D is stored in scalar var array's second cell
  //printf("Here we go!  D= %e\n", D);

  tensor r = S *(1.0/p);
  tensor temp_f1 = r - alpha;
  tensor temp_f2 = temp_f1("ij") * temp_f1("ij");
  double temp_f3 = sqrt( temp_f2.trace() );

  stresstensor n;
  if ( temp_f3 >= d_macheps() ){ 
    n = temp_f1 * (1.0 / temp_f3 );
  }
  else {
    opserr << "MDPotentialSurface01::dQods  |n_ij| = 0, divide by zero! Program exits.\n";
    exit(-1);
    //::printf(" \n\n n_ij not defined!!!! Program exits\n");
    //exit(1);
  }
  
  dQoverds =  n + I2 * (D *(1.0/3.0));

  return dQoverds;
}


//================================================================================
// tensor d2Q/dsigma_ij2: the second derivatives to the potential surface
//================================================================================

tensor MDPotentialSurface01::d2Qods2(const EPState *EPS) const 
{

  tensor d2Qoverds2;
  tensor I2("I", 2, def_dim_2);

  stresstensor stress = EPS->getStress();
  
  //double J2D = stress.Jinvariant2();
  //double q     = stress.q_deviatoric();
  double theta = stress.theta();

  tensor S = stress.deviator();
  double p = stress.p_hydrostatic();
  p = p -  Pc;
  opserr << " p = " << p << endln;
  
  tensor alpha = EPS->getTensorVar( 1 ); //the first tensor hardening var is alpha_ij
	   
  tensor r = S *(1.0/p);
  tensor temp_f1 = r - alpha;
  tensor temp_f2 = temp_f1("ij") * temp_f1("ij");
  double temp_f3 = sqrt( temp_f2.trace() );

  stresstensor n;
  if ( temp_f3 >= d_macheps() ){ 
    n = temp_f1 * (1.0 / temp_f3 );
  }
  else {
    opserr << "MDPotentialSurface01::dQods  |n_ij| = 0, divide by zero! Program exits.\n";
    exit(-1);
    //::printf(" \n\n n_ij not defined!!!! Program exits\n");
    //exit(1);
  }
  
  
  //tensor d2Qoverds2( 2, def_dim_2, 0.0); // dummy second derivatives. To be redefined.
  tensor dnds =  dnods( EPS);
  double A = 2.64;
  double Mc = 1.14;
  double kcd = 4.2;

  if (p < 0.0)
  {
    opserr << "MDPotentialSurface01::d2Qods2  p < 0, Program exits.\n";
     exit(-1);
  }
  double ec = 0.80 - 0.025 * log( p / 160 );

  double e = EPS->getScalarVar(3);
  double xi = e - ec;

  double c = 1.0;
  double cd = 0.0167;

  double dgodthetac  = dgoverdt(theta, c);
  double dgodthetacd = dgoverdt(theta, cd);
  tensor dthetaods = dthetaoverds(EPS);
  //stresstensor d = EPS->getTensorVar( 3 );   // getting  d_ij from EPState
  
  //tensor dtods = dthetaods("pq"); 

  tensor dDds1 = A*dthetaods*(dgodthetac*Mc+dgodthetacd*kcd*xi)*sqrt(2.0/3.0); 
  tensor dDds2 = A*apqdnods(EPS); 
  tensor dDds = dDds1 - dDds2;
  dDds.null_indices(); 
   
  d2Qoverds2 = dnds + I2("pq")*dDds("mn")*(0.33333333);
  d2Qoverds2.null_indices();

  return d2Qoverds2;
}

//================================================================================
// tensor dn_pq/dsigma_mn
//================================================================================
tensor MDPotentialSurface01::dnods(const EPState *EPS) const
{
  tensor dnods( 2, def_dim_2, 0.0);
  tensor I2("I", 2, def_dim_2);

  stresstensor S = EPS->getStress().deviator();
  //S.reportshort("S");
  stresstensor alpha = EPS->getTensorVar( 1 );

  double p = EPS->getStress().p_hydrostatic();
  p = p -  Pc;

  //double m = EPS->getScalarVar( 1 );
  //double fac = 1.0/(sqrt(2.0/3.0)*m);

  tensor Ipmnq = I2("pm") * I2("nq");
  tensor Imnpq = I2("mn") * I2("pq");
  tensor Apqmn = alpha("pq") * I2("mn");
  tensor X  = Ipmnq - Imnpq * (1.0/3.0) - Apqmn * (1.0/3.0);

  //Ipmnq.print("in MD_PS01", "");
  
  tensor sa = S("ij") * alpha("ij");
  double sad =sa.trace();
  tensor aa = alpha("ij") * alpha("ij");
  double aad =aa.trace();
  tensor Y = S - ( p + 0.333333333*(sad-p*aad) )*I2;
  Y.null_indices();

  tensor s_bar = S("ij") - p*alpha("ij");
  s_bar.null_indices();
  tensor norm2 = s_bar("ij") * s_bar("ij");
  double norm =  norm2.trace();

  s_bar.null_indices();
  tensor tmp = s_bar("pq")*Y("mn");
  dnods = ( X - 2.0 * tmp)*(1.0/norm);
  
  return  dnods;
}

//================================================================================
// tensor alpha_pq*dnods
//================================================================================
tensor MDPotentialSurface01::apqdnods(const EPState *EPS) const
{
  tensor ddnods( 2, def_dim_2, 0.0);
  tensor I2("I", 2, def_dim_2);

  stresstensor S = EPS->getStress().deviator();
  //S.reportshort("S");

  double p = EPS->getStress().p_hydrostatic();
  p = p -  Pc;

  //printf("Here we go!  p %f\n", p);
	      
  stresstensor alpha = EPS->getTensorVar( 1 );
  //stresstensor d = EPS->getTensorVar( 3 );   // getting  d_ij from EPState

  tensor akk = alpha("ij") * I2("ij");
  double akkd =akk.trace();

  tensor aa = alpha("ij") * alpha("ij");
  double aad =aa.trace();

  tensor aX = alpha - I2*(akkd +aad)*0.3333333;
	         
  //Ymn
  tensor sa = S("ij") * alpha("ij");
  double sad =sa.trace();
  tensor Y = S - ( p + 0.333333333*(sad-p*aad) )*I2;
  Y.null_indices();
  
  tensor s_bar = S - p*alpha;
  s_bar.null_indices();

  tensor as_bar = alpha("pq") * s_bar("pq");
  double as_bard = as_bar.trace();
  s_bar.null_indices();
  
  //Norm
  tensor norm2 = s_bar("ij") * s_bar("ij");
  double norm = norm2.trace();
  
	   	     
  ddnods = (aX - Y*2.0*as_bard)*(1.0/norm);
  return  ddnods;
}


//================================================================================
// tensor dg(theta, c)/dsigma_mn
//================================================================================
double MDPotentialSurface01::dgoverdt(double theta, double c) const
{
   //stresstensor stress = EPS->getStress();
  
   //double J2D = stress.Jinvariant2();
   //double q     = stress.q_deviatoric();
   //double theta = stress.theta();
  
   double temp = (-6*(1 - c)*c*sin(3*theta))/pow(1 + c - (1 - c)*cos(3*theta),2);  
   return temp;
}

//================================================================================
// tensor dtheta/dsigma_pq
//================================================================================
tensor MDPotentialSurface01::dthetaoverds(const EPState *EPS) const
{
   tensor ret(2, def_dim_2, 0.0);
   stresstensor s( 0.0);
   stresstensor t( 0.0);
   tensor I2("I", 2, def_dim_2);

   //double EPS = pow(d_macheps(),(1./3.));
   stresstensor stress = EPS->getStress();
  
   double J2D = stress.Jinvariant2();
   double q     = stress.q_deviatoric();
   double theta = stress.theta();

   //out    while ( theta >= 2.0*PI )
   //out      theta = theta - 2.0*PI; // if bigger than full cycle
   //out    while ( theta >= 4.0*PI/3.0 )
   //out      theta = theta - 4.0*PI/3.0; // if bigger than four thirds of half cycle
   //out    while ( theta >= 2.0*PI/3.0 )
   //out      theta = theta - 2.0*PI/3.0; // if bigger than two third of half cycle
   //out    while ( theta >= PI/3.0 )
   //out      theta = 2.0*PI/3.0 - theta; // if bigger than one third of half cycle
   //out
   //out    if ( theta < 0.0001 )
   //out      {
   //out        ::printf("theta = %.10e in Material_Model::dthetaoverds(stress)\n",
   //out                           theta);
   //out//>><<>><<>><<        return ret;
   //out      }
   
   double c3t = cos(3.0*theta);
   double s3t = sin(3.0*theta);

   double tempS = (3.0/2.0)*(c3t/(q*q*s3t));
   double tempT = (9.0/2.0)*(1.0/(q*q*q*s3t));

   s = stress.deviator();
   t = s("qk")*s("kp") - I2*(J2D*(2.0/3.0));

   s.null_indices();
   t.null_indices();

   ret = s*tempS - t*tempT;
   ret.null_indices();
   return ret;
}

OPS_Stream& operator<< (OPS_Stream& os, const MDPotentialSurface01 &PS)
{
   os << "Manzari-Dafalias Potential Surface 01( with Pc) Parameters: " << endln;
   os << "Pc = " << PS.Pc << endln;
   return os;
}



#endif

