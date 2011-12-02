
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Manzari-Dafalias  potential surface                       #
//#                    (Ref. Geotechnique v.47 No.2 255-272, 1997)               #
//# CLASS:             MDPotentialSurface                                        #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, ZHaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              
//# UPDATE HISTORY:    
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//================================================================================


#ifndef MD_PS_CPP
#define MD_PS_CPP

#include "MD_PS.h"
#include <basics.h>
#include <math.h>

#define sqrt23rd 0.8165

//================================================================================
// Normal constructor
//================================================================================

MDPotentialSurface::MDPotentialSurface( ) { } 


//================================================================================
//create a clone of itself
//================================================================================

PotentialSurface * MDPotentialSurface::newObj() {  

     PotentialSurface  *new_PS = new MDPotentialSurface();
     return new_PS;

}

//================================================================================
//  tensor dQ/dsigma_ij: the normal to the potential surface
//================================================================================

tensor MDPotentialSurface::dQods(const EPState *EPS) const { 

  tensor dQoverds( 2, def_dim_2, 0.0);
  tensor I2("I", 2, def_dim_2);

  tensor S = EPS->getStress().deviator();
  //double p = EPS->getStress().p_hydrostatic();
  
  tensor alpha = EPS->getTensorVar( 1 ); //the first tensor hardening var is alpha_ij

  //double m = EPS->getScalarVar( 1 );	 //the first scalar hardening var is m
  double D = EPS->getScalarVar( 2 );	 // D is stored in scalar var array's second cell

  stresstensor n = EPS->getTensorVar( 2 );
  //stresstensor n;

  ////-------------------------------------------------
  //stresstensor r = S *(1.0/p);
  //stresstensor r_bar = r - alpha;
  //stresstensor norm2 = r_bar("ij") * r_bar("ij");
  //double norm = sqrt( norm2.trace() );
  //if ( norm >= d_macheps() ){ 
  //  n = r_bar *(1.0 / norm );
  //}
  //else {
  //  opserr->fatal("MDYieldSurface::dQods |n_ij| = 0, divide by zero! Program exits.");
  //  exit(-1);
  //}
  //
  ////n = r_bar *(1.0 / sqrt23rd / m );
  ////-------------------------------------------------
  
  dQoverds =  n + I2 * D *(1.0/3.0);

  return dQoverds;
}


//================================================================================
// tensor d2Q/dsigma_ij2: the second derivatives to the potential surface
//================================================================================

tensor MDPotentialSurface::d2Qods2(const EPState *EPS) const 
{

  
  tensor d2Qoverds2;
  tensor I2("I", 2, def_dim_2);

  stresstensor stress = EPS->getStress();
  
  //double J2D = stress.Jinvariant2();
  //double q     = stress.q_deviatoric();
  double theta = stress.theta();

  tensor S = stress.deviator();
  double p = stress.p_hydrostatic();
  
  tensor alpha = EPS->getTensorVar( 1 ); //the first tensor hardening var is alpha_ij
	   
  /*
  tensor r = S *(1.0/p);
  stresstensor r_bar = r - alpha;
  stresstensor norm2 = r_bar("ij") * r_bar("ij");
  double norm = sqrt( norm2.trace() );

  stresstensor n;
  if ( norm >= d_macheps() ){ 
    n = r_bar *(1.0 / norm );
  }
  else {
    opserr << "MDPotentialSurface::dQods  |n_ij| = 0, divide by zero! Program exits.\n";
    exit(-1);
  }
  */
  
  //tensor d2Qoverds2( 2, def_dim_2, 0.0); // dummy second derivatives. To be redefined.
  tensor dnds =  dnods( EPS);
  double A = 2.64;
  double Mc = 1.14;
  double kcd = 4.2;
  //double m = EPS->getScalarVar( 1 );
  if (p < 0.0)
  {
     opserr << "MDPotentialSurface::d2Qods2  p < 0, Program exits.\n";
     exit(-1);
  }

  //double ec = 0.80 - 0.025 * log( p / 160 ); //p_ref = 160; (ec)ref = 0.8
  //double ec = (EPS->getec()) - (EPS->getLam()) * log( p / (EPS->getpo()) );
  //opserr << " ************" << log(2.718) << "\n";

  //double e = EPS->gete();
  double xi = EPS->getpsi();

  double c = 1.0;
  double cd = 0.0167; //0.7;

  double dgodthetac  = dgoverdt(theta, c);
  double dgodthetacd = dgoverdt(theta, cd);
  tensor dthetaods = dthetaoverds(EPS);
  //stresstensor d = EPS->getTensorVar( 3 );   // getting  d_ij from EPState
  
  //tensor dtods = dthetaods("pq"); 

  tensor dDds1 = dthetaods*A*(dgodthetac*Mc+dgodthetacd*kcd*xi)*sqrt(2.0/3.0);
  tensor dDds2 = apqdnods(EPS)*A;
  tensor dDds = dDds1 - dDds2;
  dDds.null_indices(); 
   
  d2Qoverds2 = dnds + dDds("mn")*I2("pq")*(1.0/3.0);
  d2Qoverds2.null_indices();
  

  tensor d2Qoverds2x( 4, def_dim_2, 0.0); // dummy second derivatives. To be redefined.
  return d2Qoverds2x;
  //return dnds;
}

//================================================================================
// tensor dn_pq/dsigma_mn
//================================================================================
tensor MDPotentialSurface::dnods(const EPState *EPS) const
{
  /*
  tensor xdnods;
  tensor I2("I", 2, def_dim_2);

  stresstensor S = EPS->getStress().deviator();
  //S.reportshort("S");

  double p = EPS->getStress().p_hydrostatic();
  //printf("Here we go!  p %f\n", p);
  double m = EPS->getScalarVar( 1 );

  double fac = 1.0/(sqrt(2.0/3.0)*m);

  tensor Ipmnq = I2("pm") * I2("nq");
  tensor Imnpq = I2("mn") * I2("pq");
  tensor Spqmn = S("pq") * I2("mn");
  tensor temp  = Ipmnq - Imnpq*(1.0/3.0);
  
  xdnods = (temp*p - Spqmn*(1.0/3.0))*(fac/p/p);
  return  xdnods;
  */
  
  
  tensor dnods;
  tensor I2("I", 2, def_dim_2);

  stresstensor S = EPS->getStress().deviator();
  //S.reportshort("S");
  stresstensor alpha = EPS->getTensorVar( 1 );

  double p = EPS->getStress().p_hydrostatic();
  //p = p -  Pc;

  //double m = EPS->getScalarVar( 1 );
  //double fac = 1.0/(sqrt(2.0/3.0)*m);

  tensor Ipmnq = I2("pm") * I2("nq");
  tensor Imnpq = I2("mn") * I2("pq");
  tensor Apqmn = alpha("pq") * I2("mn");
  tensor X  = Ipmnq - Imnpq*(1.0/3.0)  - Apqmn*(1.0/3.0);

  //Ipmnq.print("in MD_PS01", "");
  
  tensor s_bar = S("ij") - alpha("ij")*p;

  tensor ad = alpha - I2; //Vlado found a bug
  tensor sa = S("ij") * ad("ij");
  double sad =sa.trace();
  tensor aa = alpha("ij") * ad("ij");
  double aad =aa.trace();
  tensor Y = s_bar +I2 * (1.0/3.0) * (sad-p*aad);
  Y.null_indices();

  s_bar.null_indices();
  tensor t_norm2 = s_bar("ij") * s_bar("ij");
  double norm2 = t_norm2.trace();
  double norm =  sqrt( norm2 );

  double norm21 = 1.0/( norm2 );
  s_bar.null_indices();
  tensor tmp = s_bar("pq")*Y("mn");
  dnods = ( X - tmp*norm21)*(1.0/norm);
  
  return  dnods;
  

}

//================================================================================
// tensor alpha_pq*dnods
//================================================================================
tensor MDPotentialSurface::apqdnods(const EPState *EPS) const
{
  /* old fomulation, not correct
  tensor ddnods( 2, def_dim_2, 0.0);
  tensor I2("I", 2, def_dim_2);

  stresstensor S = EPS->getStress().deviator();
  //S.reportshort("S");

  double p = EPS->getStress().p_hydrostatic();
  //printf("Here we go!  p %f\n", p);
  double m = EPS->getScalarVar( 1 );

  double fac = 1.0/(sqrt(2.0/3.0)*m);

  stresstensor d = EPS->getTensorVar( 3 );   // getting  d_ij from EPState
  tensor dkk = d("ij") * I2("ij");
  double dkkd =dkk.trace();
  tensor dkkImn = dkkd * I2;

  tensor Spqdpq = S("pq") * d("pq");
  double Sqpdpqd =Spqdpq.trace();
  tensor dSI = Sqpdpqd*I2;
     
  tensor temp  = d  - (1.0/3.0)*dkkImn;
  
  ddnods = fac/p*(temp - 1.0/(p*3.0)*dSI);
  return  ddnods;
  */

  /*
  tensor ddnods( 2, def_dim_2, 0.0);
  tensor I2("I", 2, def_dim_2);

  stresstensor S = EPS->getStress().deviator();
  //S.reportshort("S");

  double p = EPS->getStress().p_hydrostatic();
  //p = p -  Pc;

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
  tensor Y = S - I2*( p + 0.333333333*(sad-p*aad) );
  Y.null_indices();
  
  tensor s_bar = S - alpha*p;
  s_bar.null_indices();

  tensor as_bar = alpha("pq") * s_bar("pq");
  double as_bard = as_bar.trace();
  s_bar.null_indices();
  
  //Norm
  tensor norm2 = s_bar("ij") * s_bar("ij");
  double norm = norm2.trace();
  
	   	     
  ddnods = (aX - Y*2.0*as_bard)*(1.0/norm);
  return  ddnods;
  */

  tensor xdnods = dnods(EPS);
  stresstensor alpha = EPS->getTensorVar( 1 );

  tensor adnods = alpha("pq") * xdnods("pqmn");

  return adnods;

}


//================================================================================
// tensor dg(theta, c)/dsigma_mn
//================================================================================
double MDPotentialSurface::dgoverdt(double theta, double c) const
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
tensor MDPotentialSurface::dthetaoverds(const EPState *EPS) const
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




#endif

