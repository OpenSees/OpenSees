
//################################################################################
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           OpenSees                                                  #
//# PURPOSE:           Triaxial Failure Criterion for Concrete - yield criterion #
//# CLASS:             TriFCYieldSurface                                         #
//#                                                                              #
//# VERSION:           1.0                                                       #
//# LANGUAGE:          C++ (ili tako nesto)                                      #
//# TARGET OS:                                                                   #
// DESIGNER(S):       Boris Jeremic and Zhaohui Yang [jeremic,zhyang]@ucdavis.edu| 
// PROGRAMMER(S):     Vlado Vukadin                                              |
//#                                                                              #
//#                                                                              #
//# DATE:             June 01, 2002                                              #
//# UPDATE HISTORY:    bice tako dobr da nece biti potreban update :)            #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//# SHORT EXPLANATION:                                                           #
//#                                                                              #
//# Yield surface is based on article by Menetrey, P. and William, K.J.          #
//# published in 1995 in  ACI Structural Journal pp 311-318. Purpose of the     #
//# Yield surface is to model triaxial strenght of concrete. All the necessary   #
//# explanations about functions are there                                       #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/


#ifndef Tri_a_fail_crit_YS_CPP
#define Tri_a_fail_crit_YS_CPP

#include "Tri_a_fail_crit_YS.h"


//================================================================================
// Normal constructor
//================================================================================

TriFCYieldSurface::TriFCYieldSurface (double fc, double ft, double e, double coh )
{ 

     fcomp = fc;      		// a way to fix those two factors so they are not	
     ftens = ft;	     	// in need for evolution law, do not need to be derived and could be
     el    =  e;		// suplied as a constansts (how this is done look in definition of Cam Clay model)
     c     = coh;
}
//================================================================================
// Destructor
//================================================================================

TriFCYieldSurface::~TriFCYieldSurface ( )
{

} 

//================================================================================
//create a clone of itself
//================================================================================

YieldSurface * TriFCYieldSurface ::newObj() 
{  
     YieldSurface  *new_YS = new TriFCYieldSurface (fcomp,ftens, el, c);
     return new_YS;
}

//================================================================================
//  Yield criterion evaluation function F(EPState)
//================================================================================

double TriFCYieldSurface ::f(const EPState *EPS) const 
{
  double xi    =  EPS->getStress().xi();	   //  functions to get  Haigh-Westergard   
  double rho   =  EPS->getStress().rho();	   //  stress invariants for   
  double th    =  EPS->getStress().theta();        //  explanation look Jeremic&Sture,1998  

  //double el    =  EPS->getScalarVar(1);	           // functions to select four parameters needed
  //double fcomp =  EPS->getScalarVar(3);	   // to define yield surface. All are called throug
  //double ftens =  EPS->getScalarVar(4);	   // EPS state so that scalar evolution laws could 
  //double c     =  EPS->getScalarVar(2);	           // be developed for all of them

  double a1    =  pow (1.5,(1./2.) ); 	                   // constants needed 
  double a2    =  pow (6.,(1./2.) ); 	                   // to make calculation of the yield surface 
  double a3    =  pow (3.,(1./2.) ); 	                   // more transparent 
  double a4    =  3.0*((fcomp * fcomp) - (ftens * ftens));
  double a5    =  (fcomp * ftens);
  double a6    =  (el/(1.0+ el));
  double m     =  a4*a6/a5;				   //factor   m (friction for explanation look in article)

            // William Warnke function         	
  double ww1   =  (  4.0 * (1.0 - el * el) * cos(th) * cos(th) + 
	                 (2.0 * el - 1.0) * (2.0 * el - 1.0)  );
  double ww2   =  pow ( ( (-4.0 * el) + (5.0 * el * el) + 
	                 (4.0 * (1.0 - el * el) * cos(th) * cos(th)) ), (1./2.) );
  double ww3   =  2.0 * (1.0 - el * el) * cos(th) + 
	                 (2.0 * el - 1.0) * ww2;
  double ww    =  ww1/ww3;                     	


	    // actual yield surface        	
  double yield_func =  (pow((a1*rho/fcomp),2.0) + 
	                      (m * ( (rho/(a2*fcomp)) * ww + 
																						 (xi/(a3*fcomp) ) ) ) - c);
  
  return yield_func; 
  
  		   
}

//================================================================================
// tensor dF/dsigma_ij  (this is the place wher a derivativ over yield function is coded)
//================================================================================

tensor TriFCYieldSurface ::dFods(const EPState *EPS ) const 
{
     tensor dFoverds( 2, def_dim_2, 0.0);
     
     	//double xi    =  EPS->getStress().xi();	           //  functions to get  Haigh-Westergard   
     	double rho   =  EPS->getStress().rho();	           //  stress invariants for   
     	double th    =  EPS->getStress().theta();          //  explanation look Jeremic&Sture,1998  
    
     
     	tensor DpoDs =    EPS->getStress().dpoverds();	   //  functions to get  already defined
     	tensor DqoDs =    EPS->getStress().dqoverds();	   //  drivatives of p, q and theta which are defined
     	tensor DtoDs =    EPS->getStress().dthetaoverds(); //  in stress.cpp file. With method getStress()  
     							   //  it is made sure that function is filled with curent stress 
     
     	// some parametres neccecary to define yield function and make calculation more transparent

     	double el    =  EPS->getScalarVar(1);		   // parameters called through      
     	//double c     =  EPS->getScalarVar(2);		   // definitions of 'set EPS'
     	double b1    =  pow (3.,(1./2.) ); 
     	double b2    =  pow (2.,(1./2.) ); 
     	//double b3    =  pow (1.5,(1./2.) );
     	double b4    =  pow (6.,(1./2.) ); 

     	double a4    =  3.0*((fcomp * fcomp) - (ftens * ftens));	 // variables neccecary to define 
     	double a5    =  (fcomp * ftens);				 // factor   m 
     	double a6    =  (el/(1.0+ el));				         // which is an expresion for 
     	double m     =  a4*a6/a5;					 // friction parameter
        
	
	// William warnke function
    
     	double ww1   =  4.0*(1.0-el*el)*cos(th)*cos(th) + (2.0*el-1.0)*(2.0*el-1.0);
     	double ww2   =  2.0*(1.0-el*el)*cos(th) + 
					                 (2.0*el-1)*pow((4.0*(1.0-el*el)*cos(th)*cos(th) + 
																					 5.0*el*el -4.0*el),(1./2.));
     	double ww    =  ww1/ww2; 

       	//derivation of the William Warnke fuction over theta
       	
       		// some factors enabling me to check the function. It might not be 
       		// the example of clearness but it was clear to me at the time I
       		// wrote this and to be honest function is really ugly 

     	
	double dWW1  =     (-1.0 + 2.0 * el) ;
     	double dWW2  =     (1.0 - el * el) ;
     	double dWW3  =     cos(th);
     	double dWW4  =     cos(th) * cos(th);
     	double dWW5  =     sin(th);
     	double dWW6  =     pow ( ( -4.0 * el + 5.0 * el * el + 4.0 * dWW2 * dWW4  ) , (1.0/2.0 )  );
     	double dWW7  =     -8.0 * dWW2 * dWW3 * dWW5;
     	double dWW8  =     2.0 * dWW2 * dWW3;
     	double dWW9  =     dWW1 * dWW6;
     	double dWW10 =	   dWW1 * dWW1 + 4.0 * dWW2 * dWW4;   
	double dWW11 =	   -2.0 * dWW2 * dWW5; 
	double dWW12 =	   (4.0 * dWW1 * dWW2 * dWW3 * dWW5) / dWW6;
	double dWW13 =     dWW10 * (dWW11 - dWW12);
	double dWW14 =     (2.0 * dWW2 * dWW3) + (dWW1 * dWW6);  
	 
     	//William Warnke fuction derived over theta

     	double dWWodth  =  (dWW7 / (dWW8 + dWW9)) - (dWW13 / (dWW14 * dWW14));  


     	double dFoverdp =  (  (3.0/b1) * m/(b1*fcomp) );				         // derivation of yield function over xi
     	double dFoverdq =  ((( 3.0*rho/(fcomp*fcomp)) + ( m*ww/(b4*fcomp)))* b2/b1   );    // derivation of yield function over rho
     	double dFoverdt =  (  (m*rho/(b4*fcomp) )*dWWodth    );	                 // derivation of yield function over theta
    
     dFoverds =  DpoDs   * dFoverdp +	  // factors are added to this function  to take into 
                 DqoDs  * dFoverdq +	  // account small differences between Haigh - Westergard
		 DtoDs  * dFoverdt;	  // coordinats and  defined p an q functions. Vnesel sem jih v vrste 174-175
					  // ker definica BJ tenzorja ne dopusca deljenje ce ni v tenzoski obliki
  return dFoverds;



}

//================================================================================
// double xi_s1 = dF/dS1  Derivative in terms of first scalar var  
//================================================================================

//double TriFCYieldSurface ::xi_s1(const EPState *EPS  ) const 
//{
//
//   	double rho   =  EPS->getStress().rho();	         //  stress invariants for   
//  	double th    =  EPS->getStress().theta();        //  explanation look Jeremic&Sture,1998  
// 
//  	double el    =  EPS->getScalarVar(1);	           
//   	//double c     =  EPS->getScalarVar(2);	           
// 
// 	//double a1    =  pow (1.5,(1./2.) ); 	                   // constants needed 
//   	double a2    =  pow (6.,(1./2.) ); 	                   // to make calculation of the yield surface 
//   	//double a3    =  pow (3.,(1./2.) ); 	                   // more transparent 
//   	double a4    =  3*((fcomp * fcomp) - (ftens * ftens));
//   	double a5    =  (fcomp * ftens);
//   	double a6    =  (el/(1+ el));
//   	double m     =  a4*a6/a5;	
//    	
 // 	    //derivation of the William Warnke function over el
//  // 
// 	double dWWoel1  =   (-1.0 + 2.0 * el) ;
// 	double dWWoel2  =   (1.0 - el * el) ;
// 	double dWWoel3  =   cos(th); 
// 	double dWWoel4  =   cos(th) * cos(th);
// 	double dWWoel5  =   (-4.0 + 10.0 * el - 8.0 * el * dWWoel4);  
// 	double dWWoel6  =   pow ( ( -4.0 * el + 5.0 * el * el + 4.0 * dWWoel2 * dWWoel4  ) , (1.0/2.0 )  );
// 	double dWWoel7  =   dWWoel1 + 4.0 * dWWoel2 * dWWoel4; 	
// 	double dWWoel8  =   4.0 * el * dWWoel3;
// 	double dWWoel9  =   dWWoel1 * dWWoel5/(2.0 * dWWoel6);    
// 	double dWWoel10 =   2.0 * dWWoel2 * dWWoel3 + dWWoel1 * dWWoel6;
// 	double dWWoel11 =   4.0 * dWWoel1 - 8.0 * el * dWWoel4; 
// 	
// 	double dWWoel  =    -((dWWoel7 * (dWWoel8 + dWWoel9 + 2.0 * dWWoel6))/(dWWoel10 * dWWoel10)) + (dWWoel11/dWWoel10);  
// 	
 // 	double dof1scalar =  (m *(rho/(a2*fcomp)) * dWWoel );
//   	       
//   return dof1scalar;
//   
// 
// }

//================================================================================
// double xi_s2 = dF/dS2 Derivative in terms of second scalar var 
//================================================================================

// double TriFCYieldSurface ::xi_s2( const EPState *EPS ) const 
// {
// 
//    return -1.0;
// }


//================================================================================
// double xi_t1 = dF/dt1  Derivative in terms of first tensorial var 
//================================================================================

//tensor TriFCYieldSurface ::xi_t1( ) const 
//{


//} 

//================================================================================
double TriFCYieldSurface::getfcomp() const
{       
    return fcomp;
}

double TriFCYieldSurface::getftens() const
{       
    return ftens;
}

double TriFCYieldSurface::getel() const
{       
    return el;
}

double TriFCYieldSurface::get_c() const
{       
    return c;
}


OPS_Stream& operator<< (OPS_Stream& os, const TriFCYieldSurface & YS)
{
 
   os << "Triaxial-Failure-Criterion-for-Concrete Yield Surface Parameters: " << endln;
   os << "fcomp = " << YS.fcomp << endln;
   os << "ftens = " << YS.ftens << endln;
   os << "el= " << YS.el << endln;
   os << "c= " << YS.c << endln;
   os.precision(4);
   return os;


}

#endif

