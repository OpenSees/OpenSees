
//================================================================================
// COPYRIGHT (C):     :-))                                                       |
// PROJECT:           Object Oriented Finite Element Program                     |
// PURPOSE:           General platform for elaso-plastic constitutive model      |
//                    implementation                                             |
//                                                                               |
// CLASS:             EvolutionLaw_NL_EijMD (nonlinear tensorial Evolution law)  |
//                                                                               |
//                                                                               |
// VERSION:                                                                      |
// LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )   |
// TARGET OS:         DOS || UNIX || . . .                                       |
// DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                |
// PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                                |
//                                                                               |
//                                                                               |
// DATE:              09-13-2000                                                 |
// UPDATE HISTORY:                                                               |
//                                                                               |
//                                                                               |
//                                                                               |
//                                                                               |
// SHORT EXPLANATION: This is a nonlinear evolution law for the evoltion of a    |
//                    tensorial variable alpha which depends on plastic strain   |
//                    i.e. dalpha = (3/2)^0.5*h*[2*alpha_b*de_ij_p - de_eq*	     |
//                    alpha_ij (Armstrong- Frederick Model                       |
//                                                                               |
//================================================================================

#ifndef EL_NL_EIJMD_H
#define EL_NL_EIJMD_H

#include <math.h>
#include "EL_T.h"
#include <iostream>
using std::ostream;


class EvolutionLaw_NL_EijMD : public EvolutionLaw_T
{
  // Private vars to define the MD evolution law
  private:			      
    // the exponential in G = Go (P/P_at)^a and K = Ko (P/P_at)^a
    //double eo; //Initial void ratio  ...Moved to EPState Joey 02-12-03
    //double a; 		       ...Moved to EPState Joey 02-12-03

    //Critical state parameters
    double Mc;
    double Me;
    //double Lambda; // slope of e vs. ln p ...Moved to EPState Joey 02-12-03
    //double ec_ref; // critical void ratio at reference mean effective stress p  ...Moved to EPState Joey 02-12-03
    //double p_ref;  // critical void ratio at reference mean effective stress p  ...Moved to EPState Joey 02-12-03

    //surface evolution parameters
    double kc_b;  // b stands for bounding surface
    double kc_d;  // d stands for dilatancy surface
    double ke_b;  // e extension
    double ke_d;  // c compression

    //Hardening parameters
    //double h; // Being calculated using ho and b_ij * n_ij
    double ho;
    double Cm;
    //double e;     // current void ratio ...Moved to EPState Joey 02-12-03

    //Dilatancy parameter
    double D; //also copied to EPS's second scalar var(no direct contribution to hardening)
    double Ao;

    //Get rid of the fabric tensor
    //Fabric parameters
    stresstensor F;   //Fabric tensor which will evolve like alpha_{ij}
    double Fmax;
    double Cf;
    
  public:
    //EvolutionLaw_L_EijMD( );    // default constructor---no parameters
    
    // default constructor
    EvolutionLaw_NL_EijMD( 
                           //double eod = 0.85,
    			   //double ad  = 0.5,    
    			   double Mcd = 1.14,//1.14, 
    			   double Med = 1.14,//1.14, 
    			   //double Lamdad = 0.025,
    			   //double ec_refd = 0.8, 
    			   //double p_refd = 160.0, 
    			   double kc_bd = 3.975, 
    			   double kc_dd = 4.200, 
    			   double ke_bd = 2.000, 
    			   double ke_dd = 0.07,  
    			   double hod = 1200,	 
    			   double Cmd = 0.00,
    			   double Aod = 2.64,    
    			   double Fmaxd = 100,   
    			   double Cfd = 100);
    			   //double ed  = 0.85,    
                           
    // Copy constructor
    EvolutionLaw_NL_EijMD(const EvolutionLaw_NL_EijMD &LEL );   
    
    //create a clone of itself
    //EvolutionLaw_T *newObj();
    EvolutionLaw_T *newObj();
    	  
    //double h( EPState *EPS,  double norm_dQods);     // Evaluating hardening function h
    tensor h_t( EPState *EPS, PotentialSurface *PS);    // Evaluating hardening function h
    
    //Overwrite  updateEeDm
    //Updating  E, e and D value according to current mean effective stress p and elastic strain
    int updateEeDm(EPState *EPS, double st_vol, double dLamda);
       
    //void UpdateVar( EPState *EPS, double dlamda );  // Evolve corresponding var linearly using de_eq_p
    //Moved to CDriver.cpp

    void print();
    //g++ complaining if don't have this line
    virtual ~EvolutionLaw_NL_EijMD() {}

  private:
    // some accessor functions
		         
    // set D once current stress hits the y.s. or
    // Set D value according to current EPState
    //void setD(EPState *EPS);   

    //double geta() const;
    double getMc() const;
    double getMe() const;
    //double getLambda() const;
    //double getec_ref() const;
    //double getp_ref() const; 

    double getkc_b() const;  
    double getkc_d() const;  
    double getke_b() const;  
    double getke_d() const;  
    double getho() const;
    double getCm() const;
    //double geteo() const;
    //double gete() const;

    //Dilatancy parameter
    double getAo() const;
    double getD() const;

    double getFmax() const;
    stresstensor getF() const;
    double getCf() const;

    // Interpolation function  by Agyris
    double g_A(double theta, double e);	 

    // Interpolation function  by Willan-Warkne
    double g_WW(double theta, double e); 
	       
    //================================================================================
    // Overloaded Insertion Operator	  Zhaohui Added Aug. 13, 2000
    // prints Linear EvolutionLaw's contents 
    //================================================================================
    friend OPS_Stream& operator<< (OPS_Stream& os, const EvolutionLaw_NL_EijMD & MDEL);
    //friend ostream& operator<< (ostream& os, const EvolutionLaw_NL_EijMD & MDEL);

    
};


#endif




