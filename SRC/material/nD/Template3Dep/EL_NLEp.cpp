/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           Cam clay model evolution law                    		 #
#                                                                                #
#                                                                                #
# CLASS:             EvolutionLaw_NL_Ep (on plastic volumetric strain)           #
#                                                                                #
# VERSION:                                                                       #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )    #
# TARGET OS:         DOS || UNIX || . . .                                        #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                 #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                                 #
#                                                                                #
#                                                                                #
# DATE:              Mar. 28, 2001                                               #
# UPDATE HISTORY:                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
# SHORT EXPLANATION: This is a nonlinear evolution law for the evoltion of a     #
#                    scalar variable po which depends on plastic vol. strain     #
#                    i.e. dpo = (1+eo)po/(lamda-kappa)*de_p                      #
//================================================================================
*/

#ifndef EL_NLEp_CPP
#define EL_NLEp_CPP

#include "EL_NLEp.h"
#include <basics.h>
    
//================================================================================
// Default constructor
//================================================================================
EvolutionLaw_NL_Ep::EvolutionLaw_NL_Ep( double eod, double lambdad, double kappad ) 
:eo(eod), lambda(lambdad), kappa(kappad)
{}     

//================================================================================
// Copy constructor
//================================================================================

EvolutionLaw_NL_Ep::EvolutionLaw_NL_Ep(const EvolutionLaw_NL_Ep &LE ) 
{
    this->eo     = LE.geteo();
    this->lambda = LE.getlambda();
    this->kappa  = LE.getkappa();
}


//================================================================================
//  Create a clone of itself 
//================================================================================
// Alpha machine has problem on this
//EvolutionLaw_NL_Ep * EvolutionLaw_NL_Ep::newObj() {
EvolutionLaw_S * EvolutionLaw_NL_Ep::newObj() 
{    
    EvolutionLaw_S *newEL = new EvolutionLaw_NL_Ep( *this );
    
    return newEL;
}

////================================================================================
////  Initialize some  vars in EPState				        
////  nothing						       	        
////================================================================================
//
//void EvolutionLaw_NL_Ep::InitVars(EPState  *EPS) {
//
//    // set initial E_Young corresponding to current stress state
//    //double p_atm = 100.0; //Kpa atmospheric pressure
//    //double p = EPS->getStress().p_hydrostatic();
//    //double E = EPS->getEo() * pow( (p/p_atm), geta());
//    EPS->setE( EPS->getEo() );
//      
//}   


//================================================================================
//  Set initial value of D once the current stress hit the yield surface     	
//  for MD model only    						       	
//================================================================================
//
//void EvolutionLaw_NL_Ep::setInitD(EPState  *EPS) {
//
//}   
	  	  
 
//================================================================================
// Evaluating h_s = a * pow( 2.0*Rij_dev * Rij_dev/3.0, 0.5) (For the evaluation of Kp)
//================================================================================

double EvolutionLaw_NL_Ep::h_s( EPState *EPS, PotentialSurface *PS){

    //=========================================================================
    // Getting de_eq / dLambda
    stresstensor dQods = PS->dQods( EPS );
    //dQods.reportshort("dQods");

    //Evaluate the norm of the deviator of dQods
    //temp1 =  dQods("ij")*dQods("ij");
    double dQods_p = dQods.p_hydrostatic();
    
    double de_podL = dQods_p;

    //Evaluating dSodeeq

    double po = EPS->getScalarVar( 1 );
    double dSodep = (1.0+geteo())*po/( getlambda()-getkappa() );

    double h = dSodep * de_podL;

    return h;

}


//================================================================================
//  Print vars defined in Linear Evolution Law
//================================================================================
void EvolutionLaw_NL_Ep::print()
{
    opserr << (*this);
}


//================================================================================
double EvolutionLaw_NL_Ep::geteo() const
{       
    return eo;
}

//================================================================================
double EvolutionLaw_NL_Ep::getlambda() const
{       
    return lambda;
}

//================================================================================
double EvolutionLaw_NL_Ep::getkappa() const
{       
    return kappa;
}


//================================================================================
OPS_Stream& operator<< (OPS_Stream& os, const EvolutionLaw_NL_Ep & LEL)
{
  //    os.unsetf( ios::scientific );
    os.precision(5);

    os.width(10);       
    os << endln << "Nonlinear Scalar Evolution Law(Cam Clay model)'s parameters:" << endln;
    os << "eo = " << LEL.geteo() << "; " << endln;
    os << "lambda = " << LEL.getlambda() << "; " << endln;
    os << "kappa = " << LEL.getkappa() << "; " << endln;
           
    return os;
}  

#endif

