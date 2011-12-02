/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           General platform for elaso-plastic constitutive model       #
#                    implementation                                              #
#                                                                                #
# CLASS:             EvolutionLaw_L_Eeq (on plastic equivalent strain)          #
#                                                                                #
# VERSION:                                                                       #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )    #
# TARGET OS:         DOS || UNIX || . . .                                        #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                 #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                                 #
#                                                                                #
#                                                                                #
# DATE:              09-02-2000                                                  #
# UPDATE HISTORY:                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
# SHORT EXPLANATION: This is a linear evolution law for the evoltion of a        #
#                    scalar variable k which depends on plastic equi. strain     #
#                    i.e. dk = a*de_eq_p                                         #
//================================================================================
*/

#ifndef EL_LEeq_CPP
#define EL_LEeq_CPP

#include "EL_LEeq.h"
#include <basics.h>
    
//================================================================================
// Default constructor
//================================================================================
EvolutionLaw_L_Eeq::EvolutionLaw_L_Eeq( double ad ) 
:a(ad)
{}     

//================================================================================
// Copy constructor
//================================================================================

EvolutionLaw_L_Eeq::EvolutionLaw_L_Eeq(const EvolutionLaw_L_Eeq &LE ) 
{
    this->a = LE.geta();
}


//================================================================================
//  Create a clone of itself 
//================================================================================
// Alpha machine has problem on this
//EvolutionLaw_L_Eeq * EvolutionLaw_L_Eeq::newObj() {
EvolutionLaw_S * EvolutionLaw_L_Eeq::newObj() 
{    
    EvolutionLaw_S *newEL = new EvolutionLaw_L_Eeq( *this );
    
    return newEL;
}

////================================================================================
////  Initialize some  vars in EPState				        
////  nothing						       	        
////================================================================================
//
//void EvolutionLaw_L_Eeq::InitVars(EPState  *EPS) {
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
//void EvolutionLaw_L_Eeq::setInitD(EPState  *EPS) {
//
//}   
	  	  
 
//================================================================================
// Evaluating h_s = a * pow( 2.0*Rij_dev * Rij_dev/3.0, 0.5) (For the evaluation of Kp)
//================================================================================

double EvolutionLaw_L_Eeq::h_s( EPState *EPS, PotentialSurface *PS){

    //=========================================================================
    // Getting de_eq / dLambda
    stresstensor dQods = PS->dQods( EPS );
    //dQods.reportshort("dQods");

    //Evaluate the norm of the deviator of dQods
    //temp1 =  dQods("ij")*dQods("ij");
    tensor dQods_dev = dQods.deviator();
    tensor temp1 =  dQods_dev("ij")*dQods_dev("ij");
    double norm_dQods_dev =  pow( temp1.trace(), 0.5 );
    
    double de_eqodL = pow( 2.0 / 3.0, 0.5 ) * norm_dQods_dev;

    //Evaluating dSodeeq
    double dSodeeq = geta();

    double h = dSodeeq * de_eqodL;

    return h;

}


//================================================================================
//  Print vars defined in Linear Evolution Law
//================================================================================
void EvolutionLaw_L_Eeq::print()
{
    opserr << (*this);
}


//================================================================================
double EvolutionLaw_L_Eeq::geta() const
{       
    return a;
}

//================================================================================
OPS_Stream& operator<< (OPS_Stream& os, const EvolutionLaw_L_Eeq & LEL)
{
  //    os.unsetf( ios::scientific );
    os.precision(5);

    os.width(10);       
    os << endln << "Linear Scalar Evolution Law's parameters:" << endln;
    os << "a = " << LEL.geta() << "; " << endln;
           
    return os;
}  

#endif

