/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           General platform for elaso-plastic constitutive model       #
#                    implementation                                              #
#                                                                                #
# CLASS:             EvolutionLaw_L_Eij (on plastic equivalent strain)           #
#                                                                                #
# VERSION:                                                                       #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )    #
# TARGET OS:         DOS || UNIX || . . .                                        #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                 #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                                 #
#                                                                                #
#                                                                                #
# DATE:              09-13-2000                                                  #
# UPDATE HISTORY:                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
# SHORT EXPLANATION: This is a linear evolution law for the evoltion of a        #
#                    tensorial variable k which depends on plastic strain        #
#                    i.e. dalpha = a*de_ij_p                                     #
//================================================================================
*/

#ifndef EL_LEij_CPP
#define EL_LEij_CPP

#include "EL_LEij.h"
#include <basics.h>
    
//================================================================================
// Default constructor
//================================================================================
EvolutionLaw_L_Eij::EvolutionLaw_L_Eij( double ad) 
: a(ad)
{}     

//================================================================================
// Copy constructor
//================================================================================

EvolutionLaw_L_Eij::EvolutionLaw_L_Eij(const EvolutionLaw_L_Eij &LE ) 
{
    this->a = LE.geta();
}


//================================================================================
//  Create a clone of itself 
//================================================================================
EvolutionLaw_T * EvolutionLaw_L_Eij::newObj() {
    
    EvolutionLaw_T *newEL = new EvolutionLaw_L_Eij( *this );
    
    return newEL;

}

////================================================================================
////  Initialize some  vars in EPState				        
////  nothing						       	        
////================================================================================
//
//void EvolutionLaw_L_Eij::InitVars(EPState  *EPS) {
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
//  for L model only    						       	
//							       		        
//							       	                
//================================================================================
//
//void EvolutionLaw_L_Eij::setInitD(EPState  *EPS) {
//
//}   

 
//================================================================================
// Evaluating h_s = a * dSodeij = a*Rij (For the evaluation of Kp)
//================================================================================

tensor EvolutionLaw_L_Eij::h_t( EPState *EPS, PotentialSurface *PS){

    //=========================================================================
    // Getting de_ij / dLambda
    stresstensor dQods = PS->dQods( EPS );
    //dQods.reportshort("dQods");
    
    tensor de_ijodLam = dQods;

    //Evaluating dSodeeq
    double dSodeij = geta();

    tensor h = de_ijodLam*dSodeij;

    return h;

}


//================================================================================
//  Print vars defined in Linear Evolution Law
//================================================================================
void EvolutionLaw_L_Eij::print()
{
    opserr << (*this);
}


//================================================================================
double EvolutionLaw_L_Eij::geta() const
{       
    return a;
}

//================================================================================
OPS_Stream& operator<< (OPS_Stream& os, const EvolutionLaw_L_Eij & LEL)
{
  //    os.unsetf( ios::scientific );
  os.precision(5);

  os.width(10);       
  os << endln << "Linear Tensorial Evolution Law's parameters:" << endln;
  os << "a = " << LEL.geta() << "; " << endln;
           
    return os;
}  

#endif

