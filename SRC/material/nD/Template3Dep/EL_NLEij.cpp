/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           General platform for elaso-plastic constitutive model       #
#                    implementation                                              #
#                                                                                #
# CLASS:             EvolutionLaw_NL_Eij (on plastic strain)                     #
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
# SHORT EXPLANATION: This is a nonlinear evolution law for the evoltion of a     #
#                    tensorial variable alpha which depends on plastic strain    #
#                    i.e. dalpha = 2/3*ha*dE_ij -Cr*de_eq*alpha_ij(Amstrong-     #
//                   Frederick Model                                             #
//================================================================================
*/

#ifndef EL_NLEij_CPP
#define EL_NLEij_CPP

#include "EL_NLEij.h"
#include <basics.h>
    

//================================================================================
// Copy constructor
//================================================================================

EvolutionLaw_NL_Eij::EvolutionLaw_NL_Eij(const EvolutionLaw_NL_Eij &LE ) {

    this->ha = LE.getha();
    this->Cr = LE.getCr();
}


//================================================================================
//  Create a clone of itself 
//================================================================================
EvolutionLaw_T * EvolutionLaw_NL_Eij::newObj() {
    
    EvolutionLaw_NL_Eij *newEL = new EvolutionLaw_NL_Eij( *this );
    
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

////================================================================================
////  Updating corresponding internal variable           MOVED to CDriver.cpp       
////================================================================================
//
//void EvolutionLaw_L_Eij::UpdateVar( EPState *EPS, int WhichOne) {
//   
//    //=========================================================================
//    // Updating alfa1 by dalfa1 = a* de_eq
//    
//    // Calculate  e_eq = sqrt( 2.0 * epsilon_ij * epsilon_ij / 3.0)
//    straintensor pstrain =  EPS->getdPlasticStrain(); 
//    double e_eq  = pstrain.equivalent();
//    //cout << "e_eq = " << e_eq << endln;
//
//    double dS =  e_eq * geta();
//    double S  = EPS->getScalarVar( WhichOne );
//
//    EPS->setScalarVar(WhichOne, S + dS);
//
//}


 
//================================================================================
// Evaluating h_s = pow( 2.0*Rij_dev * Rij_dev/3.0, 0.5) (For the evaluation of Kp)
//================================================================================

tensor EvolutionLaw_NL_Eij::h_t( EPState *EPS, PotentialSurface *PS){

    //=========================================================================
    // Getting de_ij / dLambda
    
    stresstensor dQods = PS->dQods( EPS );
    //dQods.reportshort("dQods");

    tensor dQods_dev = dQods.deviator();
    tensor temp1 =  dQods_dev("ij")*dQods_dev("ij");
    double norm_dQods_dev =  pow( temp1.trace(), 0.5 );
    double temp2 = pow( 2.0 / 3.0, 0.5 ) * norm_dQods_dev;
    
    double ha = getha();
    double Cr = getCr();
    tensor alpha = EPS->getTensorVar(1);
    		     
    tensor h = (2.0/3.0) * ha * dQods + Cr * temp2 * alpha;

    return h;

}


//================================================================================
//  Print vars defined in Linear Evolution Law
//================================================================================
void EvolutionLaw_NL_Eij::print()
{
    cout << (*this);
}


//================================================================================
double EvolutionLaw_NL_Eij::getha() const
{       
    return ha;
}

//================================================================================
double EvolutionLaw_NL_Eij::getCr() const
{       
    return Cr;
}

#endif

