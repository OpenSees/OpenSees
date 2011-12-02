/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           General platform for elaso-plastic constitutive model       #
#                    implementation                                              #
#                                                                                #
# CLASS:             EvolutionLaw_NL_Eeq (on plastic equivalent strain)          #
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
# SHORT EXPLANATION: This is a nonlinear evolution law for the evoltion of a     #
#                    scalar variable k which depends on plastic equi. strain     #
#                    i.e. dk = f(e_eq)*de_eq                                     #
//================================================================================
*/

#ifndef EL_NLEeq_CPP
#define EL_NLEeq_CPP

#include "EL_NLEeq.h"
#include <basics.h>
    
//================================================================================
// Copy constructor
//================================================================================

EvolutionLaw_NL_Eeq::EvolutionLaw_NL_Eeq(const EvolutionLaw_NL_Eeq &NLE ) {

    this->eeqEtaPeak  = NLE.geteeqEtaPeak();
    this->etaResidual = NLE.getetaResidual();
    this->etaStart    = NLE.getetaStart();
    this->etaPeak     = NLE.getetaPeak();
    this->e           = NLE.gete();
    this->d           = NLE.getd();
}


//================================================================================
//  Create a clone of itself 
//================================================================================
EvolutionLaw_S * EvolutionLaw_NL_Eeq::newObj() {
    
    EvolutionLaw_NL_Eeq *newEL = new EvolutionLaw_NL_Eeq( *this );
    
    return newEL;

}

////================================================================================
////  Initialize some  vars in EPState				        
////  nothing						       	        
////================================================================================
//
//void EvolutionLaw_NL_Eeq::InitVars(EPState  *EPS) {
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
//  for NL model only    						       	
//							       		        
//							       	                
//================================================================================
//
//void EvolutionLaw_NL_Eeq::setInitD(EPState  *EPS) {
//
//}   

////================================================================================
////  Updating corresponding internal variable       done in CDriver                
////================================================================================
//
//void EvolutionLaw_NL_Eeq::UpdateVar( EPState *EPS, int WhichOne) {
//   
//    //=========================================================================
//    // Updating corresponding internal var by a nonliear function f
//    
//    // Getting e_eq
//    straintensor pstrain =  EPS->getPlasticStrain();
//    double e_eq  = pstrain.equivalent();
//    double e_eq2  = e_eq * e_eq;
//    //cout << "e_eq = " << e_eq << endln;
//
//    double upper1 = getd()*getetaResidual()*e_eq2;
//    double upper2 = (gete()*getetaPeak()+2.0*getd()*geteeqEtaPeak()*(-getetaResidual()+getetaPeak()))*e_eq;
//    double upper3 = getetaStart()*getd()*pow(geteeqEtaPeak(),2.0)*(-getetaResidual()+getetaPeak())/(getetaPeak()-getetaStart());
//    double upper = upper1 + upper2 + upper3;
//    
//    double lower = getd()*e_eq2 + gete()*e_eq +getd()*pow(geteeqEtaPeak(), 2.0)*(-getetaResidual()+getetaPeak())/(getetaPeak()-getetaStart());
//    
//    double new_S =  upper/lower;
//    
//    EPS->setScalarVar(WhichOne, new_S);
//
//}


 
//================================================================================
// Evaluating h_s ( for the evaluation of Kp )
//================================================================================
//double EvolutionLaw_NL_Eeq::h( EPState *EPS, double norm_dQods_dev ) {

double EvolutionLaw_NL_Eeq::h_s( EPState *EPS, PotentialSurface *PS)
{
    double h;

    //=========================================================================
    // Getting de_eq/dLambda
    //double  de_eq = EPS->getdPlasticStrain().equivalent();
    stresstensor dQods = PS->dQods( EPS );
    tensor dQods_dev = dQods.deviator();

    //Evaluate the norm of the deviator of dQods
    //temp1 =  dQods("ij")*dQods("ij");
    tensor temp1 =  dQods_dev("ij")*dQods_dev("ij");
    double norm_dQods_dev =  pow( temp1.trace(), 0.5 );
    double de_eqodL = pow( 2.0 / 3.0, 0.5 ) * norm_dQods_dev;
    
    //Evaluating dSodeeq
    double dSodeeq;
    straintensor pstrain =  EPS->getPlasticStrain(); //bug! should be total plastic strain's deviatoric part????
    straintensor pstrain_dev =  pstrain.deviator(); 

    //Checking equivalent()
    //Using deviatoric plastic strain to evolve the scalar internal var!
    //tensor temptx  = pstrain("ij") * pstrain("ij");
    //double tempdx = temptx.trace();
    //double e_eq  = pow( 2.0 * tempdx / 3.0, 0.5 );
    ////cout << e_eq << endln;

    tensor temptx  = pstrain_dev("ij") * pstrain_dev("ij");
    double tempdx = temptx.trace();
    double e_eq  = pow( 2.0 * tempdx / 3.0, 0.5 );
    
    
    //double e_eq  = pstrain.equivalent(); ????????????why this wouldn't work?
    double e_eq2  = e_eq * e_eq;

    double upper1 = -getd()*(getetaPeak() - getetaResidual())*(getetaPeak() - getetaStart())*( e_eq- geteeqEtaPeak());
    double upper21 = gete()*(getetaPeak() - getetaStart())* (e_eq + geteeqEtaPeak()); 
    double upper22 = 2.0*getd()* geteeqEtaPeak()* (-(getetaStart()*e_eq)-getetaResidual()* geteeqEtaPeak() + getetaPeak()* (e_eq+geteeqEtaPeak()));
    double upper = upper1 * (upper21 + upper22);

    double lower1 = gete()*(getetaPeak()-getetaStart())*e_eq;
    double lower2 = getd()*(getetaPeak()*(e_eq2+pow(geteeqEtaPeak(), 2.0) )-getetaStart()*e_eq2-getetaResidual()*pow(geteeqEtaPeak(), 2.0));
    double lower = pow(lower1 + lower2, 2.0);
    
    dSodeeq = upper / lower;
    
    h = dSodeeq * de_eqodL;

    // Drucker-Prager's evolution law
    // Get the current stress's I1
    //double I1 = EPS->getStress().Iinvariant1();

    //Von Mises 
    //double k = EPS->getScalarVar(1);
    //double Kp = -2.0 * k * geta() * temp;

    return h;

}


//================================================================================
//  Print vars defined in NLinear Evolution Law
//================================================================================
void EvolutionLaw_NL_Eeq::print()
{
    cout << (*this);
}


//================================================================================
double EvolutionLaw_NL_Eeq::geteeqEtaPeak() const
{
    return eeqEtaPeak;
}

//================================================================================
double EvolutionLaw_NL_Eeq::getetaResidual() const
{
    return etaResidual;
}

//================================================================================
double EvolutionLaw_NL_Eeq::getetaStart() const
{
    return etaStart;
}

//================================================================================
double EvolutionLaw_NL_Eeq::getetaPeak() const
{
    return etaPeak;
}

//================================================================================
double EvolutionLaw_NL_Eeq::gete() const
{
    return e;
}

//================================================================================
double EvolutionLaw_NL_Eeq::getd() const
{
    return d;
}


#endif

