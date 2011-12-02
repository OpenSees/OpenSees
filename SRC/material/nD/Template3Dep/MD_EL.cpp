/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           General platform for elaso-plastic constitutive model       #
#                    implementation                                              #
# CLASS:             MDEvolutionLaw (evolution law for Manzari-Dafalias Model)   #
#                                                                                #
# VERSION:                                                                       #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )    #
# TARGET OS:         DOS || UNIX || . . .                                        #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                 #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                                 #
#                                                                                #
#                                                                                #
# DATE:              08-03-2000                                                  #
# UPDATE HISTORY:                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
# SHORT EXPLANATION: The goal is to create a platform for efficient and easy     #
#                    implemetation of any elasto-plastic constitutive model!     #
#                                                                                #
//================================================================================
*/

#ifndef MD_EL_CPP
#define MD_EL_CPP

#include "MD_EL.h"
#include <basics.h>
    

//================================================================================
// Copy constructor
//================================================================================

MDEvolutionLaw::MDEvolutionLaw(const MDEvolutionLaw &MDE ) {

    this->Mc     = MDE.getMc();
    this->Me     = MDE.getMe();
    this->Lambda = MDE.getLambda();
    this->ec_ref = MDE.getec_ref();
    this->p_ref=  MDE.getp_ref();
    this->kc_b =  MDE.getkc_b();
    this->kc_d =  MDE.getkc_d();
    this->ke_b =  MDE.getke_b();
    this->ke_d =  MDE.getke_d();
    //this->h    =  MDE.geth();
    this->ho   =  MDE.getho();
    this->Cm   =  MDE.getCm();
    this->eo   =  MDE.geteo();
    //this->D    =  MDE.D;
    this->Ao   =  MDE.getAo();
    this->Fmax =  MDE.getFmax(); 
    this->Cf   =  MDE.getCf();
}


//================================================================================
//  Create a clone of itself 
//================================================================================
MDEvolutionLaw * MDEvolutionLaw::newObj() {
    
    MDEvolutionLaw *newEL = new MDEvolutionLaw( *this );
    
    return newEL;

}

//================================================================================
//  Initialize some  vars in EPState				        
//  1. E_Young,  2. e						       	
//							       		
//  Tensor vars store | 1. alpha,|  2. F  |  3. n              	        
//  Scalar vars store | 1. m,    |  2. D  |  3. void ratio(e)  	        
//							       	        
//================================================================================

void MDEvolutionLaw::InitVars(EPState  *EPS) {

    // set initial E_Young corresponding to current stress state
    double p_atm = 100.0; //Kpa atmospheric pressure
    double p = EPS->getStress().p_hydrostatic();
    double E = EPS->getEo() * pow( (p/p_atm), geta());
    EPS->setE( E );
      
    //=========================================================================
    //set initial void ratio from hardening scalar var  (third one)
    double eo = EPS->getScalarVar( 3 );
    seteo( eo );

    ////=========================================================================
    //// Update D 
    //   
    //double m = EPS->getScalarVar(1);
    //stresstensor F = EPS->getTensorVar( 2 );   // getting  F_ij from EPState
    //tensor temp_tensor = F("ij") * n("ij");
    //double temp = temp_tensor.trace();
    //
    //if (temp < 0)   temp = 0;
    //double A = getAo() * (1.0 + temp);
    //cout << " A = " << A << endlnn;
    //
    ////Calculating the lode angle theta
    //double J2_bar = r_bar.Jinvariant2();
    //double J3_bar = r_bar.Jinvariant3();
    //double theta = acos( 3.0*pow(3.0, 0.5)/2.0*J3_bar/ pow( J2_bar, 1.5) ) / 3.0;
    //cout << " theta = " << theta << endlnn;
    //
    //double c = getMe() / getMc();
    //double cd = getke_d() / getkc_d();
    //stresstensor alpha_theta_d = n("ij") * (g_WW(theta, c) * Mc + g_WW(theta, cd) * kc_d * xi - m);
    //cout << "alpha_theta_d " << alpha_theta_d<<" g_WW(theta, c) "<< g_WW(theta, c);
    //
    //temp_tensor = A*( alpha_theta_d - alpha );
    //temp_tensor.null_indices();
    //tensor temp1 = temp_tensor("ij") * n("ij");
    //double D = temp1.trace();
    //EPS->setScalarVar(2, D);
    
}   


//================================================================================
//  Set initial value of D once the current stress hit the yield surface     	
//      						       	        
//							       		
//  Tensor vars store | 1. alpha,|  2. F  |  3. n              	        
//  Scalar vars store | 1. m,    |  2. D  |  3. void ratio(e)  	        
//							       	        
//================================================================================

void MDEvolutionLaw::setInitD(EPState  *EPS) {

    //=========================================================================
    //calculate  n_ij
    stresstensor S = EPS->getStress().deviator();
    double p = EPS->getStress().p_hydrostatic();
    stresstensor alpha = EPS->getTensorVar( 1 );  // alpha_ij

    // Find the norm of alpha
    tensor norm_alphat = alpha("ij") * alpha("ij");
    double norm_alpha = sqrt( norm_alphat.trace() );
   
    stresstensor r = S * (1.0 / p);
    //r.reportshort("r");
    stresstensor r_bar = r - alpha;
    stresstensor norm2 = r_bar("ij") * r_bar("ij");
    double norm = sqrt( norm2.trace() );
    
    stresstensor n;
    if ( norm >= d_macheps() ){ 
      n = ( r  - alpha ) *(1.0 / norm );
    }
    else {
      ::printf(" \n\n n_ij not defined!!!! Program exits\n");
      exit(1);
    }

    //Calculate the state parameters xi 
    double e = EPS->getScalarVar(3);
    double ec = getec_ref() - getLambda() * log( p/getp_ref() );
    double xi = e - ec;

    //calculating A
    double m = EPS->getScalarVar(1);
    stresstensor F = EPS->getTensorVar( 2 );   // getting  F_ij from EPState
    tensor temp_tensor = F("ij") * n("ij");
    double temp = temp_tensor.trace();
    if (temp < 0)   temp = 0;
    double A = Ao*(1.0 + temp);

    //Calculating the lode angle theta
    double J2_bar = r_bar.Jinvariant2();
    double J3_bar = r_bar.Jinvariant3();
    double tempd = 3.0*pow(3.0, 0.5)/2.0*J3_bar/ pow( J2_bar, 1.5);

    if (tempd > 1.0 ) tempd = 1.0; //bug. if tempd = 1.00000000003, acos gives nan
    if (tempd < -1.0 ) tempd = -1.0;

    double theta = acos( tempd ) / 3.0;
    //cout << "Theta = " << theta << endlnn;
    
    //=========================================================================
    //calculate the alpha_theta_b and alpha_theta_d
    double c = getMe() / getMc();

    double cd = getke_d() / getkc_d();
    double alpha_theta_dd = (g_WW(theta, c) * Mc + g_WW(theta, cd) * kc_d * xi - m);
    stresstensor alpha_theta_d = n("ij") * alpha_theta_dd * pow(2.0/3.0, 0.5);
    //cout << "alpha_theta_b " << alpha_theta_b<<" g_WW(theta, c) "<< g_WW(theta, c) << endlnn;

    stresstensor d;
    d =  alpha_theta_d - alpha;
    d.null_indices();

    tensor temp1 = d("ij") * n("ij");
    temp1.null_indices();
    double D_new = temp1.trace() * A;
    //Check the restrictions on D
    if ( (xi > 0.0) && ( D_new < 0.0) )
       D_new = 0.0;  

    EPS->setScalarVar(2, D_new);  // Updating D
    
}   

//================================================================================
//  Updating all internal vars                                          
//							       		
//  Tensor vars store | 1. alpha,|  2. F                   	        
//  Scalar vars store | 1. m,    |  2. D  |  3. void ratio(e)  	        
//							       	        
//  always use non-updated parameters to update vars, don't mix it up!!! ?? n_ij new? old?
//================================================================================

void MDEvolutionLaw::UpdateAllVars( EPState *EPS, double dlamda) {
   
    //=========================================================================
    //calculate  n_ij
    stresstensor S = EPS->getStress().deviator();
    double p = EPS->getStress().p_hydrostatic();
    stresstensor alpha = EPS->getTensorVar( 1 );  // alpha_ij

    // Find the norm of alpha
    tensor norm_alphat = alpha("ij") * alpha("ij");
    double norm_alpha = sqrt( norm_alphat.trace() );
   
    stresstensor r = S * (1.0 / p);
    //r.reportshort("r");
    stresstensor r_bar = r - alpha;
    stresstensor norm2 = r_bar("ij") * r_bar("ij");
    double norm = sqrt( norm2.trace() );
    
    stresstensor n;
    if ( norm >= d_macheps() ){ 
      n = ( r  - alpha ) *(1.0 / norm );
    }
    else {
      ::printf(" \n\n n_ij not defined!!!! Program exits\n");
      exit(1);
    }
    //EPS->setTensorVar( 3, n); //update n_ij//

    // Update E_Young corresponding to current stress state
    double p_atm = 100.0; //Kpa, atmospheric pressure
    double E = EPS->getE();  // old E_Young
    double E_new = EPS->getEo() * pow( (p/p_atm), geta() ); 
    EPS->setE( E_new );

    // Update void ratio
    
    double e = EPS->getScalarVar(3);
    double D = EPS->getScalarVar(2);
    double elastic_strain_vol = EPS->getdElasticStrain().Iinvariant1();
    double plastic_strain_vol = EPS->getdPlasticStrain().Iinvariant1();

    double de_p = -( 1.0 + e ) * plastic_strain_vol; // plastic change of void ratio ?? e or eo?
    double de_e = -( 1.0 + e ) * elastic_strain_vol; // elastic change of void ratio ????
    cout << "get dPlasticStrain-vol" << plastic_strain_vol << endlnn;
    cout << "get dElasticStrain-vol" << elastic_strain_vol << endlnn;

    cout << "^^^^^^^^^^^ de_e = " << de_e << " de_p = " << de_p << endlnn; 
    double new_e = e + de_p + de_e;

    EPS->setScalarVar( 3, new_e ); // Updating e


    //Calculate the state parameters xi 
    double ec = getec_ref() - getLambda() * log( p/getp_ref() );

    double xi = e - ec;

    // Update D 
       
    double m = EPS->getScalarVar(1);
    stresstensor F = EPS->getTensorVar( 2 );   // getting  F_ij from EPState
    tensor temp_tensor = F("ij") * n("ij");
    double temp = temp_tensor.trace();
    if (temp < 0)   temp = 0;
    double A = Ao*(1.0 + temp);

    //Calculating the lode angle theta
    double J2_bar = r_bar.Jinvariant2();
    double J3_bar = r_bar.Jinvariant3();
    double tempd = 3.0*pow(3.0, 0.5)/2.0*J3_bar/ pow( J2_bar, 1.5);

    if (tempd > 1.0 ) tempd = 1.0; //bug. if tempd = 1.00000000003, acos gives nan
    if (tempd < -1.0 ) tempd = -1.0;

    double theta = acos( tempd ) / 3.0;
    //cout << "Theta = " << theta << endlnn;
    
    //=========================================================================
    //calculate the alpha_theta_b and alpha_theta_d
    double c = getMe() / getMc();

    double cd = getke_d() / getkc_d();
    double alpha_theta_dd = (g_WW(theta, c) * Mc + g_WW(theta, cd) * kc_d * xi - m);
    stresstensor alpha_theta_d = n("ij") * alpha_theta_dd * pow(2.0/3.0, 0.5);
    //cout << "alpha_theta_b " << alpha_theta_b<<" g_WW(theta, c) "<< g_WW(theta, c) << endlnn;

    double cb = getke_b() / getkc_b();
    if ( xi > 0 ) xi = 0.0;  // < -xi >
    double alpha_theta_bd = (g_WW(theta, c) * Mc + g_WW(theta, cb) * kc_b * (-xi) - m);
    stresstensor alpha_theta_b = n("ij") *alpha_theta_bd * pow(2.0/3.0, 0.5);
    alpha_theta_b.null_indices();

    stresstensor b;
    b =  alpha_theta_b - alpha;
    b.null_indices();
    stresstensor d;
    d =  alpha_theta_d - alpha;
    d.null_indices();

    tensor temp1 = d("ij") * n("ij");
    temp1.null_indices();
    double D_new = temp1.trace() * A;
    //Check the restrictions on D
    if ( (xi > 0.0) && ( D_new < 0.0) )
       D_new = 0.0;  

    EPS->setScalarVar(2, D_new);  // Updating D
    //EPS->setScalarVar(2, 0.0);  // Updating D
    //cout << "D= " << D << endlnn;
    
    //cout << "alpha_theta_d " << alpha_theta_d<<" g_WW(theta, c) "<< g_WW(theta, c);
 
    //=========================================================================
    // Update m
    double dm = dlamda * getCm() * ( 1.0 + e ) * D;
    EPS->setScalarVar(1, m + dm); // Updating m
    cout  << endlnn << "dm = " << dm << endlnn;

    //=========================================================================
    // Update alpha

    //calculate b_ref
    double alpha_c_b = g_WW(0.0, c) * Mc + g_WW(0.0, cb) * kc_b * (-xi) - m;
    double b_ref = 2.0 * pow(2.0/3.0, 0.5) * alpha_c_b;
    
    temp1 = b("ij") * n("ij");
    double bn = temp1.trace();
    cout << "xxxxxxxxxxxxxxxxxxx  bn " << bn << endlnn;
    //cout << "alternative " << alpha_theta_bd-norm_alpha << endlnn;


    double h = getho() * fabs(bn) / ( b_ref - fabs(bn) );
    //h = h + pow(2.0/3.0, 0.5) * getCm() * ( 1.0 + geteo() ) * A * bn;

    cout << " ||b|| " << (alpha_theta_bd - norm_alpha) << endlnn;
    cout << " dlamda " << dlamda << " h = " << h << endlnn;

    stresstensor dalpha;
    dalpha = dlamda * h * b("ij");
    //dalpha.null_indices();
    cout << "delta alpha =" << dalpha << endlnn;
    
    //dalpha.reportshortpqtheta("\n dalpha ");
    alpha = alpha + dalpha;
    alpha.null_indices();
    //alpha.reportshort("Alpha");
    EPS->setTensorVar(1, alpha);

    //=========================================================================
    // Update F
    stresstensor dF;
    if ( D > 0.0 ) D = 0.0;
    dF =  dlamda * getCf() * (-D) * ( getFmax() * n("ij") + F("ij") );
    //cout << "dF" << dF;
    
    F = F - dF;
    EPS->setTensorVar(2, F);

}


//
//

//================================================================================
// calculating Kp 
//================================================================================

double MDEvolutionLaw::getKp( EPState *EPS , double dummy) {

    //cout << "el-pl EPS: " <<  *EPS ;
    
    //=========================================================================
    //calculate  n_ij
    stresstensor S = EPS->getStress().deviator();
    double p = EPS->getStress().p_hydrostatic();
    stresstensor alpha = EPS->getTensorVar( 1 );  // alpha_ij
   
    stresstensor r = S * (1.0 / p);
    //r.reportshort("r");
    stresstensor r_bar = r - alpha;
    stresstensor norm2 = r_bar("ij") * r_bar("ij");
    double norm = sqrt( norm2.trace() );
    
    stresstensor n;
    if ( norm >= d_macheps() ){ 
      n = ( r  - alpha ) *(1.0 / norm );
    }
    else {
      ::printf(" \n\n n_ij not defined!!!! Program exits\n");
      exit(1);
    }
    //cout << "nij = " << n;
    
    //=========================================================================
    //calculating b_ij

    //Calculate the state parameters xi 
    double e = EPS->getScalarVar(3);
    double ec = getec_ref() - getLambda() * log( p/getp_ref() );
    double xi = e - ec;
    //cout << "ec = " << ec << endlnn;
    //cout << "xi = " << xi << endlnn;

    //Calculating the lode angle theta
    double J2_bar = r_bar.Jinvariant2();
    double J3_bar = r_bar.Jinvariant3();

    double tempd = 3.0*pow(3.0, 0.5)/2.0*J3_bar/ pow( J2_bar, 1.5);
    if (tempd > 1.0 ) tempd = 1.0; //bug. if tempd = 1.00000000003, acos gives nan
    if (tempd < -1.0 ) tempd = -1.0;
    double theta = acos( tempd ) / 3.0;
    //cout << "theta = " << theta << endlnn;
    
    //calculate the alpha_theta_b and alpha_theta_d
    double m = EPS->getScalarVar(1);
    double c = getMe() / getMc();

    double cd = getke_d() / getkc_d();
    stresstensor alpha_theta_d = n("ij") * (g_WW(theta, c) * Mc + g_WW(theta, cd) * kc_d * xi - m) * pow(2.0/3.0, 0.5);
    //cout << "alpha_theta_d " << alpha_theta_d<<" g_WW(theta, c) "<< g_WW(theta, c) << endlnn;

    double cb = getke_b() / getkc_b();
    if ( xi > 0.0 ) xi = 0.0;  // < -xi >
    stresstensor alpha_theta_b = n("ij") * (g_WW(theta, c) * Mc - g_WW(theta, cb) * kc_b * xi - m) * pow(2.0/3.0, 0.5);
    alpha_theta_b.null_indices();

    //=========================================================================
    // calculating h
    stresstensor b;
    b =  alpha_theta_b - alpha;
    b.null_indices();
    stresstensor d;
    d =  alpha_theta_d - alpha;
    d.null_indices();

    double alpha_c_b = g_WW(0.0, c) * Mc + g_WW(0.0, cb) * kc_b * (-xi) - m;
    double b_ref = 2.0 * pow(2.0/3.0, 0.5) * alpha_c_b;
    

    tensor temp1 = b("ij") * n("ij");
    double bn = temp1.trace();

    temp1 = d("ij") * n("ij");
    double dn = temp1.trace();
    //cout << "bn =" << bn << endlnn; 
    //cout << "dn =" << dn << endlnn; 


    // Calculating A
    stresstensor F = EPS->getTensorVar( 2 );   // getting  F_ij from EPState
    temp1 = F("ij") * n("ij");
    double temp = temp1.trace();
    if (temp < 0)   temp = 0;
    double A = Ao*(1.0 + temp);

    double h = getho() * fabs(bn) / ( b_ref - fabs(bn) ); 
    cout << "ho =" << getho()  << "   h =" << h << endlnn;

    //cout << "+++++++++++p= "<< p << " +++++++++++bn= " << bn << " ++++++++++h = " << h << endlnn;
    //cout << "+++++++++++b_ref= "<< b_ref << endlnn;

    //=========================================================================

    double Kp = h * bn + pow(2.0/3.0, 0.5) * getCm() * ( 1.0 + geteo() ) * A * dn;
    //double Kp = pow(2.0/3.0, 0.5) * getCm() * ( 1.0 + geteo() ) * A * dn;
    Kp = Kp * p;

    return Kp;

}


//void UpdateAllTensorVar( EPState *EPS ) = 0;  // Evolve all tensor vars

//================================================================================
//  Interpolation function No.1  -- Agyris: g_A(theta, e) 
//================================================================================

double MDEvolutionLaw::g_A(double theta, double e) {

    double temp = 2.0 * e;
    temp  = temp / ( (1.0 + e) - (1.0 - e) * cos(3.0*theta) ); 

    return temp;
}


//================================================================================
//  Interpolation function No.1  -- Willan-Warkne: g_WW(theta, e)
//================================================================================

double MDEvolutionLaw::g_WW(double theta, double e) {


    double g1 = 4.0*( 1.0 - e*e ) * cos(theta) * cos(theta) + pow(2.0*e-1.0, 2.0);
    double d1 = 2.0*( 1.0 - e*e ) * cos( theta );
    double d2 = ( 2.0*e-1.0 ) * pow( (4.0*(1.0-e*e)*cos(theta)*cos(theta) + 5.0*e*e - 4.0*e), 0.5);
    double temp =( d1 + d2 ) / g1; 

    return temp;
}


//================================================================================
//  Print vars defined in MD Evolution Law
//================================================================================
void MDEvolutionLaw::print()
{
    cout << " Manzari-Dafalias Evolution Law's Parameters" << endlnn;
    cout << (*this);
}

//================================================================================
double MDEvolutionLaw::getMc() const 
{
    return Mc;
}

//================================================================================
double MDEvolutionLaw::getMe() const 
{
    return Me;
}


//================================================================================
double MDEvolutionLaw::getLambda() const
{
    return Lambda;
}

//================================================================================
double MDEvolutionLaw::getec_ref() const
{
    return ec_ref;
}

//================================================================================
double MDEvolutionLaw::getp_ref() const 
{
    return p_ref;
}

//================================================================================
double MDEvolutionLaw::getkc_b() const
{
    return kc_b;
}

//================================================================================
double MDEvolutionLaw::getkc_d() const  
{
    return kc_d;
}

//================================================================================
double MDEvolutionLaw::getke_b() const
{
    return ke_b;
}

//================================================================================
double MDEvolutionLaw::getke_d() const
{
    return ke_d;
}

//================================================================================
//double MDEvolutionLaw::geth() const
//{       
//    return h;
//}

//================================================================================
double MDEvolutionLaw::getho() const
{      
    return ho;
}

//================================================================================
double MDEvolutionLaw::getCm() const
{       
    return Cm;
}

//================================================================================
double MDEvolutionLaw::geteo() const
{       
    return eo;
}

//================================================================================
void MDEvolutionLaw::seteo( double eod) 
{       
    eo = eod;
}

//================================================================================
double MDEvolutionLaw::getAo() const
{       
    return Ao;
}

//================================================================================
double MDEvolutionLaw::getFmax() const
{       
    return Fmax;
}

//================================================================================
double MDEvolutionLaw::getCf() const
{       
    return Cf;
}

//================================================================================
double MDEvolutionLaw::geta() const
{       
    return a;
}


#endif

