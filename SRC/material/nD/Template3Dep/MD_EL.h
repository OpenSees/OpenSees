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

#ifndef MD_EL_H
#define MD_EL_H

#include <math.h>

#include <iostream.h>
#include <iomanip.h>

#include "EL.h"

class MDEvolutionLaw : public EvolutionLaw
{
 // private to define the evolution law for Manzari-Dafalias critical state model
 public: 
 
    //Critical state parameters
    double Mc;
    double Me;
    double Lambda; // slope of e vs. ln p
    double ec_ref;// critical void ratio at reference mean effective stress p
    double p_ref; // critical void ratio at reference mean effective stress p

    // surface evolution parameters
    double kc_b;  // b stands for bounding surface, c compression
    double kc_d;  // d stands for dilatancy surface, c compression 
    double ke_b;  // d stands for bounding surface, e compression
    double ke_d;  // d stands for dilatancy surface, e compression 

    // Parameters to calculate plastic modulus 
    //double h;  // Could be calculated using ho and b_ij * n_ij
    double ho;
    double Cm;
    double eo;	 // initial void ratio

    //Dilatancy parameter
    //double D; //Moved to EPS's second scalar var
    double Ao;

    //Parameters to define the evolution of Fabric Tensor
    double Fmax;
    double Cf;
    double a;   //exponent for elastic modulus evolution  
    
 public:
    //MDEvolutionLaw( );                   // default constructor---no parameters
    
    MDEvolutionLaw( double Mcd = 1.14,   // default constructor
                    double Med = 1.14,
                    double Lamdad = 0.025,
                    double ec_refd = 0.8,
                    double p_refd = 100, //kPa
                    double kc_bd = 3.975,
                    double kc_dd = 4.200,
                    double ke_bd = 2.000,
                    double ke_dd = 0.070,
                    double hod = 1200.0,  // old 1200
                    double Cmd = 0.0,	 // old 0.0
                    //double eod = 0.65,
		    double Aod = 2.64,	 //old 2.64
                    double Fmaxd = 100,
                    double Cfd = 100, 
                    double ad = 0.6) : 
		    Mc (Mcd), Me(Med), Lambda(Lamdad), ec_ref(ec_refd), p_ref(p_refd),
		    kc_b(kc_bd), kc_d(kc_dd), ke_b(ke_bd), ke_d(ke_dd), ho(hod), Cm(Cmd),
		    Ao(Aod), Fmax(Fmaxd), Cf(Cfd), a(ad), eo(0.0) {} 
    
    MDEvolutionLaw(const MDEvolutionLaw &MDE );   // Copy constructor
    
    MDEvolutionLaw *newObj();                     //create a colne of itself
    
    void InitVars(EPState *EPS);    // Initialize all hardening vars called only once 
                                    // after material point is formed!
    
    void setInitD(EPState  *EPS);   // set initial D once current stress hits the y.s.
    
    double getKp( EPState *EPS, double dummy );     // calculating Kp 
    
    void UpdateAllVars( EPState *EPS, double dlamda );  // Evolve all vars
    //void UpdateAllTensorVar( EPState *EPS, double dlamda );  // Evolve all tensor vars

    void print();

    // some accessor functions
    double getMc() const;
    double getMe() const;
    double getLambda() const;
    double getec_ref() const;
    double getp_ref() const; 

    double getkc_b() const;  
    double getkc_d() const;  
    double getke_b() const;  
    double getke_d() const;  
    //double geth() const;  // Could be calculated using ho and b_ij * n_ij
    double getho() const;
    double getCm() const;
    double geteo() const;
    void seteo( double eod);

    //Dilatancy parameter
    double getAo() const;

    double getFmax() const;
    double getCf() const;
    double geta() const;


    //================================================================================
    // Overloaded Insertion Operator	  Zhaohui Added Aug. 13, 2000
    // prints Manzari-Dafalia EvolutionLaw's contents 
    //================================================================================
    friend ostream& operator<< (ostream& os, const MDEvolutionLaw & MDEL)
    {
        os.unsetf( ios::scientific );
        os.precision(5);

        //os.width(10);       
        os << endln << "Manzari-Dafalias Evolution Law's parameters:" << endln;
	os << "Mc = " << MDEL.getMc() << "; ";
        //os.width(10);       
	os << "Me = "  << MDEL.getMe() << "; ";
        //os.width(10);       
	os << "Lambda = " << MDEL.getLambda() << "; ";
        //os.width(10);       
	os << "ec_ref = " << MDEL.getec_ref() << "; ";
        //os.width(10);       
	os << "p_ref = " << MDEL.getp_ref() << "kPa"  << "; " << endln;

        //os.width(10);       
        os << "kc_b = " << MDEL.getkc_b() << "; ";
        //os.width(10);       
	os << "kc_d = " << MDEL.getkc_d() << "; ";
        //os.width(10);       
	os << "ke_b = " << MDEL.getke_b() << "; ";
        //os.width(10);       
	os << "ke_d = " << MDEL.getke_d() << "; " << endln;

        //os.width(10);       
        //os << "h = " << MDEL.h << "; ";
        //os.width(10);       
	os << "ho = " << MDEL.getho() << "; ";
        //os.width(10);       
	os << "Cm = " << MDEL.getCm() << "; " << "eo = " << MDEL.geteo() << endln;

        //os.width(10);       
        //os << "D = " << MDEL.getD() << "; ";
        //os.width(10);       
	os << "Ao = " << MDEL.getAo() << "; ";
        //os.width(10);       
	os << "Fmax = " << MDEL.getFmax() << "; ";
        //os.width(10);       
	os << "Cf = " << MDEL.getCf() << "; " << endln << endln; 
               
        return os;
    }  

    double g_A(double theta, double e);	 // Interpolation function  by Agyris
    double g_WW(double theta, double e); // Interpolation function  by Willan-Warkne

    
};


#endif

// test




