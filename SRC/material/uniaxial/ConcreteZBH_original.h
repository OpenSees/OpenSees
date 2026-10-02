/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//                                                                        
// Revision: 2.0
// Date: 09/2026
// Source: /OpenSees/SRC/material/uniaxial/ConcreteZBH_original.h
//
// Written: Zignago, D., and Barbato, M. University of California - Davis
// Edited: Badal, Prakash S. IIT Madras; and Barbato, M. University of California - Davis
// Created: 09/2026
//
// Description: This file contains the class implementation for the
// ConcreteZBH_original.
//
// References:
// Zignago, D., Barbato, M., and Hu, Dan (2018). "Constitutive Model of Concrete Simultaneously Confined by FRP and Steel for Finite-Element Analysis of FRP-Confined RC Columns." 
//   Journal of Composites for Construction. 22(6), 04018064.
// Spoelstra, M. R., and Giorgio Monti. (1999). "FRP-confined concrete model." 
//   Journal of Composites for Construction. 3(3), 143-150.

#ifndef ConcreteZBH_original_h
#define ConcreteZBH_original_h

#include <UniaxialMaterial.h>

class ConcreteZBH_original : public UniaxialMaterial
{
  public:
    ConcreteZBH_original(int tag, double _fc0, double _ec0, double _Ec,
		       double _Es, double _fy, double _eults, double _s, double _As_t,
		       double _Ef, double _eultf, double _tf, double _D, double _Ds,
		       double _As_l, double _kg_f, double _ks_s, double _ks_f, double _type_reinf, double _beta_input);
    
    ConcreteZBH_original(int tag, double _fc0, double _ec0, double _Ec,
		       double _Es, double _fy, double _eults, double _s, double _As_t,
		       double _Ef, double _eultf, double _tf, double _D, double _Ds,
		       double _As_l, double _kg_f, double _ks_s, double _ks_f, double _type_reinf);
    
    ConcreteZBH_original();

    ~ConcreteZBH_original();

    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);
    double getStress(void);
    double getTangent(void);

    double getInitialTangent(void) {return Ec;};

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    UniaxialMaterial *getCopy(void);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
		 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);

  protected:

  private:
    void Conf_Pressure (double eps, double flp, double &fc, double &fl, double &el);
	void envelope(double eps, double deps, double &fl, double &sig, double &Et, double &el);
    

    // matpar : Concrete FIXED PROPERTIES
    double fc0;   // concrete compression strength           
    double ec0;   // strain at compression strength         
    double Ec;    // concrete initial stiffness    
    double Es; 
    double fy;  
    double eults;   
    double s;
    double As_t;
    double Ef;
    double eultf;
    double tf;
    double D;
    double Ds;
    double As_l;
    double kg_f;
    double ks_s;
    double ks_f;
	double type_reinf;

	double roj_f;
	double roj_s;
	double roj_sl;
	double kg_s;
	double beta;
	double eccus;

    // hstvP : Concerete HISTORY VARIABLES last committed step
    double sigp;    //  = stress at previous converged step
    double Ep;      //   stiffness modulus at last converged step
    double elp;
    double epsp;    //  = strain at previous converged step
    double eminp;
    double eunl1p;
    double eunl2p;
    double eunl3p;
    double Eunlp;
    double Eunl2p;
    double Et3p;
    double sunlp;
    double flp;
    double flunlp;
    double elunlp;
    double muunlp;
    int flaggp;

    // hstv : Concerete HISTORY VARIABLES  current step
    double sig;    //  = stress at previous converged step
    double Et;      //   stiffness modulus at last converged step
    double el;
    double eps;    //  = strain at previous converged step
    double emin;
    double eunl1;
    double eunl2;
    double eunl3;
    double Eunl;
    double Eunl2;
    double Et3;
    double sunl;
    double fl;
    double flunl;
    double elunl;
    double muunl;
    int flagg;
};
#endif