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

// $Revision: 1.1 $
// $Date: 2008/12/09 20:00:16 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/ElasticPPcpp.h,v $

#ifndef ConcreteZBH_fitted_h
#define ConcreteZBH_fitted_h

// Written: fmk
//
// Description: This file contains the class definition for
// ElasticPPcpp. ElasticPPcpp provides the abstraction
// of an elastic perfectly plastic uniaxial material,
//
// What: "@(#) ElasticPPcpp.h, revA"

#include <UniaxialMaterial.h>

class ConcreteZBH_fitted : public UniaxialMaterial
{
  public:
    ConcreteZBH_fitted(int tag, double _fc0, double _ec0, double _Ec, double _fccs, double _eccs, double _rs,
		double _e1, double _e2, double _e3, double _e4, double _e5, double _e6, double _e7, double _e8, double _e9,
		double _eps_cy, double _eps_ccuf, double _sig_ccuf, double _eps_ccus, double _sig_ccus);
    ConcreteZBH_fitted();

    ~ConcreteZBH_fitted();

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
    void envelope (double eps, double deps, double &sig, double &Et, double &el);

    // matpar : Concrete FIXED PROPERTIES
    double fc0;   // concrete compression strength          
    double ec0;   // strain at compression strength          
    double Ec;    // concrete initial stiffness    
    double fccs; 
    double eccs;  
	double rs;
    double e1;
    double e2;
    double e3;
    double e4;
    double e5;
    double e6;
    double e7;
    double e8;
	double e9;
	double eps_cy;
	double eps_ccuf;
	double sig_ccuf;
	double eps_ccus;
	double sig_ccus;

	double beta;
	double r0;

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
    double elunl;
    double muunl;
    int flagg;
};


#endif