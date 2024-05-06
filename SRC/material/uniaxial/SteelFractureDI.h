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

// $Revision: 1.0 $
// $Date: 2021-05-18 00:17:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SteelFractureDI.cpp,v $

// Written: Francisco A. Galvis and Wen-Yi Yen
// Created: 05/2021
//
// Description: This file contains the class implementation of SteelFractureDI. 
// SteelFractureDI is based on the source code for Steel02
//


#ifndef SteelFractureDI_h
#define SteelFractureDI_h

#include <UniaxialMaterial.h>
#include <vector>

class SteelFractureDI : public UniaxialMaterial
{
public:
	SteelFractureDI(int tag,
		double Fy, double Fyc, double E0, double b,
		double R0, double cR1, double cR2,
		double a1, double a2, double a3, double a4, double sigcr, double m, double sigmin, double FI_lim);
	/*SteelFracture(int tag,
	double Fy, double E0, double b);*/

	// Constructor for no isotropic hardening
	//SteelFracture(int tag,
	//	double Fy, double Fyc, double E0, double b,
	//	double R0, double cR1, double cR2);

	//// Constructor for no isotropic hardening
	//// Also provides default values for R0, cR1, and cR2
	//SteelFracture(int tag, double Fy, double Fyc, double E0, double b);

	SteelFractureDI(void);
	virtual ~SteelFractureDI();


	const char *getClassType(void) const { return "SteelFractureDI"; };

	double getInitialTangent(void);
	UniaxialMaterial *getCopy(void);

	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);
	double getDI(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel,
		FEM_ObjectBroker &theBroker);

	void Print(OPS_Stream &s, int flag = 0);

	// counting function
	void calcDI(double sigcr, double m, double sigmin, double FI_lim, int &isStart, double sig, double &sigPDI, double &DI, double &slopeP, double& sumTenP, double& sumCompP);
	int returnSign(double v);

	// override get-response function
	Response *setResponse(const char **argv, int argc, OPS_Stream &s);
	int getResponse(int responseID, Information &matInformation);

protected:

private:
	// matpar : STEEL FIXED PROPERTIES
	double Fy;  //  = matpar(1)  : yield stress
	double FyC;  //  = matpar(1)  : yield stress
	double E0;  //  = matpar(2)  : initial stiffness
	double b;   //  = matpar(3)  : hardening ratio (Esh/E0)
	double R0;  //  = matpar(4)  : exp transition elastic-plastic
	double cR1; //  = matpar(5)  : coefficient for changing R0 to R
	double cR2; //  = matpar(6)  : coefficient for changing R0 to R
	double a1;  //  = matpar(7)  : coefficient for isotropic hardening in compression
	double a2;  //  = matpar(8)  : coefficient for isotropic hardening in compression
	double a3;  //  = matpar(9)  : coefficient for isotropic hardening in tension
	double a4;  //  = matpar(10) : coefficient for isotropic hardening in tension
				//double sigini; // initial
				// hstvP : STEEL HISTORY VARIABLES
	double epsminP; //  = hstvP(1) : max eps in compression
	double epsmaxP; //  = hstvP(2) : max eps in tension
	double epsplP;  //  = hstvP(3) : plastic excursion
	double epss0P;  //  = hstvP(4) : eps at asymptotes intersection
	double sigs0P;  //  = hstvP(5) : sig at asymptotes intersection
	double sig_1P;  //  = hstvP(5) : sig at asymptotes intersection
	double epssrP;  //  = hstvP(6) : eps at last inversion point
	double sigsrP;  //  = hstvP(7) : sig at last inversion point
	int    konP;    //  = hstvP(8) : index for loading/unloading
					// hstv : STEEL HISTORY VARIABLES   
	double epsP;  //  = strain at previous converged step
	double sigP;  //  = stress at previous converged step
	double eP;    //   stiffness modulus at last converged step;



	double epsmin;
	double epsmax;
	double epspl;
	double epss0;
	double sigs0;
	double sig_1;
	double epsr;
	double sigr;
	int    kon;
	double sig;
	double e;
	double eps;   //  = strain at current step

	// ************** added for fracture ***************
	double epsContP;
	double eps_0P;
	double eps_1P;
	double eps_rP;
	int konfP;
	int konCP;

	double epsCont;
	double eps_0;
	double eps_1;
	double eps_r;
	int konf;
	int konC;

	// ************** added for fracture ***************
	// ************** added for DI ********************
	double sigcr;	
	double m;
	double FI_lim;
	double sigmin;

	double DIP;
	int isStartP;
	double sigPDIP;
	double slopePP;
	double sumTenPP;
	double sumCompPP;

	double DI;
	int isStart;
	double sigPDI;
	double slopeP;
	double sumTenP;
	double sumCompP;

	// ************** added for DI ********************
};


#endif

