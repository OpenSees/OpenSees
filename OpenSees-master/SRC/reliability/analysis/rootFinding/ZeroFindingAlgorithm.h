/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
**   Quan Gu (qgu@ucsd.edu)                                           **
**   Joel P. Conte (jpconte@ucsd.edu)                                 **
** ****************************************************************** */
                                                                        
 
//
// Written by  Quan Gu UCSD
//
 

#if !defined ZEROFINDINGAlgorithm_H 
#define  ZEROFINDINGAlgorithm_H

#include <OPS_Globals.h> 
# include <SamplingAnalysis.h>

#include "ReliabilityAnalysis.h"	// Added by ClassView
class ZeroFindingAlgorithm  
{
public:
	int getFunctionType();
	void setFunctionType(int);
	int get_ii_1();
	int get_ii_2();
	void setX1Pointer( double *);
	void setX2Pointer( double *);
	void setG1Pointer( double *);
	void setG2Pointer( double *);

	void set_ii_1(int);
	void set_ii_2(int);
	double getG2(int);
	double getX2(int);
	void saveXG2(double, double);
	void saveXG1(double, double);
	double getB();
	double getA();
	double getFb();
	double getFa();

	void setFEconvergence(bool);
	bool getFEconvergence(void);

	void setFb(double);
	void setFa(double);
	void setMaxIterNum(int);
	virtual int getIterationNumber();
	void setFunctionTolerance(double theTol){functionTol = theTol;};
	void setVariableTolerance(double theTol){variableTol = theTol;};

	double getZeroPoint();
	virtual int findZeroPoint(double)=0;
	virtual double myFunction(double)=0;

	void setA(double A);
	void setB(double B);	
	ZeroFindingAlgorithm();


	virtual ~ZeroFindingAlgorithm();

protected:
	int ii_1;  //data saved for x1 and G1
	int ii_2;  //data saved for x2 and G2
	// used to save histories
	double * x_1;
	double * x_2;
	double * G_1;
	double * G_2;

	double Fb;
	double Fa;
	int maxIterNum;
	int iterationNum;
	SamplingAnalysis * theSamplingAnalysis;

	int functionType;  // =1: samplingAnalysis. =2: upcrossing.  =3: others possible (Form ..)
	double functionTol;  // tol for function
	double variableTol;  // tol for x variables
	double zeroPoint;
	double b;
	double a;
	bool FEconvergence;
};

#endif  
