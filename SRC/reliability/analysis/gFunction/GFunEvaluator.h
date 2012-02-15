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
** ****************************************************************** */
                                                                        
// $Revision: 1.9 $
// $Date: 2010-09-13 21:39:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/GFunEvaluator.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef GFunEvaluator_h
#define GFunEvaluator_h

#include <Vector.h>
#include <ReliabilityDomain.h>
#include <Domain.h>
#include <tcl.h>
///// added buy K Fujimura /////
#include <GFunEachStepEvaluator.h>
#include <PerformanceFunctionCoefficientIter.h>

//#include <fstream>
//using std::ofstream;

class GFunEvaluator
{

public:
	GFunEvaluator();
	virtual ~GFunEvaluator();

	// Methods provided by base class
	virtual int evaluateG(const Vector &x) = 0;
	virtual double getG(void) = 0;
	virtual int runGFunAnalysis(const Vector &x) {return 0;}




	///// changed by K Fujimura /////
	//virtual int		evaluateGnoRecompute(const char* lsfExpression);
	int     initializeNumberOfEvaluations();
	int     getNumberOfEvaluations();
	

	// Methods implemented by SOME specific classes (random vibrations stuff)
	virtual void    setNsteps(int nsteps);
	virtual double  getDt();

	//////////////////////////////////////////////////////////
	//// added by K Fujimura /////////////////////////////////
	//////////////////////////////////////////////////////////
	virtual void activateSensitivty(void);
	virtual void inactivateSensitivty(void);
	virtual void setGFunEachStepEvaluator(GFunEachStepEvaluator *pGFunEachStepEvaluator);
	virtual void inactivateGFunEachStepEvaluator();
	virtual	void setThreshold(double value);
	virtual	double getThreshold();
	virtual Matrix* getEachStepResult();
	virtual Matrix* getEachStepConvFlag();
	virtual void setPerformFuncCoeffs(TaggedObjectStorage*);
	virtual void setPerformFuncCoeffIter(PerformanceFunctionCoefficientIter*);


	virtual double getG2(double g, double littleDt);


protected:
	//Tcl_Interp *theTclInterp;
	//ReliabilityDomain *theReliabilityDomain;
	//Domain *theOpenSeesDomain;

	int numberOfEvaluations;

private:
	

};

#endif

