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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/gFunction/FunctionEvaluator.h,v $

//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#ifndef FunctionEvaluator_h
#define FunctionEvaluator_h

#include <Matrix.h>

/*
///// added by K Fujimura /////
#include <GFunEachStepEvaluator.h>
#include <PerformanceFunctionCoefficientIter.h>
*/

class FunctionEvaluator
{

public:
	FunctionEvaluator();
	virtual ~FunctionEvaluator();

	// Methods provided by base class
	int     initializeNumberOfEvaluations();
	int     getNumberOfEvaluations();
    int     incrementEvaluations(void) {numberOfEvaluations++; return 0;}
	
	// pure virtual
	virtual int setVariables(void) = 0;
	virtual int setExpression(const char *expression) = 0;
	virtual int addToExpression(const char *expression) = 0;
	virtual double evaluateExpression(void) = 0;
	virtual int runAnalysis(void) = 0;
	
	// MHS hack for reliability recorders ... set value in namespace
	virtual int setResponseVariable(const char *label, int lsfTag,
					int rvTag, double value);
	virtual int setResponseVariable(const char *label, int lsfTag, 
					double value);
	// and get value from namespace
	virtual double getResponseVariable(const char *label, int lsfTag,
					   int rvTag);
	virtual double getResponseVariable(const char *label, int lsfTag);

	/*
	// Methods in base class for evaluating LSF
	int		setTclRandomVariables(const Vector &x);
	int		uParse(char *tempchar, int *node, int *dirn, char* disp, char* varName, char* arrName);
	int		nodeParse(char *tempchar, int *node, int *dirn, char* disp, char* varName, char* arrName);
	int		elementParse(char *tempchar, int *element, char* varName, char* eleArgs);
	int		nodeTclVariable(int nodeNumber, int direction, char* dispOrWhat, char* varName, char* arrName);
	int		elementTclVariable(int eleNumber, char* varName, char* restString);

	// Methods to be implemented by specific classes
	virtual int		tokenizeSpecials(TCL_Char *theExpression, Tcl_Obj *paramList) = 0;
	 */

	// Methods implemented by SOME specific classes (random vibrations stuff)
	virtual void    setNsteps(int nsteps);
	virtual double  getDt();

	//////////////////////////////////////////////////////////
	//// added by K Fujimura /////////////////////////////////
	//////////////////////////////////////////////////////////
	/*
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
	*/

protected:

private:

	int numberOfEvaluations;	

};

#endif

