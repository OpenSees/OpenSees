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
                                                                        
// $Revision: 1.5 $
// $Date: 2008-05-27 20:04:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/GradGEvaluator.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef GradGEvaluator_h
#define GradGEvaluator_h

#include <Vector.h>
#include <Matrix.h>
#include <ReliabilityDomain.h>
#include <GFunEvaluator.h>
#include <tcl.h>

class PerformanceFunctionCoefficientIter;

class GradGEvaluator
{

public:
	GradGEvaluator(ReliabilityDomain *theReliabilityDomain, GFunEvaluator *theGFunEvaluator, Tcl_Interp *theTclInterp);
	virtual ~GradGEvaluator();

	// Methods provided by the sub-classes
	virtual int		computeGradG(double gFunValue,
					     const Vector &passed_x) =0;
	virtual int		computeAllGradG(const Vector &gFunValues,
						const Vector &passed_x) =0;

	virtual Vector	getGradG() =0;
	virtual Matrix	getAllGradG() =0;

	// Methods that are provided by sub-classes, but rather specific to FE reliability
	virtual Matrix  getDgDdispl();
	virtual Matrix  getDgDpar();

    //////S added by K Fujimura ///////////
	int initializeNumberOfEvaluations();
	int getNumberOfEvaluations();
	virtual void setPerformFuncCoeffs(TaggedObjectStorage*);
	virtual void setPerformFuncCoeffIter(PerformanceFunctionCoefficientIter*);
	virtual bool getfinitedifference(){ return finitedifference; }
	virtual void setfinitedifference(bool fdf){finitedifference=fdf;}
    //////E added by K Fujimura ///////////

protected:
	int computeParameterDerivatives(double g);

	// one day these should find themselves a better home than protected
	ReliabilityDomain *theReliabilityDomain;
	GFunEvaluator *theGFunEvaluator;
	Tcl_Interp *theTclInterp;

	/////S added by K Fujimura /////
	int numberOfEvalIncSens;
	/////E added by K Fujimura /////

private:
	Matrix *DgDpar;

   /////S added by K Fujimura /////
	bool finitedifference;
	/////E added by K Fujimura ////

};

#endif
