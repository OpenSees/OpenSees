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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/OpenSeesGradGEvaluator.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef OpenSeesGradGEvaluator_h
#define OpenSeesGradGEvaluator_h

#include <GradGEvaluator.h>
#include <Vector.h>
#include <Matrix.h>
#include <ReliabilityDomain.h>
#include <Domain.h>
#include <GFunEvaluator.h>
#include <tcl.h>

class SensitivityAlgorithm;

#include <fstream>
using std::ofstream;

class OpenSeesGradGEvaluator : public GradGEvaluator
{

public:
	OpenSeesGradGEvaluator(Tcl_Interp *passedTclInterp, 
			       GFunEvaluator *passedGFunEvaluator,
				   ReliabilityDomain *passedReliabilityDomain,
				   Domain *passedOpenSeesDomain,
			       SensitivityAlgorithm *theAlgo,
			       bool doGradientCheck);
	~OpenSeesGradGEvaluator();

	int		computeGradG(double gFunValue, const Vector &passed_x);
	int		computeAllGradG(const Vector &gFunValues, const Vector &passed_x);
	double	nodeGradient(int nodeNumber, int direction, int indx, char* dispOrWhat);
	double	elementGradient(int eleNumber, int indx, char* inString);

	Vector	getGradG();
	Matrix	getAllGradG();

	Matrix  getDgDdispl();

protected:

private:
	SensitivityAlgorithm *theSensAlgo;
	Domain *theOpenSeesDomain;
	Vector *grad_g;
	Matrix *grad_g_matrix;
	Matrix *DgDdispl;
	bool doGradientCheck;

};

#endif
