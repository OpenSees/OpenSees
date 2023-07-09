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
// $Date: 2008-03-13 22:29:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/curvature/FirstPrincipalCurvature.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef FirstPrincipalCurvature_h
#define FirstPrincipalCurvature_h

#include <FindCurvatures.h>
#include <Vector.h>
#include <FORMAnalysis.h>
#include <FunctionEvaluator.h>
#include <ReliabilityDomain.h>

class FirstPrincipalCurvature : public FindCurvatures
{

public:
    FirstPrincipalCurvature(ReliabilityDomain *passedReliabilityDomain,
                            FunctionEvaluator *passedGFunEvaluator,
                            FORMAnalysis *passedFORMAnalysis);
    ~FirstPrincipalCurvature();
  
    int		computeCurvatures();
    const Vector	&getCurvatures();
    const Vector	&getPrincipalAxes();

protected:

private:	
    ReliabilityDomain *theReliabilityDomain;
    FunctionEvaluator *theFunctionEvaluator;
    FORMAnalysis *theFORMAnalysis;
    
    Vector curvatures;
    Vector principalAxes;

};

#endif
