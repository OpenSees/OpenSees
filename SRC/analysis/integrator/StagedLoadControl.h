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

// $Revision: 1.2 $
// $Date: 2003-02-14 23:00:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/StagedLoadControl.h,v $


#ifndef StagedLoadControl_h
#define StagedLoadControl_h

// File: ~/analysis/integrator/StagedLoadControl.h
//
// Written: fmk
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for StagedLoadControl.
// StagedLoadControl is an algorithmic class for performing a static analysis
// using a load control integration scheme.
//
// What: "@(#) StagedLoadControl.h, revA"

#include <LoadControl.h>
#include <StaticIntegrator.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;
class EquiSolnAlgo;
class ReliabilityDomain;

class StagedLoadControl : public LoadControl
{
public:
	StagedLoadControl();
	
    StagedLoadControl(double deltaLambda, int numIncr,
                double minLambda, double maxlambda);


    int  formTangent(int statusFlag = CURRENT_TANGENT);


};

#endif

