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
                                                                        
// $Revision: 1.4 $
// $Date: 2008-02-29 19:47:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/SensitivityIntegrator.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef SensitivityIntegrator_h
#define SensitivityIntegrator_h

#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>

class SensitivityIntegrator
{
public:

    SensitivityIntegrator();
    virtual ~SensitivityIntegrator();
    
	virtual int formSensitivityRHS(int gradNum) = 0;
	virtual int formIndependentSensitivityRHS() = 0;
	virtual int saveSensitivity   (const Vector &v, int gradNum, int numGrads) = 0;
    virtual int commitSensitivity (int gradNum, int numGrads) = 0;
	///////S added by K Fujimura ////////////////
	virtual int updateGradNumber(int passedGradNumber)=0;
	virtual int sensitivityDomainChanged(int NumGrads)=0;
	virtual bool staticSensitivity(void)=0;
	virtual bool NewSensitivity(void)=0;
	///////E added by K Fujimura ////////////////
protected:
    
private:
	AnalysisModel *theAnalysisModel;
};

#endif

