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
                                                                        
// $Revision: 1.2 $
// $Date: 2007-10-31 15:39:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/system/PCM.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef PCM_h
#define PCM_h

#include <SystemAnalysis.h>
#include <ReliabilityDomain.h>
#include <FunctionEvaluator.h>
#include <Vector.h>
#include <Matrix.h>

class PCM : public SystemAnalysis
{

public:
	PCM(ReliabilityDomain *passedReliabilityDomain,
            FunctionEvaluator *passedEvaluator,
            TCL_Char *fileName, int analysisType, 
            TCL_Char *betaName, TCL_Char *rhoName);
	~PCM();

	int		analyze(void);

protected:

private:
	double	PCMfunc(const Vector&, const Matrix&, double);
	char fileName[256];
    int analysisType;
	
};

#endif

