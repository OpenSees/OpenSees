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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/system/MVNcdf.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef MVNcdf_h
#define MVNcdf_h

#include <SystemAnalysis.h>
#include <ReliabilityDomain.h>
#include <FunctionEvaluator.h>
#include <Vector.h>
#include <Matrix.h>

class MVNcdf : public SystemAnalysis
{

public:
	MVNcdf(ReliabilityDomain*, 
           FunctionEvaluator*,
           TCL_Char*, int, TCL_Char*, TCL_Char*, 
           long int Nmax = 2e5, double errMax = 1.0e-6);
	~MVNcdf();

	int		analyze(void);

protected:

private:
	void	checkvals(long int, double);
	double	MVNcdffunc(const Vector&, const Matrix&, double);
	char fileName[256];
    int analysisType;
	long int Nmax;
	double errMax;
	
};

#endif

