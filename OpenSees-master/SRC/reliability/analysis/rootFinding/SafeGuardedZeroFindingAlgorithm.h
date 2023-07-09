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

#if !defined SAFEGUARDEDZEROFINDINGALGORITHM_H__
#define      SAFEGUARDEDZEROFINDINGALGORITHM_H__

 
#include "ZeroFindingAlgorithm.h"

class SafeGuardedZeroFindingAlgorithm : public ZeroFindingAlgorithm  
{
public:

	SafeGuardedZeroFindingAlgorithm(SamplingAnalysis *);
	int findZeroPoint(double);
	double myFunction(double);
	SafeGuardedZeroFindingAlgorithm();
	virtual ~SafeGuardedZeroFindingAlgorithm();

private:
	bool isFbValid;


};

#endif // !defined(AFX_SAFEGUARDEDZEROFINDINGALGORITHM_H__16D7E60E_8D44_4C67_B080_AE880190EA16__INCLUDED_)
