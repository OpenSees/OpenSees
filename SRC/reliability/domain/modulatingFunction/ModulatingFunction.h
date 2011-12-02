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
// $Date: 2007-02-24 01:41:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/modulatingFunction/ModulatingFunction.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef ModulatingFunction_h
#define ModulatingFunction_h

#include <Vector.h>
#include <Filter.h>
#include <ReliabilityDomainComponent.h>

class ModulatingFunction : public ReliabilityDomainComponent
{

public:
	ModulatingFunction(int tag, int classtag);
	virtual ~ModulatingFunction();

	virtual double getAmplitude(double time) = 0;
	virtual double getMaxAmplitude() = 0;
	virtual Filter *getFilter() = 0;

protected:

private:

};

#endif
